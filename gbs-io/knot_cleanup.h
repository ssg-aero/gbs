#pragma once

/**
 * @file knot_cleanup.h
 * @brief Removal, at export, of the knot spans of a B-spline shorter than a tolerance.
 *
 * Design: docs/sources/design/export_knot_cleanup.md.
 *
 * A span much shorter than the resolution of a CAD system (a knot one ulp off another, issue
 * #213) is valid for gbs, but a reader that merges knots closer than its own tolerance sees a
 * multiplicity above the degree + 1 and drops the entity. Each reader has its own settings: the
 * exporters (STEP, IGES) remove, on request, the spans whose image is shorter than a given length;
 * the algorithms of gbs keep their knots (a cut within knot_eps of a knot is made at the knot,
 * gbs::snap_to_knot, so that they do not make such spans themselves).
 *
 * A span [a, b] is removed by moving one of its knots onto the other (b onto a, or a onto b; an end
 * of the domain never moves), the poles unchanged. The basis functions whose support holds the
 * moved knot change, so the curve moves on the p + 1 spans on each side of it, by little for a
 * short span; the parameters keep their meaning, and the pcurves on a surface stay valid within the
 * deviation. At an end of the domain, the poles whose basis functions became null are dropped
 * (exact); inside, the multiplicity reached is lowered back to that of the two knots by knot removal
 * (The NURBS Book A5.8). The result is compared with the original on the window where it moved,
 * both refined to the union of their knots: the largest distance between their poles bounds the
 * distance between the two curves at the same parameter (convex hull; for a rational with positive
 * weights, the bound of deviation()). A removal beyond the tolerance is undone and the span kept.
 *
 * A surface is cleaned along u then v, all the rows of poles together. A pcurve is cleaned in the
 * parameters of its surface scaled by the speeds of the surface along it, so that a distance there
 * stands for a length on the surface. A knot vector whose ends are not of multiplicity degree + 1
 * (not clamped) is left as it is. The bounds are kept.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

#include <gbs/bscurve.h>
#include <gbs/bssurf.h>

namespace gbs::nurbs_export
{
    /// What a knot cleanup did.
    template <std::floating_point T>
    struct KnotCleanupStats
    {
        std::size_t cleaned_curves{}, cleaned_pcurves{}, cleaned_surfaces{}; ///< 3D curves, pcurves and surfaces with at least one span removed
        std::size_t removed_spans{};                                         ///< short spans removed
        std::size_t kept_spans{};                                            ///< short spans kept: their removal deviated beyond the tolerance
        T max_deviation{}; ///< largest deviation of a cleaned curve, surface or pcurve (on its surface) from the original
    };

    namespace knot_detail
    {
        template <typename T, std::size_t D>
        using Row = std::vector<std::array<T, D>>;

        /// The distinct values of a flat knot vector.
        template <typename T>
        std::vector<T> distinct(std::vector<T> k)
        {
            k.erase(std::ranges::unique(k).begin(), k.end());
            return k;
        }

        /// Whether both ends of a knot vector are of multiplicity p + 1, as the cleanup assumes.
        template <typename T>
        bool clamped(const std::vector<T> &k, std::size_t p)
        {
            return k.size() > 2 * p + 1 && std::size_t(multiplicity(k.front(), k)) == p + 1 && std::size_t(multiplicity(k.back(), k)) == p + 1;
        }

        /// Distance between two poles, Cartesian (the weight, last coordinate, divided out when rational).
        template <typename T, std::size_t D>
        T pole_distance(const std::array<T, D> &a, const std::array<T, D> &b, bool rational)
        {
            T s{};
            const std::size_t n = rational ? D - 1 : D;
            for (std::size_t i = 0; i < n; ++i)
            {
                const T d = rational ? a[i] / a[D - 1] - b[i] / b[D - 1] : a[i] - b[i];
                s += d * d;
            }
            return std::sqrt(s);
        }

        /// Diameter of the points of the rows on [a, b], sampled at the ends and quarters.
        template <typename T, std::size_t D>
        T span_size(const std::vector<T> &k, const std::vector<Row<T, D>> &rows, std::size_t p, T a, T b, bool rational)
        {
            T size{};
            for (const auto &row : rows)
            {
                std::array<std::array<T, D>, 5> q;
                for (std::size_t i = 0; i < 5; ++i)
                    q[i] = eval_value_decasteljau<T, D>(i == 4 ? b : a + (b - a) * T(i) / T(4), k, row, p);
                for (std::size_t i = 0; i < 5; ++i)
                    for (std::size_t j = i + 1; j < 5; ++j)
                        size = std::max(size, pole_distance(q[i], q[j], rational));
            }
            return size;
        }

        /// (k, rows) on [w1, w2], two of its knots: the knots and poles of the spans between them, clamped at both.
        template <typename T, std::size_t D>
        auto window(const std::vector<T> &k, const std::vector<Row<T, D>> &rows, std::size_t p, T w1, T w2)
            -> std::pair<std::vector<T>, std::vector<Row<T, D>>>
        {
            // the spans [w1, w2] use the knots from the last w1 - p to the first w2 + p, the poles from the last w1 - p
            const auto a = static_cast<std::size_t>(std::ranges::upper_bound(k, w1) - k.begin()) - 1;
            const auto b = static_cast<std::size_t>(std::ranges::lower_bound(k, w2) - k.begin());
            const std::vector<T> kw(k.begin() + std::ptrdiff_t(a - p), k.begin() + std::ptrdiff_t(b + p + 1));
            std::vector<T> kk;
            std::vector<Row<T, D>> rw;
            for (const auto &row : rows)
            {
                Row<T, D> r(row.begin() + std::ptrdiff_t(a - p), row.begin() + std::ptrdiff_t(b));
                kk = kw;
                gbs::trim(p, kk, r, w1, w2);
                rw.push_back(std::move(r));
            }
            return {std::move(kk), std::move(rw)};
        }

        /**
         * Bound of the distance at the same parameter between (k1, rows1) and (k2, rows2), both refined
         * to the union of their knots: the largest distance between their poles (convex hull). For
         * rationals with positive weights, about the centre c of the poles of the second:
         * |C1 - C2| <= (max |w1 (P1 - c) - w2 (P2 - c)| + max |P2 - c| max |w1 - w2|) / min w1, from
         * C1 - C2 = (A1 - A2) / W1 + (C2 - c) (W2 - W1) / W1 with A = sum N w (P - c). Infinity when they
         * cannot be refined to the same knots (a multiplicity above the degree to reach), or a weight
         * is not positive.
         */
        template <typename T, std::size_t D>
        T deviation(std::vector<T> k1, std::vector<Row<T, D>> rows1, std::vector<T> k2, std::vector<Row<T, D>> rows2, std::size_t p, bool rational)
        {
            auto raise = [p](std::vector<T> &k, std::vector<Row<T, D>> &rows, T u, std::size_t r) {
                std::vector<T> kk;
                for (auto &row : rows)
                {
                    kk = k;
                    insert_knots(u, p, r, kk, row);
                }
                k = std::move(kk);
            };
            auto all = k1;
            all.insert(all.end(), k2.begin(), k2.end());
            std::ranges::sort(all);
            for (auto u : distinct(all))
            {
                const auto m1 = std::size_t(multiplicity(u, k1)), m2 = std::size_t(multiplicity(u, k2));
                if (m1 < m2)
                    raise(k1, rows1, u, m2 - m1);
                else if (m2 < m1)
                    raise(k2, rows2, u, m1 - m2);
            }
            constexpr T infinity = std::numeric_limits<T>::infinity();
            if (k1 != k2 || rows1.size() != rows2.size())
                return infinity;
            for (std::size_t r = 0; r < rows1.size(); ++r)
                if (rows1[r].size() != rows2[r].size())
                    return infinity;
            T d{};
            if (!rational)
            {
                for (std::size_t r = 0; r < rows1.size(); ++r)
                    for (std::size_t i = 0; i < rows1[r].size(); ++i)
                        d = std::max(d, pole_distance(rows1[r][i], rows2[r][i], false));
                return d;
            }
            if constexpr (D >= 2)
            {
                constexpr std::size_t n = D - 1; // the weight
                std::array<T, n> lo, hi;
                lo.fill(infinity);
                hi.fill(-infinity);
                for (const auto &row : rows2)
                    for (const auto &q : row)
                    {
                        if (!(q[n] > T(0)))
                            return infinity;
                        for (std::size_t c = 0; c < n; ++c)
                            lo[c] = std::min(lo[c], q[c] / q[n]), hi[c] = std::max(hi[c], q[c] / q[n]);
                    }
                std::array<T, n> centre;
                for (std::size_t c = 0; c < n; ++c)
                    centre[c] = (lo[c] + hi[c]) / 2;
                T d_a{}, d_w{}, radius{}, w_min = infinity;
                for (std::size_t r = 0; r < rows1.size(); ++r)
                    for (std::size_t i = 0; i < rows1[r].size(); ++i)
                    {
                        const auto &q1 = rows1[r][i], &q2 = rows2[r][i];
                        if (!(q1[n] > T(0)))
                            return infinity;
                        T s_a{}, s_r{};
                        for (std::size_t c = 0; c < n; ++c)
                        {
                            const T da = (q1[c] - q1[n] * centre[c]) - (q2[c] - q2[n] * centre[c]), dr = q2[c] / q2[n] - centre[c];
                            s_a += da * da, s_r += dr * dr;
                        }
                        d_a = std::max(d_a, std::sqrt(s_a));
                        radius = std::max(radius, std::sqrt(s_r));
                        d_w = std::max(d_w, std::abs(q1[n] - q2[n]));
                        w_min = std::min(w_min, q1[n]);
                    }
                return (d_a + radius * d_w) / w_min;
            }
            return infinity;
        }

        /// Removes the knot u up to `num` times from every row, as many times as all the rows allow within tol.
        template <typename T, std::size_t D>
        std::size_t remove_from_rows(T u, std::size_t p, std::size_t num, std::vector<T> &k, std::vector<Row<T, D>> &rows, T tol)
        {
            std::size_t t = num;
            for (const auto &row : rows)
            {
                auto kk = k;
                auto P = row;
                t = std::min(t, static_cast<std::size_t>(remove_knot(u, p, t, kk, P, tol)));
                if (t == 0)
                    return 0;
            }
            std::vector<T> kk;
            for (auto &row : rows)
            {
                kk = k;
                remove_knot(u, p, t, kk, row, std::numeric_limits<T>::infinity());
            }
            k = std::move(kk);
            return t;
        }

        /**
         * The span (from, to) removed by moving the knot `from` onto `to`: the null poles dropped at an end,
         * the multiplicity lowered back inside; nullopt when it stays above what the curve had, or an end
         * is not at p + 1.
         */
        template <typename T, std::size_t D>
        auto collapse(std::vector<T> k, std::vector<Row<T, D>> rows, std::size_t p, T from, T to, T tol)
            -> std::optional<std::pair<std::vector<T>, std::vector<Row<T, D>>>>
        {
            const auto m_from = std::size_t(multiplicity(from, k)), m_to = std::size_t(multiplicity(to, k));
            const bool at_start = to == k.front(), at_end = to == k.back();
            if ((at_start || at_end) && m_to != p + 1)
                return std::nullopt; // the ends of a knot vector that is not clamped
            std::ranges::replace(k, from, to);
            const auto m = m_from + m_to;
            if (at_start || at_end)
            {
                // p + 1 + e knots at an end: the e first (last) basis functions are null on the domain
                const auto e = m - (p + 1);
                if (at_start)
                {
                    k.erase(k.begin(), k.begin() + static_cast<std::ptrdiff_t>(e));
                    for (auto &row : rows)
                        row.erase(row.begin(), row.begin() + static_cast<std::ptrdiff_t>(e));
                }
                else
                {
                    k.resize(k.size() - e);
                    for (auto &row : rows)
                        row.resize(row.size() - e);
                }
                return std::pair{std::move(k), std::move(rows)};
            }
            const auto target = std::max(m_from, m_to);
            const auto removed = remove_from_rows(to, p, m - target, k, rows, tol);
            if (m - removed > std::max(p, target))
                return std::nullopt; // a discontinuity the curve did not have
            return std::pair{std::move(k), std::move(rows)};
        }

        /**
         * Removes the spans of (k, rows) whose points are within tol of each other, while the result stays
         * within tol of the original. Returns the deviation from the original (a bound), nullopt if nothing
         * changed.
         *
         * A move changes the curve on the window of the p + 1 knots on each side of the span (the supports of
         * the basis functions on its knots): only there is it compared with the original, the rest being
         * the curve already checked, and only from there are the spans looked at again.
         */
        template <typename T, std::size_t D>
        auto clean_rows(std::vector<T> &k, std::vector<Row<T, D>> &rows, std::size_t p, T tol, bool rational, KnotCleanupStats<T> &st)
            -> std::optional<T>
        {
            if (!clamped(k, p))
                return std::nullopt;
            const auto k0 = k;
            const auto rows0 = rows;
            std::vector<std::pair<T, T>> kept;
            std::optional<T> dev;
            std::size_t first = 0; // the spans before it were found long and have not changed since
            for (;;)
            {
                const auto x = distinct(k);
                if (x.size() < 3)
                    break; // a single span: nothing to merge it with
                std::optional<std::size_t> found;
                for (std::size_t i = first; i + 1 < x.size() && !found; ++i)
                    if (std::ranges::find(kept, std::pair{x[i], x[i + 1]}) == kept.end() &&
                        span_size<T, D>(k, rows, p, x[i], x[i + 1], rational) <= tol)
                        found = i;
                if (!found)
                    break;
                const auto i = *found;
                const T a = x[i], b = x[i + 1];
                const auto fa = static_cast<std::size_t>(std::ranges::lower_bound(k, a) - k.begin());
                const auto lb = static_cast<std::size_t>(std::ranges::upper_bound(k, b) - k.begin()) - 1;
                const T w1 = k[fa > p + 1 ? fa - p - 1 : 0], w2 = k[std::min(k.size() - 1, lb + p + 1)];
                const auto [k0w, rows0w] = window<T, D>(k0, rows0, p, w1, w2);
                // an end of the domain stays; inside, b onto a first
                std::vector<std::pair<T, T>> moves;
                if (i == 0)
                    moves = {{b, a}};
                else if (i + 2 == x.size())
                    moves = {{a, b}};
                else
                    moves = {{b, a}, {a, b}};
                bool done = false;
                for (const auto &[from, to] : moves)
                {
                    auto c = collapse<T, D>(k, rows, p, from, to, tol);
                    if (!c)
                        continue;
                    const auto [kcw, rowscw] = window<T, D>(c->first, c->second, p, w1, w2);
                    const T d = deviation<T, D>(k0w, rows0w, kcw, rowscw, p, rational);
                    if (d <= tol)
                    {
                        k = std::move(c->first);
                        rows = std::move(c->second);
                        dev = std::max(dev.value_or(T(0)), d);
                        ++st.removed_spans;
                        done = true;
                        break;
                    }
                }
                if (done)
                    first = i > p + 1 ? i - p - 1 : 0;
                else
                {
                    kept.emplace_back(a, b);
                    ++st.kept_spans;
                    first = i + 1;
                }
            }
            return dev;
        }

        /// The rows of a surface along u (one per v index) or along v (one per u index), u fastest in `poles`.
        template <typename T, std::size_t D>
        std::vector<Row<T, D>> surface_rows(const std::vector<std::array<T, D>> &poles, std::size_t nu, std::size_t nv, bool along_u)
        {
            std::vector<Row<T, D>> rows(along_u ? nv : nu);
            for (std::size_t j = 0; j < nv; ++j)
                for (std::size_t i = 0; i < nu; ++i)
                    rows[along_u ? j : i].push_back(poles[i + nu * j]);
            return rows;
        }

        template <typename T, std::size_t D>
        std::vector<std::array<T, D>> surface_poles(const std::vector<Row<T, D>> &rows, bool along_u)
        {
            const std::size_t n_rows = rows.size(), n = rows.front().size();
            const std::size_t nu = along_u ? n : n_rows, nv = along_u ? n_rows : n;
            std::vector<std::array<T, D>> poles(nu * nv);
            for (std::size_t j = 0; j < nv; ++j)
                for (std::size_t i = 0; i < nu; ++i)
                    poles[i + nu * j] = along_u ? rows[j][i] : rows[i][j];
            return poles;
        }

        /// c rebuilt on the cleaned knots and poles, its bounds as they were.
        template <typename T, std::size_t dim, bool rational>
        void rebuild(BSCurveGeneral<T, dim, rational> &c, const std::vector<T> &k, const Row<T, dim + rational> &poles)
        {
            const auto [u1, u2] = c.bounds();
            using Curve = std::conditional_t<rational, BSCurveRational<T, dim>, BSCurve<T, dim>>;
            Curve cleaned(poles, k, c.degree());
            cleaned.trim(u1, u2, false);
            c = cleaned;
        }
    }

    /**
     * @brief Removes the knot spans of a B-spline curve whose image is shorter than tol, the curve
     * staying within tol of the original at the same parameter (see the file comment).
     * @return true if the curve changed
     */
    template <typename T, std::size_t dim, bool rational>
    bool clean_knots(BSCurveGeneral<T, dim, rational> &c, T tol, KnotCleanupStats<T> &st)
    {
        if (!(tol > T(0)))
            return false;
        auto k = c.knotsFlats();
        std::vector<knot_detail::Row<T, dim + rational>> rows{c.poles()};
        const auto dev = knot_detail::clean_rows(k, rows, c.degree(), tol, rational, st);
        if (!dev)
            return false;
        knot_detail::rebuild(c, k, rows.front());
        ++st.cleaned_curves;
        st.max_deviation = std::max(st.max_deviation, *dev);
        return true;
    }

    /**
     * @brief Removes the knot spans of a pcurve whose image on the surface s is shorter than tol, that
     * image staying within tol of the original at the same parameter.
     *
     * The parameters of the surface are scaled by sqrt 2 times the largest |S_u| and |S_v| along the
     * pcurve (at its knots and the quarters of its spans): a change (du, dv) moves the point of the
     * surface by at most |S_u| |du| + |S_v| |dv|, which the distance between the scaled parameters
     * bounds. The speeds being sampled, the deviation is an estimate. A pcurve along which the surface
     * does not move in one direction (on a pole), or where its derivatives cannot be evaluated, is
     * left as it is.
     * @return true if the pcurve changed
     */
    template <typename T, bool rational>
    bool clean_knots(BSCurveGeneral<T, 2, rational> &pc, const Surface<T, 3> &s, T tol, KnotCleanupStats<T> &st)
    {
        if (!(tol > T(0)))
            return false;
        auto k = pc.knotsFlats();
        std::vector<knot_detail::Row<T, 2 + rational>> rows{pc.poles()};
        const auto [u1, u2, v1, v2] = s.bounds();
        T su{}, sv{};
        const auto x = knot_detail::distinct(k);
        try
        {
            for (std::size_t i = 0; i + 1 < x.size(); ++i)
                for (std::size_t q = 0; q < 5; ++q)
                {
                    const auto h = eval_value_decasteljau<T, 2 + rational>(x[i] + (x[i + 1] - x[i]) * T(q) / T(4), k, rows.front(), pc.degree());
                    const T w = rational ? h[2] : T(1);
                    const T u = std::clamp(h[0] / w, u1, u2), v = std::clamp(h[1] / w, v1, v2);
                    su = std::max(su, norm(s.value(u, v, 1, 0)));
                    sv = std::max(sv, norm(s.value(u, v, 0, 1)));
                }
        }
        catch (const std::exception &)
        {
            return false; // a surface without derivatives there: the pcurve as it is
        }
        if (!(su > T(0)) || !(sv > T(0)) || !std::isfinite(su) || !std::isfinite(sv))
            return false;
        const T fu = std::numbers::sqrt2_v<T> * su, fv = std::numbers::sqrt2_v<T> * sv;
        for (auto &h : rows.front())
            h[0] *= fu, h[1] *= fv;
        const auto dev = knot_detail::clean_rows(k, rows, pc.degree(), tol, rational, st);
        if (!dev)
            return false;
        for (auto &h : rows.front())
            h[0] /= fu, h[1] /= fv;
        knot_detail::rebuild(pc, k, rows.front());
        ++st.cleaned_pcurves;
        st.max_deviation = std::max(st.max_deviation, *dev);
        return true;
    }

    /**
     * @brief Removes the knot spans of a B-spline surface, along u then v, whose image is shorter
     * than tol, the surface staying within tol of the original at the same parameters.
     * @return true if the surface changed
     */
    template <typename T, std::size_t dim, bool rational>
    bool clean_knots(BSSurfaceGeneral<T, dim, rational> &s, T tol, KnotCleanupStats<T> &st)
    {
        if (!(tol > T(0)))
            return false;
        constexpr std::size_t D = dim + rational;
        auto ku = s.knotsFlatsU(), kv = s.knotsFlatsV();
        auto poles = s.poles();
        std::size_t nu = s.nPolesU(), nv = s.nPolesV();
        // the tolerance is shared by the two directions: their deviations add up
        T dev{};
        bool changed = false;
        for (bool along_u : {true, false})
        {
            auto rows = knot_detail::surface_rows<T, D>(poles, nu, nv, along_u);
            auto &k = along_u ? ku : kv;
            const auto d = knot_detail::clean_rows(k, rows, along_u ? s.degreeU() : s.degreeV(), tol - dev, rational, st);
            if (!d)
                continue;
            poles = knot_detail::surface_poles<T, D>(rows, along_u);
            (along_u ? nu : nv) = rows.front().size();
            dev += *d;
            changed = true;
        }
        if (!changed)
            return false;
        using Surface = std::conditional_t<rational, BSSurfaceRational<T, dim>, BSSurface<T, dim>>;
        static_cast<BSSurfaceGeneral<T, dim, rational> &>(s) = Surface(poles, ku, kv, s.degreeU(), s.degreeV());
        ++st.cleaned_surfaces;
        st.max_deviation = std::max(st.max_deviation, dev);
        return true;
    }
}
