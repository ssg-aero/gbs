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
 * exporters (STEP, IGES) remove, on request, the spans whose image is shorter than a given length,
 * the algorithms of gbs keep their knots.
 *
 * A span [a, b] is removed by moving one of its knots onto the other (b onto a, or a onto b; an end
 * of the domain never moves), the poles unchanged: the parameters outside [a, b] keep their points.
 * At an end of the domain, the poles whose basis functions became null are dropped (exact); inside,
 * the multiplicity reached is lowered back to that of the two knots by knot removal (The NURBS Book
 * A5.8). The result is compared with the original by refining both to the union of their knots: the
 * largest distance between their poles bounds the distance between the two curves at the same
 * parameter (convex hull; for a rational, the distance between the Cartesian poles, an estimate). A
 * removal beyond the tolerance is undone and the span kept.
 *
 * A surface is cleaned along u then v, all the rows of poles together. The bounds are kept.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
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
        std::size_t cleaned_curves{}, cleaned_surfaces{}; ///< curves and surfaces with at least one span removed
        std::size_t removed_spans{};                      ///< short spans removed
        std::size_t kept_spans{};                         ///< short spans kept: their removal deviated beyond the tolerance
        T max_deviation{};                                ///< largest deviation of a cleaned curve or surface from the original
    };

    namespace knot_detail
    {
        template <typename T, std::size_t D>
        using Row = std::vector<std::array<T, D>>;

        template <typename T>
        std::size_t mult(const std::vector<T> &k, T u)
        {
            return static_cast<std::size_t>(std::upper_bound(k.begin(), k.end(), u) - std::lower_bound(k.begin(), k.end(), u));
        }

        template <typename T>
        std::vector<T> distinct(const std::vector<T> &k)
        {
            std::vector<T> x;
            for (auto u : k)
                if (x.empty() || u != x.back())
                    x.push_back(u);
            return x;
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

        /**
         * Largest distance between the poles of (k1, rows1) and (k2, rows2), both refined to the union
         * of their knots; infinity when they cannot be (a multiplicity above the degree to reach).
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
                const auto m1 = mult(k1, u), m2 = mult(k2, u);
                if (m1 < m2)
                    raise(k1, rows1, u, m2 - m1);
                else if (m2 < m1)
                    raise(k2, rows2, u, m1 - m2);
            }
            if (k1 != k2 || rows1.size() != rows2.size())
                return std::numeric_limits<T>::infinity();
            T d{};
            for (std::size_t r = 0; r < rows1.size(); ++r)
            {
                if (rows1[r].size() != rows2[r].size())
                    return std::numeric_limits<T>::infinity();
                for (std::size_t i = 0; i < rows1[r].size(); ++i)
                    d = std::max(d, pole_distance(rows1[r][i], rows2[r][i], rational));
            }
            return d;
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
         * the multiplicity lowered back inside; nullopt when it stays above what the curve had.
         */
        template <typename T, std::size_t D>
        auto collapse(std::vector<T> k, std::vector<Row<T, D>> rows, std::size_t p, T from, T to, T tol)
            -> std::optional<std::pair<std::vector<T>, std::vector<Row<T, D>>>>
        {
            const auto m_from = mult(k, from), m_to = mult(k, to);
            const bool at_start = to == k.front(), at_end = to == k.back();
            std::ranges::replace(k, from, to);
            const auto m = mult(k, to);
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
         * within tol of the original. Returns the deviation from the original, nullopt if nothing changed.
         */
        template <typename T, std::size_t D>
        auto clean_rows(std::vector<T> &k, std::vector<Row<T, D>> &rows, std::size_t p, T tol, bool rational, KnotCleanupStats<T> &st)
            -> std::optional<T>
        {
            const auto k0 = k;
            const auto rows0 = rows;
            std::vector<std::pair<T, T>> kept;
            std::optional<T> dev;
            for (;;)
            {
                const auto x = distinct(k);
                if (x.size() < 3)
                    break; // a single span: nothing to merge it with
                std::optional<std::size_t> found;
                for (std::size_t i = 0; i + 1 < x.size() && !found; ++i)
                    if (std::ranges::find(kept, std::pair{x[i], x[i + 1]}) == kept.end() &&
                        span_size<T, D>(k, rows, p, x[i], x[i + 1], rational) <= tol)
                        found = i;
                if (!found)
                    break;
                const auto i = *found;
                const T a = x[i], b = x[i + 1];
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
                    const T d = deviation<T, D>(k0, rows0, c->first, c->second, p, rational);
                    if (d <= tol)
                    {
                        k = std::move(c->first);
                        rows = std::move(c->second);
                        dev = d;
                        ++st.removed_spans;
                        done = true;
                        break;
                    }
                }
                if (!done)
                {
                    kept.emplace_back(a, b);
                    ++st.kept_spans;
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
        const auto [u1, u2] = c.bounds();
        using Curve = std::conditional_t<rational, BSCurveRational<T, dim>, BSCurve<T, dim>>;
        Curve cleaned(rows.front(), k, c.degree());
        cleaned.trim(u1, u2, false); // the bounds as they were
        static_cast<BSCurveGeneral<T, dim, rational> &>(c) = cleaned;
        ++st.cleaned_curves;
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
