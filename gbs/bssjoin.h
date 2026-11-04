#pragma once
/**
 * @file bssjoin.h
 * @brief Two B-spline surfaces that share a side joined into one, as smooth as their junction
 * allows (`docs/sources/design/brep4_pr01ter_join_bsplines.md`).
 *
 * `join_surfaces(a, b, tol)` looks for a side of a and a side of b that are one curve within tol
 * (u = u_min, u_max, v = v_min or v_max of each, either way along; a side collapsed to a point, a
 * pole or an apex, is no side to join along), then builds one surface; of the pairs of sides that
 * join, the smoothest (the most knots removed, then the least deviation). `join_surfaces(a, b,
 * side_a, side_b, tol)` joins along the sides given only.
 *
 * - The parameters of a are kept, exactly: a point (u, v) of a is the same point (u, v) of the
 *   joined surface, to the deviation reported, and the bounds of the joined surface are those of a
 *   beyond the side, bit for bit. b is put beyond the side, its parameters mapped by an affine map
 *   (a swap, reversals, a scale and a shift) given with the result; its bounds map onto the others
 *   to the rounding of the map (a few ulps, either way: an image of b's bounds is put back onto
 *   them before the joined surface is evaluated there, as `unify_faces` does).
 * - Along the side, the degrees are raised to the larger one and each surface gets the knots of
 *   the other (degree elevation and knot insertion, exact); its parameter on b is mapped onto
 *   that of a.
 * - Across the side, the degrees are raised to the larger one and b is scaled so that the length
 *   of its first derivative across meets that of a (least squares along the side): the condition
 *   for the junction to be more than C0. The common row of poles is the mean of the two.
 * - The knot of the junction, of multiplicity p (C0), is then removed as many times as every row
 *   allows (the same number of times on every row): once for C1, p times when the two surfaces are
 *   pieces of one polynomial. `knots_removed` says how many.
 * - `deviation` bounds the distance between the joined surface and a or b at the same parameters:
 *   half the gap of the common rows plus the change of the poles by the removal (the knot inserted
 *   back is exact, and a B-spline is within the hull of its poles). The removal is measured so, and
 *   kept within what tol leaves: `deviation` <= tol.
 *
 * Only non-rational surfaces, clamped at all their bounds, with outward directions that agree
 * (the map of b has a positive determinant, so that the normal of b is that of the joined
 * surface). Nothing otherwise.
 *
 * `join_curves(a, b, tol)` is the same for two curves end to end (b starting where a ends, within
 * tol), both clamped: `join(a, b)` (b shifted to start where a ends, at its own speed: the
 * parameters of both are kept, an edge and its pcurves stay parametrized alike; the common pole
 * their mean), then the knot of the junction removed as far as tol allows.
 */
#include <gbs/bsctools.h>
#include <gbs/bssurf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import knots_functions;
    import vecop;
#endif

namespace gbs
{
    /// A side of a surface, at one of its bounds.
    enum class SurfaceSide
    {
        u_min,
        u_max,
        v_min,
        v_max
    };

    /// What `join_surfaces` gives.
    template <std::floating_point T, std::size_t dim>
    struct SurfaceJoin
    {
        BSSurface<T, dim> surface;   ///< on the parameters of a, b beyond the side they share
        std::array<T, 6> map{};      ///< (u, v) of b -> (m0 u + m1 v + m2, m3 u + m4 v + m5) on surface
        std::size_t knots_removed{}; ///< times the knot of the junction was removed: 0 for C0, 1 for C1, ...
        T deviation{};               ///< bound of the distance between surface and a or b at the same parameters (<= tol)
        SurfaceSide side_a{};        ///< the side of a joined
        SurfaceSide side_b{};        ///< the side of b joined
    };

    /// What `join_curves` gives.
    template <std::floating_point T, std::size_t dim>
    struct CurveJoin
    {
        BSCurve<T, dim> curve;       ///< a, then b shifted to start where a ends
        std::size_t knots_removed{}; ///< times the knot of the junction was removed: 0 for C0, 1 for C1, ...
        T deviation{};               ///< bound of the distance between curve and a or b at the same parameters (<= tol)
    };

    namespace bssjoin_detail
    {
        /// Whether a flat knot vector of degree p is clamped at both ends.
        template <std::floating_point T>
        auto clamped(const std::vector<T> &k, std::size_t p) -> bool
        {
            return k.size() >= 2 * (p + 1) && multiplicity(k, k.front()) == p + 1 && multiplicity(k, k.back()) == p + 1;
        }

        /// The knot u of a clamped curve removed as many times as possible, at most m, the change of
        /// the curve within tol: the times, and the change, measured as a bound (the knot inserted
        /// back, exact, against the poles before). A5.8 (NURBS book) on every row against tol tells
        /// how many times to try; each pass of A5.8 is checked against tol alone and the passes add
        /// up, so the measure decides, one time fewer while it exceeds tol.
        template <std::floating_point T, std::size_t dim>
        auto remove_within(std::vector<T> &k, std::vector<std::vector<std::array<T, dim>>> &rows, T u, std::size_t p, std::size_t m, T tol)
            -> std::pair<std::size_t, T>
        {
            std::size_t times = m;
            for (const auto &r : rows)
            {
                auto kt = k;
                auto q = r;
                times = std::min(times, remove_knot(u, p, m, kt, q, tol));
            }
            for (; times > 0; --times)
            {
                auto removed = rows;
                std::vector<T> kn;
                for (auto &r : removed)
                {
                    kn = k;
                    remove_knot(u, times, kn, r, p);
                }
                T dev{};
                for (std::size_t i{}; i < rows.size(); ++i)
                {
                    auto kb = kn;
                    auto back = removed[i];
                    insert_knots(u, p, times, kb, back);
                    for (std::size_t j{}; j < back.size(); ++j)
                        dev = std::max(dev, norm(back[j] - rows[i][j]));
                }
                if (dev <= tol)
                {
                    k = std::move(kn);
                    rows = std::move(removed);
                    return {times, dev};
                }
            }
            return {0, T(0)};
        }

        /// (u, v) -> (a[0] u + a[1] v + a[2], a[3] u + a[4] v + a[5])
        template <std::floating_point T>
        struct Map
        {
            std::array<T, 6> a{T(1), T(0), T(0), T(0), T(1), T(0)};

            /// This map after `first`.
            auto after(const Map &first) const -> Map
            {
                const auto &f = first.a;
                return {{a[0] * f[0] + a[1] * f[3], a[0] * f[1] + a[1] * f[4], a[0] * f[2] + a[1] * f[5] + a[2],
                         a[3] * f[0] + a[4] * f[3], a[3] * f[1] + a[4] * f[4], a[3] * f[2] + a[4] * f[5] + a[5]}};
            }
        };

        /// A surface as rows of poles along u (rows[j][i], j along v), its knots and degrees.
        ///
        /// The turns (a swap, reversals as u -> -u) are exact on the knots and on the map: undone,
        /// they give the knots back bit for bit.
        template <std::floating_point T, std::size_t dim>
        struct Grid
        {
            std::vector<std::vector<std::array<T, dim>>> rows;
            std::vector<T> ku, kv;
            std::size_t p{}, q{};
            Map<T> map; ///< the parameters of the surface it came from -> its own

            explicit Grid(const BSSurface<T, dim> &s) : ku{s.knotsFlatsU()}, kv{s.knotsFlatsV()}, p{s.degreeU()}, q{s.degreeV()}
            {
                for (std::size_t j{}; j < s.nPolesV(); ++j)
                    rows.push_back(s.polesU(j));
            }

            auto nu() const -> std::size_t { return rows.front().size(); }
            auto nv() const -> std::size_t { return rows.size(); }

            auto surface() const -> BSSurface<T, dim>
            {
                std::vector<std::array<T, dim>> poles;
                poles.reserve(nu() * nv());
                for (const auto &r : rows)
                    poles.insert(poles.end(), r.begin(), r.end());
                return {poles, ku, kv, p, q};
            }

            /// u and v swapped.
            void transpose()
            {
                std::vector<std::vector<std::array<T, dim>>> t(nu(), std::vector<std::array<T, dim>>(nv()));
                for (std::size_t j{}; j < nv(); ++j)
                    for (std::size_t i{}; i < nu(); ++i)
                        t[i][j] = rows[j][i];
                rows = std::move(t);
                std::swap(ku, kv);
                std::swap(p, q);
                map = Map<T>{{T(0), T(1), T(0), T(1), T(0), T(0)}}.after(map);
            }

            /// u -> -u: the surface run the other way along u, exactly (twice: as it was).
            void reverse_u()
            {
                for (auto &r : rows)
                    std::ranges::reverse(r);
                std::ranges::reverse(ku);
                for (auto &k : ku)
                    k = -k;
                map = Map<T>{{T(-1), T(0), T(0), T(0), T(1), T(0)}}.after(map);
            }

            void reverse_v()
            {
                transpose();
                reverse_u();
                transpose();
            }

            /// u -> s u + d on the knots.
            void scale_u(T s, T d)
            {
                for (auto &k : ku)
                    k = s * k + d;
                map = Map<T>{{s, T(0), d, T(0), T(1), T(0)}}.after(map);
            }

            void raise_u(std::size_t to)
            {
                if (to <= p)
                    return;
                std::vector<T> k;
                for (auto &r : rows)
                {
                    k = ku;
                    increase_degree(k, r, p, to - p);
                }
                ku = std::move(k);
                p = to;
            }

            /// The knot u inserted until its multiplicity is m (exact).
            void insert_u(T u, std::size_t m)
            {
                const auto have = multiplicity(ku, u);
                if (have >= m)
                    return;
                std::vector<T> k;
                for (auto &r : rows)
                {
                    k = ku;
                    insert_knots(u, p, m - have, k, r);
                }
                ku = std::move(k);
            }
        };

        /// The sign of the determinant of the map of `turned(s, side, last, flip)`: -1 for each
        /// swap and reversal.
        inline auto turn_sign(int side, bool last, bool flip) -> int
        {
            return (side >= 2 ? -1 : 1) * ((side % 2 == 1) != last ? -1 : 1) * (flip ? -1 : 1);
        }

        /// The grid of s turned so that its side `side` (SurfaceSide: 0 u_min, 1 u_max, 2 v_min, 3
        /// v_max) is u = u_max (`last`) or u = u_min, its v reversed with `flip`.
        template <std::floating_point T, std::size_t dim>
        auto turned(const BSSurface<T, dim> &s, int side, bool last, bool flip) -> Grid<T, dim>
        {
            Grid<T, dim> g{s};
            if (side >= 2)
                g.transpose();
            if ((side % 2 == 1) != last)
                g.reverse_u();
            if (flip)
                g.reverse_v();
            return g;
        }

        /// The poles of the side `side` of a clamped surface, along it (v for a side at a bound of
        /// u, u for one at a bound of v): the side itself.
        template <std::floating_point T, std::size_t dim>
        auto side_poles(const BSSurface<T, dim> &s, int side) -> std::vector<std::array<T, dim>>
        {
            const std::size_t nu = s.nPolesU(), nv = s.nPolesV();
            const auto &poles = s.poles();
            std::vector<std::array<T, dim>> r;
            if (side < 2)
                for (std::size_t j{}; j < nv; ++j)
                    r.push_back(poles[j * nu + (side == 0 ? 0 : nu - 1)]);
            else
                for (std::size_t i{}; i < nu; ++i)
                    r.push_back(poles[(side == 2 ? 0 : nv - 1) * nu + i]);
            return r;
        }

        /// Whether the poles of a side are all within tol of its first one: a side collapsed to a
        /// point (a pole, an apex).
        template <std::floating_point T, std::size_t dim>
        auto collapsed(const std::vector<std::array<T, dim>> &side, T tol) -> bool
        {
            return std::ranges::all_of(side, [&](const auto &x) { return norm(x - side.front()) <= tol; });
        }

        /// Both grids given the same degree and the same knots along v, the knots of b along v first
        /// mapped onto the range of a; false if their ranges cannot be mapped.
        template <std::floating_point T, std::size_t dim>
        auto same_along_v(Grid<T, dim> &a, Grid<T, dim> &b) -> bool
        {
            a.transpose();
            b.transpose();
            const T a0 = a.ku.front(), a1 = a.ku.back(), b0 = b.ku.front(), b1 = b.ku.back();
            bool ok = b1 > b0 && a1 > a0;
            if (ok)
            {
                b.scale_u((a1 - a0) / (b1 - b0), a0 - b0 * (a1 - a0) / (b1 - b0));
                const std::size_t deg = std::max(a.p, b.p);
                a.raise_u(deg);
                b.raise_u(deg);
                // the knots of each in the other, a knot of b within eps of one of a taken as it
                const T eps = T(1e-10) * (a1 - a0);
                for (auto &k : b.ku)
                    for (const auto ka : a.ku)
                        if (std::abs(k - ka) <= eps)
                            k = ka;
                const auto [ua, ma] = knots_and_mults(a.ku);
                const auto [ub, mb] = knots_and_mults(b.ku);
                for (std::size_t i{}; i < ub.size(); ++i)
                    a.insert_u(ub[i], std::min(mb[i], deg));
                for (std::size_t i{}; i < ua.size(); ++i)
                    b.insert_u(ua[i], std::min(ma[i], deg));
                ok = a.ku.size() == b.ku.size() &&
                     std::ranges::equal(a.ku, b.ku, [&](T x, T y) { return std::abs(x - y) <= eps; });
                b.ku = a.ku; // (the same values)
            }
            a.transpose();
            b.transpose();
            return ok;
        }

        /// The largest distance between two lists of points of the same size (not empty).
        template <std::floating_point T, std::size_t dim>
        auto gap(const std::vector<std::array<T, dim>> &x, const std::vector<std::array<T, dim>> &y) -> T
        {
            return std::ranges::max(std::views::zip(x, y) | std::views::transform([](const auto &pq) { return norm(std::get<0>(pq) - std::get<1>(pq)); }));
        }

        template <std::floating_point T, std::size_t dim>
        struct Joined
        {
            Grid<T, dim> grid; ///< its map: that of b, from its parameters to the joined ones
            std::size_t times{};
            T dev{};
        };

        /// a and b turned (a's side last, b's first, both clamped, their corners within tol),
        /// joined; nothing if their sides differ.
        template <std::floating_point T, std::size_t dim>
        auto joined(Grid<T, dim> a, Grid<T, dim> b, T tol) -> std::optional<Joined<T, dim>>
        {
            if (!same_along_v(a, b))
                return std::nullopt;
            const std::size_t p = std::max(a.p, b.p);
            a.raise_u(p);
            b.raise_u(p);
            // the common rows
            std::vector<std::array<T, dim>> side_a, side_b;
            for (std::size_t j{}; j < a.nv(); ++j)
            {
                side_a.push_back(a.rows[j].back());
                side_b.push_back(b.rows[j].front());
            }
            const T common = gap(side_a, side_b);
            if (!(common <= tol))
                return std::nullopt;
            // b's u scaled by s, so that its derivative across, db / s, meets that of a: s minimises
            // the sum of (|db| - s |da|)² along the side
            const T ua = a.ku.back(), ub = b.ku.front();
            const auto sa = a.surface(), sb = b.surface();
            const auto [v0, v1] = std::array{a.kv.front(), a.kv.back()};
            T num{}, den{};
            const std::size_t n = 4 * (a.nv() + 1);
            for (std::size_t k{}; k <= n; ++k)
            {
                const T v = std::min(v0 + (v1 - v0) * T(k) / T(n), v1); // (the rounding past v1 refused)
                const auto da = sa.value(ua, v, 1, 0), db = sb.value(ub, v, 1, 0);
                num += norm(da) * norm(db); // (lengths: a kink has derivatives at an angle)
                den += da * da;
            }
            const T s = num / den;
            if (!(s > T(0)) || !std::isfinite(s))
                return std::nullopt;
            b.scale_u(s, ua - s * ub);
            std::ranges::fill_n(b.ku.begin(), std::ptrdiff_t(p + 1), ua); // (exactly)
            // one grid: a's rows, then b's after their first pole, the common pole their mean
            Grid<T, dim> j = a;
            for (std::size_t r{}; r < j.nv(); ++r)
            {
                j.rows[r].back() = T(0.5) * (a.rows[r].back() + b.rows[r].front());
                j.rows[r].insert(j.rows[r].end(), std::next(b.rows[r].begin()), b.rows[r].end());
            }
            j.ku.pop_back();
            j.ku.insert(j.ku.end(), std::next(b.ku.begin(), std::ptrdiff_t(p + 1)), b.ku.end());
            // the knot of the junction removed as many times as every row allows, within what the
            // mean of the common rows leaves of tol
            const auto [times, removal] = remove_within(j.ku, j.rows, ua, p, p, tol - common / T(2));
            j.map = b.map;
            return Joined<T, dim>{std::move(j), times, common / T(2) + removal}; // (both bounds add up)
        }

        /// a and b joined along the sides given (any of them with `any`), the best join kept: the
        /// most knots removed, then the least deviation.
        template <std::floating_point T, std::size_t dim>
        auto join_along(const BSSurface<T, dim> &a, const BSSurface<T, dim> &b, std::optional<std::array<int, 2>> sides, T tol)
            -> std::optional<SurfaceJoin<T, dim>>
        {
            if (!clamped(a.knotsFlatsU(), a.degreeU()) || !clamped(a.knotsFlatsV(), a.degreeV()) ||
                !clamped(b.knotsFlatsU(), b.degreeU()) || !clamped(b.knotsFlatsV(), b.degreeV()))
                return std::nullopt;
            std::optional<SurfaceJoin<T, dim>> best;
            for (int sa = 0; sa < 4; ++sa)
            {
                if (sides && (*sides)[0] != sa)
                    continue;
                const auto pa = side_poles(a, sa);
                if (collapsed(pa, tol))
                    continue;
                std::optional<Grid<T, dim>> ga; // (turned once, when a side of b may match)
                for (int sb = 0; sb < 4; ++sb)
                {
                    if (sides && (*sides)[1] != sb)
                        continue;
                    const auto pb = side_poles(b, sb);
                    if (collapsed(pb, tol))
                        continue;
                    for (const bool flip : {false, true})
                    {
                        if (turn_sign(sa, true, false) * turn_sign(sb, false, flip) < 0)
                            continue; // b would face the other way
                        // the corners first, on the surfaces as they are: a cheap refusal
                        const auto &b0 = flip ? pb.back() : pb.front(), &b1 = flip ? pb.front() : pb.back();
                        if (!(norm(pa.front() - b0) <= tol && norm(pa.back() - b1) <= tol))
                            continue;
                        if (!ga)
                            ga = turned(a, sa, true, false);
                        auto j = joined(*ga, turned(b, sb, false, flip), tol);
                        if (!j)
                            continue;
                        if (best && (j->times < best->knots_removed || (j->times == best->knots_removed && !(j->dev < best->deviation))))
                            continue;
                        // back to the parameters of a: its turns undone in reverse order, each its
                        // own inverse (`turned`: a transpose for a side along v, then a reversal of
                        // u for a side at u_min); b's map follows
                        auto &g = j->grid;
                        if (sa % 2 == 0)
                            g.reverse_u();
                        if (sa >= 2)
                            g.transpose();
                        best = SurfaceJoin<T, dim>{g.surface(), g.map.a, j->times, j->dev, SurfaceSide(sa), SurfaceSide(sb)};
                    }
                }
            }
            return best;
        }
    }

    /**
     * @brief a and b joined into one surface along a side they share within tol, the smoothest
     * join of all such sides (see the file comment); nothing if they share none, are not clamped,
     * or face opposite ways.
     */
    template <std::floating_point T, std::size_t dim>
    auto join_surfaces(const BSSurface<T, dim> &a, const BSSurface<T, dim> &b, T tol) -> std::optional<SurfaceJoin<T, dim>>
    {
        return bssjoin_detail::join_along(a, b, std::nullopt, tol);
    }

    /**
     * @brief a and b joined into one surface along the side `side_a` of a and `side_b` of b, which
     * are one curve within tol, either way along (see the file comment); nothing if they are not,
     * a side is collapsed to a point, a surface is not clamped, or they face opposite ways.
     */
    template <std::floating_point T, std::size_t dim>
    auto join_surfaces(const BSSurface<T, dim> &a, const BSSurface<T, dim> &b, SurfaceSide side_a, SurfaceSide side_b, T tol)
        -> std::optional<SurfaceJoin<T, dim>>
    {
        return bssjoin_detail::join_along(a, b, std::array{int(side_a), int(side_b)}, tol);
    }

    /**
     * @brief a then b, end to end (b's start within tol of a's end), as one curve (see the file
     * comment); nothing if they do not meet, or are not clamped.
     */
    template <std::floating_point T, std::size_t dim>
    auto join_curves(const BSCurve<T, dim> &a, const BSCurve<T, dim> &b, T tol) -> std::optional<CurveJoin<T, dim>>
    {
        using namespace bssjoin_detail;
        if (!clamped(a.knotsFlats(), a.degree()) || !clamped(b.knotsFlats(), b.degree()))
            return std::nullopt;
        const T gap = norm(a.poles().back() - b.poles().front()); // (their ends: clamped)
        if (!(gap <= tol))
            return std::nullopt;
        const auto c = join(a, b);
        auto k = c.knotsFlats();
        std::vector<std::vector<std::array<T, dim>>> rows{c.poles()};
        const auto [times, removal] = remove_within(k, rows, a.knotsFlats().back(), c.degree(), c.degree(), tol - gap / T(2));
        return CurveJoin<T, dim>{BSCurve<T, dim>{rows.front(), k, c.degree()}, times, gap / T(2) + removal};
    }
}
