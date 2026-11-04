#pragma once
/**
 * @file bssjoin.h
 * @brief Two B-spline surfaces that share a side joined into one, as smooth as their junction
 * allows (`docs/sources/design/brep4_pr01ter_join_bsplines.md`).
 *
 * `join_surfaces(a, b, tol)` looks for a side of a and a side of b that are one curve within tol
 * (u = u_min, u_max, v = v_min or v_max of each, either way along), then builds one surface:
 *
 * - The parameters of a are kept: a point (u, v) of a is the same point (u, v) of the joined
 *   surface, to the deviation reported. b is put beyond the side, its parameters mapped by an
 *   affine map (a swap, reversals, a scale and a shift) given with the result.
 * - Along the side, the degrees are raised to the larger one and each surface gets the knots of
 *   the other (degree elevation and knot insertion, exact); its parameter on b is mapped onto
 *   that of a.
 * - Across the side, the degrees are raised to the larger one and b is scaled so that the length
 *   of its first derivative across meets that of a (least squares along the side): the condition
 *   for the junction to be more than C0. The common row of poles is the mean of the two.
 * - The knot of the junction, of multiplicity p (C0), is then removed as many times as every row
 *   allows within tol (NURBS book A5.8, the same number of times on every row): once for C1, p
 *   times when the two surfaces are pieces of one polynomial. `knots_removed` says how many.
 * - `deviation` bounds the distance between the joined surface and a or b at the same
 *   parameters: half the gap of the common rows plus the change of the poles by the removal (the
 *   knot inserted back is exact, and a B-spline is within the hull of its poles).
 *
 * Only non-rational surfaces, clamped on the side they share, with outward directions that agree
 * (the map of b has a positive determinant, so that the normal of b is that of the joined
 * surface). Nothing otherwise.
 *
 * `join_curves(a, b, tol)` is the same for two curves end to end (b starting where a ends, within
 * tol): b shifted to start where a ends, at its own speed (the parameters of both are kept: an
 * edge and its pcurves stay parametrized alike), the common pole their mean, the knot of the
 * junction removed as far as tol allows.
 */
#include <gbs/bssurf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import knots_functions;
    import vecop;
#endif

namespace gbs
{
    /// What `join_surfaces` gives.
    template <std::floating_point T, std::size_t dim>
    struct SurfaceJoin
    {
        BSSurface<T, dim> surface;   ///< on the parameters of a, b beyond the side they share
        std::array<T, 6> map{};      ///< (u, v) of b -> (m0 u + m1 v + m2, m3 u + m4 v + m5) on surface
        std::size_t knots_removed{}; ///< times the knot of the junction was removed: 0 for C0, 1 for C1, ...
        T deviation{};               ///< bound of the distance between surface and a or b at the same parameters
    };

    /// What `join_curves` gives.
    template <std::floating_point T, std::size_t dim>
    struct CurveJoin
    {
        BSCurve<T, dim> curve;       ///< a, then b shifted to start where a ends
        std::size_t knots_removed{}; ///< times the knot of the junction was removed: 0 for C0, 1 for C1, ...
        T deviation{};               ///< bound of the distance between curve and a or b at the same parameters
    };

    namespace bssjoin_detail
    {
        /// The knot u of a clamped curve removed as many times as tol allows (at most m): the times,
        /// and a bound of the change of the curve (the knot inserted back, exact, against the poles
        /// before).
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
            if (times == 0)
                return {0, T(0)};
            const auto before = rows;
            std::vector<T> kn;
            for (auto &r : rows)
            {
                kn = k;
                remove_knot(u, times, kn, r, p);
            }
            k = std::move(kn);
            T dev{};
            for (std::size_t i{}; i < rows.size(); ++i)
            {
                auto kb = k;
                auto back = rows[i];
                insert_knots(u, p, times, kb, back);
                for (std::size_t j{}; j < back.size(); ++j)
                    dev = std::max(dev, norm(back[j] - before[i][j]));
            }
            return {times, dev};
        }

        /// (u, v) -> (a[0] u + a[1] v + a[2], a[3] u + a[4] v + a[5])
        template <std::floating_point T>
        struct Map
        {
            std::array<T, 6> a{T(1), T(0), T(0), T(0), T(1), T(0)};

            auto operator()(T u, T v) const -> std::array<T, 2> { return {a[0] * u + a[1] * v + a[2], a[3] * u + a[4] * v + a[5]}; }
            /// This map after `first`.
            auto after(const Map &first) const -> Map
            {
                const auto &f = first.a;
                return {{a[0] * f[0] + a[1] * f[3], a[0] * f[1] + a[1] * f[4], a[0] * f[2] + a[1] * f[5] + a[2],
                         a[3] * f[0] + a[4] * f[3], a[3] * f[1] + a[4] * f[4], a[3] * f[2] + a[4] * f[5] + a[5]}};
            }
            auto det() const -> T { return a[0] * a[4] - a[1] * a[3]; }
        };

        /// A surface as rows of poles along u (rows[j][i], j along v), its knots and degrees.
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

            /// u -> c - u, by default c = u_min + u_max: the same range, run the other way.
            void reverse_u(std::optional<T> c_ = std::nullopt)
            {
                const T c = c_.value_or(ku.front() + ku.back());
                for (auto &r : rows)
                    std::ranges::reverse(r);
                std::ranges::reverse(ku);
                for (auto &k : ku)
                    k = c - k;
                map = Map<T>{{T(-1), T(0), c, T(0), T(1), T(0)}}.after(map);
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

            auto clamped_u() const -> bool { return multiplicity(ku, ku.front()) == p + 1 && multiplicity(ku, ku.back()) == p + 1; }
        };

        /// The grid of s turned so that its side `side` (0: u = u_min, 1: u = u_max, 2: v = v_min, 3:
        /// v = v_max) is u = u_max (`last`) or u = u_min, its v reversed with `flip`.
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

        /// The largest distance between two lists of points of the same size.
        template <std::floating_point T, std::size_t dim>
        auto gap(const std::vector<std::array<T, dim>> &x, const std::vector<std::array<T, dim>> &y) -> T
        {
            T d{};
            for (std::size_t k{}; k < x.size(); ++k)
                d = std::max(d, norm(x[k] - y[k]));
            return d;
        }

        template <std::floating_point T, std::size_t dim>
        struct Joined
        {
            Grid<T, dim> grid; ///< its map: that of b, from its parameters to the joined ones
            std::size_t times{};
            T dev{};
        };

        /// a and b turned (a's side last, b's first), joined; nothing if their sides differ.
        template <std::floating_point T, std::size_t dim>
        auto joined(Grid<T, dim> a, Grid<T, dim> b, T tol) -> std::optional<Joined<T, dim>>
        {
            if (!a.clamped_u() || !b.clamped_u())
                return std::nullopt;
            // the corners first: a cheap refusal
            const auto corner_gap = [&](std::size_t j_a, std::size_t j_b) { return norm(a.rows[j_a].back() - b.rows[j_b].front()); };
            if (!(corner_gap(0, 0) <= tol && corner_gap(a.nv() - 1, b.nv() - 1) <= tol))
                return std::nullopt;
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
                const T v = v0 + (v1 - v0) * T(k) / T(n);
                const auto da = sa.value(ua, v, 1, 0), db = sb.value(ub, v, 1, 0);
                num += norm(da) * norm(db); // (lengths: a kink has derivatives at an angle)
                den += da * da;
            }
            const T s = num / den;
            if (!(s > T(0)) || !std::isfinite(s))
                return std::nullopt;
            b.scale_u(s, ua - s * ub);
            b.ku.front() = ua; // (exactly)
            for (std::size_t k = 1; k <= p; ++k)
                b.ku[k] = ua;
            // one grid: a's rows, then b's after their first pole, the common pole their mean
            Grid<T, dim> j = a;
            for (std::size_t r{}; r < j.nv(); ++r)
            {
                j.rows[r].back() = T(0.5) * (a.rows[r].back() + b.rows[r].front());
                j.rows[r].insert(j.rows[r].end(), std::next(b.rows[r].begin()), b.rows[r].end());
            }
            j.ku.pop_back();
            j.ku.insert(j.ku.end(), std::next(b.ku.begin(), std::ptrdiff_t(p + 1)), b.ku.end());
            T dev = common / T(2);
            // the knot of the junction removed as many times as every row allows within tol
            const auto [times, removal] = remove_within(j.ku, j.rows, ua, p, p, tol);
            dev += removal; // (both bounds add up)
            j.map = b.map;
            return Joined<T, dim>{std::move(j), times, dev};
        }
    }

    /**
     * @brief a and b joined into one surface along a side they share within tol (see the file
     * comment); nothing if they share none, are not clamped there, or face opposite ways.
     */
    template <std::floating_point T, std::size_t dim>
    auto join_surfaces(const BSSurface<T, dim> &a, const BSSurface<T, dim> &b, T tol) -> std::optional<SurfaceJoin<T, dim>>
    {
        using namespace bssjoin_detail;
        for (int sa = 0; sa < 4; ++sa)
            for (int sb = 0; sb < 4; ++sb)
                for (const bool flip : {false, true})
                {
                    auto ga = turned(a, sa, true, false);
                    auto gb = turned(b, sb, false, flip);
                    if (ga.map.det() * gb.map.det() < T(0))
                        continue; // b would face the other way
                    auto j = joined(ga, gb, tol);
                    if (!j)
                        continue;
                    // back to the parameters of a: its turns undone in reverse order, each its own
                    // inverse with the same constant (`turned`: a transpose for a side along v, then
                    // a reversal of u for a side at u_min); b's map follows
                    auto &g = j->grid;
                    if (sa % 2 == 0)
                        g.reverse_u(ga.ku.front() + ga.ku.back());
                    if (sa >= 2)
                        g.transpose();
                    return SurfaceJoin<T, dim>{g.surface(), g.map.a, j->times, j->dev};
                }
        return std::nullopt;
    }

    /**
     * @brief a then b, end to end (b's start within tol of a's end), as one curve (see the file
     * comment); nothing if they do not meet, or are not clamped where they do.
     */
    template <std::floating_point T, std::size_t dim>
    auto join_curves(const BSCurve<T, dim> &a, const BSCurve<T, dim> &b, T tol) -> std::optional<CurveJoin<T, dim>>
    {
        auto ka = a.knotsFlats(), kb = b.knotsFlats();
        auto pa = a.poles(), pb = b.poles();
        std::size_t p = a.degree(), q = b.degree();
        if (multiplicity(ka, ka.back()) != p + 1 || multiplicity(kb, kb.front()) != q + 1)
            return std::nullopt;
        if (q < p)
            increase_degree(kb, pb, q, p - q);
        if (p < q)
            increase_degree(ka, pa, p, q - p);
        p = std::max(p, q);
        const T gap = norm(pa.back() - pb.front());
        if (!(gap <= tol))
            return std::nullopt;
        const T u = ka.back(), d = u - kb.front();
        for (auto &k : kb)
            k += d;
        for (std::size_t i{}; i <= p; ++i)
            kb[i] = u; // (exactly)
        pa.back() = T(0.5) * (pa.back() + pb.front());
        pa.insert(pa.end(), std::next(pb.begin()), pb.end());
        ka.pop_back();
        ka.insert(ka.end(), std::next(kb.begin(), std::ptrdiff_t(p + 1)), kb.end());
        std::vector<std::vector<std::array<T, dim>>> rows{std::move(pa)};
        const auto [times, removal] = bssjoin_detail::remove_within(ka, rows, u, p, p, tol);
        return CurveJoin<T, dim>{BSCurve<T, dim>{rows.front(), ka, p}, times, gap / T(2) + removal};
    }
}
