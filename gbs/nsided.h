#pragma once
/**
 * @file nsided.h
 * @brief A surface filling an n-sided hole, tangent to the surfaces round it (stage 6, PR 2,
 * `docs/sources/design/brep_palier6_vertex_blends.md`, § 4.3 to § 4.6; note
 * `brep6_pr02_nsided.md`).
 *
 * Each side of the hole is a **ribbon**: a curve P(t) and a cross direction X(t), tangent to the
 * surface next to that side and pointing into the hole. Its direction sets the tangent plane; its
 * length, as it varies along the side, shapes the cross derivatives (step 2, b linear against X as
 * given): the derivatives of the surface next to the side across it, as they are, make the best
 * surface (a surface of the second degree found again), and a constant factor on a side does not
 * count. The sides follow one another round the hole (the end of one at the start of the next),
 * turning counter-clockwise seen from the side the normals P' x X point to.
 *
 * 1. **Domain.** The regular polygon of m sides inscribed in the unit circle, its side i from the
 *    vertex V_i to V_{i+1}, counter-clockwise, the image of [t0, t1] of ribbon i (linearly).
 * 2. **Cross derivatives** compatible at the corners: along side i, the derivative of the surface
 *    along the inner normal of the side in the domain, a(s) P'(s) + b(s) X(s), a and b linear in
 *    s; at each corner, the derivatives along both sides of the domain are the tangents of the two
 *    curves meeting there, which sets a and b at the ends.
 * 3. **B-spline.** A bicubic B-spline on the square [-1, 1]^2 holding the polygon, N x N uniform
 *    spans, of least bending energy (|S_uu|^2 + 2 |S_uv|^2 + |S_vv|^2 integrated over the polygon
 *    by a rule of degree 4 on small triangles; over the rest of the square, a thousand times
 *    smaller, to hold its corners), by weighted least squares, each row weighted by the square
 *    root of its share of a side or of the area: its points on the sides on the curves and its
 *    derivatives there normal to the ribbon (n . S_u = n . S_v = 0, n = P' x X normalized;
 *    `side_weight`), its cross derivatives near those of step 2 (`cross_weight`, fading to none at
 *    the corners, where the twists the two sides ask for seldom agree). The cross derivatives
 *    hold the parameters of the surface along the sides: with the normals alone, the bending
 *    energy slides them and the shape does not settle as the spans are doubled. A surface of the
 *    second degree whose sides are images of the polygon by an affine map is found again (it is
 *    biharmonic). A Kato surface on the ribbons, tried as a guide, sagged by a third of the depth
 *    of a pentagonal hole in a paraboloid: it is not used.
 * 4. **Refinement.** N is doubled from 4 until the surface is within `approx_tol` of the curves
 *    and its normal within `angular_tol` of the ribbons' along every side (sampled finely), and
 *    inside, along its normal, within `shape_tol` of where the doublings lead (sampled on a grid:
 *    the move from the last doubling, or what is left to come where the moves fall
 *    geometrically, the move times rho / (1 - rho), rho the ratio of the last two moves), up to
 *    `max_spans`. Every weight is relative to the size of the hole: the same surface to scale.
 *    The normal equations of the least squares are summed as the rows come, in a band round each
 *    pole: a fit takes about 70 MB at 64 spans, whatever the number of sides, a quarter of it at
 *    32, 300 MB and 1.6 s at 128.
 * 5. **Guide sphere** (`NSidedOptions::sphere`, none by default). Each point of the surface pulled
 *    toward the foot on the sphere of the surface before (along the direction from its centre;
 *    twice at each doubling, three fits in all), weighted against the bending energy and fading
 *    as the square of the distance to the sides (the normals there are the ribbons'): the hole of
 *    a vertex blend whose fillets nearly agree on a sphere, filled near that sphere. Each
 *    coordinate on its own: a row along the radius alone would couple the three at every point
 *    inside (up to 128 spans, 860 MB and 26 s instead of 340 MB and 7 s).
 *
 * The surface is meant to be trimmed by the polygon: its sides are segments of the domain.
 * Refused (`std::invalid_argument`): fewer than three ribbons; tolerances or weights out of range,
 * fewer than 8 spans allowed; a ribbon without curve, cross direction or interval; a curve whose
 * derivative vanishes; a cross direction null or along its curve; sides that do not follow one
 * another, or whose gap at a corner is above twice `approx_tol`; ribbons whose normals disagree at
 * a corner by more than `angular_tol` allows (the sine of half their angle above it); a corner
 * where the hole is not convex (180 degrees or more), or a cross direction that does not point
 * into it there. (`std::runtime_error`) tolerances that `max_spans` cannot reach; the message gives
 * the three criteria at the last doubling.
 */
#ifdef GBS_USE_MODULES
import basis_functions;
import vecop;
#endif
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <format>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace gbs
{
    /// A side of an n-sided hole: the curve on [t0, t1] and the cross direction X(t), tangent to the surface next to it, into the hole (best its derivative across the side; a constant factor does not count).
    template <std::floating_point T>
    struct Ribbon
    {
        std::shared_ptr<Curve<T, 3>> curve;
        std::shared_ptr<Curve<T, 3>> cross;
        T t0{}, t1{1};
    };

    /// A sphere the surface leans toward inside the hole (`NSidedOptions::sphere`): its centre, its radius, and the weight of the gap along its radius against the bending energy.
    template <std::floating_point T>
    struct GuideSphere
    {
        std::array<T, 3> centre{};
        T radius{};
        T weight{};
    };

    template <std::floating_point T>
    struct NSidedOptions
    {
        T approx_tol{T(1e-5)};      ///< largest distance of the surface to the curves along the sides
        T angular_tol{T(1e-5)};     ///< largest sine of the angle of its normal to the ribbons' along the sides
        std::size_t max_spans{64};  ///< largest number of spans of the B-spline in each direction (8 at least: the shape is compared from one doubling to the next)
        T side_weight{T(1e5)};       ///< weight of the positions and normals along the sides against the bending energy (positive)
        T cross_weight{T(1e3)};      ///< weight of the cross derivatives along the sides, none at the corners (they hold the parameters there; 0: none)
        T shape_tol{T(0)};           ///< largest gap of the surface inside, along its normal, to where the doublings of the spans lead: the move from the last doubling, or the move left to come where the moves fall geometrically, whichever is smaller (0: 2e-3 of the size of the hole, 10 approx_tol at least)
        std::size_t side_samples{0}; ///< samples on each side for the checks (0: 64 per span, 256 at least)
        /// a sphere the surface leans toward inside (the hole of a vertex blend near a sphere of
        /// its fillets): none, the surface of least bending energy
        std::optional<GuideSphere<T>> sphere{};
    };

    /// The surface filling an n-sided hole: the B-spline, the polygon of its domain (side i of ribbon i), and its measured gaps.
    template <std::floating_point T>
    struct NSidedSurface
    {
        std::shared_ptr<BSSurface<T, 3>> surface;
        std::vector<std::array<T, 2>> polygon;
        T deviation{}; ///< largest distance to the curves along the sides
        T angle{};     ///< largest sine of the angle of the normal to the ribbons' along the sides
        std::size_t spans{};
    };

    namespace nsided_detail
    {
        template <std::floating_point T>
        auto unit(const std::array<T, 3> &v) -> std::array<T, 3>
        {
            const T n = norm(v);
            return {v[0] / n, v[1] / n, v[2] / n};
        }

        /// The regular polygon of m sides inscribed in the unit circle, its first side at the bottom, counter-clockwise.
        template <std::floating_point T>
        auto regular_polygon(std::size_t m) -> std::vector<std::array<T, 2>>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            std::vector<std::array<T, 2>> v;
            for (std::size_t k = 0; k < m; ++k)
            {
                const T a = -pi / T(2) - pi / T(m) + T(2) * pi * T(k) / T(m);
                v.push_back({std::cos(a), std::sin(a)});
            }
            return v;
        }

        /// The distances of p to the lines of the sides of the polygon (positive inside).
        template <std::floating_point T>
        auto distances(const std::vector<std::array<T, 2>> &poly, const std::array<T, 2> &p) -> std::vector<T>
        {
            const std::size_t m = poly.size();
            std::vector<T> d(m);
            for (std::size_t i = 0; i < m; ++i)
            {
                const auto &a = poly[i], &b = poly[(i + 1) % m];
                const T ex = b[0] - a[0], ey = b[1] - a[1], l = std::hypot(ex, ey);
                d[i] = (ex * (p[1] - a[1]) - ey * (p[0] - a[0])) / l;
            }
            return d;
        }

        /// The point of the domain at the parameter s of the side i.
        template <std::floating_point T>
        auto on_side(const std::vector<std::array<T, 2>> &poly, std::size_t i, T s) -> std::array<T, 2>
        {
            const auto &a = poly[i], &b = poly[(i + 1) % poly.size()];
            return {a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])};
        }

        /// The ribbons evaluated: the curve at s in [0, 1] of the side, its derivative in s, the unit cross direction.
        template <std::floating_point T>
        struct Side
        {
            const Ribbon<T> *r;
            auto t(T s) const { return r->t0 + s * (r->t1 - r->t0); }
            auto p(T s) const { return r->curve->value(t(s)); }
            auto dp(T s) const
            {
                const auto d = r->curve->value(t(s), 1);
                const T k = r->t1 - r->t0;
                return std::array<T, 3>{k * d[0], k * d[1], k * d[2]};
            }
            auto x(T s) const { return unit<T>(r->cross->value(t(s))); }
            auto x_given(T s) const { return r->cross->value(t(s)); }
            auto normal(T s) const { return unit<T>(cross(dp(s), x(s))); }
        };

        /// The basis functions of a uniform cubic B-spline on [-1, 1] of n spans at u: the first pole index and the values of order d.
        template <std::floating_point T>
        auto basis(const std::vector<T> &k, std::size_t n_poles, T u, std::size_t d) -> std::pair<std::size_t, std::array<T, 4>>
        {
            const std::size_t span = std::size_t(find_span(n_poles, 3, u, k) - k.begin());
            std::array<T, 4> nd{};
            basis_ders(span, 3, d, u, k, nd.data());
            return {span - 3, nd};
        }

        /// The surface measured along the sides: largest distance to the curves and sine of the angle to the ribbons' normals.
        template <std::floating_point T>
        auto measure(const BSSurface<T, 3> &srf, const std::vector<Side<T>> &sides, const std::vector<std::array<T, 2>> &poly, std::size_t samples)
            -> std::pair<T, T>
        {
            T dev{}, ang{};
            for (std::size_t i = 0; i < sides.size(); ++i)
                for (std::size_t k = 0; k <= samples; ++k)
                {
                    const T s = T(k) / T(samples);
                    const auto uv = on_side(poly, i, s);
                    const auto q = srf.value(uv[0], uv[1]);
                    dev = std::max(dev, norm(q - sides[i].p(s)));
                    const auto su = srf.value(uv[0], uv[1], 1, 0), sv = srf.value(uv[0], uv[1], 0, 1);
                    const auto ns = cross(su, sv);
                    const T l = norm(ns);
                    const auto n = sides[i].normal(s);
                    const T c = (ns * n) / l;
                    ang = std::max(ang, c > T(0) ? norm(cross(ns, n)) / l : T(1));
                }
            return {dev, ang};
        }

        /**
         * The cross derivatives of the sides, compatible at the corners: along side i, the
         * derivative of the surface along the inner normal of the side in the domain, a(s) P'(s) +
         * b(s) X(s), a and b linear in s, their values at the ends drawn from the tangents of the
         * two curves meeting at each corner (the derivatives of the surface along both sides of
         * the domain there). Refused where b <= 0 at a corner: the hole not convex there (the
         * previous side coming back behind X), or X not into it.
         */
        template <std::floating_point T>
        auto cross_derivatives(const std::vector<Side<T>> &sides, const std::vector<std::array<T, 2>> &poly) -> std::vector<std::array<T, 4>>
        {
            const std::size_t m = sides.size();
            auto dir = [&](std::size_t from, std::size_t to) {
                const T x = poly[to][0] - poly[from][0], y = poly[to][1] - poly[from][1], l = std::hypot(x, y);
                return std::array<T, 3>{x / l, y / l, l};
            };
            // the derivative S_nu = alpha S_e + beta S_f in the plane of P' and X: (a, b)
            auto split = [](const std::array<T, 3> &v, const std::array<T, 3> &dp, const std::array<T, 3> &x) {
                const T g11 = dp * dp, g12 = dp * x, g22 = x * x, r1 = v * dp, r2 = v * x, det = g11 * g22 - g12 * g12;
                return std::array<T, 2>{(r1 * g22 - r2 * g12) / det, (g11 * r2 - g12 * r1) / det};
            };
            std::vector<std::array<T, 4>> ab(m);
            for (std::size_t i = 0; i < m; ++i)
            {
                const std::size_t ip = (i + m - 1) % m, in = (i + 1) % m, inn = (i + 2) % m;
                const auto e = dir(i, in);
                const T L = e[2];
                const std::array<T, 2> nu{-e[1], e[0]};
                // at V_i: along e (side i) and f (side i - 1, backward)
                {
                    const auto f = dir(i, ip);
                    const T det = e[0] * f[1] - e[1] * f[0];
                    const T al = (nu[0] * f[1] - nu[1] * f[0]) / det, be = (e[0] * nu[1] - e[1] * nu[0]) / det;
                    const auto se = sides[i].dp(T(0)), sf = sides[ip].dp(T(1));
                    std::array<T, 3> v;
                    for (std::size_t c = 0; c < 3; ++c)
                        v[c] = (al * se[c] - be * sf[c]) / L;
                    const auto q = split(v, sides[i].dp(T(0)), sides[i].x_given(T(0)));
                    ab[i][0] = q[0], ab[i][1] = q[1];
                }
                // at V_{i+1}: along -e (side i, backward) and g (side i + 1)
                {
                    const auto g = dir(in, inn);
                    const std::array<T, 2> me{-e[0], -e[1]};
                    const T det = me[0] * g[1] - me[1] * g[0];
                    const T al = (nu[0] * g[1] - nu[1] * g[0]) / det, be = (me[0] * nu[1] - me[1] * nu[0]) / det;
                    const auto se = sides[i].dp(T(1)), sg = sides[in].dp(T(0));
                    std::array<T, 3> v;
                    for (std::size_t c = 0; c < 3; ++c)
                        v[c] = (-al * se[c] + be * sg[c]) / L;
                    const auto q = split(v, sides[i].dp(T(1)), sides[i].x_given(T(1)));
                    ab[i][2] = q[0], ab[i][3] = q[1];
                }
                if (!(ab[i][1] > T(0)) || !(ab[i][3] > T(0)))
                    throw std::invalid_argument("fill_n_sided: a corner where the hole is not convex, or a cross direction that does not point into it (side " +
                                                std::to_string(i) + ")");
            }
            return ab;
        }

        /**
         * The normal equations of the weighted least squares of `fit`, summed as the rows come
         * rather than from the rows stored (the bending energy alone makes millions of them at 64
         * spans): A^T A in a band of 7 x 7 poles round each pole (a row of a bicubic touches 4 x 4
         * of them), the same for the three coordinates (the rows on one coordinate) and between
         * them (the rows on the normal), and A^T b.
         */
        template <std::floating_point T>
        struct NormalEquations
        {
            std::size_t np, P;
            std::vector<T> same, mixed, atb;
            explicit NormalEquations(std::size_t n_poles) : np{n_poles}, P{n_poles * n_poles}, same(P * 49), mixed(P * 49 * 9), atb(3 * P) {}
            /// The place in the band of the pole b of a row from that of the pole a (4 x 4, u first).
            static auto band(std::size_t a, std::size_t b) { return (b % 4 + 3 - a % 4) * 7 + (b / 4 + 3 - a / 4); }
            auto pole(std::size_t iu, std::size_t iv, std::size_t a) const { return iu + a % 4 + np * (iv + a / 4); }
            /// A row on each coordinate c alike: z . (the poles from (iu, iv)) = value[c].
            void add_each(std::size_t iu, std::size_t iv, const std::array<T, 16> &z, const std::array<T, 3> &value)
            {
                for (std::size_t a = 0; a < 16; ++a)
                {
                    if (z[a] == T(0))
                        continue;
                    const std::size_t i = pole(iu, iv, a);
                    T *row = &same[i * 49];
                    for (std::size_t b = 0; b < 16; ++b)
                        row[band(a, b)] += z[a] * z[b];
                    for (std::size_t c = 0; c < 3; ++c)
                        atb[c * P + i] += z[a] * value[c];
                }
            }
            /// A row on the three coordinates: the sum over c of coef[c] z . (the poles from (iu, iv) of the coordinate c) = 0.
            void add_mixed(std::size_t iu, std::size_t iv, const std::array<T, 16> &z, const std::array<T, 3> &coef)
            {
                for (std::size_t a = 0; a < 16; ++a)
                {
                    if (z[a] == T(0))
                        continue;
                    T *row = &mixed[pole(iu, iv, a) * 49 * 9];
                    for (std::size_t b = 0; b < 16; ++b)
                        for (std::size_t c = 0; c < 3; ++c)
                            for (std::size_t d = 0; d < 3; ++d)
                                row[band(a, b) * 9 + c * 3 + d] += coef[c] * z[a] * coef[d] * z[b];
                }
            }
            /// A^T A, its lower triangle (the one `SimplicialLDLT` reads).
            auto matrix() const -> Eigen::SparseMatrix<T>
            {
                std::vector<Eigen::Triplet<T>> t;
                for (std::size_t i = 0; i < P; ++i)
                {
                    const std::size_t pu = i % np, pv = i / np;
                    for (std::size_t o = 0; o < 49; ++o)
                    {
                        const std::size_t du = o / 7, dv = o % 7;
                        if (pu + du < 3 || pu + du >= np + 3 || pv + dv < 3 || pv + dv >= np + 3)
                            continue;
                        const std::size_t j = pu + du - 3 + np * (pv + dv - 3);
                        for (std::size_t c = 0; c < 3; ++c)
                            for (std::size_t d = 0; d < 3; ++d)
                            {
                                const T v = (c == d ? same[i * 49 + o] : T(0)) + mixed[(i * 49 + o) * 9 + c * 3 + d];
                                if (v != T(0) && c * P + i >= d * P + j)
                                    t.emplace_back(int(c * P + i), int(d * P + j), v);
                            }
                    }
                }
                Eigen::SparseMatrix<T> M(Eigen::Index(3 * P), Eigen::Index(3 * P));
                M.setFromTriplets(t.begin(), t.end());
                return M;
            }
        };

        /**
         * The bicubic B-spline of n x n spans on [-1, 1]^2 of least bending energy along the sides
         * (file header, step 3): its points on the curves and its normal that of the ribbons
         * (weight `side_weight`), its cross derivatives near those given (weight `cross_weight`).
         */
        template <std::floating_point T>
        auto fit(const std::vector<Side<T>> &sides, const std::vector<std::array<T, 2>> &poly, const std::vector<std::array<T, 4>> &ab, T h, std::size_t n,
                 T side_weight, T cross_weight, const GuideSphere<T> *sphere = nullptr, const BSSurface<T, 3> *prev = nullptr) -> BSSurface<T, 3>
        {
            const std::size_t np = n + 3, P = np * np;
            std::vector<T> k(n + 7);
            for (std::size_t i = 0; i < k.size(); ++i)
                k[i] = i < 4 ? T(-1) : i >= n + 3 ? T(1) : T(-1) + T(2) * T(i - 3) / T(n);
            NormalEquations<T> ne(np);
            // the weights of the 4 x 4 poles from (iu, iv) in the derivative of the surface along
            // (a, b) in the domain at (u, v) (order 1; a = b = 0: of the orders (du, dv))
            auto weights = [&](T u, T v, std::size_t du, std::size_t dv, T a, T b) {
                const bool along = a != T(0) || b != T(0);
                const auto [iu, bu] = basis(k, np, u, along ? 0 : du);
                const auto [iv, bv] = basis(k, np, v, along ? 0 : dv);
                std::array<T, 4> bu1{}, bv1{};
                if (along)
                    bu1 = basis(k, np, u, 1).second, bv1 = basis(k, np, v, 1).second;
                std::array<T, 16> z;
                for (std::size_t x = 0; x < 4; ++x)
                    for (std::size_t y = 0; y < 4; ++y)
                        z[x + 4 * y] = along ? a * bu1[x] * bv[y] + b * bu[x] * bv1[y] : bu[x] * bv[y];
                return std::tuple{iu, iv, z};
            };
            // a row on each coordinate, weighted by w
            auto each = [&](T u, T v, std::size_t du, std::size_t dv, T a, T b, T w, const std::array<T, 3> &value) {
                auto [iu, iv, z] = weights(u, v, du, dv, a, b);
                for (auto &e : z)
                    e *= w;
                ne.add_each(iu, iv, z, {w * value[0], w * value[1], w * value[2]});
            };
            // the sides (each row weighted by the square root of its share of the length of its
            // side, of the area for the energy: sums that do not depend on how finely they are
            // sampled)
            const std::size_t per_side = std::max<std::size_t>(24, 6 * n);
            const T ds = std::sqrt(T(1) / T(per_side));
            const T wp = side_weight * ds / h, wn = side_weight * ds / h, wc = cross_weight * ds / h;
            for (std::size_t i = 0; i < sides.size(); ++i)
            {
                const auto &a = poly[i], &b = poly[(i + 1) % poly.size()];
                const T l = std::hypot(b[0] - a[0], b[1] - a[1]);
                const T nx = -(b[1] - a[1]) / l, ny = (b[0] - a[0]) / l;
                for (std::size_t q = 0; q <= per_side; ++q)
                {
                    const T s = T(q) / T(per_side);
                    const auto uv = on_side(poly, i, s);
                    each(uv[0], uv[1], 0, 0, T(0), T(0), wp, sides[i].p(s));
                    const auto nr = sides[i].normal(s);
                    const std::array<T, 3> coef{wn * nr[0], wn * nr[1], wn * nr[2]};
                    for (const auto &[da, db] : {std::pair{T(1), T(0)}, std::pair{T(0), T(1)}})
                    {
                        const auto [iu, iv, z] = weights(uv[0], uv[1], 0, 0, da, db);
                        ne.add_mixed(iu, iv, z, coef);
                    }
                    // (none at the corners, where the cross derivatives of both sides would ask
                    // for twists that seldom agree)
                    if (const T fade = T(16) * s * s * (T(1) - s) * (T(1) - s); wc > T(0) && fade > T(0))
                    {
                        const T ca = ab[i][0] + s * (ab[i][2] - ab[i][0]), cb = ab[i][1] + s * (ab[i][3] - ab[i][1]);
                        const auto dp = sides[i].dp(s), x = sides[i].x_given(s);
                        each(uv[0], uv[1], 0, 0, nx, ny, wc * fade, {ca * dp[0] + cb * x[0], ca * dp[1] + cb * x[1], ca * dp[2] + cb * x[2]});
                    }
                }
            }
            // the bending energy: inside the polygon, integrated by the 6 points rule of degree 4
            // (Dunavant) on the triangles of a fan from the centre, each cut in (2 n)^2 (a rule
            // over the triangles inside only made the energy that of a smaller polygon: the
            // surface 1e-3 off a paraboloid it should find again; a few 1e-6 so); over the
            // square, a much smaller one on a grid, to hold the corners out of the polygon
            auto energy = [&](T u, T v, T w) {
                each(u, v, 2, 0, T(0), T(0), w, {});
                each(u, v, 1, 1, T(0), T(0), w * std::sqrt(T(2)), {});
                each(u, v, 0, 2, T(0), T(0), w, {});
            };
            static constexpr std::array<std::array<double, 3>, 6> rule{{{0.445948490915965, 0.445948490915965, 0.223381589678011},
                                                                        {0.445948490915965, 0.108103018168070, 0.223381589678011},
                                                                        {0.108103018168070, 0.445948490915965, 0.223381589678011},
                                                                        {0.091576213509771, 0.091576213509771, 0.109951743655322},
                                                                        {0.091576213509771, 0.816847572980459, 0.109951743655322},
                                                                        {0.816847572980459, 0.091576213509771, 0.109951743655322}}};
            const std::size_t nt = 2 * n;
            const T inradius = std::cos(std::numbers::pi_v<T> / T(poly.size()));
            for (std::size_t i = 0; i < poly.size(); ++i)
            {
                const auto &a = poly[i], &b = poly[(i + 1) % poly.size()];
                // the triangle (0, a, b) in nt^2 small ones: (j, l) up and down
                auto at = [&](T x, T y) { return std::array<T, 2>{x * a[0] + y * b[0], x * a[1] + y * b[1]}; };
                const T area = (a[0] * b[1] - a[1] * b[0]) / (T(2) * T(nt * nt));
                for (std::size_t j = 0; j < nt; ++j)
                    for (std::size_t l = 0; j + l < nt; ++l)
                        for (int up = 0; up < 2; ++up)
                        {
                            if (up && j + l + 1 >= nt)
                                continue;
                            // the corners of the small triangle in (x, y), x + y <= 1
                            std::array<std::array<T, 2>, 3> c;
                            if (!up)
                                c = {{{T(j), T(l)}, {T(j + 1), T(l)}, {T(j), T(l + 1)}}};
                            else
                                c = {{{T(j + 1), T(l + 1)}, {T(j), T(l + 1)}, {T(j + 1), T(l)}}};
                            for (const auto &q : rule)
                            {
                                const T r0 = T(q[0]), r1 = T(q[1]), r2 = T(1) - r0 - r1;
                                const T x = (r0 * c[0][0] + r1 * c[1][0] + r2 * c[2][0]) / T(nt);
                                const T y = (r0 * c[0][1] + r1 * c[1][1] + r2 * c[2][1]) / T(nt);
                                const auto p = at(x, y);
                                energy(p[0], p[1], std::sqrt(area * T(q[2])) / h);
                                // the guide sphere: the point pulled toward the foot c + r d on the
                                // sphere of the surface before (d the direction from the centre),
                                // each coordinate on its own (a row along d alone would couple
                                // the three coordinates at every point inside: up to 128 spans,
                                // 860 MB and 26 s instead of 340 MB and 7 s), fading to none at
                                // the sides (their normals are those of the ribbons)
                                if (sphere && prev)
                                {
                                    const auto d = unit<T>(prev->value(p[0], p[1]) - sphere->centre);
                                    const T in = std::ranges::min(distances(poly, p)) / inradius;
                                    const T w = sphere->weight * in * in * std::sqrt(area * T(q[2])) / h;
                                    each(p[0], p[1], 0, 0, T(0), T(0), w, sphere->centre + sphere->radius * d);
                                }
                            }
                        }
            }
            const std::size_t g = 2 * n + 2;
            const T wout = T(2) / T(g) * T(1e-3) / h;
            for (std::size_t a = 0; a <= g; ++a)
                for (std::size_t b = 0; b <= g; ++b)
                    energy(T(-1) + T(2) * T(a) / T(g), T(-1) + T(2) * T(b) / T(g), wout);
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver(ne.matrix());
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("fill_n_sided: the least squares system cannot be factorized");
            const Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> atb(ne.atb.data(), Eigen::Index(3 * P));
            const Eigen::Matrix<T, Eigen::Dynamic, 1> x = solver.solve(atb);
            std::vector<std::array<T, 3>> poles(P);
            for (std::size_t i = 0; i < P; ++i)
                poles[i] = {x[Eigen::Index(i)], x[Eigen::Index(P + i)], x[Eigen::Index(2 * P + i)]};
            return BSSurface<T, 3>(poles, k, k, 3, 3);
        }
    }

    /**
     * @brief The surface filling the n-sided hole bounded by `ribbons` (file header): a bicubic
     * B-spline on [-1, 1]^2 to be trimmed by the regular polygon of n sides, within
     * `approx_tol` of the curves and `angular_tol` of the ribbons' normals along them.
     */
    template <std::floating_point T>
    auto fill_n_sided(const std::vector<Ribbon<T>> &ribbons, const NSidedOptions<T> &o = {}) -> NSidedSurface<T>
    {
        using namespace nsided_detail;
        const std::size_t m = ribbons.size();
        if (m < 3)
            throw std::invalid_argument("fill_n_sided: fewer than three ribbons");
        auto finite = [](T v) { return std::isfinite(v); };
        if (!(o.approx_tol > T(0)) || !(o.angular_tol > T(0)) || !(o.side_weight > T(0)) || !(o.cross_weight >= T(0)) || !(o.shape_tol >= T(0)) ||
            !finite(o.approx_tol) || !finite(o.angular_tol) || !finite(o.side_weight) || !finite(o.cross_weight) || !finite(o.shape_tol))
            throw std::invalid_argument("fill_n_sided: tolerances or side_weight not positive, cross_weight or shape_tol negative, or one of them not finite");
        if (o.sphere && (!(o.sphere->radius > T(0)) || !(o.sphere->weight >= T(0)) || !finite(o.sphere->radius) || !finite(o.sphere->weight) ||
                         !finite(o.sphere->centre[0]) || !finite(o.sphere->centre[1]) || !finite(o.sphere->centre[2])))
            throw std::invalid_argument("fill_n_sided: a guide sphere whose radius is not positive, weight negative, or a value not finite");
        // (the shape is compared from one doubling of the spans to the next: two at least)
        if (o.max_spans < 8)
            throw std::invalid_argument("fill_n_sided: fewer than 8 spans allowed");
        std::vector<Side<T>> sides;
        T perimeter{};
        for (const auto &r : ribbons)
        {
            if (!r.curve || !r.cross || !(r.t1 != r.t0))
                throw std::invalid_argument("fill_n_sided: a ribbon without curve, cross direction or interval");
            sides.push_back(Side<T>{&r});
        }
        // each side (sampled), then each corner
        std::vector<T> len(m);
        for (std::size_t i = 0; i < m; ++i)
        {
            const auto &a = sides[i];
            auto prev = a.p(T(0));
            for (int q = 0; q <= 32; ++q)
            {
                const T s = T(q) / T(32);
                const auto cur = a.p(s);
                len[i] += norm(cur - prev);
                prev = cur;
                const auto d = a.dp(s);
                if (const T l = norm(d); !(l > T(0)) || !finite(l))
                    throw std::invalid_argument("fill_n_sided: a curve whose derivative vanishes (side " + std::to_string(i) + ")");
                if (!(norm(cross(unit<T>(d), a.x(s))) >= T(1e-6)))
                    throw std::invalid_argument("fill_n_sided: a cross direction null or along its curve (side " + std::to_string(i) + ")");
            }
            perimeter += len[i];
        }
        for (std::size_t i = 0; i < m; ++i)
        {
            const auto &a = sides[i], &b = sides[(i + 1) % m];
            // (beyond twice approx_tol, no point of the surface at the corner is near both)
            if (const T gap = norm(a.p(T(1)) - b.p(T(0))); !(gap <= std::min(T(1e-6) * (T(1) + len[i]), T(2) * o.approx_tol)))
                throw std::invalid_argument(std::format("fill_n_sided: the sides do not follow one another (side {}: {:.3g} from its end to the start of the next)", i, gap));
            // (|n_a - n_b| = 2 sin(angle / 2): the least the larger angle of the normal of the surface
            // to both can be)
            if (const T gap = norm(a.normal(T(1)) - b.normal(T(0))) / T(2); !(gap <= o.angular_tol))
                throw std::invalid_argument(std::format("fill_n_sided: the ribbons' normals disagree at the end of side {} (sine of half their angle {:.3g}, "
                                                        "above angular_tol)",
                                                        i, gap));
        }
        const auto poly = regular_polygon<T>(m);
        const auto ab = cross_derivatives(sides, poly);
        const T h = perimeter / (T(m) * norm(std::array<T, 2>{poly[1][0] - poly[0][0], poly[1][1] - poly[0][1]}));
        // the surface inside, sampled on a grid of the polygon, to tell when doubling the spans
        // no longer changes it
        std::vector<std::array<T, 2>> probes;
        for (int a = 0; a <= 24; ++a)
            for (int b = 0; b <= 24; ++b)
                if (const std::array<T, 2> p{T(-1) + T(a) / T(12), T(-1) + T(b) / T(12)}; std::ranges::min(distances(poly, p)) > T(0))
                    probes.push_back(p);
        // (the size of the hole: the largest distance between two of its corners)
        T size{};
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = i + 1; j < m; ++j)
                size = std::max(size, norm(sides[i].p(T(0)) - sides[j].p(T(0))));
        const T shape_tol = o.shape_tol > T(0) ? o.shape_tol : std::max(T(10) * o.approx_tol, T(2e-3) * size);
        T dev{}, ang{}, change{}, last{std::numeric_limits<T>::infinity()}, error{};
        std::size_t spans{};
        std::vector<std::array<T, 3>> before;
        for (std::size_t n = 4; n <= o.max_spans; n *= 2)
        {
            auto srf = fit(sides, poly, ab, h, n, o.side_weight, o.cross_weight);
            // leaning toward the guide sphere: its radial gap linearized about the surface before,
            // twice
            if (o.sphere)
                for (int it = 0; it < 2; ++it)
                    srf = fit(sides, poly, ab, h, n, o.side_weight, o.cross_weight, &*o.sphere, &srf);
            const std::size_t samples = o.side_samples ? o.side_samples : std::max<std::size_t>(256, 64 * n);
            std::tie(dev, ang) = measure(srf, sides, poly, samples);
            // (along the normal: a slide of the parameters along the surface does not change it)
            std::vector<std::array<T, 3>> now;
            change = before.empty() ? std::numeric_limits<T>::infinity() : T(0);
            for (std::size_t q = 0; q < probes.size(); ++q)
            {
                const auto &p = probes[q];
                now.push_back(srf.value(p[0], p[1]));
                if (!before.empty())
                {
                    const auto nr = unit<T>(cross(srf.value(p[0], p[1], 1, 0), srf.value(p[0], p[1], 0, 1)));
                    change = std::max(change, std::abs((now.back() - before[q]) * nr));
                }
            }
            before = std::move(now);
            spans = n;

            // the error of this surface inside: from the moves of the last two doublings, falling
            // geometrically by rho, what is left to come, change rho / (1 - rho) (the move itself
            // from the first doubling, or where they do not fall)
            const T rho = std::isfinite(last) && last > T(0) ? change / last : T(1);
            error = rho < T(1) ? change * rho / (T(1) - rho) : change;
            last = change;
            if (dev <= o.approx_tol && ang <= o.angular_tol && std::min(change, error) <= shape_tol)
                return {std::make_shared<BSSurface<T, 3>>(std::move(srf)), poly, dev, ang, n};
        }
        throw std::runtime_error(std::format("fill_n_sided: the tolerances are not reached with {} spans ({:.3g} in position for {:.3g}, {:.3g} in angle "
                                             "for {:.3g}, the shape moved by {:.3g} for {:.3g} from {} spans)",
                                             spans, dev, o.approx_tol, ang, o.angular_tol, change, shape_tol, spans / 2));
    }

}
