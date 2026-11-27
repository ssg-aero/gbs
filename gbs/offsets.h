#pragma once
/**
 * @file offsets.h
 * @brief Offset surfaces: S + d n with exact derivatives, the exact offset of a primitive, the
 * fold of an offset (d beyond a radius of curvature), and the B-spline approximation of an offset.
 *
 * Notation. S(u, v) the basis, N = S_u x S_v, w = |N|, n = N / w; the offset is P = S + d n.
 * Writing X^(a,b) for the derivative of order a in u and b in v and C(a, i) for the binomial
 * coefficients, Leibniz's rule on the products N = S_u x S_v, w² = N . N and N = w n gives, by
 * increasing a + b,
 *
 *   N^(a,b) = sum_{i<=a, j<=b} C(a,i) C(b,j) S^(i+1,j) x S^(a-i,b-j+1)
 *   w^(a,b) = ( (N . N)^(a,b) - sum' C(a,i) C(b,j) w^(i,j) w^(a-i,b-j) ) / 2w   (sum' without (0,0), (a,b))
 *   n^(a,b) = ( N^(a,b) - sum_{(i,j) != (0,0)} C(a,i) C(b,j) w^(i,j) n^(a-i,b-j) ) / w
 *   P^(a,b) = S^(a,b) + d n^(a,b)
 *
 * so that P^(a,b) takes the derivatives of S up to (a + 1, b + 1) and no others (each (a, b)
 * needs only the orders below it in both variables), and is exact (the recurrence of the
 * derivatives of a quotient, as for the rational B-splines: L. Piegl, W. Tiller, *The NURBS
 * Book*, 2nd ed., § 4.5; the offset surfaces themselves: T. Maekawa, *An overview of offset
 * curves and surfaces*, CAD 31(3), 1999). The recurrence is that of the offset curves
 * (`curveoffset.h`, `offset_detail::unit_derivatives`), in two variables.
 *
 * The offset is regular where 1 - d k > 0 for both principal curvatures k of S (k counted
 * positive where S bends toward n, II = S_uu . n: -1/R on a sphere of outward normal), since
 * P_u x P_v = (1 - d k1)(1 - d k2) N; it folds beyond (Maekawa, § 3).
 */
#include <gbs/bssapprox.h>
#include <gbs/bssanalysis.h>
#include <gbs/bssurf.h>
#include <gbs/curveoffset.h>
#include <gbs/execution.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    template <typename T>
    auto offset_points(const BSSfunction<T> &f_ep, const Surface<T, 3> &cl_srf, const std::vector<T> &u, const std::vector<T> &v)
    {
        const size_t nu = u.size();
        const size_t nv = v.size();
        // Flatten the (u,v) grid in the existing U-first, V-outer order, then run
        // the offset kernel (value + du + dv + normal + offset) as one size-gated
        // parallel transform writing each node to its own slot (#89). Each node is
        // an independent surface eval, so the result is bit-identical to the old
        // nested push_back loop — just no longer serialized by a growing container.
        std::vector<std::array<T, 2>> uv(nu * nv);
        for (size_t j = 0; j < nv; ++j)
            for (size_t i = 0; i < nu; ++i)
                uv[j * nu + i] = {u[i], v[j]};

        points_vector<T, 3> pts(nu * nv);
        transform_threshold(
            uv.begin(), uv.end(), pts.begin(),
            [&f_ep, &cl_srf](const std::array<T, 2> &uv_) -> point<T, 3>
            {
                const T u_ = uv_[0], v_ = uv_[1];
                auto pt = cl_srf(u_, v_);
                auto tu = cl_srf(u_, v_, 1, 0);
                auto tv = cl_srf(u_, v_, 0, 1);
                auto n = tu ^ tv;
                n = n / norm(n);
                return pt + n * f_ep(u_, v_);
            });
        return pts;
    }

    template <typename T, bool rational>
    auto offset_approx(const BSSfunction<T> &f_ep, const BSSurfaceGeneral<T, 3, rational> &cl_srf,size_t nu , size_t nv )
    {
        auto [u1, u2, v1, v2] = cl_srf.bounds();
        auto u = make_range<T>(u1, u2, nu);
        auto v = make_range<T>(v1, v2, nv);

        auto pts = offset_points(f_ep,cl_srf,u,v);

        auto ku_flat = cl_srf.knotsFlatsU();
        auto kv_flat = cl_srf.knotsFlatsV();
        auto p = cl_srf.degreeU();
        auto q = cl_srf.degreeV();
        auto poles = approx(pts, ku_flat, kv_flat, u, v, p, q);

        using bs_type = typename std::conditional<rational, BSSurfaceRational<T, 3>, BSSurface<T, 3>>::type;

        return bs_type{
            poles,
            ku_flat,
            kv_flat,
            p,
            q};
    }

    namespace offset_detail
    {
        /// The lengths of the bounds of s (1 where infinite).
        template <std::floating_point T>
        auto lengths(const Surface<T, 3> &s) -> std::array<T, 2>
        {
            const auto b = s.bounds();
            auto len = [](T x, T y) { return std::isfinite(y - x) && y > x ? y - x : T(1); };
            return {len(b[0], b[1]), len(b[2], b[3])};
        }

        /**
         * Is (S_u, S_v) a frame of the tangent plane? Compared as S_u lu and S_v lv (lu, lv the
         * lengths of the bounds), neither is shorter than 1e-10 of the other (at a pole, along a
         * collapsed row of poles, the shorter one is the rounding of a null derivative) and the sine
         * of their angle is above 1e-10.
         */
        template <std::floating_point T>
        bool regular(const point<T, 3> &su, const point<T, 3> &sv, T lu, T lv)
        {
            const T a = norm(su) * lu, b = norm(sv) * lv, l = norm(cross(su, sv)) * lu * lv;
            return std::isfinite(l) && std::min(a, b) > T(1e-10) * std::max(a, b) && l > T(1e-10) * a * b;
        }

        /**
         * `regular`, and away from a degenerate point: neither S_u nor S_v would vanish within 1e-4
         * of the bounds (|S_u| > 1e-4 |S_uv| lv, |S_v| > 1e-4 |S_uv| lu). The normal and the
         * curvatures divide by |S_u x S_v|: at a distance e of a pole (a fraction of the bounds) the
         * derivatives of order k of n lose about k e digits of e (rounding over e^k), so that they
         * are taken there by extrapolation instead.
         */
        template <std::floating_point T>
        bool clear(const point<T, 3> &su, const point<T, 3> &sv, const point<T, 3> &suv, T lu, T lv)
        {
            return regular(su, sv, lu, lv) && norm(su) > T(1e-4) * norm(suv) * lv && norm(sv) > T(1e-4) * norm(suv) * lu;
        }

        /// A buffer of n values: on the stack up to Cap, on the heap beyond.
        template <typename X, std::size_t Cap>
        class Scratch
        {
            std::array<X, Cap> a_{};
            std::vector<X> v_;
            std::span<X> s_;

        public:
            explicit Scratch(std::size_t n)
            {
                if (n <= Cap)
                    s_ = std::span<X>{a_.data(), n};
                else
                {
                    v_.resize(n);
                    s_ = std::span<X>{v_};
                }
            }
            Scratch(const Scratch &) = delete;
            Scratch &operator=(const Scratch &) = delete;
            auto span() -> std::span<X> { return s_; }
            auto data() -> X * { return s_.data(); }
            auto operator[](std::size_t i) -> X & { return s_[i]; }
        };

        /**
         * Exact derivatives of S + d n on the rectangle a <= A, b <= B at (u, v), into P (row-major,
         * B + 1 per row), from those of S up to (A + 1, B + 1) in one pass; false near a
         * degenerate point (not `clear`).
         */
        template <std::floating_point T>
        bool offset_derivatives(const Surface<T, 3> &s, T d, T u, T v, std::size_t A, std::size_t B, std::span<point<T, 3>> P)
        {
            const std::size_t m = B + 2, q = B + 1;
            Scratch<point<T, 3>, 25> S((A + 2) * m);
            s.derivatives(u, v, A + 1, B + 1, S.data());
            auto at = [&](std::size_t a, std::size_t b) -> const point<T, 3> & { return S[a * m + b]; };
            const auto [lu, lv] = lengths(s);
            if (!clear(at(1, 0), at(0, 1), at(1, 1), lu, lv))
                return false;
            Scratch<point<T, 3>, 16> N((A + 1) * q), n((A + 1) * q);
            Scratch<T, 16> w((A + 1) * q);
            for (std::size_t a = 0; a <= A; ++a)
                for (std::size_t b = 0; b <= B; ++b)
                {
                    point<T, 3> x{};
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                            x = x + (binomial<T>(a, i) * binomial<T>(b, j)) * cross(at(i + 1, j), at(a - i, b - j + 1));
                    N[a * q + b] = x;
                }
            unit_derivatives<T, 3>(N.span(), A, B, w.span(), n.span());
            for (std::size_t a = 0; a <= A; ++a)
                for (std::size_t b = 0; b <= B; ++b)
                    P[a * q + b] = at(a, b) + d * n[a * q + b];
            return true;
        }
    }

    /**
     * @brief The offset S + d n of a surface S, n = S_u x S_v / |S_u x S_v| (a negative d goes the
     * other way), on the parameters of S.
     *
     * The derivatives are exact, P^(a,b) from those of S up to (a + 1, b + 1) (file header): the
     * basis must give them when it has second derivatives; a basis without second derivatives
     * gives the first ones by central differences, and none beyond. At a degenerate point (a pole,
     * a collapsed row of poles, an apex) n is its limit (`normal_direction`); there and next to
     * it (not `offset_detail::clear`: within 1e-4 of the bounds, where the recurrences lose their
     * digits), the derivatives are extrapolated (cubic) from four points 5e-4 of the bounds apart
     * on the way to the centre of the bounds, or to a corner where that way stays near the
     * degenerate line (then 5e-3, 5e-2 apart): within about 1e-7 of their limits.
     */
    template <std::floating_point T>
    class OffsetSurface : public Surface<T, 3>
    {
        std::shared_ptr<Surface<T, 3>> s_;
        T d_;

        /// The derivatives on the rectangle (A, B) at (u, v), extrapolated near a degenerate point.
        void exact(T u, T v, std::size_t A, std::size_t B, std::span<point<T, 3>> P) const
        {
            if (offset_detail::offset_derivatives(*s_, d_, u, v, A, B, P))
                return;
            const auto b = bounds();
            const auto [lu, lv] = offset_detail::lengths(*s_);
            const std::array<std::array<T, 2>, 5> targets{{{(b[0] + b[1]) / T(2), (b[2] + b[3]) / T(2)},
                                                           {b[0], b[2]}, {b[1], b[2]}, {b[0], b[3]}, {b[1], b[3]}}};
            const std::size_t n = P.size();
            offset_detail::Scratch<point<T, 3>, 64> Q(4 * n);
            for (T h : {T(5e-4), T(5e-3), T(5e-2)})
                for (const auto &[tu, tv] : targets)
                {
                    const T su = std::clamp(tu - u, -h * lu, h * lu), sv = std::clamp(tv - v, -h * lv, h * lv);
                    if (!std::isfinite(su) || !std::isfinite(sv) || (su == T(0) && sv == T(0)))
                        continue;
                    bool ok = true;
                    for (std::size_t k = 0; k < 4 && ok; ++k)
                        ok = offset_detail::offset_derivatives(*s_, d_, u + T(k + 1) * su, v + T(k + 1) * sv, A, B, Q.span().subspan(k * n, n));
                    if (!ok)
                        continue;
                    // the cubic through the four points, at 0
                    for (std::size_t i = 0; i < n; ++i)
                        P[i] = T(4) * Q[i] - T(6) * Q[n + i] + T(4) * Q[2 * n + i] - Q[3 * n + i];
                    return;
                }
            throw std::domain_error("OffsetSurface: no normal near the point");
        }

        auto central(T u, T v, std::size_t du) const -> point<T, 3>
        {
            const auto b = bounds();
            const std::size_t i = du == 1 ? 0 : 2;
            const T h = T(1e-6) * offset_detail::lengths(*s_)[i / 2];
            T x = du == 1 ? u : v;
            const T a = std::max(b[i], x - h), c = std::min(b[i + 1], x + h);
            return du == 1 ? (value(c, v) - value(a, v)) / (c - a) : (value(u, c) - value(u, a)) / (c - a);
        }

    public:
        OffsetSurface(std::shared_ptr<Surface<T, 3>> s, T d) : s_{std::move(s)}, d_{d}
        {
            if (!s_)
                throw std::invalid_argument("OffsetSurface: null basis");
        }
        [[nodiscard]] auto basis() const -> const std::shared_ptr<Surface<T, 3>> & { return s_; }
        [[nodiscard]] T distance() const { return d_; }
        auto bounds() const -> std::array<T, 4> override { return s_->bounds(); }
        auto has_second_derivatives() const -> bool override { return s_->has_second_derivatives(); }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            if (du + dv == 0)
                return s_->value(u, v) + d_ * normal_direction(*s_, u, v);
            if (!s_->has_second_derivatives())
            {
                if (du + dv > 1)
                    throw std::runtime_error("OffsetSurface: only first derivatives on a basis without second ones");
                return central(u, v, du);
            }
            offset_detail::Scratch<point<T, 3>, 16> P((du + 1) * (dv + 1));
            exact(u, v, du, dv, P.span());
            return P[du * (dv + 1) + dv];
        }
        auto derivatives(T u, T v, size_t nu, size_t nv, std::array<T, 3> *SKL) const -> void override
        {
            if (!s_->has_second_derivatives() || nu + nv == 0)
            {
                Surface<T, 3>::derivatives(u, v, nu, nv, SKL);
                return;
            }
            exact(u, v, nu, nv, std::span<point<T, 3>>{SKL, (nu + 1) * (nv + 1)});
            SKL[0] = value(u, v);
        }
    };

    /// The smallest of 1 - d k over the principal curvatures k of a surface on a box, and where.
    template <std::floating_point T>
    struct OffsetMargin
    {
        T margin{std::numeric_limits<T>::infinity()};
        std::array<T, 2> uv{};
    };

    namespace offset_detail
    {
        /// 1 - d k for the principal curvature k that brings it lowest at (u, v); nullopt near a
        /// degenerate point (not `clear`), where the curvatures lose their digits.
        template <std::floating_point T>
        auto margin_at(const Surface<T, 3> &s, T d, T u, T v) -> std::optional<T>
        {
            point<T, 3> su, sv, suu, suv, svv;
            if (s.has_second_derivatives())
            {
                std::array<point<T, 3>, 9> S;
                s.derivatives(u, v, 2, 2, S.data());
                su = S[3], sv = S[1], suu = S[6], suv = S[4], svv = S[2];
            }
            else
            {
                const auto b = s.bounds();
                const T hu = T(1e-5) * (b[1] - b[0]), hv = T(1e-5) * (b[3] - b[2]);
                su = s.value(u, v, 1, 0), sv = s.value(u, v, 0, 1);
                const T u1 = std::max(b[0], u - hu), u2 = std::min(b[1], u + hu), v1 = std::max(b[2], v - hv), v2 = std::min(b[3], v + hv);
                suu = (s.value(u2, v, 1, 0) - s.value(u1, v, 1, 0)) / (u2 - u1);
                suv = (s.value(u, v2, 1, 0) - s.value(u, v1, 1, 0)) / (v2 - v1);
                svv = (s.value(u, v2, 0, 1) - s.value(u, v1, 0, 1)) / (v2 - v1);
            }
            const auto [lu, lv] = lengths(s);
            if (!clear(su, sv, suv, lu, lv))
                return std::nullopt;
            const auto N = cross(su, sv);
            const T w = norm(N);
            const auto n = N / w;
            // the second fundamental form on an orthonormal basis of the tangent plane, e1 along
            // S_u: S_u = a e1, S_v = b e1 + c e2; symmetric, its eigenvalues are the principal
            // curvatures, their half difference a hypot (no cancellation at an umbilic)
            const T a = norm(su), b = (sv * su) / a, c = w / a;
            const T L = suu * n, M = suv * n, Nn = svv * n, t = b / a;
            const T p = L / (a * a), q = (M - t * L) / (a * c), r = (Nn - T(2) * t * M + t * t * L) / (c * c);
            const T H = (p + r) / T(2), h = std::hypot((p - r) / T(2), q);
            return T(1) - d * (d > T(0) ? H + h : H - h);
        }

        /// Sample parameters on [a, b]: 32 even intervals, and 4 per knot span of a B-spline.
        template <std::floating_point T>
        auto samples(T a, T b, const std::vector<T> *knots) -> std::vector<T>
        {
            std::vector<T> x;
            for (int i = 0; i <= 32; ++i)
                x.push_back(a + (b - a) * T(i) / T(32));
            if (knots)
            {
                std::vector<T> k;
                for (T t : *knots)
                    if (t >= a && t <= b)
                        k.push_back(t);
                k.push_back(a);
                k.push_back(b);
                std::ranges::sort(k);
                const auto [first, last] = std::ranges::unique(k);
                k.erase(first, last);
                for (std::size_t i = 0; i + 1 < k.size(); ++i)
                    for (int j = 0; j < 4; ++j)
                        x.push_back(k[i] + (k[i + 1] - k[i]) * T(j) / T(4));
            }
            std::ranges::sort(x);
            const auto [first, last] = std::ranges::unique(x);
            x.erase(first, last);
            return x;
        }
    }

    /**
     * @brief The flat knots in u and in v of a B-spline surface (rational or not), or of the basis of an offset
     * surface (its parameters), where the surface may not be smooth; nullopt in both directions for any other
     * surface.
     */
    template <std::floating_point T>
    [[nodiscard]] auto surface_flat_knots(const Surface<T, 3> &s) -> std::array<std::optional<std::vector<T>>, 2>
    {
        if (const auto *o = dynamic_cast<const OffsetSurface<T> *>(&s))
            return surface_flat_knots(*o->basis());
        if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(&s))
            return {b->knotsFlatsU(), b->knotsFlatsV()};
        if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
            return {b->knotsFlatsU(), b->knotsFlatsV()};
        return {};
    }

    /**
     * @brief Where the offset of s at the distance d comes nearest to fold on `box`: the smallest
     * of 1 - d k over the principal curvatures k (file header), regular while it is positive.
     *
     * Sampled (33 even values per direction, and 4 per knot span of a B-spline or of the basis of
     * an offset surface), then refined
     * around the 4 lowest samples by a pattern search down to 1e-6 of a sample step. The points
     * within 1e-4 of the bounds of a degenerate point (a pole, an apex) are skipped, the curvatures
     * being continuous there. A fold narrower than a sample step around a
     * sample that is not among the lowest may be missed.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_margin(const Surface<T, 3> &s, T d, const std::array<T, 4> &box) -> OffsetMargin<T>
    {
        using namespace offset_detail;
        const auto k = surface_flat_knots(s);
        const auto us = samples(box[0], box[1], k[0] ? &*k[0] : nullptr), vs = samples(box[2], box[3], k[1] ? &*k[1] : nullptr);
        std::vector<std::pair<T, std::array<T, 2>>> best;
        for (T u : us)
            for (T v : vs)
                if (const auto m = margin_at(s, d, u, v))
                    best.push_back({*m, {u, v}});
        OffsetMargin<T> r;
        if (best.empty())
            return r;
        const auto n = std::min<std::size_t>(4, best.size());
        std::ranges::partial_sort(best, best.begin() + std::ptrdiff_t(n), {}, &std::pair<T, std::array<T, 2>>::first);
        const T hu0 = (box[1] - box[0]) / T(32), hv0 = (box[3] - box[2]) / T(32);
        for (std::size_t i = 0; i < n; ++i)
        {
            auto [m, uv] = best[i];
            for (T h = T(1); h > T(1e-6); h /= T(2))
            {
                bool moved = true;
                while (moved)
                {
                    moved = false;
                    for (const auto &[du, dv] : {std::pair{1, 0}, std::pair{-1, 0}, std::pair{0, 1}, std::pair{0, -1}})
                    {
                        const T u = std::clamp(uv[0] + T(du) * h * hu0, box[0], box[1]), v = std::clamp(uv[1] + T(dv) * h * hv0, box[2], box[3]);
                        if (const auto x = margin_at(s, d, u, v); x && *x < m)
                        {
                            m = *x;
                            uv = {u, v};
                            moved = true;
                        }
                    }
                }
            }
            if (m < r.margin)
                r = {m, uv};
        }
        return r;
    }

    /// An offset surface and the shift of its v: S'(u, v + dv) = S(u, v) + d n (a cone moves its v).
    template <std::floating_point T>
    struct ExactOffset
    {
        std::shared_ptr<Surface<T, 3>> surface;
        T dv{};
    };

    namespace offset_detail
    {
        /// The least and the largest cosine over [a, b].
        template <std::floating_point T>
        auto cos_range(T a, T b) -> std::array<T, 2>
        {
            constexpr T tp = T(2) * std::numbers::pi_v<T>, pi = std::numbers::pi_v<T>;
            std::array<T, 2> r{std::min(std::cos(a), std::cos(b)), std::max(std::cos(a), std::cos(b))};
            if (!(b - a < tp))
                return {T(-1), T(1)};
            if (std::ceil(a / tp) * tp <= b)
                r[1] = T(1);
            if (std::ceil((a - pi) / tp) * tp + pi <= b)
                r[0] = T(-1);
            return r;
        }
    }

    /**
     * @brief The regular part of the offset at the distance d of a torus of radii R and r (r + d > 0,
     * the normal cos v X(u) + sin v Z where R + r cos v > 0): the v where cos v is above the value
     * returned, max(-R / r, -R / (r + d)), where the circles of centres of the torus and of its
     * offset, R + r cos v and R + (r + d) cos v, stay off the axis; below -1, every v.
     *
     * primitive_offset requires a box strictly inside it; the STEP reader admits its bounds, where
     * the offset closes onto the axis.
     */
    template <std::floating_point T>
    [[nodiscard]] auto torus_offset_cos_min(T R, T r, T d) -> T
    {
        return std::max(-R / r, -R / (r + d));
    }

    /**
     * @brief The apexes of a cone of radius R at v = 0 and half-angle a (radius R + v tan a, v along
     * the axis) and of its offset at the distance d at the same v (radius R + d cos a + v tan a):
     * the v where each radius vanishes, {-R / tan a, -(R + d cos a) / tan a}. The offset is regular
     * on the side of both, v beyond the larger for tan a > 0, below the smaller otherwise.
     *
     * primitive_offset admits the apex of the cone on its box, not that of the offset; the STEP
     * reader admits both.
     */
    template <std::floating_point T>
    [[nodiscard]] auto cone_offset_apexes(T R, T a, T d) -> std::array<T, 2>
    {
        const T t = std::tan(a);
        return {-R / t, -(R + d * std::cos(a)) / t};
    }

    /**
     * @brief The exact offset of a primitive: plane moved, cylinder, sphere and cone of another
     * radius, torus of another minor radius, with the shift of v of the cone; checked against
     * S + d n at 4 points within tol. nullopt where the offset folds on `box` (the bounds of s by
     * default; the surface returned keeps the bounds of s): a cylinder or a sphere beyond its
     * radius, a torus whose tube or whose circle of centres R + r cos v crosses the axis there, a
     * cone whose radius crosses its apex there (or a box across the apex, two nappes).
     */
    template <std::floating_point T>
    [[nodiscard]] auto primitive_offset(const ElementarySurface<T> &s, T d, T tol, std::optional<std::array<std::type_identity_t<T>, 4>> box = std::nullopt)
        -> std::optional<ExactOffset<T>>
    {
        const auto &f = s.frame();
        const auto sb = s.bounds();
        const auto bx = box.value_or(sb);
        std::vector<std::pair<std::shared_ptr<ElementarySurface<T>>, T>> cands;
        if (dynamic_cast<const Plane<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                auto g = f;
                g.o = f.o + sg * d * f.z;
                cands.push_back({std::make_shared<Plane<T>>(g, sb), T(0)});
            }
        else if (const auto *c = dynamic_cast<const CylindricalSurface<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                if (c->radius() + sg * d > T(0))
                    cands.push_back({std::make_shared<CylindricalSurface<T>>(c->radius() + sg * d, f, sb), T(0)});
            }
        else if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                if (sp->radius() + sg * d > T(0))
                    cands.push_back({std::make_shared<SphericalSurface<T>>(sp->radius() + sg * d, f, sb), T(0)});
            }
        else if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(&s))
        {
            // the circle of centres R + r cos v stays off the axis on the box, for the torus and its offset
            const T R = to->major_radius(), r = to->minor_radius();
            const auto [c1, c2] = offset_detail::cos_range(bx[2], bx[3]);
            for (T sg : {T(1), T(-1)})
            {
                const T r2 = r + sg * d;
                if (r2 > T(0) && c1 > torus_offset_cos_min(R, r, sg * d))
                    cands.push_back({std::make_shared<ToroidalSurface<T>>(R, r2, f, sb), T(0)});
            }
        }
        else if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(&s))
        {
            // the radius R + v tan a on the box, and R' + v' tan a on the box shifted by dv, stay on
            // one side of the apex (the cone at its apex: 0); the reference radius of the offset at
            // v' = 0 may be negative: the frame then moves along the axis to the start of the box
            const T a = co->semi_angle(), t = std::tan(a), R = co->radius();
            for (T sg : {T(1), T(-1)})
            {
                // the near end of the box (the radii grow away from the apexes): on the apex of the
                // cone, not on that of the offset
                const auto [v0, v0d] = cone_offset_apexes(R, a, sg * d);
                const T vn = t > T(0) ? bx[2] : bx[3], side = t > T(0) ? T(1) : T(-1);
                if (!(side * (vn - v0) >= T(0)) || !(side * (vn - v0d) > T(0)))
                    continue;
                const T dv = -sg * d * std::sin(a), R2 = R + sg * d / std::cos(a);
                const T h = R2 >= T(0) ? T(0) : bx[2] + dv; // the new origin, along the axis
                auto g = f;
                g.o = f.o + h * f.z;
                cands.push_back({std::make_shared<ConicalSurface<T>>(R2 + h * t, a, g, std::array<T, 4>{sb[0], sb[1], sb[2] + dv - h, sb[3] + dv - h}), dv - h});
            }
        }
        // the candidate where S'(u, v + dv) = S(u, v) + d n
        for (const auto &[c, dv] : cands)
        {
            bool ok = true;
            for (T fu : {T(0.3), T(0.7)})
                for (T fv : {T(0.4), T(0.6)})
                {
                    const T u = bx[0] + fu * (bx[1] - bx[0]), v = bx[2] + fv * (bx[3] - bx[2]);
                    const auto n = cross(s.value(u, v, 1, 0), s.value(u, v, 0, 1));
                    const auto target = s.value(u, v) + (d / norm(n)) * n;
                    ok = ok && norm(c->value(u, v + dv) - target) <= tol;
                }
            if (ok)
                return ExactOffset<T>{c, dv};
        }
        return std::nullopt;
    }

    /// The margin at or under which `offset_surface` takes the offset of a surface that is not a primitive to fold.
    template <std::floating_point T>
    inline constexpr T offset_fold = T(1e-9);

    /**
     * @brief The offset of s at the distance d on `box` (its bounds by default), exact when it
     * can be: a primitive gives a primitive (`primitive_offset`); an offset surface the offset of
     * its basis at the sum of the distances, while its own offset does not fold (the normals then
     * agree), exact for a primitive basis; any other surface an `OffsetSurface`. nullopt where the
     * offset folds on the box: a primitive beyond its radius of curvature, another surface where
     * `offset_margin` is at most `fold`. `found`, if given, receives the margin of that other
     * surface (the search made, so that a caller refused reads the fold on it); it is left as it
     * is for a primitive and an offset surface.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_surface(const std::shared_ptr<Surface<std::type_identity_t<T>, 3>> &s, T d, std::type_identity_t<T> tol,
                                      std::optional<std::array<std::type_identity_t<T>, 4>> box = std::nullopt,
                                      std::type_identity_t<T> fold = offset_fold<T>, OffsetMargin<std::type_identity_t<T>> *found = nullptr)
        -> std::optional<ExactOffset<T>>
    {
        if (!s)
            return std::nullopt;
        if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(s.get()))
            return primitive_offset(*e, d, tol, box);
        if (const auto *o = dynamic_cast<const OffsetSurface<T> *>(s.get()))
        {
            if (!offset_surface<T>(o->basis(), o->distance(), tol, box, fold))
                return std::nullopt;
            return offset_surface<T>(o->basis(), o->distance() + d, tol, box, fold);
        }
        const auto margin = offset_margin(*s, d, box.value_or(s->bounds()));
        if (found)
            *found = margin;
        if (!(margin.margin > fold))
            return std::nullopt;
        return ExactOffset<T>{std::make_shared<OffsetSurface<T>>(s, d), T(0)};
    }
}
