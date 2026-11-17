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
 * so that a derivative of order k of P takes the derivatives of S up to k + 1, and is exact (the
 * recurrence of the derivatives of a quotient, as for the rational B-splines: L. Piegl, W. Tiller,
 * *The NURBS Book*, 2nd ed., § 4.5; the offset surfaces themselves: T. Maekawa, *An overview of
 * offset curves and surfaces*, CAD 31(3), 1999).
 *
 * The offset is regular where 1 - d k > 0 for both principal curvatures k of S (k counted
 * positive where S bends away from n: -1/R on a sphere of outward normal), since
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
#include <optional>
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
        /**
         * Derivatives of n = N / |N| up to the order k, from those of N (N[a (k + 1) + b] = N^(a,b)
         * for a + b <= k, the others unused), by the recurrences of the file header; nullopt where
         * |N| <= small.
         */
        template <std::floating_point T, std::size_t dim>
        auto unit_derivatives(const std::vector<point<T, dim>> &N, std::size_t k, T small) -> std::optional<std::vector<point<T, dim>>>
        {
            const std::size_t q = k + 1;
            std::vector<T> w(q * q, T(0));
            std::vector<point<T, dim>> n(q * q, point<T, dim>{});
            const T w0 = norm(N[0]);
            if (!(w0 > small) || !std::isfinite(w0))
                return std::nullopt;
            for (std::size_t o = 0; o <= k; ++o)
                for (std::size_t a = 0; a <= o; ++a)
                {
                    const std::size_t b = o - a;
                    if (o == 0)
                    {
                        w[0] = w0;
                        n[0] = N[0] / w0;
                        continue;
                    }
                    T g{}, ww{};
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                        {
                            const T c = binomial<T>(a, i) * binomial<T>(b, j);
                            g += c * (N[i * q + j] * N[(a - i) * q + (b - j)]);
                            if (!(i == 0 && j == 0) && !(i == a && j == b))
                                ww += c * w[i * q + j] * w[(a - i) * q + (b - j)];
                        }
                    w[a * q + b] = (g - ww) / (T(2) * w0);
                    auto r = N[a * q + b];
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                            if (!(i == 0 && j == 0))
                                r = r - (binomial<T>(a, i) * binomial<T>(b, j) * w[i * q + j]) * n[(a - i) * q + (b - j)];
                    n[a * q + b] = r / w0;
                }
            return n;
        }

        /// Exact derivatives of S + d n up to the order k at (u, v) (indexed as unit_derivatives),
        /// nullopt where the normal of S vanishes (|S_u x S_v| <= 1e-12 |S_u| |S_v|).
        template <std::floating_point T>
        auto offset_derivatives(const Surface<T, 3> &s, T d, T u, T v, std::size_t k) -> std::optional<std::vector<point<T, 3>>>
        {
            const std::size_t m = k + 2, q = k + 1;
            std::vector<point<T, 3>> S(m * m);
            s.derivatives(u, v, k + 1, k + 1, S.data());
            auto at = [&](std::size_t a, std::size_t b) -> const point<T, 3> & { return S[a * m + b]; };
            std::vector<point<T, 3>> N(q * q, point<T, 3>{});
            for (std::size_t a = 0; a <= k; ++a)
                for (std::size_t b = 0; a + b <= k; ++b)
                {
                    point<T, 3> x{};
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                            x = x + (binomial<T>(a, i) * binomial<T>(b, j)) * cross(at(i + 1, j), at(a - i, b - j + 1));
                    N[a * q + b] = x;
                }
            const auto n = unit_derivatives<T, 3>(N, k, T(1e-12) * norm(at(1, 0)) * norm(at(0, 1)));
            if (!n)
                return std::nullopt;
            std::vector<point<T, 3>> P(q * q, point<T, 3>{});
            for (std::size_t a = 0; a <= k; ++a)
                for (std::size_t b = 0; a + b <= k; ++b)
                    P[a * q + b] = at(a, b) + d * (*n)[a * q + b];
            return P;
        }
    }

    /**
     * @brief The offset S + d n of a surface S, n = S_u x S_v / |S_u x S_v| (a negative d goes the
     * other way), on the parameters of S.
     *
     * The derivatives are exact, of any order k when S has those of order k + 1 (file header);
     * a basis without second derivatives gives the first ones by central differences, and none
     * beyond. Where the normal of S vanishes (a pole, a collapsed row of poles), n is its limit
     * (`normal_direction`) and the derivatives are those of a point moved toward the centre of
     * the bounds by 1e-7 of their size (then 1e-6, 1e-5): limits within that step.
     */
    template <std::floating_point T>
    class OffsetSurface : public Surface<T, 3>
    {
        std::shared_ptr<Surface<T, 3>> s_;
        T d_;

        /// The derivatives up to k at (u, v), or at a point moved toward the centre at a pole.
        auto exact(T u, T v, std::size_t k) const -> std::vector<point<T, 3>>
        {
            if (auto P = offset_detail::offset_derivatives(*s_, d_, u, v, k))
                return std::move(*P);
            const auto b = bounds();
            const T cu = (b[0] + b[1]) / T(2), cv = (b[2] + b[3]) / T(2);
            for (T h : {T(1e-7), T(1e-6), T(1e-5)})
            {
                const T lu = std::isfinite(b[1] - b[0]) ? b[1] - b[0] : T(1), lv = std::isfinite(b[3] - b[2]) ? b[3] - b[2] : T(1);
                const T mu = u + std::clamp(cu - u, -h * lu, h * lu), mv = v + std::clamp(cv - v, -h * lv, h * lv);
                if (auto P = offset_detail::offset_derivatives(*s_, d_, mu, mv, k))
                    return std::move(*P);
            }
            throw std::domain_error("OffsetSurface: no normal near the point");
        }

        auto central(T u, T v, std::size_t du) const -> point<T, 3>
        {
            const auto b = bounds();
            const std::size_t i = du == 1 ? 0 : 2;
            const T h = T(1e-6) * (b[i + 1] - b[i]);
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
            const std::size_t k = du + dv;
            return exact(u, v, k)[du * (k + 1) + dv];
        }
        auto derivatives(T u, T v, size_t nu, size_t nv, std::array<T, 3> *SKL) const -> void override
        {
            if (!s_->has_second_derivatives() || nu + nv == 0)
            {
                Surface<T, 3>::derivatives(u, v, nu, nv, SKL);
                return;
            }
            const std::size_t k = nu + nv;
            const auto P = exact(u, v, k);
            for (std::size_t a = 0; a <= nu; ++a)
                for (std::size_t b = 0; b <= nv; ++b)
                    SKL[a * (nv + 1) + b] = a + b == 0 ? value(u, v) : P[a * (k + 1) + b];
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
        /// 1 - d k for the principal curvature k that brings it lowest at (u, v); nullopt where the
        /// normal of s vanishes.
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
            const auto N = cross(su, sv);
            const T w = norm(N);
            if (!(w > T(1e-12) * norm(su) * norm(sv)))
                return std::nullopt;
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

        template <std::floating_point T>
        auto knots(const Surface<T, 3> &s) -> std::array<std::optional<std::vector<T>>, 2>
        {
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(&s))
                return {b->knotsFlatsU(), b->knotsFlatsV()};
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                return {b->knotsFlatsU(), b->knotsFlatsV()};
            return {};
        }
    }

    /**
     * @brief Where the offset of s at the distance d comes nearest to fold on `box`: the smallest
     * of 1 - d k over the principal curvatures k (file header), regular while it is positive.
     *
     * Sampled (33 even values per direction, and 4 per knot span of a B-spline), then refined
     * around the 4 lowest samples by a pattern search down to 1e-6 of a sample step. The points
     * where the normal of s vanishes are skipped. A fold narrower than a sample step around a
     * sample that is not among the lowest may be missed.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_margin(const Surface<T, 3> &s, T d, const std::array<T, 4> &box) -> OffsetMargin<T>
    {
        using namespace offset_detail;
        const auto k = knots(s);
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

    /**
     * @brief The exact offset of a primitive: plane moved, cylinder, sphere and cone of another
     * radius, torus of another minor radius, with the shift of v of the cone; checked against
     * S + d n at 4 points within tol. nullopt beyond the radius of curvature (no primitive of the
     * same kind fits: the offset folds or reduces to a point or a line).
     */
    template <std::floating_point T>
    [[nodiscard]] auto primitive_offset(const ElementarySurface<T> &s, T d, T tol) -> std::optional<ExactOffset<T>>
    {
        const auto &f = s.frame();
        const auto box = s.bounds();
        std::vector<std::pair<std::shared_ptr<ElementarySurface<T>>, T>> cands;
        if (dynamic_cast<const Plane<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                auto g = f;
                g.o = f.o + sg * d * f.z;
                cands.push_back({std::make_shared<Plane<T>>(g, box), T(0)});
            }
        else if (const auto *c = dynamic_cast<const CylindricalSurface<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                if (c->radius() + sg * d > T(0))
                    cands.push_back({std::make_shared<CylindricalSurface<T>>(c->radius() + sg * d, f, box), T(0)});
            }
        else if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                if (sp->radius() + sg * d > T(0))
                    cands.push_back({std::make_shared<SphericalSurface<T>>(sp->radius() + sg * d, f, box), T(0)});
            }
        else if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(&s))
            for (T sg : {T(1), T(-1)})
            {
                if (to->minor_radius() + sg * d > T(0))
                    cands.push_back({std::make_shared<ToroidalSurface<T>>(to->major_radius(), to->minor_radius() + sg * d, f, box), T(0)});
            }
        else if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(&s))
        {
            const T a = co->semi_angle();
            for (T sg : {T(1), T(-1)})
            {
                const T dv = -sg * d * std::sin(a), R = co->radius() + sg * d / std::cos(a);
                if (R >= T(0))
                    cands.push_back({std::make_shared<ConicalSurface<T>>(R, a, f, std::array<T, 4>{box[0], box[1], box[2] + dv, box[3] + dv}), dv});
            }
        }
        // the candidate where S'(u, v + dv) = S(u, v) + d n
        for (const auto &[c, dv] : cands)
        {
            bool ok = true;
            for (T fu : {T(0.3), T(0.7)})
                for (T fv : {T(0.4), T(0.6)})
                {
                    const T u = box[0] + fu * (box[1] - box[0]), v = box[2] + fv * (box[3] - box[2]);
                    const auto n = cross(s.value(u, v, 1, 0), s.value(u, v, 0, 1));
                    const auto target = s.value(u, v) + (d / norm(n)) * n;
                    ok = ok && norm(c->value(u, v + dv) - target) <= tol;
                }
            if (ok)
                return ExactOffset<T>{c, dv};
        }
        return std::nullopt;
    }

    /**
     * @brief The offset of s at the distance d on `box` (its bounds by default), exact when it
     * can be: a primitive gives a primitive (`primitive_offset`), an offset surface the offset of
     * its basis at the sum of the distances (same normal while it does not fold), any other
     * surface an `OffsetSurface`. nullopt where the offset folds on the box: a primitive beyond
     * its radius of curvature, another surface where `offset_margin` is at most `fold`.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_surface(const std::shared_ptr<Surface<std::type_identity_t<T>, 3>> &s, T d, std::type_identity_t<T> tol,
                                      std::optional<std::array<std::type_identity_t<T>, 4>> box = std::nullopt,
                                      std::type_identity_t<T> fold = T(1e-9)) -> std::optional<ExactOffset<T>>
    {
        if (!s)
            return std::nullopt;
        if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(s.get()))
            return primitive_offset(*e, d, tol);
        if (!(offset_margin(*s, d, box.value_or(s->bounds())).margin > fold))
            return std::nullopt;
        if (const auto *o = dynamic_cast<const OffsetSurface<T> *>(s.get()))
            return ExactOffset<T>{std::make_shared<OffsetSurface<T>>(o->basis(), o->distance() + d), T(0)};
        return ExactOffset<T>{std::make_shared<OffsetSurface<T>>(s, d), T(0)};
    }
}
