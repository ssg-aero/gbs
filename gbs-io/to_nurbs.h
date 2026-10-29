#pragma once

/**
 * @file to_nurbs.h
 * @brief Conversion of any gbs curve or surface to a NURBS for the exporters
 * (IGES, STEP), with the parameter preserved so that pcurves stay valid.
 *
 * A B-spline is copied exactly (restricted to the requested range for a
 * curve); any other curve or surface is interpolated at its own parameters by
 * a cubic B-spline, refined where it deviates: a span (of a curve) or a row or
 * column of the grid (of a surface) is halved while the deviation measured at
 * its middle and quarters exceeds tol. An approximation that cannot reach tol
 * within the limits (16384 spans for a curve, 512 per direction for a surface)
 * is counted in `ApproxStats::over_tolerance`, its deviation in `max_deviation`.
 * Shared by gbs-io/iges_brep.h and the STEP writer (design
 * docs/sources/design/step_writer.md § 3.3).
 */

#include <algorithm>
#include <cmath>
#include <variant>
#include <vector>

#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/curveonsurface.h>
#include <gbs/primitives.h>

namespace gbs::nurbs_export
{
    /// How much an export had to approximate.
    template <std::floating_point T>
    struct ApproxStats
    {
        std::size_t approximated_curves{}, approximated_surfaces{};
        std::size_t over_tolerance{}; ///< approximations whose deviation stayed above the tolerance
        T max_deviation{};            ///< largest deviation of an approximated curve or surface
    };

    template <std::floating_point T, std::size_t dim>
    using AnyBSCurve = std::variant<BSCurve<T, dim>, BSCurveRational<T, dim>>;
    template <std::floating_point T>
    using AnyBSSurface = std::variant<BSSurface<T, 3>, BSSurfaceRational<T, 3>>;

    namespace detail
    {
        /// Largest gap over [t1, t2]: at its middle and quarters, and at the top of the parabola
        /// through them when it is concave (the error of an interpolation peaks off the middle).
        template <std::floating_point T, typename Gap>
        T span_gap(const Gap &gap, T t1, T t2)
        {
            const T h = (t2 - t1) / 4, tm = (t1 + t2) / 2;
            const T d1 = gap(tm - h), d2 = gap(tm), d3 = gap(tm + h), c = d1 - 2 * d2 + d3;
            T d = std::max({d1, d2, d3});
            if (c < T(0))
                d = std::max(d, gap(std::clamp(tm + h * (d1 - d3) / (2 * c), t1, t2)));
            return d;
        }

        /// Halves the spans of x marked in `bad` (wider than `min_width`), up to `max_spans` spans; false if none.
        template <std::floating_point T>
        bool refine(std::vector<T> &x, const std::vector<char> &bad, std::size_t max_spans, T min_width)
        {
            std::vector<T> y{x.front()};
            std::size_t spans = x.size() - 1;
            bool any = false;
            for (std::size_t i = 0; i + 1 < x.size(); ++i)
            {
                if (bad[i] && spans < max_spans && x[i + 1] - x[i] > min_width)
                {
                    y.push_back((x[i] + x[i + 1]) / 2);
                    ++spans, any = true;
                }
                y.push_back(x[i + 1]);
            }
            x = std::move(y);
            return any;
        }

        template <std::floating_point T>
        std::vector<T> uniform(T a, T b, std::size_t n)
        {
            std::vector<T> x(n);
            for (std::size_t i = 0; i < n; ++i)
                x[i] = a + (b - a) * T(i) / T(n - 1);
            x.back() = b;
            return x;
        }

        /// Bicubic interpolation of the grid q (u fastest) at the parameters u, v: two banded collocation solves.
        template <std::floating_point T>
        auto interpolate_grid(const points_vector<T, 3> &q, const std::vector<T> &u, const std::vector<T> &v, std::size_t p) -> BSSurface<T, 3>
        {
            const auto nu = Eigen::Index(u.size()), nv = Eigen::Index(v.size());
            auto ku = build_simple_mult_flat_knots<T>(u, p), kv = build_simple_mult_flat_knots<T>(v, p);
            collocation_solver<T> su, sv;
            su.prepare(ku, u, p, 1, nu);
            sv.prepare(kv, v, p, 1, nv);
            MatrixX<T> b(nu, nv * 3);
            for (Eigen::Index j = 0; j < nv; ++j)
                for (Eigen::Index i = 0; i < nu; ++i)
                    for (Eigen::Index d = 0; d < 3; ++d)
                        b(i, 3 * j + d) = q[std::size_t(i + nu * j)][std::size_t(d)];
            const MatrixX<T> a = su.solve(b); // Mu^-1 Q
            MatrixX<T> c(nv, nu * 3);
            for (Eigen::Index j = 0; j < nv; ++j)
                for (Eigen::Index i = 0; i < nu; ++i)
                    for (Eigen::Index d = 0; d < 3; ++d)
                        c(j, 3 * i + d) = a(i, 3 * j + d);
            const MatrixX<T> x = sv.solve(c); // (Mu^-1 Q Mv^-T)^T
            points_vector<T, 3> poles(q.size());
            for (Eigen::Index j = 0; j < nv; ++j)
                for (Eigen::Index i = 0; i < nu; ++i)
                    for (Eigen::Index d = 0; d < 3; ++d)
                        poles[std::size_t(i + nu * j)][std::size_t(d)] = x(j, 3 * i + d);
            return {poles, ku, kv, p, p};
        }
    }

    /**
     * Cubic B-spline through f(t) at parameters of [a, b], parameter preserved: a span is halved
     * while the deviation gap(bs, s) at its middle and quarters exceeds tol, up to `max_spans`
     * spans. Counted in `rep` as an approximated curve.
     */
    template <std::floating_point T, std::size_t dim, typename F, typename Gap>
    auto interpolate_within(const F &f, T a, T b, const Gap &gap, T tol, ApproxStats<T> &rep, std::size_t max_spans = 16384) -> BSCurve<T, dim>
    {
        ++rep.approximated_curves;
        auto t = detail::uniform(a, b, 9);
        for (;;)
        {
            std::vector<constrType<T, dim, 1>> q(t.size());
            for (std::size_t i{}; i < t.size(); ++i)
                q[i] = {f(t[i])};
            auto bs = interpolate<T, dim>(q, t, std::size_t{3});
            std::vector<char> bad(t.size() - 1);
            T dev{};
            for (std::size_t i{}; i + 1 < t.size(); ++i)
            {
                const T d = detail::span_gap([&](T s) { return gap(bs, s); }, t[i], t[i + 1]);
                dev = std::max(dev, d);
                bad[i] = d > tol;
            }
            if (dev <= tol || !detail::refine(t, bad, max_spans, T(1e-9) * (b - a)))
            {
                if (dev > tol)
                    ++rep.over_tolerance;
                rep.max_deviation = std::max(rep.max_deviation, dev);
                return bs;
            }
        }
    }

    /**
     * The iso line that a curve on a surface is, exactly, when its pcurve is an iso segment parametrized
     * like it (the moving parameter is the curve parameter): an ElementaryCurve on [a, b] for a primitive,
     * the B-spline iso of a B-spline surface; nullptr otherwise.
     */
    template <std::floating_point T>
    auto iso_curve(const CurveOnSurface<T, 3> &cos, T a, T b) -> std::shared_ptr<Curve<T, 3>>
    {
        const auto pc = std::dynamic_pointer_cast<BSCurve<T, 2>>(cos.p_basisCurve());
        if (!pc || pc->degree() != 1 || pc->poles().size() != 2)
            return nullptr;
        const auto p0 = pc->poles()[0], p1 = pc->poles()[1];
        const auto [k1, k2] = pc->bounds();
        for (std::size_t k = 0; k < 2; ++k)
        {
            const std::size_t o = 1 - k; // the coordinate that moves
            if (std::abs(p1[k] - p0[k]) > T(1e-14) * (1 + std::abs(p0[k])))
                continue;
            // the moving coordinate must equal the curve parameter
            if (std::abs((p1[o] - p0[o]) - (k2 - k1)) > T(1e-12) * (1 + std::abs(k2 - k1)) || std::abs(p0[o] - k1) > T(1e-12) * (1 + std::abs(k1)))
                return nullptr;
            const auto srf = cos.p_basisSurface();
            if (auto e = std::dynamic_pointer_cast<ElementarySurface<T>>(srf))
                return k == 0 ? e->iso_u(p0[0], a, b) : e->iso_v(p0[1], a, b);
            if (auto r = std::dynamic_pointer_cast<BSSurfaceRational<T, 3>>(srf))
                return std::make_shared<BSCurveRational<T, 3>>(k == 0 ? r->isoU(p0[0]) : r->isoV(p0[1]));
            if (auto q = std::dynamic_pointer_cast<BSSurface<T, 3>>(srf))
                return std::make_shared<BSCurve<T, 3>>(k == 0 ? q->isoU(p0[0]) : q->isoV(p0[1]));
        }
        return nullptr;
    }

    /**
     * NURBS of `c` on [a, b] with its parameter preserved: an exact copy
     * (trimmed) for a B-spline, otherwise a cubic interpolation refined until within tol.
     */
    template <std::floating_point T, std::size_t dim>
    auto to_nurbs(const Curve<T, dim> &c, T a, T b, T tol, ApproxStats<T> &rep) -> AnyBSCurve<T, dim>
    {
        auto trimmed = [&](auto copy) -> AnyBSCurve<T, dim> {
            const auto [c1, c2] = copy.bounds();
            if (a > c1 + knot_eps<T> || b < c2 - knot_eps<T>)
                copy.trim(a, b, true);
            return copy;
        };
        if (auto bs = dynamic_cast<const BSCurve<T, dim> *>(&c))
            return trimmed(*bs);
        if (auto bsr = dynamic_cast<const BSCurveRational<T, dim> *>(&c))
            return trimmed(*bsr);
        return interpolate_within<T, dim>([&](T t) { return c.value(t); }, a, b,
                                          [&](const BSCurve<T, dim> &bs, T s) { return distance(bs.value(s), c.value(s)); }, tol, rep);
    }

    /**
     * NURBS of a surface on its whole parametric rectangle, parameters preserved: an exact copy
     * for a B-spline, otherwise a bicubic interpolation on a grid whose rows and columns are
     * halved where the deviation at the middles and quarters of a cell exceeds tol.
     */
    template <std::floating_point T>
    auto to_nurbs(const Surface<T, 3> &s, T tol, ApproxStats<T> &rep) -> AnyBSSurface<T>
    {
        if (auto bs = dynamic_cast<const BSSurface<T, 3> *>(&s))
            return *bs;
        if (auto bsr = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
            return *bsr;

        ++rep.approximated_surfaces;
        const auto [u1, u2, v1, v2] = s.bounds();
        auto u = detail::uniform(u1, u2, 9), v = detail::uniform(v1, v2, 9);
        constexpr std::size_t max_spans = 512;
        for (;;)
        {
            points_vector<T, 3> q(u.size() * v.size());
            for (std::size_t j{}; j < v.size(); ++j)
                for (std::size_t i{}; i < u.size(); ++i)
                    q[i + u.size() * j] = s.value(u[i], v[j]);
            auto bs = detail::interpolate_grid(q, u, v, 3);
            std::vector<char> bad_u(u.size() - 1), bad_v(v.size() - 1);
            T dev{};
            for (std::size_t j{}; j + 1 < v.size(); ++j)
                for (std::size_t i{}; i + 1 < u.size(); ++i)
                {
                    T d{};
                    for (int k = 1; k < 4; ++k)
                        for (int l = 1; l < 4; ++l)
                        {
                            const T uu = u[i] + (u[i + 1] - u[i]) * k / 4, vv = v[j] + (v[j + 1] - v[j]) * l / 4;
                            d = std::max(d, distance(bs.value(uu, vv), s.value(uu, vv)));
                        }
                    dev = std::max(dev, d);
                    if (d > tol)
                        bad_u[i] = bad_v[j] = 1;
                }
            bool refined = false;
            if (dev > tol)
            {
                const bool ru = detail::refine(u, bad_u, max_spans, T(1e-9) * (u2 - u1));
                const bool rv = detail::refine(v, bad_v, max_spans, T(1e-9) * (v2 - v1));
                refined = ru || rv;
            }
            if (!refined)
            {
                if (dev > tol)
                    ++rep.over_tolerance;
                rep.max_deviation = std::max(rep.max_deviation, dev);
                return bs;
            }
        }
    }

}
