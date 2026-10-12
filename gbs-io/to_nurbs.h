#pragma once

/**
 * @file to_nurbs.h
 * @brief Conversion of any gbs curve or surface to a NURBS for the exporters
 * (IGES, STEP), with the parameter preserved so that pcurves stay valid.
 *
 * A B-spline is copied exactly (restricted to the requested range for a
 * curve); any other curve or surface is interpolated at its own parameters by
 * a cubic B-spline, refined until the deviation at the midpoints is <= tol.
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

namespace gbs::nurbs_export
{
    /// How much an export had to approximate.
    template <std::floating_point T>
    struct ApproxStats
    {
        std::size_t approximated_curves{}, approximated_surfaces{};
        T max_deviation{}; ///< largest deviation of an approximated curve or surface
    };

    template <std::floating_point T, std::size_t dim>
    using AnyBSCurve = std::variant<BSCurve<T, dim>, BSCurveRational<T, dim>>;
    template <std::floating_point T>
    using AnyBSSurface = std::variant<BSSurface<T, 3>, BSSurfaceRational<T, 3>>;

    /**
     * NURBS of `c` on [a, b] with its parameter preserved: an exact copy
     * (trimmed) for a B-spline, otherwise a cubic interpolation at uniform
     * parameters, refined until the deviation at the midpoints is <= tol.
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

        ++rep.approximated_curves;
        std::size_t n = 9;
        for (;;)
        {
            std::vector<T> t(n);
            std::vector<constrType<T, dim, 1>> q(n);
            for (std::size_t i{}; i < n; ++i)
            {
                t[i] = a + (b - a) * T(i) / T(n - 1);
                q[i] = {c.value(t[i])};
            }
            auto bs = interpolate<T, dim>(q, t, std::min<std::size_t>(3, n - 1));
            T dev{};
            for (std::size_t i{}; i + 1 < n; ++i)
            {
                const T tm = (t[i] + t[i + 1]) / 2;
                dev = std::max(dev, distance(bs.value(tm), c.value(tm)));
            }
            if (dev <= tol || n >= 1025)
            {
                rep.max_deviation = std::max(rep.max_deviation, dev);
                return bs;
            }
            n = 2 * n - 1;
        }
    }

    /// NURBS of a surface on its whole parametric rectangle, parameters preserved.
    template <std::floating_point T>
    auto to_nurbs(const Surface<T, 3> &s, T tol, ApproxStats<T> &rep) -> AnyBSSurface<T>
    {
        if (auto bs = dynamic_cast<const BSSurface<T, 3> *>(&s))
            return *bs;
        if (auto bsr = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
            return *bsr;

        ++rep.approximated_surfaces;
        const auto [u1, u2, v1, v2] = s.bounds();
        std::size_t n = 9;
        for (;;)
        {
            std::vector<T> u(n), v(n);
            for (std::size_t i{}; i < n; ++i)
            {
                u[i] = u1 + (u2 - u1) * T(i) / T(n - 1);
                v[i] = v1 + (v2 - v1) * T(i) / T(n - 1);
            }
            points_vector<T, 3> q(n * n);
            for (std::size_t j{}; j < n; ++j)
                for (std::size_t i{}; i < n; ++i)
                    q[i + n * j] = s.value(u[i], v[j]);
            const std::size_t p = std::min<std::size_t>(3, n - 1);
            auto ku = build_simple_mult_flat_knots<T>(u, p);
            auto kv = build_simple_mult_flat_knots<T>(v, p);
            BSSurface<T, 3> bs{build_poles(q, ku, kv, u, v, p, p), ku, kv, p, p};
            T dev{};
            for (std::size_t j{}; j + 1 < n; ++j)
                for (std::size_t i{}; i + 1 < n; ++i)
                {
                    const T um = (u[i] + u[i + 1]) / 2, vm = (v[j] + v[j + 1]) / 2;
                    dev = std::max(dev, distance(bs.value(um, vm), s.value(um, vm)));
                }
            if (dev <= tol || n >= 129)
            {
                rep.max_deviation = std::max(rep.max_deviation, dev);
                return bs;
            }
            n = 2 * n - 1;
        }
    }

}
