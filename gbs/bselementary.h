#pragma once
/**
 * @file bselementary.h
 * @brief Exact rational NURBS of conic arcs and elementary surfaces.
 *
 * Arcs of circle and ellipse between two angles, rational revolution of a NURBS around an axis,
 * linear extrusion of a NURBS, and the elementary surfaces built from them (cylinder, cone,
 * sphere, torus), with the parametrization conventions of ISO 10303-42 (STEP).
 *
 * An arc of angle Δ is split into n = ceil(Δ / 90°) rational quadratic segments of equal angle
 * (The NURBS Book, § 7.5). Its knots are the segment boundaries **in radians**: the curve
 * parameter equals the angle at every knot, and in between ellipse_arc_parameter() and
 * ellipse_arc_angle() convert one into the other exactly.
 */
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include "vecop.ixx"
#endif

namespace gbs
{
    namespace elementary_detail
    {
        // Pole of a unit arc: point = O + c X + s Y, weight w.
        template <std::floating_point T>
        struct UnitArcPole
        {
            T c, s, w;
        };

        template <std::floating_point T>
        auto arc_segments(T theta1, T theta2) -> size_t
        {
            const T delta = theta2 - theta1;
            if (!(delta > T(0)) || !std::isfinite(delta))
                throw std::invalid_argument("arc: theta2 must be greater than theta1");
            if (delta > T(2) * std::numbers::pi_v<T> * (T(1) + T(1e-12)))
                throw std::invalid_argument("arc: angle span larger than 2 pi");
            const auto n = static_cast<size_t>(std::ceil(delta / (T(0.5) * std::numbers::pi_v<T>) - T(1e-9)));
            return std::max<size_t>(n, 1);
        }

        // Poles (2n + 1) and flat knots (2n + 4) of the unit arc from theta1 to theta2.
        template <std::floating_point T>
        auto unit_arc(T theta1, T theta2) -> std::pair<std::vector<UnitArcPole<T>>, std::vector<T>>
        {
            const auto n = arc_segments(theta1, theta2);
            const T d = (theta2 - theta1) / T(n);
            const T w = std::cos(T(0.5) * d);
            std::vector<UnitArcPole<T>> poles;
            poles.reserve(2 * n + 1);
            std::vector<T> knots{theta1, theta1, theta1};
            poles.push_back({std::cos(theta1), std::sin(theta1), T(1)});
            for (size_t k = 0; k < n; ++k)
            {
                const T a = theta1 + T(k) * d;
                const T b = k + 1 == n ? theta2 : a + d;
                const T m = T(0.5) * (a + b);
                // middle pole at the intersection of the end tangents, distance 1 / cos(d/2) from the center
                poles.push_back({std::cos(m) / w, std::sin(m) / w, w});
                poles.push_back({std::cos(b), std::sin(b), T(1)});
                if (k + 1 < n)
                    knots.insert(knots.end(), {b, b});
            }
            knots.insert(knots.end(), {theta2, theta2, theta2});
            return {std::move(poles), std::move(knots)};
        }

        template <typename T, size_t dim>
        auto homogeneous(const point<T, dim> &p, T w) -> std::array<T, dim + 1>
        {
            std::array<T, dim + 1> h;
            for (size_t i = 0; i < dim; ++i)
                h[i] = p[i] * w;
            h[dim] = w;
            return h;
        }

        // Unit vector of a 3D frame direction; throws on a null vector.
        template <typename T>
        auto unit(const point<T, 3> &v, const char *what) -> point<T, 3>
        {
            const T n = norm(v);
            if (!(n > T(0)) || !std::isfinite(n))
                throw std::invalid_argument(what);
            return v / n;
        }

        // Orthonormal (X, Y, Z) of an ax2 {origin, axis, ref_direction}, as STEP axis2_placement_3d:
        // X is the ref direction made orthogonal to the axis, Y = Z x X.
        template <typename T>
        auto frame(const ax2<T, 3> &ax) -> std::array<point<T, 3>, 3>
        {
            const auto Z = unit(ax[1], "frame: null axis");
            const auto X = unit(ax[2] - (ax[2] * Z) * Z, "frame: reference direction parallel to the axis");
            return {X, cross(Z, X), Z};
        }

        template <typename T>
        auto check_positive(T r, const char *what) -> void
        {
            if (!(r > T(0)) || !std::isfinite(r))
                throw std::invalid_argument(what);
        }
    }

    /**
     * @brief Parameter of the arc point at angle theta, for an arc built from theta1 to theta2
     * by build_ellipse_arc() or build_circle_arc(). Exact; equals theta at the knots.
     */
    template <std::floating_point T>
    auto ellipse_arc_parameter(T theta1, T theta2, T theta) -> T
    {
        const auto n = elementary_detail::arc_segments(theta1, theta2);
        const T d = (theta2 - theta1) / T(n);
        const auto k = std::min<size_t>(static_cast<size_t>(std::max(T(0), std::floor((theta - theta1) / d))), n - 1);
        const T a = theta1 + T(k) * d;
        // rational quadratic arc of angle d with weights (1, cos(d/2), 1): tan((theta - mid)/2) = s tan(d/4), s in [-1, 1]
        const T s = std::tan(T(0.5) * (theta - a - T(0.5) * d)) / std::tan(T(0.25) * d);
        return a + T(0.5) * (s + T(1)) * d;
    }

    /**
     * @brief Angle of the arc point at parameter u, inverse of ellipse_arc_parameter().
     */
    template <std::floating_point T>
    auto ellipse_arc_angle(T theta1, T theta2, T u) -> T
    {
        const auto n = elementary_detail::arc_segments(theta1, theta2);
        const T d = (theta2 - theta1) / T(n);
        const auto k = std::min<size_t>(static_cast<size_t>(std::max(T(0), std::floor((u - theta1) / d))), n - 1);
        const T a = theta1 + T(k) * d;
        const T s = T(2) * (u - a) / d - T(1);
        return a + T(0.5) * d + T(2) * std::atan(s * std::tan(T(0.25) * d));
    }

    /**
     * @brief Exact arc of ellipse C(θ) = O + r1 cos θ X + r2 sin θ Y, θ from theta1 to theta2,
     * in any dimension. X and Y are not normalized: any affine image of a circle is exact.
     *
     * The parameter range is [theta1, theta2]; the parameter is the angle at the knots
     * (see ellipse_arc_parameter()). Degree 2, 2n + 1 poles for n = ceil(Δθ / 90°).
     */
    template <typename T, size_t dim>
    auto build_ellipse_arc(T r1, T r2, T theta1, T theta2,
                           const point<T, dim> &O, const point<T, dim> &X, const point<T, dim> &Y) -> BSCurveRational<T, dim>
    {
        auto [arc, knots] = elementary_detail::unit_arc(theta1, theta2);
        std::vector<std::array<T, dim + 1>> poles;
        poles.reserve(arc.size());
        for (const auto &a : arc)
            poles.push_back(elementary_detail::homogeneous<T, dim>(O + (r1 * a.c) * X + (r2 * a.s) * Y, a.w));
        return BSCurveRational<T, dim>{poles, knots, 2};
    }

    /**
     * @brief Exact arc of ellipse in the frame ax = {center, axis, major direction} (STEP ellipse):
     * C(θ) = O + r1 cos θ X + r2 sin θ Y, Y = Z × X.
     */
    template <typename T>
    auto build_ellipse_arc(T r1, T r2, T theta1, T theta2, const ax2<T, 3> &ax) -> BSCurveRational<T, 3>
    {
        elementary_detail::check_positive(r1, "ellipse: radius must be positive");
        elementary_detail::check_positive(r2, "ellipse: radius must be positive");
        const auto [X, Y, Z] = elementary_detail::frame(ax);
        return build_ellipse_arc<T, 3>(r1, r2, theta1, theta2, ax[0], X, Y);
    }

    /**
     * @brief Exact arc of circle in the frame ax = {center, axis, start direction} (STEP circle).
     */
    template <typename T>
    auto build_circle_arc(T r, T theta1, T theta2, const ax2<T, 3> &ax) -> BSCurveRational<T, 3>
    {
        return build_ellipse_arc<T>(r, r, theta1, theta2, ax);
    }

    /**
     * @brief Exact 2D arc of circle of center O, angles measured from the x axis.
     */
    template <typename T>
    auto build_circle_arc(T r, T theta1, T theta2, const point<T, 2> &O = {}) -> BSCurveRational<T, 2>
    {
        elementary_detail::check_positive(r, "circle: radius must be positive");
        return build_ellipse_arc<T, 2>(r, r, theta1, theta2, O, {T(1), T(0)}, {T(0), T(1)});
    }

    /**
     * @brief Exact rational surface of revolution of a NURBS generatrix around an axis (STEP
     * surface_of_revolution): S(u, v) = generatrix(v) rotated by the angle u around the axis,
     * right-handed, u from theta1 to theta2.
     *
     * u is the rotation (parameter equal to the angle at the knots, as for the arcs) and v the
     * generatrix parameter, unchanged. A generatrix pole on the axis gives a degenerate row (pole).
     */
    template <typename T, bool rational>
    auto build_revolution(const BSCurveGeneral<T, 3, rational> &generatrix, const ax1<T, 3> &axis,
                          T theta1 = T(0), T theta2 = T(2) * std::numbers::pi_v<T>) -> BSSurfaceRational<T, 3>
    {
        const auto Z = elementary_detail::unit(axis[1], "revolution: null axis");
        auto [arc, knots_u] = elementary_detail::unit_arc(theta1, theta2);
        const auto &gp = generatrix.poles();
        std::vector<std::array<T, 4>> poles;
        poles.reserve(arc.size() * gp.size());
        for (const auto &g : gp) // v rows, u (rotation) fastest
        {
            point<T, 3> P;
            T wg{1};
            if constexpr (rational)
            {
                wg = g[3];
                P = point<T, 3>{g[0], g[1], g[2]} / wg;
            }
            else
                P = g;
            const auto O = axis[0] + ((P - axis[0]) * Z) * Z;
            const auto X = P - O;
            const auto Y = cross(Z, X);
            for (const auto &a : arc)
                poles.push_back(elementary_detail::homogeneous<T, 3>(O + a.c * X + a.s * Y, a.w * wg));
        }
        return BSSurfaceRational<T, 3>{poles, knots_u, generatrix.knotsFlats(), 2, generatrix.degree()};
    }

    /**
     * @brief Exact linear extrusion of a NURBS (STEP surface_of_linear_extrusion):
     * S(u, v) = C(u) + v V, v from v1 to v2, degree 1 in v. Rational if the curve is.
     */
    template <typename T, size_t dim, bool rational>
    auto build_extrusion(const BSCurveGeneral<T, dim, rational> &crv, const point<T, dim> &V, T v1 = T(0), T v2 = T(1))
    {
        if (!(v2 > v1))
            throw std::invalid_argument("extrusion: v2 must be greater than v1");
        const auto &cp = crv.poles();
        std::vector<std::array<T, dim + rational>> poles;
        poles.reserve(2 * cp.size());
        for (T v : {v1, v2})
            for (auto p : cp)
            {
                T w{1};
                if constexpr (rational)
                    w = p[dim];
                for (size_t i = 0; i < dim; ++i)
                    p[i] += w * v * V[i]; // homogeneous coordinates when rational
                poles.push_back(p);
            }
        using Surf = std::conditional_t<rational, BSSurfaceRational<T, dim>, BSSurface<T, dim>>;
        return Surf{poles, crv.knotsFlats(), std::vector<T>{v1, v1, v2, v2}, crv.degree(), 1};
    }

    /**
     * @brief Exact cylinder (STEP cylindrical_surface):
     * S(u, v) = O + R (cos u X + sin u Y) + v Z, u from theta1 to theta2, v from v1 to v2.
     */
    template <typename T>
    auto build_cylinder(T R, const ax2<T, 3> &ax, T v1, T v2,
                        T theta1 = T(0), T theta2 = T(2) * std::numbers::pi_v<T>) -> BSSurfaceRational<T, 3>
    {
        const auto [X, Y, Z] = elementary_detail::frame(ax);
        return build_extrusion(build_circle_arc<T>(R, theta1, theta2, ax), Z, v1, v2);
    }

    /**
     * @brief Exact cone (STEP conical_surface, semi_angle in radians):
     * S(u, v) = O + (R + v tan α)(cos u X + sin u Y) + v Z, v from v1 to v2.
     * A radius of zero at v1 or v2 (the apex) gives a degenerate side.
     */
    template <typename T>
    auto build_cone(T R, T semi_angle, const ax2<T, 3> &ax, T v1, T v2,
                    T theta1 = T(0), T theta2 = T(2) * std::numbers::pi_v<T>) -> BSSurfaceRational<T, 3>
    {
        if (!(v2 > v1))
            throw std::invalid_argument("cone: v2 must be greater than v1");
        if (!(std::abs(semi_angle) < T(0.5) * std::numbers::pi_v<T>))
            throw std::invalid_argument("cone: semi angle must be in ]-pi/2, pi/2[");
        const auto [X, Y, Z] = elementary_detail::frame(ax);
        const T t = std::tan(semi_angle);
        auto line = [&](T v) { return ax[0] + (R + v * t) * X + v * Z; };
        const BSCurve<T, 3> generatrix{points_vector<T, 3>{line(v1), line(v2)}, std::vector<T>{v1, v1, v2, v2}, 1};
        return build_revolution(generatrix, ax1<T, 3>{ax[0], Z}, theta1, theta2);
    }

    /**
     * @brief Exact sphere (STEP spherical_surface):
     * S(u, v) = O + R cos v (cos u X + sin u Y) + R sin v Z, v (latitude) from v1 to v2 in
     * [-pi/2, pi/2]; the parameter is the latitude at the knots in v, as for the arcs.
     */
    template <typename T>
    auto build_sphere(T R, const ax2<T, 3> &ax,
                      T theta1 = T(0), T theta2 = T(2) * std::numbers::pi_v<T>,
                      T v1 = -T(0.5) * std::numbers::pi_v<T>, T v2 = T(0.5) * std::numbers::pi_v<T>) -> BSSurfaceRational<T, 3>
    {
        elementary_detail::check_positive(R, "sphere: radius must be positive");
        if (v1 < -T(0.5) * std::numbers::pi_v<T> * (T(1) + T(1e-12)) || v2 > T(0.5) * std::numbers::pi_v<T> * (T(1) + T(1e-12)))
            throw std::invalid_argument("sphere: latitude out of [-pi/2, pi/2]");
        const auto [X, Y, Z] = elementary_detail::frame(ax);
        const auto meridian = build_ellipse_arc<T, 3>(R, R, v1, v2, ax[0], X, Z);
        return build_revolution(meridian, ax1<T, 3>{ax[0], Z}, theta1, theta2);
    }

    /**
     * @brief Exact torus (STEP toroidal_surface, major radius R, minor radius r):
     * S(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z.
     */
    template <typename T>
    auto build_torus(T R, T r, const ax2<T, 3> &ax,
                     T theta1 = T(0), T theta2 = T(2) * std::numbers::pi_v<T>,
                     T v1 = T(0), T v2 = T(2) * std::numbers::pi_v<T>) -> BSSurfaceRational<T, 3>
    {
        elementary_detail::check_positive(r, "torus: minor radius must be positive");
        const auto [X, Y, Z] = elementary_detail::frame(ax);
        const auto section = build_ellipse_arc<T, 3>(r, r, v1, v2, ax[0] + R * X, X, Z);
        return build_revolution(section, ax1<T, 3>{ax[0], Z}, theta1, theta2);
    }
}
