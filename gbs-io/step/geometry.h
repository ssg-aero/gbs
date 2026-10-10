#pragma once

/**
 * @file geometry.h
 * @brief Translation of STEP geometry (ISO 10303-42) to gbs geometry.
 *
 * Design: docs/sources/design/step_reader.md, sections 3.2, 3.3 and 4;
 * architecture note docs/sources/design/step_pr04_geometry.md.
 *
 * Two stages:
 *  1. GeometryReader reads an instance into a *definition* (CurveDef, SurfaceDef)
 *     that keeps the STEP parametrization, already converted to the target units
 *     (lengths) and to radians (angles): lines and planes are infinite, circles
 *     are parametrized by the angle. A definition evaluates points, inverts
 *     points to parameters (analytically for lines, conics and elementary
 *     surfaces) and knows the kind of each parameter direction.
 *  2. to_nurbs() converts a definition restricted to a parameter range into an
 *     exact NURBS (gbs/bselementary.h for conics and elementary surfaces). The
 *     range of an edge comes from its vertices, the box of a face from its
 *     boundary (parameter_box()).
 */

#include <gbs-io/step/p21.h>
#include <gbs-io/step/units.h>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
#include <gbs/bselementary.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <optional>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <unordered_map>
#include <variant>
#include <vector>

namespace gbs::step
{
    using Real = double;
    template <std::size_t dim>
    using Point = point<Real, dim>;

    /// What a parameter measures, which decides its unit conversion and its periodicity.
    enum class ParamKind : std::uint8_t
    {
        None,   ///< dimensionless (B-spline, line with a vector of given magnitude)
        Length, ///< a length (plane, height of a cylinder or cone)
        Angle,  ///< an angle in radians, periodic on a full turn
    };

    // =========================================================================
    // Curves
    // =========================================================================

    template <std::size_t dim>
    struct CurveDef;

    /// C(t) = origin + t V
    template <std::size_t dim>
    struct LineDef
    {
        Point<dim> origin, vector;
    };

    /// C(t) = center + r1 cos t X + r2 sin t Y, X and Y orthonormal (circle: r1 = r2)
    template <std::size_t dim>
    struct ConicDef
    {
        Point<dim> center, x, y;
        Real r1, r2;
    };

    /// A bounded NURBS in its own parameter (B-spline, Bezier, polyline, composite curve).
    template <std::size_t dim>
    struct SplineDef
    {
        std::shared_ptr<const Curve<Real, dim>> curve;
    };

    /// Basis curve restricted to [t1, t2]; reversed when the curve runs from t2 to t1.
    template <std::size_t dim>
    struct TrimmedDef
    {
        std::shared_ptr<const CurveDef<dim>> basis;
        Real t1, t2;
        bool reversed;
    };

    /// NURBS of an edge with the NURBS parameters of the requested STEP range.
    template <std::size_t dim>
    struct NurbsCurve
    {
        std::shared_ptr<Curve<Real, dim>> curve;
        Real u1, u2;
    };

    namespace detail
    {
        constexpr Real two_pi = 2. * std::numbers::pi;
        constexpr Real inf = std::numeric_limits<Real>::infinity();

        template <std::size_t dim>
        Real closest_parameter(const Curve<Real, dim> &c, const Point<dim> &p)
        {
            const auto [a, b] = c.bounds();
            constexpr int n = 200;
            Real best = a, dmin = inf;
            for (int i = 0; i <= n; ++i)
            {
                const Real t = a + (b - a) * i / n;
                const Real d = sq_norm(c.value(t) - p);
                if (d < dmin)
                    dmin = d, best = t;
            }
            for (int it = 0; it < 30; ++it) // Newton on (C(t) - p) . C'(t) = 0
            {
                const auto r = c.value(best) - p;
                const auto d1 = c.value(best, 1), d2 = c.value(best, 2);
                const Real g = r * d1, h = d1 * d1 + r * d2;
                if (!(std::abs(h) > 0.))
                    break;
                const Real next = std::clamp(best - g / h, a, b);
                if (std::abs(next - best) < 1e-14 * std::max<Real>(1., std::abs(b - a)))
                {
                    best = next;
                    break;
                }
                best = next;
            }
            return best;
        }

        inline Real wrap_to(Real t, Real lo) // t + 2 pi k in [lo, lo + 2 pi)
        {
            return t - two_pi * std::floor((t - lo) / two_pi);
        }
    }

    /**
     * @brief A STEP curve in its own parametrization (see the file comment).
     */
    template <std::size_t dim>
    struct CurveDef
    {
        std::variant<LineDef<dim>, ConicDef<dim>, SplineDef<dim>, TrimmedDef<dim>> def;
        std::uint64_t id{};

        [[nodiscard]] Point<dim> value(Real t) const
        {
            return std::visit([t](const auto &d) -> Point<dim> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return d.origin + t * d.vector;
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                    return d.center + (d.r1 * std::cos(t)) * d.x + (d.r2 * std::sin(t)) * d.y;
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return d.curve->value(t);
                else
                    return d.basis->value(t);
            }, def);
        }

        /// Parameter of the curve point nearest to p (exact for a point on a line or conic).
        [[nodiscard]] Real parameter(const Point<dim> &p) const
        {
            return std::visit([&p](const auto &d) -> Real {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return ((p - d.origin) * d.vector) / sq_norm(d.vector);
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                    return detail::wrap_to(std::atan2(((p - d.center) * d.y) / d.r2, ((p - d.center) * d.x) / d.r1), 0.);
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return detail::closest_parameter(*d.curve, p);
                else
                {
                    const Real t = d.basis->parameter(p);
                    if (!d.basis->periodic())
                        return t;
                    // the turn of t closest to [t1, t2]
                    const Real m = 0.5 * (d.t1 + d.t2);
                    return detail::wrap_to(t, m - std::numbers::pi);
                }
            }, def);
        }

        /// Parameter domain: infinite for a line, [0, 2 pi] for a full conic.
        [[nodiscard]] std::array<Real, 2> domain() const
        {
            return std::visit([](const auto &d) -> std::array<Real, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return {-detail::inf, detail::inf};
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                    return {0., detail::two_pi};
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return d.curve->bounds();
                else
                    return {d.t1, d.t2};
            }, def);
        }

        /// True for a full circle or ellipse: any range of at most 2 pi is valid.
        [[nodiscard]] bool periodic() const { return std::holds_alternative<ConicDef<dim>>(def); }

        /// True for a trimmed curve running against its basis parameter.
        [[nodiscard]] bool reversed() const
        {
            auto t = std::get_if<TrimmedDef<dim>>(&def);
            return t && t->reversed;
        }

        [[nodiscard]] ParamKind kind() const
        {
            if (std::holds_alternative<ConicDef<dim>>(def))
                return ParamKind::Angle;
            if (auto t = std::get_if<TrimmedDef<dim>>(&def))
                return t->basis->kind();
            return ParamKind::None;
        }

        /// NURBS parameter of the point at parameter t, for the NURBS built by to_nurbs(t1, t2).
        [[nodiscard]] Real nurbs_parameter(Real t1, Real t2, Real t) const
        {
            if (std::holds_alternative<ConicDef<dim>>(def))
                return ellipse_arc_parameter(t1, t2, t);
            if (auto tr = std::get_if<TrimmedDef<dim>>(&def))
                return tr->basis->nurbs_parameter(t1, t2, t);
            return t;
        }

        /**
         * @brief Exact NURBS of the curve between t1 < t2 (STEP parameters).
         *
         * A line gives a segment parametrized by t, a conic an arc whose ends are
         * at the parameters t1 and t2 (2 pi at most), a B-spline itself with u1 = t1, u2 = t2.
         */
        [[nodiscard]] NurbsCurve<dim> to_nurbs(Real t1, Real t2) const
        {
            if (!(t2 > t1) || !std::isfinite(t1) || !std::isfinite(t2))
                throw StepError(id, "invalid curve range");
            return std::visit([&](const auto &d) -> NurbsCurve<dim> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                {
                    points_vector<Real, dim> poles{value(t1), value(t2)};
                    return {std::make_shared<BSCurve<Real, dim>>(poles, std::vector<Real>{t1, t1, t2, t2}, 1), t1, t2};
                }
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                {
                    if (t2 - t1 > detail::two_pi * (1. + 1e-12))
                        throw StepError(id, "conic range larger than a turn");
                    auto c = build_ellipse_arc<Real, dim>(d.r1, d.r2, t1, std::min(t2, t1 + detail::two_pi), d.center, d.x, d.y);
                    return {std::make_shared<BSCurveRational<Real, dim>>(std::move(c)), t1, t2};
                }
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return {std::const_pointer_cast<Curve<Real, dim>>(d.curve), t1, t2};
                else
                    return d.basis->to_nurbs(t1, t2);
            }, def);
        }
    };

    // =========================================================================
    // Surfaces
    // =========================================================================

    /// Orthonormal frame of an axis2_placement_3d.
    struct Frame
    {
        Point<3> origin, x, y, z;
        [[nodiscard]] ax2<Real, 3> ax() const { return {origin, z, x}; }
    };

    struct PlaneDef { Frame f; };                          ///< S = O + u X + v Y
    struct CylinderDef { Frame f; Real r; };               ///< S = O + R (cos u X + sin u Y) + v Z
    struct ConeDef { Frame f; Real r, semi_angle; };       ///< S = O + (R + v tan a)(cos u X + sin u Y) + v Z
    struct SphereDef { Frame f; Real r; };                 ///< S = O + R cos v (cos u X + sin u Y) + R sin v Z
    struct TorusDef { Frame f; Real major, minor; };       ///< S = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
    /// S(u, v) = generatrix(v) rotated by u around (origin, z); x: radial direction of the generatrix (u = 0)
    struct RevolutionDef { std::shared_ptr<const CurveDef<3>> generatrix; Point<3> origin, z, x; };
    /// S(u, v) = C(u) + v V
    struct ExtrusionDef { std::shared_ptr<const CurveDef<3>> curve; Point<3> v; };
    /// B-spline surface in its own parameters
    struct SplineSurfaceDef { std::shared_ptr<const Surface<Real, 3>> surface; };
    struct SurfaceDef;
    /// Basis surface restricted to [u1, u2] x [v1, v2]
    struct TrimmedSurfaceDef { std::shared_ptr<const SurfaceDef> basis; Real u1, u2, v1, v2; };

    /// Parameter box of a surface: [u1, u2] x [v1, v2].
    struct UVBox
    {
        Real u1, u2, v1, v2;
    };

    /**
     * @brief A STEP surface in its own parametrization (see the file comment).
     */
    struct SurfaceDef
    {
        std::variant<PlaneDef, CylinderDef, ConeDef, SphereDef, TorusDef, RevolutionDef, ExtrusionDef, SplineSurfaceDef, TrimmedSurfaceDef> def;
        std::uint64_t id{};

        [[nodiscard]] Point<3> value(Real u, Real v) const
        {
            return std::visit([&](const auto &d) -> Point<3> {
                using D = std::decay_t<decltype(d)>;
                auto radial = [u](const Frame &f) { return std::cos(u) * f.x + std::sin(u) * f.y; };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return d.f.origin + u * d.f.x + v * d.f.y;
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return d.f.origin + d.r * radial(d.f) + v * d.f.z;
                else if constexpr (std::is_same_v<D, ConeDef>)
                    return d.f.origin + (d.r + v * std::tan(d.semi_angle)) * radial(d.f) + v * d.f.z;
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return d.f.origin + (d.r * std::cos(v)) * radial(d.f) + (d.r * std::sin(v)) * d.f.z;
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return d.f.origin + (d.major + d.minor * std::cos(v)) * radial(d.f) + (d.minor * std::sin(v)) * d.f.z;
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    const auto p = d.generatrix->value(v) - d.origin;
                    const auto h = (p * d.z) * d.z;
                    const auto r = p - h;
                    return d.origin + h + std::cos(u) * r + std::sin(u) * cross(d.z, r);
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return d.curve->value(u) + v * d.v;
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return d.surface->value(u, v);
                else
                    return d.basis->value(u, v);
            }, def);
        }

        /// Kinds of the u and v parameters.
        [[nodiscard]] std::array<ParamKind, 2> kinds() const
        {
            using K = ParamKind;
            return std::visit([](const auto &d) -> std::array<K, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {K::Length, K::Length};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return {K::Angle, K::Length};
                else if constexpr (std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                    return {K::Angle, K::Angle};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {K::Angle, d.generatrix->kind()};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->kind(), K::None};
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return {K::None, K::None};
                else
                    return d.basis->kinds();
            }, def);
        }

        /// Natural parameter domain (infinite in the open directions, [-pi/2, pi/2] for a sphere latitude…).
        [[nodiscard]] UVBox domain() const
        {
            using detail::inf;
            constexpr Real tp = detail::two_pi, hp = 0.5 * std::numbers::pi;
            return std::visit([](const auto &d) -> UVBox {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {-inf, inf, -inf, inf};
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return {0., tp, -inf, inf};
                else if constexpr (std::is_same_v<D, ConeDef>)
                {
                    // the half cone of the reference radius, up to the apex
                    const Real t = std::tan(d.semi_angle);
                    if (t == 0.)
                        return {0., tp, -inf, inf};
                    const Real apex = -d.r / t;
                    return t > 0. ? UVBox{0., tp, apex, inf} : UVBox{0., tp, -inf, apex};
                }
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return {0., tp, -hp, hp};
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return {0., tp, 0., tp};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    auto [a, b] = d.generatrix->domain();
                    return {0., tp, a, b};
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    auto [a, b] = d.curve->domain();
                    return {a, b, -inf, inf};
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                {
                    auto b = d.surface->bounds();
                    return {b[0], b[1], b[2], b[3]};
                }
                else
                    return {d.u1, d.u2, d.v1, d.v2};
            }, def);
        }

        /// Periodicity of u and v (full turn, not trimmed).
        [[nodiscard]] std::array<bool, 2> periodic() const
        {
            return std::visit([](const auto &d) -> std::array<bool, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef> || std::is_same_v<D, SphereDef>)
                    return {true, false};
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return {true, true};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {true, d.generatrix->periodic()};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->periodic(), false};
                else
                    return {false, false};
            }, def);
        }

        /**
         * @brief Parameters (u, v) of a point on the surface, analytic except for
         * B-spline surfaces (for which this returns nullopt: use the natural domain).
         * u is in [0, 2 pi) in a periodic direction. A point on the axis or at a pole has u = 0.
         */
        [[nodiscard]] std::optional<std::array<Real, 2>> parameters(const Point<3> &p) const
        {
            return std::visit([&p](const auto &d) -> std::optional<std::array<Real, 2>> {
                using D = std::decay_t<decltype(d)>;
                auto angle = [&p](const Frame &f) {
                    const auto q = p - f.origin;
                    return detail::wrap_to(std::atan2(q * f.y, q * f.x), 0.);
                };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return std::array{(p - d.f.origin) * d.f.x, (p - d.f.origin) * d.f.y};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return std::array{angle(d.f), (p - d.f.origin) * d.f.z};
                else if constexpr (std::is_same_v<D, SphereDef>)
                {
                    const auto q = p - d.f.origin;
                    const Real h = q * d.f.z;
                    const Real rho = norm(q - h * d.f.z);
                    return std::array{angle(d.f), std::atan2(h, rho)};
                }
                else if constexpr (std::is_same_v<D, TorusDef>)
                {
                    const auto q = p - d.f.origin;
                    const Real h = q * d.f.z;
                    const Real rho = norm(q - h * d.f.z);
                    return std::array{angle(d.f), detail::wrap_to(std::atan2(h, rho - d.major), 0.)};
                }
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    const Frame f{d.origin, d.x, cross(d.z, d.x), d.z};
                    const Real u = angle(f);
                    // rotate p back into the half-plane of the generatrix
                    const auto q = p - d.origin;
                    const auto h = (q * d.z) * d.z;
                    const auto r = q - h;
                    const auto back = d.origin + h + std::cos(u) * r - std::sin(u) * cross(d.z, r);
                    return std::array{u, d.generatrix->parameter(back)};
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    const Real vv = sq_norm(d.v);
                    Real u;
                    if (auto l = std::get_if<LineDef<3>>(&d.curve->def))
                    {
                        // p - O = u L + v V, least squares
                        const auto w = p - l->origin;
                        const Real a = sq_norm(l->vector), b = l->vector * d.v, det = a * vv - b * b;
                        u = (vv * (w * l->vector) - b * (w * d.v)) / det;
                    }
                    else if (auto c = std::get_if<ConicDef<3>>(&d.curve->def))
                    {
                        // slide p along V onto the plane of the conic
                        const auto n = cross(c->x, c->y);
                        const Real vn = d.v * n;
                        const auto q = std::abs(vn) > 1e-12 * norm(d.v) ? p - (((p - c->center) * n) / vn) * d.v : p;
                        u = d.curve->parameter(q);
                    }
                    else
                        u = d.curve->parameter(p - ((p - d.curve->value(d.curve->domain()[0])) * d.v / vv) * d.v);
                    return std::array{u, ((p - d.curve->value(u)) * d.v) / vv};
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return std::nullopt;
                else
                    return d.basis->parameters(p);
            }, def);
        }

        /// NURBS parameters of the point at (u, v), for the NURBS built by to_nurbs(box).
        [[nodiscard]] std::array<Real, 2> nurbs_parameters(const UVBox &b, Real u, Real v) const
        {
            return std::visit([&](const auto &d) -> std::array<Real, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {u, v};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), v};
                else if constexpr (std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), ellipse_arc_parameter(b.v1, b.v2, v)};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), d.generatrix->nurbs_parameter(b.v1, b.v2, v)};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->nurbs_parameter(b.u1, b.u2, u), v};
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return {u, v};
                else
                    return d.basis->nurbs_parameters(b, u, v);
            }, def);
        }

        /// Exact NURBS of the surface over the box (STEP parameters; a B-spline surface is returned whole).
        [[nodiscard]] std::shared_ptr<Surface<Real, 3>> to_nurbs(const UVBox &b) const
        {
            if (!(b.u2 > b.u1) || !(b.v2 > b.v1) || !std::isfinite(b.u1 + b.u2 + b.v1 + b.v2))
                throw StepError(id, "invalid surface parameter box");
            return std::visit([&](const auto &d) -> std::shared_ptr<Surface<Real, 3>> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                {
                    points_vector<Real, 3> poles{value(b.u1, b.v1), value(b.u2, b.v1), value(b.u1, b.v2), value(b.u2, b.v2)};
                    return std::make_shared<BSSurface<Real, 3>>(poles, std::vector<Real>{b.u1, b.u1, b.u2, b.u2},
                                                                std::vector<Real>{b.v1, b.v1, b.v2, b.v2}, 1, 1);
                }
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_cylinder<Real>(d.r, d.f.ax(), b.v1, b.v2, b.u1, b.u2));
                else if constexpr (std::is_same_v<D, ConeDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_cone<Real>(d.r, d.semi_angle, d.f.ax(), b.v1, b.v2, b.u1, b.u2));
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_sphere<Real>(d.r, d.f.ax(), b.u1, b.u2, b.v1, b.v2));
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_torus<Real>(d.major, d.minor, d.f.ax(), b.u1, b.u2, b.v1, b.v2));
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    auto g = d.generatrix->to_nurbs(b.v1, b.v2).curve;
                    const ax1<Real, 3> axis{d.origin, d.z};
                    if (auto r = std::dynamic_pointer_cast<BSCurveRational<Real, 3>>(g))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_revolution(*r, axis, b.u1, b.u2));
                    if (auto p = std::dynamic_pointer_cast<BSCurve<Real, 3>>(g))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_revolution(*p, axis, b.u1, b.u2));
                    throw StepError(id, "generatrix is not a NURBS");
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    auto c = d.curve->to_nurbs(b.u1, b.u2).curve;
                    if (auto r = std::dynamic_pointer_cast<BSCurveRational<Real, 3>>(c))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_extrusion(*r, d.v, b.v1, b.v2));
                    if (auto p = std::dynamic_pointer_cast<BSCurve<Real, 3>>(c))
                        return std::make_shared<BSSurface<Real, 3>>(build_extrusion(*p, d.v, b.v1, b.v2));
                    throw StepError(id, "swept curve is not a NURBS");
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return std::const_pointer_cast<Surface<Real, 3>>(d.surface);
                else
                    return d.basis->to_nurbs(b);
            }, def);
        }
    };

    /**
     * @brief Parameter box of the face of a surface, from points sampled on its boundary loops
     * (in the target units), enlarged by `margin` of its size (design question 4).
     *
     * - A B-spline or trimmed surface gives its own domain.
     * - In a periodic direction, a loop that winds around the axis makes the face cover a full
     *   turn: the box is [0, 2 pi], the seam of the STEP surface. Otherwise the box is the angular
     *   extent of the loops, which may cross 0 (e.g. [-0.3, 0.4]): such a face does not meet the seam.
     * - In other directions, the extent of the points, clamped to the natural domain (sphere
     *   latitude, cone apex).
     * Points on the axis or at a pole carry no angle and are ignored in u.
     */
    inline UVBox parameter_box(const SurfaceDef &s, std::span<const std::vector<Point<3>>> loops, Real margin = 0.1)
    {
        const auto dom = s.domain();
        if (std::holds_alternative<SplineSurfaceDef>(s.def) || std::holds_alternative<TrimmedSurfaceDef>(s.def))
            return dom;
        const auto per = s.periodic();
        // axis of rotation, to drop points without angle
        std::optional<std::pair<Point<3>, Point<3>>> axis;
        std::visit([&](const auto &d) {
            using D = std::decay_t<decltype(d)>;
            if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef> || std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                axis = std::pair{d.f.origin, d.f.z};
            else if constexpr (std::is_same_v<D, RevolutionDef>)
                axis = std::pair{d.origin, d.z};
        }, s.def);

        Real scale = 0.;
        for (const auto &l : loops)
            for (const auto &p : l)
                scale = std::max(scale, norm(p));
        const Real eps = 1e-9 * std::max<Real>(1., scale);

        std::array<Real, 2> lo{detail::inf, detail::inf}, hi{-detail::inf, -detail::inf};
        std::array<bool, 2> full{false, false};
        std::array<std::optional<Real>, 2> ref; // centre of the first loop in a periodic direction
        for (const auto &l : loops)
        {
            std::array<Real, 2> llo{detail::inf, detail::inf}, lhi{-detail::inf, -detail::inf};
            std::array<std::optional<Real>, 2> first, prev;
            std::array<Real, 2> unwrapped{};
            for (const auto &p : l)
            {
                auto uv = s.parameters(p);
                if (!uv)
                    return dom;
                for (int k = 0; k < 2; ++k)
                {
                    if (k == 0 && per[0] && axis)
                    {
                        const auto q = p - axis->first;
                        if (norm(q - (q * axis->second) * axis->second) < eps)
                            continue; // on the axis: no angle
                    }
                    Real t = (*uv)[k];
                    if (per[k])
                    {
                        if (prev[k])
                        {
                            Real d = t - *prev[k];
                            d -= detail::two_pi * std::round(d / detail::two_pi);
                            unwrapped[k] += d;
                        }
                        else
                            unwrapped[k] = t, first[k] = t;
                        prev[k] = t;
                        t = unwrapped[k];
                    }
                    llo[k] = std::min(llo[k], t);
                    lhi[k] = std::max(lhi[k], t);
                }
            }
            for (int k = 0; k < 2; ++k)
            {
                if (!(lhi[k] >= llo[k]))
                    continue;
                if (per[k])
                {
                    // closing step back to the first point: a winding loop has turned by 2 pi
                    Real d = *first[k] - *prev[k];
                    d -= detail::two_pi * std::round(d / detail::two_pi);
                    if (std::abs(unwrapped[k] + d - *first[k]) > std::numbers::pi)
                        full[k] = true;
                    // bring the loop next to the first one
                    const Real c = 0.5 * (llo[k] + lhi[k]);
                    if (ref[k])
                    {
                        const Real shift = detail::two_pi * std::round((*ref[k] - c) / detail::two_pi);
                        llo[k] += shift, lhi[k] += shift;
                    }
                    else
                        ref[k] = c;
                }
                lo[k] = std::min(lo[k], llo[k]);
                hi[k] = std::max(hi[k], lhi[k]);
            }
        }

        std::array<Real, 4> box{dom.u1, dom.u2, dom.v1, dom.v2};
        const Real size = std::max({hi[0] - lo[0], hi[1] - lo[1], Real(0)});
        for (int k = 0; k < 2; ++k)
        {
            Real &b1 = box[2 * k], &b2 = box[2 * k + 1];
            if (!(hi[k] >= lo[k])) // no usable point (a face reduced to a pole): keep the domain if finite
            {
                if (!std::isfinite(b1) || !std::isfinite(b2))
                    throw StepError(s.id, "cannot bound the surface: no boundary point");
                continue;
            }
            if (per[k])
            {
                // the turn whose centre is in [0, 2 pi): [-20, 30] degrees rather than [340, 390]
                const Real shift = detail::two_pi * std::floor(0.5 * (lo[k] + hi[k]) / detail::two_pi);
                lo[k] -= shift, hi[k] -= shift;
                const Real span = hi[k] - lo[k];
                const Real m = std::min(margin * span, 0.5 * (detail::two_pi - span));
                if (full[k] || m <= 0.)
                    b1 = 0., b2 = detail::two_pi;
                else
                    b1 = lo[k] - m, b2 = hi[k] + m;
                continue;
            }
            Real m = margin * (hi[k] - lo[k]);
            if (!(m > 0.))
                m = margin * size > 0. ? margin * size : Real(1);
            b1 = std::max(b1, lo[k] - m);
            b2 = std::min(b2, hi[k] + m);
        }
        return {box[0], box[1], box[2], box[3]};
    }

    // =========================================================================
    // Reader
    // =========================================================================

    namespace detail
    {
        // Insert the knot u once into a clamped or unclamped knot vector, for every row of poles (Boehm).
        template <typename P>
        void insert_knot(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p, Real u)
        {
            const auto span = static_cast<std::size_t>(std::upper_bound(k.begin(), k.end(), u) - k.begin()) - 1;
            const auto s = static_cast<std::size_t>(std::count(k.begin(), k.end(), u));
            for (auto &poles : rows)
            {
                std::vector<P> q;
                q.reserve(poles.size() + 1);
                for (std::size_t i = 0; i <= poles.size(); ++i)
                {
                    if (i + p <= span)
                        q.push_back(poles[i]);
                    else if (i + s > span)
                        q.push_back(poles[i - 1]);
                    else
                    {
                        const Real a = (u - k[i]) / (k[i + p] - k[i]);
                        q.push_back(a * poles[i] + (1. - a) * poles[i - 1]);
                    }
                }
                poles = std::move(q);
            }
            k.insert(k.begin() + static_cast<std::ptrdiff_t>(span) + 1, u);
        }

        // Clamp the start of a knot vector at its domain start k[p] (multiplicity p + 1).
        template <typename P>
        void clamp_start(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p)
        {
            const Real a = k[p];
            while (static_cast<std::size_t>(std::count(k.begin(), k.end(), a)) < p)
                insert_knot(k, rows, p, a);
            const auto j = static_cast<std::size_t>(std::lower_bound(k.begin(), k.end(), a) - k.begin());
            if (j == 0)
                return; // already clamped
            k.erase(k.begin(), k.begin() + static_cast<std::ptrdiff_t>(j));
            k.insert(k.begin(), a);
            for (auto &poles : rows)
                poles.erase(poles.begin(), poles.begin() + static_cast<std::ptrdiff_t>(j - 1));
        }

        // Clamp both ends of a (possibly unclamped, e.g. periodic) knot vector.
        template <typename P>
        void clamp(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p)
        {
            clamp_start(k, rows, p);
            std::ranges::reverse(k);
            for (auto &x : k)
                x = -x;
            for (auto &poles : rows)
                std::ranges::reverse(poles);
            clamp_start(k, rows, p);
            std::ranges::reverse(k);
            for (auto &x : k)
                x = -x;
            for (auto &poles : rows)
                std::ranges::reverse(poles);
        }

        // Flat knots of the B-spline subtypes defined by their degree and number of poles.
        inline std::vector<Real> implicit_knots(std::string_view kind, std::size_t p, std::size_t n, std::uint64_t id)
        {
            std::vector<Real> k;
            if (kind == "BEZIER")
            {
                // n = p * segments + 1 ; knots 0..segments, ends p + 1 times, inner p times
                if (p == 0 || (n - 1) % p != 0)
                    throw StepError(id, "Bezier poles do not match the degree");
                const std::size_t segs = (n - 1) / p;
                for (std::size_t i = 0; i <= segs; ++i)
                    k.insert(k.end(), i == 0 || i == segs ? p + 1 : p, Real(i));
            }
            else if (kind == "QUASI_UNIFORM")
            {
                const std::size_t segs = n - p;
                for (std::size_t i = 0; i <= segs; ++i)
                    k.insert(k.end(), i == 0 || i == segs ? p + 1 : 1, Real(i));
            }
            else // UNIFORM: unclamped, k_i = i - p, clamped afterwards
                for (std::size_t i = 0; i < n + p + 1; ++i)
                    k.push_back(Real(i) - Real(p));
            return k;
        }

        inline std::vector<Real> flat_knots_of(ListView mults, ListView knots, std::uint64_t id)
        {
            if (mults.size() != knots.size())
                throw StepError(id, "knot multiplicities and knots differ in size");
            std::vector<Real> k;
            for (std::size_t i = 0; i < knots.size(); ++i)
                k.insert(k.end(), static_cast<std::size_t>(mults[i].as_int()), knots[i].as_real());
            return k;
        }

        inline bool clamped(const std::vector<Real> &k, std::size_t p)
        {
            return k.size() > 2 * p + 1 && k[0] == k[p] && k[k.size() - 1] == k[k.size() - 1 - p];
        }
    }

    /**
     * @brief Reads STEP geometric entities of a file into definitions, in target units.
     *
     * Results are cached by #id. Errors: StepUnsupported for a valid but unsupported entity,
     * StepError for invalid content, P21AccessError for arguments of the wrong kind.
     * 2D entities (pcurves of definitional representations) are read in the parameter space of
     * their STEP surface, without unit conversion.
     */
    class GeometryReader
    {
        const P21File &f_;
        Units u_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const CurveDef<2>>> curves2_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const CurveDef<3>>> curves3_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const SurfaceDef>> surfaces_;

        template <std::size_t dim>
        auto &cache()
        {
            if constexpr (dim == 2)
                return curves2_;
            else
                return curves3_;
        }

        // Own attributes of entity `type` in a simple instance of a subtype (after `skip`
        // inherited ones) or in the partial entity of a complex instance.
        static std::vector<ValueView> own(const InstanceView &in, std::string_view type, std::size_t skip, std::size_t count)
        {
            ListView a = in.is_complex() ? in.args(type) : in.args();
            const std::size_t first = in.is_complex() ? 0 : skip;
            if (a.size() < first + count)
                throw StepError(in.id(), "missing arguments of " + std::string(type));
            std::vector<ValueView> r;
            for (std::size_t i = 0; i < count; ++i)
                r.push_back(a[first + i]);
            return r;
        }

        // A simple instance of one of the types, or a complex one with one of them as partial type
        // (a simple instance does not list its supertypes).
        static bool is_one_of(const InstanceView &in, std::initializer_list<std::string_view> types)
        {
            return std::ranges::any_of(types, [&](std::string_view t) { return in.has_type(t); });
        }

        template <std::size_t dim>
        Point<dim> coords(std::uint64_t id, Real scale)
        {
            const auto in = f_.instance(id);
            const auto l = in.args()[1].as_list();
            if (l.size() < dim)
                throw StepError(id, "expected " + std::to_string(dim) + " coordinates");
            Point<dim> p{};
            for (std::size_t i = 0; i < dim; ++i)
                p[i] = l[i].as_real() * scale;
            return p;
        }

    public:
        GeometryReader(const P21File &f, Units u) : f_{f}, u_{std::move(u)} {}

        [[nodiscard]] const Units &units() const noexcept { return u_; }
        [[nodiscard]] const P21File &file() const noexcept { return f_; }

        /// Scale of a length in 3D, none in 2D (parameter space).
        template <std::size_t dim>
        [[nodiscard]] Real length_scale() const { return dim == 3 ? u_.length : 1.; }

        template <std::size_t dim>
        Point<dim> point(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            if (!in.has_type("CARTESIAN_POINT"))
                throw StepUnsupported(id, in.type());
            return coords<dim>(id, length_scale<dim>());
        }

        /// Unit direction (DIRECTION).
        template <std::size_t dim>
        Point<dim> direction(std::uint64_t id)
        {
            auto d = coords<dim>(id, 1.);
            const Real n = norm(d);
            if (!(n > 0.))
                throw StepError(id, "null direction");
            return d / n;
        }

        /// VECTOR(orientation, magnitude) as a vector in target units.
        template <std::size_t dim>
        Point<dim> vector(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            return (a[2].as_real() * length_scale<dim>()) * direction<dim>(a[1].as_ref());
        }

        /// AXIS2_PLACEMENT_3D: orthonormal frame, defaults as ISO 10303-42 (axis z, reference x).
        Frame frame(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            if (!in.has_type("AXIS2_PLACEMENT_3D"))
                throw StepUnsupported(id, in.type());
            const auto a = in.args();
            const auto o = point<3>(a[1].as_ref());
            const Point<3> z = a[2].is_null() ? Point<3>{0., 0., 1.} : direction<3>(a[2].as_ref());
            Point<3> r = a[3].is_null() ? Point<3>{1., 0., 0.} : direction<3>(a[3].as_ref());
            if (a[3].is_null() && norm(cross(z, r)) < 1e-12)
                r = {0., 0., 1.};
            auto x = r - (r * z) * z;
            const Real n = norm(x);
            if (!(n > 1e-12))
                throw StepError(id, "reference direction parallel to the axis");
            x = x / n;
            return {o, x, cross(z, x), z};
        }

        /// AXIS2_PLACEMENT_2D: origin and orthonormal (x, y).
        std::array<Point<2>, 3> frame2d(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            const auto o = point<2>(a[1].as_ref());
            const Point<2> x = a[2].is_null() ? Point<2>{1., 0.} : direction<2>(a[2].as_ref());
            return {o, x, Point<2>{-x[1], x[0]}};
        }

        /// AXIS1_PLACEMENT: origin and unit axis.
        std::pair<Point<3>, Point<3>> axis1(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            return {point<3>(a[1].as_ref()), a[2].is_null() ? Point<3>{0., 0., 1.} : direction<3>(a[2].as_ref())};
        }

        /// Curve of any supported type, 2D or 3D.
        template <std::size_t dim>
        std::shared_ptr<const CurveDef<dim>> curve(std::uint64_t id)
        {
            auto &c = cache<dim>();
            if (auto it = c.find(id); it != c.end())
                return it->second;
            auto r = std::make_shared<const CurveDef<dim>>(read_curve<dim>(id));
            c.emplace(id, r);
            return r;
        }

        /// Surface of any supported type.
        std::shared_ptr<const SurfaceDef> surface(std::uint64_t id)
        {
            if (auto it = surfaces_.find(id); it != surfaces_.end())
                return it->second;
            auto r = std::make_shared<const SurfaceDef>(read_surface(id));
            surfaces_.emplace(id, r);
            return r;
        }

    private:
        template <std::size_t dim>
        CurveDef<dim> read_curve(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            const auto a = in.args();
            const auto t = in.type();
            CurveDef<dim> c{.id = id};
            if (t == "LINE")
                c.def = LineDef<dim>{point<dim>(a[1].as_ref()), vector<dim>(a[2].as_ref())};
            else if (t == "CIRCLE" || t == "ELLIPSE")
            {
                const Real s = length_scale<dim>();
                const Real r1 = a[2].as_real() * s, r2 = (t == "CIRCLE" ? a[2] : a[3]).as_real() * s;
                if (!(r1 > 0.) || !(r2 > 0.))
                    throw StepError(id, "radius must be positive");
                if constexpr (dim == 3)
                {
                    const auto fr = frame(a[1].as_ref());
                    c.def = ConicDef<3>{fr.origin, fr.x, fr.y, r1, r2};
                }
                else
                {
                    const auto [o, x, y] = frame2d(a[1].as_ref());
                    c.def = ConicDef<2>{o, x, y, r1, r2};
                }
            }
            else if (is_one_of(in, {"B_SPLINE_CURVE", "B_SPLINE_CURVE_WITH_KNOTS", "BEZIER_CURVE", "QUASI_UNIFORM_CURVE", "UNIFORM_CURVE"}))
                c.def = SplineDef<dim>{read_bspline_curve<dim>(in)};
            else if (t == "POLYLINE")
            {
                points_vector<Real, dim> poles;
                for (auto p : a[1].as_list())
                    poles.push_back(point<dim>(p.as_ref()));
                if (poles.size() < 2)
                    throw StepError(id, "polyline with less than two points");
                std::vector<Real> k{0.};
                for (std::size_t i = 0; i < poles.size(); ++i)
                    k.push_back(Real(i));
                k.push_back(Real(poles.size() - 1));
                c.def = SplineDef<dim>{std::make_shared<BSCurve<Real, dim>>(poles, k, 1)};
            }
            else if (t == "TRIMMED_CURVE")
                c.def = read_trimmed<dim>(in);
            else if (t == "COMPOSITE_CURVE")
                c.def = SplineDef<dim>{read_composite<dim>(in)};
            else if (t == "SURFACE_CURVE" || t == "SEAM_CURVE" || t == "INTERSECTION_CURVE")
                return CurveDef<dim>{curve<dim>(a[1].as_ref())->def, id};
            else if (t == "PCURVE" && dim == 2)
                throw StepUnsupported(id, t); // read through its definitional representation
            else
                throw StepUnsupported(id, t);
            return c;
        }

        template <std::size_t dim>
        std::shared_ptr<const Curve<Real, dim>> read_bspline_curve(const InstanceView &in)
        {
            const auto id = in.id();
            const auto b = own(in, "B_SPLINE_CURVE", 1, 5);
            const auto p = static_cast<std::size_t>(b[0].as_int());
            std::vector<Point<dim>> cps;
            for (auto r : b[1].as_list())
                cps.push_back(point<dim>(r.as_ref()));
            const auto n = cps.size();
            std::vector<Real> k;
            if (in.has_type("B_SPLINE_CURVE_WITH_KNOTS"))
            {
                const auto kn = own(in, "B_SPLINE_CURVE_WITH_KNOTS", 6, 3);
                k = detail::flat_knots_of(kn[0].as_list(), kn[1].as_list(), id);
            }
            else if (in.has_type("BEZIER_CURVE"))
                k = detail::implicit_knots("BEZIER", p, n, id);
            else if (in.has_type("QUASI_UNIFORM_CURVE"))
                k = detail::implicit_knots("QUASI_UNIFORM", p, n, id);
            else if (in.has_type("UNIFORM_CURVE"))
                k = detail::implicit_knots("UNIFORM", p, n, id);
            else
                throw StepUnsupported(id, in.type());
            if (k.size() != n + p + 1)
                throw StepError(id, "B-spline poles, knots and degree do not match");

            std::vector<Real> w;
            if (auto rp = in.part("RATIONAL_B_SPLINE_CURVE"))
                for (auto x : rp->args()[0].as_list())
                    w.push_back(x.as_real());
            if (!w.empty() && w.size() != n)
                throw StepError(id, "weights and poles differ in size");
            if (!w.empty())
            {
                std::vector<std::vector<std::array<Real, dim + 1>>> rows(1);
                for (std::size_t i = 0; i < n; ++i)
                {
                    if (!(w[i] > 0.))
                        throw StepError(id, "weights must be positive");
                    std::array<Real, dim + 1> h;
                    for (std::size_t j = 0; j < dim; ++j)
                        h[j] = cps[i][j] * w[i]; // homogeneous coordinates
                    h[dim] = w[i];
                    rows[0].push_back(h);
                }
                if (!detail::clamped(k, p))
                    detail::clamp(k, rows, p);
                return std::make_shared<BSCurveRational<Real, dim>>(rows[0], k, p);
            }
            std::vector<std::vector<Point<dim>>> rows{cps};
            if (!detail::clamped(k, p))
                detail::clamp(k, rows, p);
            return std::make_shared<BSCurve<Real, dim>>(rows[0], k, p);
        }

        template <std::size_t dim>
        TrimmedDef<dim> read_trimmed(const InstanceView &in)
        {
            const auto id = in.id();
            const auto a = in.args();
            auto basis = curve<dim>(a[1].as_ref());
            if (std::holds_alternative<TrimmedDef<dim>>(basis->def))
                throw StepUnsupported(id, "TRIMMED_CURVE of a trimmed curve");
            const bool sense = a[4].as_bool();
            const bool prefer_point = a.size() > 5 && a[5].is(ValueKind::Enum) && a[5].as_enum() == "CARTESIAN";
            auto trim = [&](ValueView v) -> Real {
                std::optional<Real> param, by_point;
                for (auto s : v.as_list())
                {
                    if (s.is(ValueKind::Reference))
                        by_point = basis->parameter(point<dim>(s.as_ref()));
                    else if (s.is(ValueKind::Typed))
                    {
                        Real t = s.as_measure();
                        if (basis->kind() == ParamKind::Angle)
                            t *= u_.angle;
                        param = t;
                    }
                }
                if (by_point && (prefer_point || !param))
                    return *by_point;
                if (param)
                    return *param;
                throw StepError(id, "trimming select without parameter nor point");
            };
            Real t1 = trim(a[2]), t2 = trim(a[3]);
            if (basis->periodic())
            {
                // from t1 to t2 in the sense of the curve, at most a turn (equal trims: a full turn)
                if (sense)
                    t2 = detail::wrap_to(t2, t1 + 1e-12);
                else
                    t1 = detail::wrap_to(t1, t2 + 1e-12);
            }
            if (!sense)
                std::swap(t1, t2);
            if (!(t2 > t1))
                throw StepError(id, "empty trimmed curve");
            return {basis, t1, t2, !sense};
        }

        template <std::size_t dim>
        std::shared_ptr<const Curve<Real, dim>> read_composite(const InstanceView &in)
        {
            std::optional<BSCurveRational<Real, dim>> joined;
            for (auto sref : in.args()[1].as_list())
            {
                const auto seg = f_.instance(sref.as_ref());
                if (!seg.has_type("COMPOSITE_CURVE_SEGMENT"))
                    throw StepUnsupported(seg.id(), seg.type());
                const auto sa = seg.args("COMPOSITE_CURVE_SEGMENT");
                const bool same = sa[1].as_bool();
                auto c = curve<dim>(sa[2].as_ref());
                const auto [t1, t2] = c->domain();
                if (!std::isfinite(t1) || !std::isfinite(t2))
                    throw StepError(seg.id(), "unbounded composite curve segment");
                auto nc = c->to_nurbs(t1, t2);
                BSCurveRational<Real, dim> r = [&] {
                    if (auto rr = std::dynamic_pointer_cast<BSCurveRational<Real, dim>>(nc.curve))
                        return *rr;
                    return BSCurveRational<Real, dim>{*std::dynamic_pointer_cast<BSCurve<Real, dim>>(nc.curve)};
                }();
                if (nc.u1 > r.bounds()[0] || nc.u2 < r.bounds()[1])
                    r.trim(nc.u1, nc.u2);
                if (same == c->reversed())
                    r.reverse();
                if (!joined)
                    joined = r;
                else
                {
                    if (norm(joined->end() - r.begin()) > 1e-6 * std::max<Real>(1., norm(r.begin())))
                        throw StepError(in.id(), "composite curve segments are not connected");
                    joined = BSCurveRational<Real, dim>{join(*joined, r)};
                }
            }
            if (!joined)
                throw StepError(in.id(), "empty composite curve");
            return std::make_shared<BSCurveRational<Real, dim>>(*joined);
        }

        SurfaceDef read_surface(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            const auto a = in.args();
            const auto t = in.type();
            const Real L = u_.length;
            SurfaceDef s{.id = id};
            auto positive = [&](Real r) {
                if (!(r > 0.))
                    throw StepError(id, "radius must be positive");
                return r;
            };
            if (t == "PLANE")
                s.def = PlaneDef{frame(a[1].as_ref())};
            else if (t == "CYLINDRICAL_SURFACE")
                s.def = CylinderDef{frame(a[1].as_ref()), positive(a[2].as_real() * L)};
            else if (t == "CONICAL_SURFACE")
            {
                const Real r = a[2].as_real() * L, alpha = a[3].as_measure() * u_.angle;
                if (r < 0. || !(std::abs(alpha) < 0.5 * std::numbers::pi) || alpha == 0.)
                    throw StepError(id, "invalid cone radius or semi angle");
                s.def = ConeDef{frame(a[1].as_ref()), r, alpha};
            }
            else if (t == "SPHERICAL_SURFACE")
                s.def = SphereDef{frame(a[1].as_ref()), positive(a[2].as_real() * L)};
            else if (t == "TOROIDAL_SURFACE" || t == "DEGENERATE_TOROIDAL_SURFACE")
                s.def = TorusDef{frame(a[1].as_ref()), a[2].as_real() * L, positive(a[3].as_real() * L)};
            else if (is_one_of(in, {"B_SPLINE_SURFACE", "B_SPLINE_SURFACE_WITH_KNOTS", "BEZIER_SURFACE", "QUASI_UNIFORM_SURFACE", "UNIFORM_SURFACE"}))
                s.def = SplineSurfaceDef{read_bspline_surface(in)};
            else if (t == "SURFACE_OF_REVOLUTION")
            {
                auto g = curve<3>(a[1].as_ref());
                const auto [o, z] = axis1(a[2].as_ref());
                // angle 0: radial direction of the first generatrix point off the axis
                const auto [g1, g2] = g->domain();
                const Real lo = std::isfinite(g1) ? g1 : (std::isfinite(g2) ? g2 - 1. : 0.);
                const Real hi = std::isfinite(g2) ? g2 : lo + 1.;
                std::optional<Point<3>> x;
                for (int i = 0; i <= 16 && !x; ++i)
                {
                    const auto q = g->value(lo + (hi - lo) * i / 16.) - o;
                    const auto r = q - (q * z) * z;
                    if (norm(r) > 1e-9 * std::max<Real>(1., norm(q)))
                        x = r / norm(r);
                }
                if (!x)
                    throw StepError(id, "generatrix on the axis of revolution");
                s.def = RevolutionDef{g, o, z, *x};
            }
            else if (t == "SURFACE_OF_LINEAR_EXTRUSION")
            {
                const auto v = vector<3>(a[2].as_ref());
                if (!(norm(v) > 0.))
                    throw StepError(id, "null extrusion vector");
                s.def = ExtrusionDef{curve<3>(a[1].as_ref()), v};
            }
            else if (t == "RECTANGULAR_TRIMMED_SURFACE")
            {
                auto basis = surface(a[1].as_ref());
                const auto k = basis->kinds();
                auto conv = [&](ValueView v, ParamKind kind) {
                    const Real x = v.as_measure();
                    return kind == ParamKind::Angle ? x * u_.angle : kind == ParamKind::Length ? x * L : x;
                };
                Real u1 = conv(a[2], k[0]), u2 = conv(a[3], k[0]), v1 = conv(a[4], k[1]), v2 = conv(a[5], k[1]);
                if (u1 > u2)
                    std::swap(u1, u2);
                if (v1 > v2)
                    std::swap(v1, v2);
                s.def = TrimmedSurfaceDef{basis, u1, u2, v1, v2};
            }
            else if (t == "OFFSET_SURFACE")
                s.def = read_offset(in);
            else
                throw StepUnsupported(id, t);
            return s;
        }

        // Offset of an elementary surface: the elementary surface of offset radius.
        decltype(SurfaceDef::def) read_offset(const InstanceView &in)
        {
            const auto a = in.args();
            const auto basis = surface(a[1].as_ref());
            const Real d = a[2].as_measure() * u_.length;
            return std::visit([&](const auto &b) -> decltype(SurfaceDef::def) {
                using D = std::decay_t<decltype(b)>;
                auto positive = [&](Real r) {
                    if (!(r > 0.))
                        throw StepError(in.id(), "offset gives a null radius");
                    return r;
                };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return PlaneDef{{b.f.origin + d * b.f.z, b.f.x, b.f.y, b.f.z}};
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return CylinderDef{b.f, positive(b.r + d)};
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return SphereDef{b.f, positive(b.r + d)};
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return TorusDef{b.f, b.major, positive(b.minor + d)};
                else if constexpr (std::is_same_v<D, ConeDef>)
                {
                    // normal (cos a, -sin a) in the (radial, axis) plane: same half-angle and origin,
                    // radius R + d / cos a (the v parameter shifts by -d sin a)
                    const Real r = b.r + d / std::cos(b.semi_angle);
                    if (r < 0.)
                        throw StepError(in.id(), "offset cone beyond its apex");
                    return ConeDef{b.f, r, b.semi_angle};
                }
                else
                    throw StepUnsupported(in.id(), "OFFSET_SURFACE of a non elementary surface");
            }, basis->def);
        }

        std::shared_ptr<const Surface<Real, 3>> read_bspline_surface(const InstanceView &in)
        {
            const auto id = in.id();
            const auto b = own(in, "B_SPLINE_SURFACE", 1, 7);
            const auto pu = static_cast<std::size_t>(b[0].as_int()), pv = static_cast<std::size_t>(b[1].as_int());
            const auto grid = b[2].as_list();
            const std::size_t nu = grid.size();
            const std::size_t nv = nu ? grid[0].as_list().size() : 0;
            if (nu == 0 || nv == 0)
                throw StepError(id, "empty control point grid");
            std::vector<Real> ku, kv;
            if (in.has_type("B_SPLINE_SURFACE_WITH_KNOTS"))
            {
                const auto kn = own(in, "B_SPLINE_SURFACE_WITH_KNOTS", 8, 5);
                ku = detail::flat_knots_of(kn[0].as_list(), kn[2].as_list(), id);
                kv = detail::flat_knots_of(kn[1].as_list(), kn[3].as_list(), id);
            }
            else
            {
                std::string_view kind = in.has_type("BEZIER_SURFACE")          ? "BEZIER"
                                        : in.has_type("QUASI_UNIFORM_SURFACE") ? "QUASI_UNIFORM"
                                        : in.has_type("UNIFORM_SURFACE")       ? "UNIFORM"
                                                                               : "";
                if (kind.empty())
                    throw StepUnsupported(id, in.type());
                ku = detail::implicit_knots(kind, pu, nu, id);
                kv = detail::implicit_knots(kind, pv, nv, id);
            }
            if (ku.size() != nu + pu + 1 || kv.size() != nv + pv + 1)
                throw StepError(id, "B-spline surface poles, knots and degrees do not match");
            std::optional<ListView> weights;
            if (auto rp = in.part("RATIONAL_B_SPLINE_SURFACE"))
                weights = rp->args()[0].as_list();

            // rows along u (one per v), homogeneous 4D poles; STEP stores [u][v]
            std::vector<std::vector<std::array<Real, 4>>> rows(nv, std::vector<std::array<Real, 4>>(nu));
            for (std::size_t i = 0; i < nu; ++i)
            {
                const auto col = grid[i].as_list();
                if (col.size() != nv)
                    throw StepError(id, "ragged control point grid");
                for (std::size_t j = 0; j < nv; ++j)
                {
                    const auto p = point<3>(col[j].as_ref());
                    const Real w = weights ? (*weights)[i].as_list()[j].as_real() : 1.;
                    if (!(w > 0.))
                        throw StepError(id, "weights must be positive");
                    rows[j][i] = {p[0] * w, p[1] * w, p[2] * w, w};
                }
            }
            if (!detail::clamped(ku, pu))
                detail::clamp(ku, rows, pu);
            if (!detail::clamped(kv, pv))
            {
                // transpose: rows along v (one per u)
                std::vector<std::vector<std::array<Real, 4>>> cols(rows[0].size(), std::vector<std::array<Real, 4>>(rows.size()));
                for (std::size_t j = 0; j < rows.size(); ++j)
                    for (std::size_t i = 0; i < rows[j].size(); ++i)
                        cols[i][j] = rows[j][i];
                detail::clamp(kv, cols, pv);
                rows.assign(cols[0].size(), std::vector<std::array<Real, 4>>(cols.size()));
                for (std::size_t i = 0; i < cols.size(); ++i)
                    for (std::size_t j = 0; j < cols[i].size(); ++j)
                        rows[j][i] = cols[i][j];
            }
            std::vector<std::array<Real, 4>> hp; // u fastest
            for (const auto &r : rows)
                hp.insert(hp.end(), r.begin(), r.end());
            if (weights)
                return std::make_shared<BSSurfaceRational<Real, 3>>(hp, ku, kv, pu, pv);
            points_vector<Real, 3> poles;
            poles.reserve(hp.size());
            for (const auto &h : hp)
                poles.push_back({h[0], h[1], h[2]});
            return std::make_shared<BSSurface<Real, 3>>(poles, ku, kv, pu, pv);
        }
    };
}
