#pragma once

/**
 * @file geometry_writer.h
 * @brief STEP writer, geometry layer: gbs curves and surfaces to STEP entities.
 *
 * Design: docs/sources/design/step_writer.md, sections 3 and 6; architecture
 * note docs/sources/design/step_wpr02_geometry.md.
 *
 * - Primitives (gbs/primitives.h) become the analytic entities (line, circle,
 *   ellipse, plane, cylindrical, conical, spherical, toroidal surfaces), with
 *   the same parametrization: pcurves stay valid as they are.
 * - B-splines become b_spline_curve_with_knots / b_spline_surface_with_knots
 *   (complex instances when rational, Cartesian poles and weights), a curve
 *   restricted to the requested range; 2D degree-1 segments become lines.
 * - A curve on a surface along an iso line of a primitive or of a B-spline
 *   becomes that iso line, exactly.
 * - Anything else is approximated within approx_tol, parameter preserved
 *   (gbs-io/to_nurbs.h), and counted.
 *
 * Coordinates are written as they are, in the model unit described by the
 * representation context (context()).
 */

#include <gbs-io/step/p21_writer.h>
#include <gbs-io/to_nurbs.h>
#include <gbs/primitives.h>
#include <gbs/curveonsurface.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace gbs::step
{
    class GeometryWriter
    {
        using T = double;
        using V = P21Value;

        P21Writer &w_;
        T approx_tol_;
        nurbs_export::ApproxStats<T> stats_;
        std::uint64_t uncertainty_{}, length_{}; // uncertainty measure and length unit of the context

        static void knots_and_mults(const std::vector<T> &flat, std::vector<std::int64_t> &mults, std::vector<T> &knots)
        {
            std::vector<std::size_t> m;
            unflat_knots(flat, m, knots);
            mults.assign(m.begin(), m.end());
        }

        template <std::size_t dim>
        std::uint64_t point_dim(const point<T, dim> &p)
        {
            std::vector<T> c(p.begin(), p.end());
            return w_.add("CARTESIAN_POINT", {V::string(""), V::reals(c)});
        }

        template <std::size_t dim>
        std::uint64_t direction_dim(const point<T, dim> &d)
        {
            const T n = norm(d);
            if (!(n > T(0)))
                throw std::invalid_argument("STEP direction: null vector");
            std::vector<T> c;
            for (auto x : d)
                c.push_back(x / n);
            return w_.add("DIRECTION", {V::string(""), V::reals(c)});
        }

        template <std::size_t dim>
        std::uint64_t line_dim(const point<T, dim> &origin, const point<T, dim> &vector)
        {
            const auto vec = w_.add("VECTOR", {V::string(""), V::ref(direction_dim<dim>(vector)), V::real(norm(vector))});
            return w_.add("LINE", {V::string(""), V::ref(point_dim<dim>(origin)), V::ref(vec)});
        }

        // B-spline curve (homogeneous poles when rational) as STEP entities.
        template <std::size_t dim, bool rational>
        std::uint64_t bspline_curve(const BSCurveGeneral<T, dim, rational> &c)
        {
            std::vector<std::uint64_t> cps;
            std::vector<T> weights;
            for (const auto &h : c.poles())
            {
                point<T, dim> p;
                const T wgt = rational ? h[dim] : T(1);
                for (std::size_t i = 0; i < dim; ++i)
                    p[i] = h[i] / wgt; // Cartesian poles, the weights apart
                cps.push_back(point_dim<dim>(p));
                weights.push_back(wgt);
            }
            std::vector<std::int64_t> mults;
            std::vector<T> knots;
            knots_and_mults(c.knotsFlats(), mults, knots);
            const auto deg = V::integer(static_cast<std::int64_t>(c.degree()));
            if constexpr (!rational)
                return w_.add("B_SPLINE_CURVE_WITH_KNOTS", {V::string(""), deg, V::refs(cps), V::enumeration("UNSPECIFIED"), V::boolean(false),
                                                            V::boolean(false), V::integers(mults), V::reals(knots), V::enumeration("UNSPECIFIED")});
            else
                return w_.add_complex({{"BOUNDED_CURVE", {}},
                                       {"B_SPLINE_CURVE", {deg, V::refs(cps), V::enumeration("UNSPECIFIED"), V::boolean(false), V::boolean(false)}},
                                       {"B_SPLINE_CURVE_WITH_KNOTS", {V::integers(mults), V::reals(knots), V::enumeration("UNSPECIFIED")}},
                                       {"CURVE", {}},
                                       {"GEOMETRIC_REPRESENTATION_ITEM", {}},
                                       {"RATIONAL_B_SPLINE_CURVE", {V::reals(weights)}},
                                       {"REPRESENTATION_ITEM", {V::string("")}}});
        }

        template <std::size_t dim>
        std::uint64_t any_bspline_curve(const nurbs_export::AnyBSCurve<T, dim> &c)
        {
            return std::visit([this](const auto &b) { return bspline_curve(b); }, c);
        }

        template <bool rational>
        std::uint64_t bspline_surface(const BSSurfaceGeneral<T, 3, rational> &s)
        {
            const auto nu = s.nPolesU(), nv = s.nPolesV();
            const auto &poles = s.poles();
            std::vector<P21Value> grid, wgrid; // STEP [u][v], gbs u fastest
            for (std::size_t i = 0; i < nu; ++i)
            {
                std::vector<std::uint64_t> row;
                std::vector<T> wrow;
                for (std::size_t j = 0; j < nv; ++j)
                {
                    const auto &h = poles[i + nu * j];
                    const T wgt = rational ? h[3] : T(1);
                    row.push_back(point_dim<3>({h[0] / wgt, h[1] / wgt, h[2] / wgt}));
                    wrow.push_back(wgt);
                }
                grid.push_back(V::refs(row));
                wgrid.push_back(V::reals(wrow));
            }
            std::vector<std::int64_t> mu, mv;
            std::vector<T> ku, kv;
            knots_and_mults(s.knotsFlatsU(), mu, ku);
            knots_and_mults(s.knotsFlatsV(), mv, kv);
            const auto du = V::integer(static_cast<std::int64_t>(s.degreeU())), dv = V::integer(static_cast<std::int64_t>(s.degreeV()));
            const auto f = V::boolean(false);
            if constexpr (!rational)
                return w_.add("B_SPLINE_SURFACE_WITH_KNOTS", {V::string(""), du, dv, V::list(grid), V::enumeration("UNSPECIFIED"), f, f, f,
                                                              V::integers(mu), V::integers(mv), V::reals(ku), V::reals(kv), V::enumeration("UNSPECIFIED")});
            else
                return w_.add_complex({{"BOUNDED_SURFACE", {}},
                                       {"B_SPLINE_SURFACE", {du, dv, V::list(grid), V::enumeration("UNSPECIFIED"), f, f, f}},
                                       {"B_SPLINE_SURFACE_WITH_KNOTS", {V::integers(mu), V::integers(mv), V::reals(ku), V::reals(kv), V::enumeration("UNSPECIFIED")}},
                                       {"GEOMETRIC_REPRESENTATION_ITEM", {}},
                                       {"RATIONAL_B_SPLINE_SURFACE", {V::list(wgrid)}},
                                       {"REPRESENTATION_ITEM", {V::string("")}},
                                       {"SURFACE", {}}});
        }

        std::uint64_t elementary_curve(const ElementaryCurve<T> &c)
        {
            if (auto l = dynamic_cast<const LineCurve<T> *>(&c))
                return line_dim<3>(l->origin(), l->vector());
            const auto &e = dynamic_cast<const EllipseCurve<T> &>(c);
            const auto ax = axis2(e.frame());
            if (e.is_circle())
                return w_.add("CIRCLE", {V::string(""), V::ref(ax), V::real(e.radius1())});
            return w_.add("ELLIPSE", {V::string(""), V::ref(ax), V::real(e.radius1()), V::real(e.radius2())});
        }

        // The iso line of a curve on a surface whose pcurve is an iso segment parametrized like the line.
        std::optional<std::uint64_t> iso_line(const CurveOnSurface<T, 3> &cos, T a, T b)
        {
            const auto iso = nurbs_export::iso_curve(cos, a, b);
            if (!iso)
                return std::nullopt;
            if (auto e = std::dynamic_pointer_cast<ElementaryCurve<T>>(iso))
                return elementary_curve(*e);
            return any_bspline_curve<3>(nurbs_export::to_nurbs<T, 3>(*iso, a, b, approx_tol_, stats_));
        }

        std::vector<P21Value> uncertainty_args(T uncertainty) const
        {
            return {V::typed("LENGTH_MEASURE", V::real(uncertainty)), V::ref(length_), V::string("distance_accuracy_value"), V::string("confusion accuracy")};
        }

    public:
        GeometryWriter(P21Writer &w, T approx_tol = T(1e-5)) : w_{w}, approx_tol_{approx_tol} {}

        [[nodiscard]] const nurbs_export::ApproxStats<T> &stats() const noexcept { return stats_; }
        [[nodiscard]] P21Writer &writer() noexcept { return w_; }

        /**
         * @brief Geometric representation context in the model unit (mm per unit), radians,
         * steradians, with the distance uncertainty `uncertainty` (model units).
         */
        std::uint64_t context(T unit_mm, T uncertainty)
        {
            if (!(unit_mm > 0.) || !std::isfinite(unit_mm) || !(uncertainty > 0.))
                throw std::invalid_argument("STEP context: invalid unit or uncertainty");
            auto si = [&](std::optional<std::string_view> prefix, std::string_view name, std::string_view kind) {
                return w_.add_complex({{std::string(kind), {}},
                                       {"NAMED_UNIT", {V::derived()}},
                                       {"SI_UNIT", {prefix ? V::enumeration(*prefix) : V::omitted(), V::enumeration(name)}}});
            };
            std::uint64_t length;
            const auto near = [unit_mm](T x) { return std::abs(unit_mm - x) <= 1e-12 * x; };
            if (near(1.))
                length = si("MILLI", "METRE", "LENGTH_UNIT");
            else if (near(10.))
                length = si("CENTI", "METRE", "LENGTH_UNIT");
            else if (near(1000.))
                length = si(std::nullopt, "METRE", "LENGTH_UNIT");
            else if (near(1e-3))
                length = si("MICRO", "METRE", "LENGTH_UNIT");
            else
            {
                // conversion based unit, defined in millimetres
                const auto mm = si("MILLI", "METRE", "LENGTH_UNIT");
                const auto measure = w_.add("LENGTH_MEASURE_WITH_UNIT", {V::typed("LENGTH_MEASURE", V::real(unit_mm)), V::ref(mm)});
                const auto dims = w_.add("DIMENSIONAL_EXPONENTS", {V::real(1.), V::real(0.), V::real(0.), V::real(0.), V::real(0.), V::real(0.), V::real(0.)});
                const std::string name = near(25.4) ? "INCH" : near(304.8) ? "FOOT" : "MODEL UNIT";
                length = w_.add_complex({{"CONVERSION_BASED_UNIT", {V::string(name), V::ref(measure)}},
                                         {"LENGTH_UNIT", {}},
                                         {"NAMED_UNIT", {V::ref(dims)}}});
            }
            const auto angle = si(std::nullopt, "RADIAN", "PLANE_ANGLE_UNIT");
            const auto solid = si(std::nullopt, "STERADIAN", "SOLID_ANGLE_UNIT");
            length_ = length;
            uncertainty_ = w_.add("UNCERTAINTY_MEASURE_WITH_UNIT", uncertainty_args(uncertainty));
            const auto unc = uncertainty_;
            return w_.add_complex({{"GEOMETRIC_REPRESENTATION_CONTEXT", {V::integer(3)}},
                                   {"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT", {V::refs({unc})}},
                                   {"GLOBAL_UNIT_ASSIGNED_CONTEXT", {V::refs({length, angle, solid})}},
                                   {"REPRESENTATION_CONTEXT", {V::string("Context #1"), V::string("3D Context with UNIT and UNCERTAINTY")}}});
        }

        /// Changes the uncertainty of the context written by context() (model units).
        void set_uncertainty(T uncertainty)
        {
            if (!uncertainty_ || !(uncertainty > 0.) || !std::isfinite(uncertainty))
                throw std::invalid_argument("STEP context: invalid uncertainty");
            w_.replace(uncertainty_, "UNCERTAINTY_MEASURE_WITH_UNIT", uncertainty_args(uncertainty));
        }

        /// Context of the 2D parameter space of pcurves.
        std::uint64_t context2d()
        {
            return w_.add_complex({{"GEOMETRIC_REPRESENTATION_CONTEXT", {V::integer(2)}},
                                   {"PARAMETRIC_REPRESENTATION_CONTEXT", {}},
                                   {"REPRESENTATION_CONTEXT", {V::string("2D SPACE"), V::string("")}}});
        }

        std::uint64_t cartesian_point(const gbs::point<T, 3> &p) { return point_dim<3>(p); }
        std::uint64_t cartesian_point2d(const gbs::point<T, 2> &p) { return point_dim<2>(p); }
        std::uint64_t direction(const gbs::point<T, 3> &d) { return direction_dim<3>(d); }

        /// AXIS2_PLACEMENT_3D of a frame.
        std::uint64_t axis2(const Frame3<T> &f)
        {
            return w_.add("AXIS2_PLACEMENT_3D", {V::string(""), V::ref(cartesian_point(f.o)), V::ref(direction(f.z)), V::ref(direction(f.x))});
        }

        /**
         * @brief 3D curve of an edge on [a, b]: the primitive, the B-spline restricted to [a, b],
         * the exact iso line of a curve on surface, or an approximation.
         */
        std::uint64_t curve(const Curve<T, 3> &c, T a, T b)
        {
            if (auto e = dynamic_cast<const ElementaryCurve<T> *>(&c))
                return elementary_curve(*e);
            if (auto cs = dynamic_cast<const CurveOnSurface<T, 3> *>(&c))
                if (auto iso = iso_line(*cs, a, b))
                    return *iso;
            return any_bspline_curve<3>(nurbs_export::to_nurbs<T, 3>(c, a, b, approx_tol_, stats_));
        }

        /// 2D curve (pcurve) on [a, b]: a line for a degree-1 segment, a B-spline otherwise.
        std::uint64_t curve2d(const Curve<T, 2> &c, T a, T b)
        {
            if (auto l = dynamic_cast<const BSCurve<T, 2> *>(&c); l && l->degree() == 1 && l->poles().size() == 2)
            {
                const auto [k1, k2] = l->bounds();
                const auto v = (1. / (k2 - k1)) * (l->poles()[1] - l->poles()[0]);
                if (norm(v) > 0.)
                    return line_dim<2>(l->poles()[0] - k1 * v, v); // C(t) = P0 + (t - k1) V
            }
            return any_bspline_curve<2>(nurbs_export::to_nurbs<T, 2>(c, a, b, approx_tol_, stats_));
        }

        /// Surface of a face: the primitive (unbounded STEP entity), the B-spline, or an approximation.
        std::uint64_t surface(const Surface<T, 3> &s)
        {
            if (auto p = dynamic_cast<const Plane<T> *>(&s))
                return w_.add("PLANE", {V::string(""), V::ref(axis2(p->frame()))});
            if (auto c = dynamic_cast<const CylindricalSurface<T> *>(&s))
                return w_.add("CYLINDRICAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(c->radius())});
            if (auto c = dynamic_cast<const ConicalSurface<T> *>(&s))
                return w_.add("CONICAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(c->radius()), V::real(c->semi_angle())});
            if (auto c = dynamic_cast<const SphericalSurface<T> *>(&s))
                return w_.add("SPHERICAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(c->radius())});
            if (auto c = dynamic_cast<const ToroidalSurface<T> *>(&s))
            {
                const T R = c->major_radius(), r = c->minor_radius();
                if (R >= r)
                    return w_.add("TOROIDAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(R), V::real(r)});
                // self-intersecting (r > R): select_outer chooses the apple, v in [-t, t], or the lemon,
                // v in [t, 2 pi - t] (cos t = -R / r), as the v range of the face's surface lies; a range
                // over both is not one STEP surface: approximated
                const auto b = c->bounds();
                const T t = std::acos(-R / r), mid = (b[2] + b[3]) / 2, eps = T(1e-9) * (1 + std::abs(mid));
                const T sa = 2 * std::numbers::pi * std::round(mid / (2 * std::numbers::pi)), sl = 2 * std::numbers::pi * std::floor(mid / (2 * std::numbers::pi));
                const bool apple = b[2] - sa >= -t - eps && b[3] - sa <= t + eps;
                const bool lemon = b[2] - sl >= t - eps && b[3] - sl <= 2 * std::numbers::pi - t + eps;
                if (apple || lemon)
                    return w_.add("DEGENERATE_TOROIDAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(R), V::real(r), V::boolean(apple)});
            }
            return std::visit([this](const auto &b) { return bspline_surface(b); }, nurbs_export::to_nurbs<T>(s, approx_tol_, stats_));
        }
    };
}
