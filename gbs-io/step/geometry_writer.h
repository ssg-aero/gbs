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
 *   the same parametrization: pcurves stay valid as they are. A cone of
 *   negative semi-angle, which ISO 10303-42 does not allow, is written as the
 *   same cone of axis -Z and semi-angle |a|, parameters (-u, -v): its pcurves,
 *   and those of an offset of it, are written negated (flipped()).
 * - B-splines become b_spline_curve_with_knots / b_spline_surface_with_knots
 *   (complex instances when rational, Cartesian poles and weights), a curve
 *   restricted to the requested range; 2D degree-1 segments become lines.
 * - An offset surface (gbs/offsets.h) whose basis is written so and is C1 by
 *   its knots becomes an offset_surface of it (S + d N on the parameters of S,
 *   as in gbs); on a basis with a knot repeated degree times (the NURBS of a
 *   circle), it is approximated: OCCT refuses such a basis.
 * - A curve on a surface along an iso line of a primitive or of a B-spline
 *   becomes that iso line, exactly.
 * - Anything else is approximated within approx_tol, parameter preserved
 *   (gbs-io/to_nurbs.h), and counted.
 * - With a knot tolerance, the B-splines written (3D curves, surfaces, and the
 *   pcurves given their surface) lose their spans shorter than it
 *   (gbs-io/knot_cleanup.h).
 *
 * Coordinates are written as they are, in the model unit described by the
 * representation context (context()).
 */

#include <gbs-io/step/p21_writer.h>
#include <gbs-io/to_nurbs.h>
#include <gbs-io/knot_cleanup.h>
#include <gbs/primitives.h>
#include <gbs/offsets.h>
#include <gbs/curveonsurface.h>

#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace gbs::step
{
    class GeometryWriter
    {
        using T = double;
        using V = P21Value;

        P21Writer &w_;
        T approx_tol_;
        T knot_tol_;
        nurbs_export::ApproxStats<T> stats_;
        nurbs_export::KnotCleanupStats<T> knot_stats_;
        std::uint64_t uncertainty_{}, length_{}; // uncertainty measure and length unit of the context

        // Distinct knots and multiplicities, by exact equality: two knots an ulp apart stay two knots (grouped
        // within knot_eps, as unflat_knots does, they made a multiplicity above degree + 1: an invalid entity).
        static void knots_and_mults(const std::vector<T> &flat, std::vector<std::int64_t> &mults, std::vector<T> &knots)
        {
            knots.clear();
            mults.clear();
            for (const auto same : flat | std::views::chunk_by(std::ranges::equal_to{}))
            {
                knots.push_back(same.front());
                mults.push_back(std::ranges::distance(same));
            }
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

        // B-spline curve (homogeneous poles when rational) as STEP entities; sign: -1 writes its
        // poles negated (a pcurve on a flipped cone).
        template <std::size_t dim, bool rational>
        std::uint64_t bspline_curve(const BSCurveGeneral<T, dim, rational> &c, T sign = T(1))
        {
            std::vector<std::uint64_t> cps;
            std::vector<T> weights;
            for (const auto &h : c.poles())
            {
                point<T, dim> p;
                const T wgt = rational ? h[dim] : T(1);
                for (std::size_t i = 0; i < dim; ++i)
                    p[i] = sign * h[i] / wgt; // Cartesian poles, the weights apart
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

        /// on: the surface of a pcurve (dim 2), along which its spans are measured; none, the pcurve as it is.
        /// A pcurve on a flipped surface (flipped()) is written negated, after its spans are cleaned along `on`.
        template <std::size_t dim>
        std::uint64_t any_bspline_curve(nurbs_export::AnyBSCurve<T, dim> c, const Surface<T, 3> *on = nullptr)
        {
            return std::visit([&](auto &b) {
                if constexpr (dim == 3)
                    nurbs_export::clean_knots(b, knot_tol_, knot_stats_);
                else if (on)
                    nurbs_export::clean_knots(b, *on, knot_tol_, knot_stats_);
                return bspline_curve(b, dim == 2 && on && flipped(*on) ? T(-1) : T(1));
            },
                              c);
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
        /// knot_tol: length below which the spans of the B-splines written are removed (0: none)
        GeometryWriter(P21Writer &w, T approx_tol = T(1e-5), T knot_tol = T(0)) : w_{w}, approx_tol_{approx_tol}, knot_tol_{knot_tol} {}

        [[nodiscard]] const nurbs_export::ApproxStats<T> &stats() const noexcept { return stats_; }
        [[nodiscard]] const nurbs_export::KnotCleanupStats<T> &knot_stats() const noexcept { return knot_stats_; }
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

        /// 2D curve (pcurve) on [a, b]: a line for a degree-1 segment, a B-spline otherwise, cleaned along
        /// its surface `on` when given; on a surface written flipped (a cone of negative semi-angle, or an
        /// offset of one), mapped by (u, v) -> (-u, -v), the parameters of the STEP surface.
        std::uint64_t curve2d(const Curve<T, 2> &c, T a, T b, const Surface<T, 3> *on = nullptr)
        {
            if (auto l = dynamic_cast<const BSCurve<T, 2> *>(&c); l && l->degree() == 1 && l->poles().size() == 2)
            {
                const auto [k1, k2] = l->bounds();
                const T sign = on && flipped(*on) ? T(-1) : T(1);
                const auto v = (sign / (k2 - k1)) * (l->poles()[1] - l->poles()[0]);
                if (norm(v) > 0.)
                    return line_dim<2>(sign * l->poles()[0] - k1 * v, v); // C(t) = P0 + (t - k1) V, negated when flipped
            }
            return any_bspline_curve<2>(nurbs_export::to_nurbs<T, 2>(c, a, b, approx_tol_, stats_), on);
        }

        /// Entities of the surfaces written, by surface (those of a part: its faces and the bases of their offsets).
        using SurfaceIds = std::unordered_map<const Surface<T, 3> *, std::uint64_t>;

        /**
         * @brief Surface of a face: the primitive (unbounded STEP entity), the B-spline, the
         * OFFSET_SURFACE of an offset surface whose basis is one of them (or an offset of one) and
         * is C1 by its knots, or an approximation.
         *
         * With `ids`, a surface already written there, a face or the basis of an offset, is not
         * written again; the surfaces must outlive `ids`.
         */
        std::uint64_t surface(const Surface<T, 3> &s, SurfaceIds *ids = nullptr)
        {
            if (ids)
                if (auto it = ids->find(&s); it != ids->end())
                    return it->second;
            std::uint64_t id;
            if (auto e = exact_surface(s, ids))
                id = *e;
            else
                id = std::visit([this](auto b) {
                    nurbs_export::clean_knots(b, knot_tol_, knot_stats_);
                    return bspline_surface(b);
                },
                                nurbs_export::to_nurbs<T>(s, approx_tol_, stats_));
            if (ids)
                ids->emplace(&s, id);
            return id;
        }

        /**
         * @brief True for a surface written with parameters (-u, -v): a cone of negative semi-angle,
         * which ISO 10303-42 does not allow (positive_plane_angle_measure), written as the same cone
         * of axis -Z, reference direction X and semi-angle |a| (Y' = -Y), and an offset written exactly
         * of a flipped surface (its normal S_u x S_v is unchanged). Its pcurves are written negated.
         */
        static bool flipped(const Surface<T, 3> &s)
        {
            if (auto c = dynamic_cast<const ConicalSurface<T> *>(&s))
                return c->semi_angle() < T(0);
            if (auto o = dynamic_cast<const OffsetSurface<T> *>(&s))
                return exact(s) && flipped(*o->basis());
            return false;
        }

    private:
        /// True for a surface that exact_surface() writes.
        static bool exact(const Surface<T, 3> &s)
        {
            if (dynamic_cast<const Plane<T> *>(&s) || dynamic_cast<const CylindricalSurface<T> *>(&s) || dynamic_cast<const ConicalSurface<T> *>(&s) ||
                dynamic_cast<const SphericalSurface<T> *>(&s) || dynamic_cast<const BSSurface<T, 3> *>(&s) || dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                return true;
            if (auto c = dynamic_cast<const ToroidalSurface<T> *>(&s))
                return c->major_radius() >= c->minor_radius() || torus_part(*c).has_value();
            if (auto o = dynamic_cast<const OffsetSurface<T> *>(&s))
                return exact(*o->basis()) && c1(*o->basis()) && !lemon(*o->basis());
            return false;
        }

        // The inner part of a self-intersecting torus (select_outer .F.): its offset is not read (by
        // gbs: the normal turns towards the inside of the tube there), it is approximated.
        static bool lemon(const Surface<T, 3> &s)
        {
            const auto *c = dynamic_cast<const ToroidalSurface<T> *>(&s);
            return c && c->major_radius() < c->minor_radius() && torus_part(*c) == false;
        }

        // C1 by its knots: no inner knot of a B-spline repeated degree times or more (the NURBS of a
        // circle, with double knots, is G1 but not C1). The basis of an offset written as such: the
        // offset of a surface that is not C1 is not defined at its knots, and OCCT refuses it
        // (Geom_OffsetSurface: "Offset with no C1 Surface", when removing the knots fails).
        static bool c1(const Surface<T, 3> &s)
        {
            auto knots_c1 = [](const std::vector<T> &flat, std::size_t p) {
                if (flat.size() < 2 * (p + 1))
                    return true;
                const auto inner = std::span{flat}.subspan(p + 1, flat.size() - 2 * (p + 1));
                for (const auto same : inner | std::views::chunk_by(std::ranges::equal_to{}))
                    if (static_cast<std::size_t>(std::ranges::distance(same)) >= p)
                        return false;
                return true;
            };
            auto bspline = [&](const auto &b) { return knots_c1(b.knotsFlatsU(), b.degreeU()) && knots_c1(b.knotsFlatsV(), b.degreeV()); };
            if (auto b = dynamic_cast<const BSSurface<T, 3> *>(&s))
                return bspline(*b);
            if (auto b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                return bspline(*b);
            if (auto o = dynamic_cast<const OffsetSurface<T> *>(&s))
                return c1(*o->basis());
            return true; // a primitive
        }

        // A self-intersecting torus (r > R): select_outer chooses the apple, v in [-t, t], or the lemon,
        // v in [t, 2 pi - t] (cos t = -R / r), as the v range of the face's surface lies (true: the
        // apple); a range over both is not one STEP surface.
        static std::optional<bool> torus_part(const ToroidalSurface<T> &c)
        {
            const T R = c.major_radius(), r = c.minor_radius();
            const auto b = c.bounds();
            const T t = std::acos(-R / r), mid = (b[2] + b[3]) / 2, eps = T(1e-9) * (1 + std::abs(mid));
            const T sa = 2 * std::numbers::pi * std::round(mid / (2 * std::numbers::pi)), sl = 2 * std::numbers::pi * std::floor(mid / (2 * std::numbers::pi));
            if (b[2] - sa >= -t - eps && b[3] - sa <= t + eps)
                return true;
            if (b[2] - sl >= t - eps && b[3] - sl <= 2 * std::numbers::pi - t + eps)
                return false;
            return std::nullopt;
        }

        /// The STEP entity of a surface written as it is, with its parameters; nullopt, and nothing
        /// written, for another one.
        std::optional<std::uint64_t> exact_surface(const Surface<T, 3> &s, SurfaceIds *ids)
        {
            if (!exact(s))
                return std::nullopt;
            if (auto o = dynamic_cast<const OffsetSurface<T> *>(&s))
            {
                // S + d N, N = S_u x S_v / |S_u x S_v|: the normal of ISO 10303-42, on the parameters
                // of S; S written once per part (an offset face and a face on its basis)
                const auto basis = surface(*o->basis(), ids);
                return w_.add("OFFSET_SURFACE", {V::string(""), V::ref(basis), V::real(o->distance()), V::logical(false)});
            }
            if (auto p = dynamic_cast<const Plane<T> *>(&s))
                return w_.add("PLANE", {V::string(""), V::ref(axis2(p->frame()))});
            if (auto c = dynamic_cast<const CylindricalSurface<T> *>(&s))
                return w_.add("CYLINDRICAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(c->radius())});
            if (auto c = dynamic_cast<const ConicalSurface<T> *>(&s))
            {
                // S = O + (R + v tan a)(cos u X + sin u Y) + v Z; a < 0: the cone of axis -Z, reference
                // direction X, semi-angle -a, of parameters u' = -u, v' = -v (see flipped())
                auto f = c->frame();
                if (c->semi_angle() < T(0))
                    f.z = T(-1) * f.z, f.y = T(-1) * f.y;
                return w_.add("CONICAL_SURFACE", {V::string(""), V::ref(axis2(f)), V::real(c->radius()), V::real(std::abs(c->semi_angle()))});
            }
            if (auto c = dynamic_cast<const SphericalSurface<T> *>(&s))
                return w_.add("SPHERICAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(c->radius())});
            if (auto c = dynamic_cast<const ToroidalSurface<T> *>(&s))
            {
                const T R = c->major_radius(), r = c->minor_radius();
                if (R >= r)
                    return w_.add("TOROIDAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(R), V::real(r)});
                return w_.add("DEGENERATE_TOROIDAL_SURFACE", {V::string(""), V::ref(axis2(c->frame())), V::real(R), V::real(r), V::boolean(*torus_part(*c))});
            }
            // a B-spline, copied
            return std::visit([this](auto b) {
                nurbs_export::clean_knots(b, knot_tol_, knot_stats_);
                return bspline_surface(b);
            },
                              nurbs_export::to_nurbs<T>(s, approx_tol_, stats_));
        }
    };
}
