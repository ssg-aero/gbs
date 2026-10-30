#pragma once

/**
 * @file iges_brep.h
 * @brief IGES export of native BREP shapes (gbs-brep) through libIGES.
 *
 * Design: docs/sources/design/brep_core.md, section 7; architecture note
 * docs/sources/design/brep_pr09_iges.md.
 *
 * Every face becomes a trimmed parametric surface (entity 144) over a NURBS
 * surface (128); every wire of the face becomes a curve on parametric surface
 * (142) whose parameter-space and model-space representations are composite
 * curves (102) of NURBS curves (126), one segment per co-edge, in the co-edge
 * order and sense. Edges of the shape not bounding an exported face are
 * written as 126 curves. Shells, solids and compounds are written as their
 * faces: the BREP entities 186/514/510/508/504/502 are not exposed by the
 * libIGES DLL API (design question 10).
 *
 * NURBS geometry is written exactly (rational poles converted from the gbs
 * homogeneous form to the Cartesian points and weights libIGES expects; the
 * OCCT round trip of gbs-occt/tests checks it). Primitives too (gbs/primitives.h):
 * lines, circles and ellipses as their exact NURBS, planes, cylinders, cones,
 * spheres and tori as their exact NURBS, on which the pcurves are carried to the
 * NURBS parameters (`ElementarySurface::nurbs_parameters`) within `approx_tol`.
 * Any other curve or surface is converted to a cubic NURBS by interpolation at its
 * own parameters (the pcurves stay valid), refined until `approx_tol`.
 *
 * Coordinates are written in the unit of the IGES file (its units flag, millimetre
 * by default; `IgesWriter::model().SetUnitsFlag()` changes it), converted from the
 * unit of the model (`Model::unitScale`); the minimum resolution of the file is
 * raised to the tolerances of the shape.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <unordered_set>
#include <variant>
#include <vector>

#include <gbs-io/iges.h>
#include <gbs-io/to_nurbs.h>
#include <gbs-brep/brep>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/bsctools.h>
#include <gbs/bselementary.h>
#include <gbs/primitives.h>

namespace gbs
{
    template <std::floating_point T>
    struct IgesExportOptions
    {
        T approx_tol = T(1e-5); ///< max deviation of a non-NURBS curve or surface converted to NURBS
        T scale = T(1);         ///< factor applied to model-space coordinates, after the conversion to the unit of the file (not to parameters)
    };

    /// What the export wrote and how much it had to approximate (approximated_curves,
    /// approximated_surfaces, max_deviation: see nurbs_export::ApproxStats).
    template <std::floating_point T>
    struct IgesExportReport : nurbs_export::ApproxStats<T>
    {
        std::size_t faces{}, wires{}, free_edges{};
    };

    namespace iges_detail
    {
        // to_nurbs (exact copy of a B-spline, interpolation otherwise) lives in gbs-io/to_nurbs.h.
        using nurbs_export::AnyBSCurve;
        using nurbs_export::AnyBSSurface;
        using nurbs_export::to_nurbs;

        // set_126 / set_128 (NURBS writers, Cartesian poles + weights) live in gbs-io/iges.h.

        /// Millimetres per unit of an IGES units flag (0 for a unit given by name only).
        inline double unit_mm(IGES_UNIT u)
        {
            switch (u)
            {
            case UNIT_INCH: return 25.4;
            case UNIT_MILLIMETER: return 1.;
            case UNIT_FOOT: return 304.8;
            case UNIT_MILE: return 1609344.;
            case UNIT_METER: return 1000.;
            case UNIT_KILOMETER: return 1e6;
            case UNIT_MIL: return 0.0254;
            case UNIT_MICRON: return 1e-3;
            case UNIT_CENTIMETER: return 10.;
            case UNIT_MICROINCH: return 2.54e-5;
            default: return 0.;
            }
        }

        /// Exact NURBS of a line, a circle, an ellipse or an iso line of a surface on [a, b]; nullopt for another curve.
        template <std::floating_point T>
        auto exact_nurbs(const Curve<T, 3> &c, T a, T b) -> std::optional<AnyBSCurve<T, 3>>
        {
            if (auto l = dynamic_cast<const LineCurve<T> *>(&c))
                return BSCurve<T, 3>{points_vector<T, 3>{l->value(a), l->value(b)}, std::vector<T>{a, a, b, b}, 1};
            if (auto e = dynamic_cast<const EllipseCurve<T> *>(&c); e && b - a <= 2 * std::numbers::pi_v<T> * (1 + T(1e-12)))
            {
                const auto &f = e->frame();
                return build_ellipse_arc<T>(e->radius1(), e->radius2(), a, std::min(b, a + 2 * std::numbers::pi_v<T>), ax2<T, 3>{f.o, f.z, f.x});
            }
            // an iso line of a curve on a surface (edges of natural faces): the primitive, or the B-spline iso
            if (auto cs = dynamic_cast<const CurveOnSurface<T, 3> *>(&c))
                if (const auto iso = nurbs_export::iso_curve(*cs, a, b))
                {
                    if (dynamic_cast<const ElementaryCurve<T> *>(iso.get()))
                        return exact_nurbs(*iso, a, b);
                    nurbs_export::ApproxStats<T> exact; // a B-spline is copied: no approximation
                    return to_nurbs(*iso, a, b, T(0), exact);
                }
            return std::nullopt;
        }

        /// 126 entity for a curve on [a, b], reversed if requested: a primitive exactly (3D), otherwise to_nurbs.
        template <std::floating_point T, std::size_t dim>
        auto curve_126(DLL_IGES &model, const Curve<T, dim> &c, T a, T b, bool reversed, T tol, T scale, IgesExportReport<T> &rep)
            -> std::unique_ptr<DLL_IGES_ENTITY_126>
        {
            auto e = std::make_unique<DLL_IGES_ENTITY_126>(model, true);
            std::optional<AnyBSCurve<T, dim>> nurbs;
            if constexpr (dim == 3)
                nurbs = exact_nurbs(c, a, b);
            if (!nurbs)
                nurbs = to_nurbs(c, a, b, tol, rep);
            std::visit([&](auto &bs) {
                if (reversed)
                    bs.reverse();
                set_126(*e, bs, scale);
            },
                       *nurbs);
            return e;
        }

        /**
         * 126 entity for the pcurve of a co-edge on a primitive written as its exact NURBS `nurbs`:
         * the pcurve carried to the NURBS parameters, refined until the NURBS point is within tol of
         * the point of the primitive.
         */
        template <std::floating_point T>
        auto pcurve_126(DLL_IGES &model, const Curve<T, 2> &pc, const ElementarySurface<T> &es, const Surface<T, 3> &nurbs, T a, T b, bool reversed, T tol,
                        IgesExportReport<T> &rep) -> std::unique_ptr<DLL_IGES_ENTITY_126>
        {
            const auto [s1, s2, w1, w2] = nurbs.bounds();
            auto f = [&](T t) {
                const auto uv = pc.value(t);
                const auto w = es.nurbs_parameters(uv[0], uv[1]);
                return point<T, 2>{w[0], w[1]};
            };
            auto gap = [&](const BSCurve<T, 2> &bs, T t) {
                const auto w = bs.value(t);
                const auto uv = pc.value(t);
                return distance(nurbs.value(std::clamp(w[0], s1, s2), std::clamp(w[1], w1, w2)), es.value(uv[0], uv[1]));
            };
            auto bs = nurbs_export::interpolate_within<T, 2>(f, a, b, gap, tol, rep);
            if (reversed)
                bs.reverse();
            auto e = std::make_unique<DLL_IGES_ENTITY_126>(model, true);
            set_126(*e, bs, T(1));
            return e;
        }

        /// 126 entity for a degenerate edge: a degree-1 curve collapsed on the vertex.
        template <std::floating_point T>
        auto point_126(DLL_IGES &model, const point<T, 3> &p, T a, T b, T scale) -> std::unique_ptr<DLL_IGES_ENTITY_126>
        {
            auto e = std::make_unique<DLL_IGES_ENTITY_126>(model, true);
            set_126(*e, BSCurve<T, 3>{points_vector<T, 3>{p, p}, std::vector<T>{a, a, b, b}, 1}, scale);
            return e;
        }

        /**
         * Entity 142 of a face wire: composite curves of the pcurves and of the 3D curves, co-edge order and
         * sense. On a primitive written as its exact NURBS `exact` (`es` the primitive), the pcurves are
         * carried to the NURBS parameters; `scale` converts the coordinates.
         */
        template <std::floating_point T>
        auto wire_142(DLL_IGES &model, const brep::Model<T> &m, brep::WireId wid, DLL_IGES_ENTITY_128 &srf, const ElementarySurface<T> *es,
                      const Surface<T, 3> *exact, T scale, const IgesExportOptions<T> &opts, IgesExportReport<T> &rep) -> std::unique_ptr<DLL_IGES_ENTITY_142>
        {
            DLL_IGES_ENTITY_102 bptr(model, true), cptr(model, true);
            for (const auto &ce : m.wire(wid).coedges)
            {
                const auto &e = m.edge(ce.edge);
                const bool rev = ce.orient == brep::Orientation::Reversed;
                auto p2 = es && exact ? pcurve_126<T>(model, *ce.pcurve, *es, *exact, e.u1, e.u2, rev, opts.approx_tol, rep)
                                      : curve_126<T, 2>(model, *ce.pcurve, e.u1, e.u2, rev, opts.approx_tol, T(1), rep);
                bptr.AddSegment(*p2);
                auto p3 = e.degenerate || !e.curve ? point_126(model, m.vertex(e.v1).pnt, e.u1, e.u2, scale)
                                                   : curve_126<T, 3>(model, *e.curve, e.u1, e.u2, rev, opts.approx_tol, scale, rep);
                cptr.AddSegment(*p3);
            }
            auto c = std::make_unique<DLL_IGES_ENTITY_142>(model, true);
            c->SetSurface(srf);
            c->SetParameterSpaceBound(bptr);
            c->SetModelSpaceBound(cptr);
            c->SetCurvePreference(BOUND_PREF_ANY);
            c->SetCurveCreationFlag(CURVE_CREATE_UNSPECIFIED);
            ++rep.wires;
            return c;
        }

        /// Names an entity: its label (8 columns of its directory entry, which libIGES truncates to) and, for a
        /// longer name, the whole of it as a Name property (entity 406, form 15).
        inline void set_name(DLL_IGES &model, DLL_IGES_ENTITY &e, const std::string &name)
        {
            e.SetLabel(name.c_str());
            if (name.size() <= 8)
                return;
            DLL_IGES_ENTITY_406 prop(model, true);
            if (prop.SetEntityForm(15) && prop.SetProperty_Name(name.c_str()))
                e.AddOptionalEntity(prop);
        }
    } // namespace iges_detail

    /**
     * @brief Writes the faces (and the free edges) of `shape` into `model`.
     * Each face is named with its own name (Model::name) when it has one, otherwise
     * with `name` followed by its rank when `name` is not empty; a name longer than the 8 columns of an IGES
     * label is also written whole as a Name property (entity 406, form 15).
     */
    template <std::floating_point T>
    auto add_brep(DLL_IGES &model, const brep::Model<T> &m, const brep::ShapeId &shape, const std::string &name = "",
                  IgesExportOptions<T> opts = {}) -> IgesExportReport<T>
    {
        using namespace brep;
        IgesExportReport<T> rep;
        // model units to the units of the file
        IGES_UNIT unit{UNIT_MILLIMETER};
        model.GetUnitsFlag(unit);
        const double unit_mm = iges_detail::unit_mm(unit);
        const T scale = opts.scale * (unit_mm > 0. ? T(m.unitScale() / unit_mm) : T(1));
        std::unordered_set<EdgeId> face_edges;
        T max_tol{};
        std::size_t rank = 0;
        for (auto fid : explore<FaceId>(m, shape))
        {
            const auto &f = m.face(fid);
            if (!f.surface || f.wires.empty())
                continue;
            max_tol = std::max(max_tol, f.tol);
            // a primitive as its exact NURBS; a plane keeps its parameters, the others get their pcurves carried
            const auto *es = dynamic_cast<const ElementarySurface<T> *>(f.surface.get());
            const auto exact = es ? es->to_nurbs() : nullptr;
            const auto *mapped = es && !dynamic_cast<const Plane<T> *>(es) ? es : nullptr;
            DLL_IGES_ENTITY_128 srf(model, true);
            if (auto r = std::dynamic_pointer_cast<BSSurfaceRational<T, 3>>(exact))
                iges_detail::set_128(srf, *r, scale);
            else if (auto b = std::dynamic_pointer_cast<BSSurface<T, 3>>(exact))
                iges_detail::set_128(srf, *b, scale);
            else
                std::visit([&](const auto &bs) { iges_detail::set_128(srf, bs, scale); }, iges_detail::to_nurbs(*f.surface, opts.approx_tol, rep));

            DLL_IGES_ENTITY_144 tf(model, true);
            tf.SetSurface(srf);
            for (std::size_t i{}; i < f.wires.size(); ++i)
            {
                auto c = iges_detail::wire_142(model, m, f.wires[i], srf, mapped, mapped ? exact.get() : nullptr, scale, opts, rep);
                if (i == 0)
                    tf.SetBoundCurve(*c);
                else
                    tf.AddCutout(*c);
                for (const auto &ce : m.wire(f.wires[i]).coedges)
                {
                    face_edges.insert(ce.edge);
                    max_tol = std::max({max_tol, m.edge(ce.edge).tol, m.vertex(m.edge(ce.edge).v1).tol});
                }
            }
            ++rank;
            if (const auto own = m.name(fid); !own.empty())
                iges_detail::set_name(model, tf, own); // the face's own name (e.g. read from STEP)
            else if (!name.empty())
                iges_detail::set_name(model, tf, name + std::to_string(rank));
            ++rep.faces;
        }
        for (auto eid : explore<EdgeId>(m, shape))
        {
            if (face_edges.contains(eid))
                continue;
            const auto &e = m.edge(eid);
            if (e.degenerate || !e.curve)
                continue;
            max_tol = std::max(max_tol, e.tol);
            auto c = iges_detail::curve_126<T, 3>(model, *e.curve, e.u1, e.u2, false, opts.approx_tol, scale, rep);
            if (!name.empty())
                iges_detail::set_name(model, *c, name);
            ++rep.free_edges;
        }
        // the resolution of the file covers the tolerances of the shape and the approximations
        if (double res{}; model.GetMinResolution(res))
            model.SetMinResolution(std::max(res, double(std::max(max_tol, rep.max_deviation) * scale)));
        return rep;
    }

    /// Same, through an `IgesWriter`.
    template <std::floating_point T>
    auto add_brep(IgesWriter<T> &w, const brep::Model<T> &m, const brep::ShapeId &shape, const std::string &name = "",
                  IgesExportOptions<T> opts = {}) -> IgesExportReport<T>
    {
        return add_brep(w.model(), m, shape, name, opts);
    }

} // namespace gbs
