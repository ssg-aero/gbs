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
 * OCCT round trip of gbs-occt/tests checks it); any other curve or surface is converted
 * to a cubic NURBS by interpolation at its own parameters (the pcurves stay
 * valid), refined until `approx_tol`.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <gbs-io/iges.h>
#include <gbs-brep/brep>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/bsctools.h>

namespace gbs
{
    template <std::floating_point T>
    struct IgesExportOptions
    {
        T approx_tol = T(1e-5); ///< max deviation of a non-NURBS curve or surface converted to NURBS
        T scale = T(1);         ///< factor applied to model-space coordinates (not to parameters)
    };

    /// What the export wrote and how much it had to approximate.
    template <std::floating_point T>
    struct IgesExportReport
    {
        std::size_t faces{}, wires{}, free_edges{};
        std::size_t approximated_curves{}, approximated_surfaces{};
        T max_deviation{}; ///< largest deviation of an approximated curve or surface
    };

    namespace iges_detail
    {
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
        auto to_nurbs(const Curve<T, dim> &c, T a, T b, T tol, IgesExportReport<T> &rep) -> AnyBSCurve<T, dim>
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
        auto to_nurbs(const Surface<T, 3> &s, T tol, IgesExportReport<T> &rep) -> AnyBSSurface<T>
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

        // set_126 / set_128 (NURBS writers, Cartesian poles + weights) live in gbs-io/iges.h.

        /// 126 entity for a curve on [a, b], reversed if requested.
        template <std::floating_point T, std::size_t dim>
        auto curve_126(DLL_IGES &model, const Curve<T, dim> &c, T a, T b, bool reversed, T tol, T scale, IgesExportReport<T> &rep)
            -> std::unique_ptr<DLL_IGES_ENTITY_126>
        {
            auto e = std::make_unique<DLL_IGES_ENTITY_126>(model, true);
            auto nurbs = to_nurbs(c, a, b, tol, rep);
            std::visit([&](auto &bs) {
                if (reversed)
                    bs.reverse();
                set_126(*e, bs, scale);
            },
                       nurbs);
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

        /// Entity 142 of a face wire: composite curves of the pcurves and of the 3D curves, co-edge order and sense.
        template <std::floating_point T>
        auto wire_142(DLL_IGES &model, const brep::Model<T> &m, brep::WireId wid, DLL_IGES_ENTITY_128 &srf,
                      const IgesExportOptions<T> &opts, IgesExportReport<T> &rep) -> std::unique_ptr<DLL_IGES_ENTITY_142>
        {
            DLL_IGES_ENTITY_102 bptr(model, true), cptr(model, true);
            for (const auto &ce : m.wire(wid).coedges)
            {
                const auto &e = m.edge(ce.edge);
                const bool rev = ce.orient == brep::Orientation::Reversed;
                auto p2 = curve_126<T, 2>(model, *ce.pcurve, e.u1, e.u2, rev, opts.approx_tol, T(1), rep);
                bptr.AddSegment(*p2);
                auto p3 = e.degenerate || !e.curve ? point_126(model, m.vertex(e.v1).pnt, e.u1, e.u2, opts.scale)
                                                   : curve_126<T, 3>(model, *e.curve, e.u1, e.u2, rev, opts.approx_tol, opts.scale, rep);
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
    } // namespace iges_detail

    /**
     * @brief Writes the faces (and the free edges) of `shape` into `model`.
     * Each face is labelled with its own name (Model::name) when it has one, otherwise
     * with `name` followed by its rank when `name` is not empty.
     */
    template <std::floating_point T>
    auto add_brep(DLL_IGES &model, const brep::Model<T> &m, const brep::ShapeId &shape, const std::string &name = "",
                  IgesExportOptions<T> opts = {}) -> IgesExportReport<T>
    {
        using namespace brep;
        IgesExportReport<T> rep;
        std::vector<EdgeId> face_edges;
        std::size_t rank = 0;
        for (auto fid : explore<FaceId>(m, shape))
        {
            const auto &f = m.face(fid);
            if (!f.surface || f.wires.empty())
                continue;
            DLL_IGES_ENTITY_128 srf(model, true);
            std::visit([&](const auto &bs) { iges_detail::set_128(srf, bs, opts.scale); }, iges_detail::to_nurbs(*f.surface, opts.approx_tol, rep));

            DLL_IGES_ENTITY_144 tf(model, true);
            tf.SetSurface(srf);
            for (std::size_t i{}; i < f.wires.size(); ++i)
            {
                auto c = iges_detail::wire_142(model, m, f.wires[i], srf, opts, rep);
                if (i == 0)
                    tf.SetBoundCurve(*c);
                else
                    tf.AddCutout(*c);
                for (const auto &ce : m.wire(f.wires[i]).coedges)
                    face_edges.push_back(ce.edge);
            }
            ++rank;
            if (const auto own = m.name(fid); !own.empty())
                tf.SetLabel(std::string(own).c_str()); // the face's own name (e.g. read from STEP)
            else if (!name.empty())
                tf.SetLabel((name + std::to_string(rank)).c_str());
            ++rep.faces;
        }
        for (auto eid : explore<EdgeId>(m, shape))
        {
            if (std::ranges::find(face_edges, eid) != face_edges.end())
                continue;
            const auto &e = m.edge(eid);
            if (e.degenerate || !e.curve)
                continue;
            auto c = iges_detail::curve_126<T, 3>(model, *e.curve, e.u1, e.u2, false, opts.approx_tol, opts.scale, rep);
            if (!name.empty())
                c->SetLabel(name.c_str());
            ++rep.free_edges;
        }
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
