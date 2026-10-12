#pragma once

/**
 * @file edit.h
 * @brief Native BREP core: local edits of an existing model (first operation
 * of stage 2, needed by the STEP reader to insert seams).
 *
 * Architecture note: docs/sources/design/step_pr07_seam.md.
 */

#include <gbs-brep/builders.h>
#include <gbs-brep/placement.h>
#include <gbs/transform.h>
#include <gbs/primitives.h>

#include <array>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gbs::brep
{
    /// Result of split_edge.
    struct SplitEdge
    {
        EdgeId first;   ///< the original edge, now from v1 to the new vertex
        EdgeId second;  ///< new edge from the new vertex to the former v2
        VertexId vertex; ///< new vertex at the split point
    };

    /**
     * @brief Splits the edge `e` at its parameter `u` (u1 < u < u2).
     *
     * The edge keeps its id and becomes [u1, u]; a new edge covers [u, u2];
     * both share the 3D curve and meet at a new vertex (tolerance of the edge).
     * Every wire of the model using `e` is updated: a Forward co-edge becomes
     * (first, second), a Reversed one (second, first), both keeping the
     * co-edge's pcurve, which is parametrized like the edge. Name and external
     * id are copied to the new edge. A closed edge becomes two edges between
     * its vertex and the new one.
     * On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split_edge(Model<T> &m, EdgeId e, std::type_identity_t<T> u) -> BuildResult<SplitEdge>
    {
        if (!m.alive(e))
            return detail::fail(BuildErrc::InvalidId, "split_edge: dead or invalid edge", {e});
        const auto ed = m.edge(e);
        if (ed.degenerate || !ed.curve)
            return detail::fail(BuildErrc::DegenerateEdge, "split_edge: degenerate edge", {e});
        const T eps = T(1e-9) * (ed.u2 - ed.u1);
        if (!std::isfinite(u) || !(u > ed.u1 + eps) || !(u < ed.u2 - eps))
            return detail::fail(BuildErrc::InvalidBounds, "split_edge: parameter not strictly inside the edge", {e});
        const auto p = ed.curve->value(u);
        if (!detail::finite(p))
            return detail::fail(BuildErrc::InvalidPoint, "split_edge: curve evaluates to a non finite point", {e});

        // ---- write
        const auto v = m.add(Vertex<T>{p, ed.tol});
        auto &first = m.edge(e);
        first.u2 = u;
        first.v2 = v;
        const auto second = m.add(Edge<T>{ed.curve, u, ed.u2, v, ed.v2, ed.tol});
        if (auto n = m.name(e); !n.empty())
            m.setName(second, std::string(n));
        m.setExternalId(second, m.externalId(e));

        for (auto wid : m.template ids<WireId>())
        {
            auto &cos = m.wire(wid).coedges;
            for (std::size_t i = 0; i < cos.size(); ++i)
            {
                if (cos[i].edge != e)
                    continue;
                const auto ce = cos[i];
                const CoEdge<T> half{second, ce.orient, ce.pcurve};
                // Forward: e then second ; Reversed: second then e
                cos.insert(cos.begin() + static_cast<std::ptrdiff_t>(ce.orient == Orientation::Forward ? i + 1 : i), half);
                ++i;
            }
        }
        return SplitEdge{e, second, v};
    }

    // =========================================================================
    // Transformed copy
    // =========================================================================

    namespace detail
    {
        template <std::floating_point T>
        class ShapeCopier
        {
            Model<T> &m_;
            Placement<T> p_;
            std::unordered_map<VertexId, VertexId> v_;
            std::unordered_map<EdgeId, EdgeId> e_;
            std::unordered_map<WireId, WireId> w_;
            std::unordered_map<FaceId, FaceId> f_;
            std::unordered_map<ShellId, ShellId> sh_;
            std::unordered_map<SolidId, SolidId> so_;
            std::unordered_map<CompoundId, CompoundId> c_;
            std::unordered_map<const Surface<T, 3> *, std::shared_ptr<Surface<T, 3>>> srf_;
            std::unordered_map<const Curve<T, 3> *, std::shared_ptr<Curve<T, 3>>> crv_;

            template <typename Id>
            void attributes(Id from, Id to)
            {
                if (auto n = m_.name(from); !n.empty())
                    m_.setName(to, std::string(n));
                m_.setExternalId(to, m_.externalId(from));
                if (auto inst = m_.instance(from)) // the copy is an instance of the same part, moved
                    m_.setInstance(to, InstanceOf<T>{inst->part, p_ * inst->place});
            }

            std::shared_ptr<Surface<T, 3>> surface(const std::shared_ptr<Surface<T, 3>> &s)
            {
                if (!s)
                    return nullptr;
                if (auto it = srf_.find(s.get()); it != srf_.end())
                    return it->second;
                auto trf = [this](const point<T, 3> &x) { return p_.apply(x); };
                std::shared_ptr<Surface<T, 3>> r;
                if (auto e = std::dynamic_pointer_cast<ElementarySurface<T>>(s))
                    r = e->moved(trf);
                else if (auto q = std::dynamic_pointer_cast<BSSurfaceRational<T, 3>>(s))
                {
                    auto c = std::make_shared<BSSurfaceRational<T, 3>>(*q);
                    gbs::transform(*c, trf);
                    r = c;
                }
                else if (auto q = std::dynamic_pointer_cast<BSSurface<T, 3>>(s))
                {
                    auto c = std::make_shared<BSSurface<T, 3>>(*q);
                    gbs::transform(*c, trf);
                    r = c;
                }
                else
                    throw BRepError("copy_shape: only primitive and B-spline surfaces can be transformed");
                srf_.emplace(s.get(), r);
                return r;
            }

            std::shared_ptr<Curve<T, 3>> curve(const std::shared_ptr<Curve<T, 3>> &c)
            {
                if (!c)
                    return nullptr;
                if (auto it = crv_.find(c.get()); it != crv_.end())
                    return it->second;
                auto trf = [this](const point<T, 3> &x) { return p_.apply(x); };
                std::shared_ptr<Curve<T, 3>> r;
                if (auto e = std::dynamic_pointer_cast<ElementaryCurve<T>>(c))
                    r = e->moved(trf);
                else if (auto q = std::dynamic_pointer_cast<BSCurveRational<T, 3>>(c))
                {
                    auto n = std::make_shared<BSCurveRational<T, 3>>(*q);
                    gbs::transform(*n, trf);
                    r = n;
                }
                else if (auto q = std::dynamic_pointer_cast<BSCurve<T, 3>>(c))
                {
                    auto n = std::make_shared<BSCurve<T, 3>>(*q);
                    gbs::transform(*n, trf);
                    r = n;
                }
                else if (auto q = std::dynamic_pointer_cast<CurveOnSurface<T, 3>>(c)) // natural face sides
                    r = std::make_shared<CurveOnSurface<T, 3>>(q->p_basisCurve(), surface(q->p_basisSurface()));
                else
                    throw BRepError("copy_shape: only primitive and B-spline curves, and curves on such surfaces, can be transformed");
                crv_.emplace(c.get(), r);
                return r;
            }

        public:
            ShapeCopier(Model<T> &m, const Placement<T> &p) : m_{m}, p_{p} {}

            VertexId copy(VertexId id)
            {
                if (auto it = v_.find(id); it != v_.end())
                    return it->second;
                auto v = m_.vertex(id);
                v.pnt = p_.apply(v.pnt);
                const auto r = m_.add(std::move(v));
                attributes(id, r);
                return v_[id] = r;
            }
            EdgeId copy(EdgeId id)
            {
                if (auto it = e_.find(id); it != e_.end())
                    return it->second;
                auto e = m_.edge(id);
                e.curve = curve(e.curve);
                e.v1 = copy(e.v1);
                e.v2 = copy(e.v2);
                const auto r = m_.add(std::move(e));
                attributes(id, r);
                return e_[id] = r;
            }
            WireId copy(WireId id)
            {
                if (auto it = w_.find(id); it != w_.end())
                    return it->second;
                auto w = m_.wire(id);
                for (auto &ce : w.coedges)
                    ce.edge = copy(ce.edge); // pcurves unchanged: the parametrization is kept
                const auto r = m_.add(std::move(w));
                attributes(id, r);
                return w_[id] = r;
            }
            FaceId copy(FaceId id)
            {
                if (auto it = f_.find(id); it != f_.end())
                    return it->second;
                auto f = m_.face(id);
                f.surface = surface(f.surface);
                for (auto &w : f.wires)
                    w = copy(w);
                const auto r = m_.add(std::move(f));
                attributes(id, r);
                return f_[id] = r;
            }
            ShellId copy(ShellId id)
            {
                if (auto it = sh_.find(id); it != sh_.end())
                    return it->second;
                auto s = m_.shell(id);
                for (auto &fu : s.faces)
                    fu.face = copy(fu.face);
                const auto r = m_.add(std::move(s));
                attributes(id, r);
                return sh_[id] = r;
            }
            SolidId copy(SolidId id)
            {
                if (auto it = so_.find(id); it != so_.end())
                    return it->second;
                auto s = m_.solid(id);
                s.outer = copy(s.outer);
                for (auto &v : s.voids)
                    v = copy(v);
                const auto r = m_.add(std::move(s));
                attributes(id, r);
                return so_[id] = r;
            }
            CompoundId copy(CompoundId id)
            {
                if (auto it = c_.find(id); it != c_.end())
                    return it->second;
                auto c = m_.compound(id);
                for (auto &s : c.shapes)
                    s = copy(s);
                const auto r = m_.add(std::move(c));
                attributes(id, r);
                return c_[id] = r;
            }
            ShapeId copy(const ShapeId &id)
            {
                return std::visit([this](auto h) -> ShapeId { return copy(h); }, id);
            }
        };
    }

    /**
     * @brief Copy of a shape and of all its sub-shapes, geometry moved by the rigid placement `p`.
     *
     * Primitives (gbs/primitives.h) move their frame; B-spline curves and surfaces are
     * transformed (rational ones through their Cartesian poles); curves on surfaces follow
     * their transformed surface; pcurves are shared, the
     * parametrization being unchanged. Entities shared inside the shape stay shared in the copy.
     * Names and external ids are copied; an instance of a part (Model::instance) stays an
     * instance of that part, its placement composed with `p`. Throws BRepError on another
     * geometry type.
     */
    template <std::floating_point T>
    [[nodiscard]] auto copy_shape(Model<T> &m, const ShapeId &shape, const Placement<T> &p) -> ShapeId
    {
        if (!m.alive(shape))
            throw BRepError("copy_shape: dead or invalid shape");
        detail::ShapeCopier<T> c{m, p};
        return c.copy(shape);
    }
}
