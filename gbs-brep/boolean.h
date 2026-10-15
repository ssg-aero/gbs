#pragma once
/**
 * @file boolean.h
 * @brief Boolean operations on solids: fuse, cut, common.
 *
 * A simplified general fuse, on copies of the operands (which are not modified):
 * 1. Section of the copies (`section.h`): the section edges are shared by the faces of both.
 * 2. Split of the faces of both copies along the section (`split.h`).
 * 3. Classification of every face against the other operand, by an interior point
 *    (`classify.h`): inside or outside.
 * 4. Selection: fuse keeps A outside B and B outside A; common keeps A inside B and B inside A;
 *    cut keeps A outside B and B inside A, reversed.
 * 5. Assembly: the kept faces already share their edges; shells are their connected components,
 *    solids their closed shells, cavities the shells of negative volume, given to the smallest
 *    solid around them.
 *
 * Generic configurations only: faces meeting along a tangency, coincident faces, an edge lying on
 * a face of the other operand or a singular intersection give `BuildErrc::NonGeneric`.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/section.h>
#include <gbs-brep/split.h>

#include <algorithm>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <vector>

namespace gbs::brep
{
    enum class BooleanOp
    {
        Fuse,
        Cut,
        Common
    };

    template <std::floating_point T>
    struct BooleanOptions
    {
        T tol{brep_default_tolerance<T>};
        IntersectOptions<T> intersect{};
    };

    /// Operand a face of the result comes from.
    enum class FaceOrigin
    {
        A,
        B
    };

    template <std::floating_point T>
    struct BooleanResult
    {
        ShapeId shape;                                   ///< a solid, or a compound of solids (empty compound for an empty result)
        std::vector<SolidId> solids;
        std::vector<std::pair<FaceId, FaceOrigin>> faces; ///< faces of the result and their operand
    };

    namespace boolean_detail
    {
        template <std::floating_point T>
        auto fail_non_generic(const std::string &why) -> std::unexpected<BuildError>
        {
            return detail::fail(BuildErrc::NonGeneric, "boolean: " + why);
        }

        /// Erases the faces, wires, shells, solids and compounds reachable from a shape.
        template <std::floating_point T>
        void erase_shape(Model<T> &m, const ShapeId &s, const std::set<std::uint32_t> &keep_faces = {})
        {
            for (auto c : explore<CompoundId>(m, s))
                m.erase(c);
            for (auto so : explore<SolidId>(m, s))
                m.erase(so);
            for (auto sh : explore<ShellId>(m, s))
                m.erase(sh);
            for (auto f : explore<FaceId>(m, s))
                if (!keep_faces.contains(f.index) && m.alive(f))
                {
                    for (auto w : m.face(f).wires)
                        m.erase(w);
                    m.erase(f);
                }
        }

        /// Position of a face against the solids of an operand: inside if inside one of them.
        template <std::floating_point T>
        auto locate(const Model<T> &m, FaceId f, const std::vector<SolidId> &solids, T tol) -> std::optional<Location>
        {
            Location r = Location::Out;
            for (auto s : solids)
            {
                const auto l = classify(m, f, s, tol);
                if (!l)
                    return std::nullopt;
                if (*l == Location::On)
                    return Location::On;
                if (*l == Location::In)
                    r = Location::In;
            }
            return r;
        }
    }

    /**
     * @brief Boolean operation `op` between the solids of `a` and those of `b` (a solid or a
     * compound of solids). The operands are not modified; the result is built from copies of
     * their faces.
     */
    template <std::floating_point T>
    [[nodiscard]] auto boolean(Model<T> &m, BooleanOp op, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
        -> BuildResult<BooleanResult<T>>
    {
        using namespace boolean_detail;
        if (!m.alive(a) || !m.alive(b))
            return detail::fail(BuildErrc::InvalidId, "boolean: dead or invalid shape");
        if (explore<SolidId>(m, a).empty() || explore<SolidId>(m, b).empty())
            return detail::fail(BuildErrc::InvalidId, "boolean: both operands must contain solids");
        const T tol = opts.tol;

        // ---- copies, so that the operands are not modified
        const auto ca = copy_shape(m, a, Placement<T>{}), cb = copy_shape(m, b, Placement<T>{});
        auto abort = [&](BuildError e) -> BuildResult<BooleanResult<T>> {
            erase_shape(m, ca);
            erase_shape(m, cb);
            return std::unexpected(std::move(e));
        };
        const auto solids_a = explore<SolidId>(m, ca), solids_b = explore<SolidId>(m, cb);

        // ---- 1. section
        SectionOptions<T> so{tol, opts.intersect};
        auto sec = section(m, ca, cb, so);
        if (!sec)
            return abort(sec.error());
        if (sec->singular || !sec->complete)
            return abort(fail_non_generic<T>("singular or incomplete intersection").error());
        if (!sec->coincident.empty())
            return abort(fail_non_generic<T>("coincident faces").error());
        if (!sec->boundary_edges.empty() || !sec->points.empty())
            return abort(fail_non_generic<T>("an edge or a vertex lies on a face of the other operand").error());
        for (const auto &se : sec->edges)
            if (se.tangent)
                return abort(fail_non_generic<T>("faces tangent along a section edge").error());

        // ---- 2. split of the faces of both copies by the same edges
        std::map<std::uint32_t, std::vector<FaceCut<T>>> cuts;
        for (const auto &se : sec->edges)
        {
            cuts[se.face_a.index].push_back({se.edge, se.pcurve_a});
            cuts[se.face_b.index].push_back({se.edge, se.pcurve_b});
        }
        for (auto &[index, c] : cuts)
        {
            auto s = split_face(m, FaceId{index}, c, tol);
            if (!s)
                return abort(s.error());
            if (!s->dangling.empty())
                return abort(fail_non_generic<T>("a section edge does not split a face").error());
        }
        // the section wires and compound are not part of the result (their edges are)
        for (auto w : sec->wires)
            m.erase(w);
        m.erase(sec->compound);

        // ---- 3. classification and 4. selection
        struct Kept
        {
            FaceId face;
            Orientation orient;
            FaceOrigin origin;
        };
        std::vector<Kept> kept;
        for (auto [shape, other, origin] : {std::tuple{ca, &solids_b, FaceOrigin::A}, std::tuple{cb, &solids_a, FaceOrigin::B}})
            for (auto sh : explore<ShellId>(m, shape))
                for (const auto &fu : m.shell(sh).faces)
                {
                    const auto where = locate(m, fu.face, *other, tol);
                    if (!where || *where == Location::On)
                        return abort(fail_non_generic<T>("a face lies on the boundary of the other operand").error());
                    const bool in = *where == Location::In;
                    bool keep{}, flip{};
                    switch (op)
                    {
                    case BooleanOp::Fuse: keep = !in; break;
                    case BooleanOp::Common: keep = in; break;
                    case BooleanOp::Cut:
                        keep = origin == FaceOrigin::A ? !in : in;
                        flip = origin == FaceOrigin::B;
                        break;
                    }
                    if (keep)
                        kept.push_back({fu.face, flip ? reverse(fu.orient) : fu.orient, origin});
                }

        // ---- 5. assembly: shells by shared edges, solids by closed shells
        std::map<std::uint32_t, std::vector<std::size_t>> by_edge;
        for (std::size_t k = 0; k < kept.size(); ++k)
            for (auto e : explore<EdgeId>(m, ShapeId{kept[k].face}))
                by_edge[e.index].push_back(k);
        std::vector<std::size_t> parent(kept.size());
        for (std::size_t k = 0; k < kept.size(); ++k)
            parent[k] = k;
        std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = root(parent[x]); };
        for (const auto &[e, ks] : by_edge)
            for (std::size_t i = 1; i < ks.size(); ++i)
                parent[root(ks[i])] = root(ks[0]);
        std::map<std::size_t, Shell> shells;
        for (std::size_t k = 0; k < kept.size(); ++k)
            shells[root(k)].faces.push_back(FaceUse{kept[k].face, kept[k].orient});
        std::vector<std::pair<ShellId, T>> outers, voids;
        std::set<std::uint32_t> kept_faces;
        for (const auto &k : kept)
            kept_faces.insert(k.face.index);
        for (auto &[r, sh] : shells)
        {
            const auto id = m.add(sh);
            if (!is_closed(m, id))
            {
                erase_shape(m, ca, kept_faces);
                erase_shape(m, cb, kept_faces);
                return detail::fail(BuildErrc::ShellNotClosed, "boolean: a shell of the result is open", {id});
            }
            m.shell(id).closed = true;
            const T v = signed_volume(m, id, 64);
            (v > T(0) ? outers : voids).emplace_back(id, v);
        }
        // each cavity in the smallest solid whose box contains it
        std::map<std::uint32_t, std::vector<ShellId>> cavities;
        for (const auto &[vid, vv] : voids)
        {
            const auto vb = bounding_box(m, ShapeId{vid});
            std::optional<std::pair<ShellId, T>> owner;
            for (const auto &[oid, ov] : outers)
            {
                const auto ob = bounding_box(m, ShapeId{oid});
                if (ob.contains(vb.min) && ob.contains(vb.max) && (!owner || ov < owner->second))
                    owner = std::pair{oid, ov};
            }
            if (!owner)
            {
                erase_shape(m, ca, kept_faces);
                erase_shape(m, cb, kept_faces);
                return detail::fail(BuildErrc::HoleOutsideOuter, "boolean: a cavity of the result is outside every solid", {vid});
            }
            cavities[owner->first.index].push_back(vid);
        }
        // the copies give their faces to the result
        erase_shape(m, ca, kept_faces);
        erase_shape(m, cb, kept_faces);

        BooleanResult<T> res;
        for (const auto &[oid, ov] : outers)
        {
            auto s = make_solid(m, oid, cavities[oid.index]);
            if (!s)
                return std::unexpected(s.error());
            res.solids.push_back(*s);
        }
        for (const auto &k : kept)
            res.faces.emplace_back(k.face, k.origin);
        if (res.solids.size() == 1)
            res.shape = res.solids.front();
        else
        {
            std::vector<ShapeId> shapes(res.solids.begin(), res.solids.end());
            res.shape = unwrap(make_compound(m, shapes));
        }
        return res;
    }

    template <std::floating_point T>
    [[nodiscard]] auto fuse(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Fuse, a, b, opts);
    }

    template <std::floating_point T>
    [[nodiscard]] auto cut(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Cut, a, b, opts);
    }

    template <std::floating_point T>
    [[nodiscard]] auto common(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Common, a, b, opts);
    }
}
