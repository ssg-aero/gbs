#pragma once

/**
 * @file writer.h
 * @brief STEP writer: shapes of a gbs::brep::Model to a STEP file.
 *
 * Design: docs/sources/design/step_writer.md, sections 4, 5, 7 and 8;
 * architecture note docs/sources/design/step_wpr03_topology.md.
 *
 * - Solids, shells and faces become parts (a product with its shape
 *   representation), compounds become assemblies whose members are
 *   components placed by an item_defined_transformation.
 * - The instances of a same part (Model::instance, recorded by the STEP
 *   reader) are written once, in the coordinates of the part, and every
 *   occurrence carries its placement; option `flatten` writes everything as
 *   one part instead.
 * - Topology as in the model: vertices and edges shared, edges with their
 *   pcurves (surface_curve, seam_curve for an edge used twice by a face),
 *   loops in the order of the wires, degenerate edges left out (STEP has
 *   none), face sense from the face use (same_sense, bound orientation).
 * - Geometry through GeometryWriter: primitives as analytic entities,
 *   B-splines exactly, the rest approximated and reported.
 */

#include <gbs-io/step/p21_writer.h>
#include <gbs-io/step/geometry_writer.h>
#include <gbs-brep/brep>

#include <filesystem>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace gbs::step
{
    struct StepWriteOptions
    {
        std::string schema{"AP214"};        ///< "AP214" or "AP242"
        double approx_tol{1e-5};            ///< deviation allowed when approximating non-NURBS geometry
        bool pcurves{true};                 ///< write the pcurves of the edges
        bool flatten{false};                ///< one part with every shape instead of assemblies
        std::optional<double> uncertainty;  ///< distance uncertainty written; default brep_default_tolerance
        std::string product_name;           ///< name of a part or assembly whose shape has no name
        std::string author, organization;
    };

    struct StepWriteReport
    {
        std::size_t entities{};                       ///< instances written
        std::size_t parts{}, assemblies{}, occurrences{};
        std::size_t solids{}, shells{}, faces{}, edges{}, vertices{};
        std::size_t approximated_curves{}, approximated_surfaces{};
        double max_deviation{};                       ///< of the approximations
        std::vector<brep::ShapeId> skipped;           ///< shapes without BREP equivalent (wires, edges, vertices)
    };

    namespace detail
    {
        class TopologyWriter
        {
            using T = double;
            using V = P21Value;
            using Model = brep::Model<T>;
            using brep_place = brep::Placement<T>;

            const Model &m0_;
            const StepWriteOptions &o_;
            P21Writer w_;
            GeometryWriter g_;
            std::uint64_t ctx_{}, ctx2d_{}, app_{}, pc_{}, pdc_{};
            std::vector<std::uint64_t> products_;
            StepWriteReport rep_;
            std::size_t unnamed_{};

            // a copy of the model in which moved instances are brought back to the part coordinates
            std::optional<Model> scratch_;
            brep::IdRemap scratch_map_;

            struct Written // a product and the representation of its shape
            {
                std::uint64_t pd, rep, origin;
                std::array<std::size_t, 3> signature; // faces, edges, vertices, to check that instances agree
            };
            std::unordered_map<std::uint64_t, Written> parts_by_key_;

            // ---------------------------------------------------------------- products

            std::string name_of(const Model &m, const brep::ShapeId &s, const char *kind)
            {
                if (auto n = m.name(s); !n.empty())
                    return std::string(n);
                if (!o_.product_name.empty())
                    return o_.product_name;
                return std::string(kind) + std::to_string(++unnamed_);
            }

            std::uint64_t product(const std::string &name)
            {
                const auto p = w_.add("PRODUCT", {V::string(name), V::string(name), V::string(""), V::refs({pc_})});
                products_.push_back(p);
                const auto pdf = w_.add("PRODUCT_DEFINITION_FORMATION", {V::string(""), V::string(""), V::ref(p)});
                return w_.add("PRODUCT_DEFINITION", {V::string("design"), V::string(""), V::ref(pdf), V::ref(pdc_)});
            }

            void shape_definition(std::uint64_t pd, std::uint64_t rep)
            {
                const auto pds = w_.add("PRODUCT_DEFINITION_SHAPE", {V::string(""), V::string(""), V::ref(pd)});
                w_.add("SHAPE_DEFINITION_REPRESENTATION", {V::ref(pds), V::ref(rep)});
            }

            std::uint64_t origin_axis() { return g_.axis2(Frame3<T>{}); }

            std::uint64_t axis_of(const brep_place &p)
            {
                Frame3<T> f;
                f.o = p.t;
                for (std::size_t i = 0; i < 3; ++i)
                    f.x[i] = p.r[i][0], f.y[i] = p.r[i][1], f.z[i] = p.r[i][2];
                return g_.axis2(f);
            }

            // ---------------------------------------------------------------- topology of a part

            struct PartCache
            {
                std::unordered_map<brep::VertexId, std::uint64_t> v;
                std::unordered_map<brep::EdgeId, std::uint64_t> e;
                std::unordered_map<const Surface<T, 3> *, std::uint64_t> srf;
                std::map<std::pair<std::uint32_t, int>, std::uint64_t> faces; // (face index, sense)
                std::unordered_map<brep::EdgeId, std::vector<std::pair<brep::FaceId, const brep::CoEdge<T> *>>> uses;
            };

            std::uint64_t surface_of(const Model &m, PartCache &c, brep::FaceId f)
            {
                const auto &s = m.face(f).surface;
                if (!s)
                    throw BRepError("face without surface");
                if (auto it = c.srf.find(s.get()); it != c.srf.end())
                    return it->second;
                return c.srf[s.get()] = g_.surface(*s);
            }

            std::uint64_t vertex(const Model &m, PartCache &c, brep::VertexId v)
            {
                if (auto it = c.v.find(v); it != c.v.end())
                    return it->second;
                ++rep_.vertices;
                return c.v[v] = w_.add("VERTEX_POINT", {V::string(""), V::ref(g_.cartesian_point(m.vertex(v).pnt))});
            }

            std::uint64_t edge(const Model &m, PartCache &c, brep::EdgeId eid)
            {
                if (auto it = c.e.find(eid); it != c.e.end())
                    return it->second;
                const auto &e = m.edge(eid);
                const auto v1 = vertex(m, c, e.v1), v2 = vertex(m, c, e.v2);
                std::uint64_t geom = g_.curve(*e.curve, e.u1, e.u2);
                if (o_.pcurves)
                {
                    // pcurves of its co-edges, the Forward use first for a seam
                    auto uses = c.uses[eid];
                    std::ranges::stable_sort(uses, {}, [](const auto &u) { return u.second->orient == brep::Orientation::Forward ? 0 : 1; });
                    std::vector<std::uint64_t> pcs;
                    bool seam = false;
                    for (std::size_t i = 0; i < uses.size() && pcs.size() < 2; ++i)
                    {
                        const auto &[fid, ce] = uses[i];
                        if (!ce->pcurve)
                            continue;
                        seam = seam || std::ranges::any_of(uses, [&](const auto &u) { return u.first == fid && u.second != ce; });
                        const auto c2 = g_.curve2d(*ce->pcurve, e.u1, e.u2);
                        const auto dr = w_.add("DEFINITIONAL_REPRESENTATION", {V::string(""), V::refs({c2}), V::ref(ctx2d_)});
                        pcs.push_back(w_.add("PCURVE", {V::string(""), V::ref(surface_of(m, c, fid)), V::ref(dr)}));
                    }
                    if (!pcs.empty())
                        geom = w_.add(seam ? "SEAM_CURVE" : "SURFACE_CURVE", {V::string(""), V::ref(geom), V::refs(pcs), V::enumeration("PCURVE_S1")});
                }
                ++rep_.edges;
                const auto id = w_.add("EDGE_CURVE", {V::string(std::string(m.name(eid))), V::ref(v1), V::ref(v2), V::ref(geom), V::boolean(true)});
                return c.e[eid] = id;
            }

            std::uint64_t face(const Model &m, PartCache &c, brep::FaceId fid, brep::Orientation sense)
            {
                const int key = sense == brep::Orientation::Forward ? 0 : 1;
                if (auto it = c.faces.find({fid.index, key}); it != c.faces.end())
                    return it->second;
                // the other use already written: an oriented face of it
                if (auto it = c.faces.find({fid.index, 1 - key}); it != c.faces.end())
                    return c.faces[{fid.index, key}] = w_.add("ORIENTED_FACE", {V::string(""), V::derived(), V::ref(it->second), V::boolean(false)});
                const auto &f = m.face(fid);
                const bool same = sense == brep::Orientation::Forward;
                const auto srf = surface_of(m, c, fid);
                std::vector<std::uint64_t> bounds;
                for (std::size_t wi = 0; wi < f.wires.size(); ++wi)
                {
                    std::vector<std::uint64_t> oes;
                    for (const auto &ce : m.wire(f.wires[wi]).coedges)
                    {
                        if (m.edge(ce.edge).degenerate)
                            continue; // STEP has no degenerate edge: readers insert them back at poles
                        oes.push_back(w_.add("ORIENTED_EDGE", {V::string(""), V::derived(), V::derived(), V::ref(edge(m, c, ce.edge)),
                                                               V::boolean(ce.orient == brep::Orientation::Forward)}));
                    }
                    if (oes.empty())
                        continue;
                    const auto loop = w_.add("EDGE_LOOP", {V::string(""), V::refs(oes)});
                    // loops are counter-clockwise around the surface normal; a reversed use reverses them
                    bounds.push_back(w_.add(wi == 0 ? "FACE_OUTER_BOUND" : "FACE_BOUND", {V::string(""), V::ref(loop), V::boolean(same)}));
                }
                ++rep_.faces;
                return c.faces[{fid.index, key}] = w_.add("ADVANCED_FACE", {V::string(std::string(m.name(fid))), V::refs(bounds), V::ref(srf), V::boolean(same)});
            }

            std::uint64_t shell(const Model &m, PartCache &c, brep::ShellId sid, bool flip = false)
            {
                std::vector<std::uint64_t> faces;
                for (const auto &fu : m.shell(sid).faces)
                    faces.push_back(face(m, c, fu.face, flip ? brep::reverse(fu.orient) : fu.orient));
                ++rep_.shells;
                return w_.add(brep::is_closed(m, sid) ? "CLOSED_SHELL" : "OPEN_SHELL", {V::string(std::string(m.name(sid))), V::refs(faces)});
            }

            /// Writes the shapes of a part; returns the representation of its shape and its origin axis.
            std::pair<std::uint64_t, std::uint64_t> part_geometry(const Model &m, const std::vector<brep::ShapeId> &shapes)
            {
                PartCache c;
                // edge uses (face, co-edge) for the pcurves, and the surfaces first
                for (const auto &s : shapes)
                    for (auto fid : brep::explore<brep::FaceId>(m, s))
                    {
                        surface_of(m, c, fid);
                        for (auto wid : m.face(fid).wires)
                            for (const auto &ce : m.wire(wid).coedges)
                                c.uses[ce.edge].emplace_back(fid, &ce);
                    }
                std::vector<std::uint64_t> solids, surfaces;
                for (const auto &s : shapes)
                {
                    if (auto so = std::get_if<brep::SolidId>(&s))
                    {
                        const auto &sd = m.solid(*so);
                        const auto outer = shell(m, c, sd.outer);
                        ++rep_.solids;
                        if (sd.voids.empty())
                            solids.push_back(w_.add("MANIFOLD_SOLID_BREP", {V::string(std::string(m.name(s))), V::ref(outer)}));
                        else
                        {
                            std::vector<std::uint64_t> voids;
                            for (auto v : sd.voids) // cavities point inwards: written outwards, used reversed
                                voids.push_back(w_.add("ORIENTED_CLOSED_SHELL", {V::string(""), V::derived(), V::ref(shell(m, c, v, true)), V::boolean(false)}));
                            solids.push_back(w_.add("BREP_WITH_VOIDS", {V::string(std::string(m.name(s))), V::ref(outer), V::refs(voids)}));
                        }
                    }
                    else if (auto sh = std::get_if<brep::ShellId>(&s))
                        surfaces.push_back(w_.add("SHELL_BASED_SURFACE_MODEL", {V::string(std::string(m.name(s))), V::refs({shell(m, c, *sh)})}));
                    else if (auto f = std::get_if<brep::FaceId>(&s))
                    {
                        const auto open = w_.add("OPEN_SHELL", {V::string(""), V::refs({face(m, c, *f, brep::Orientation::Forward)})});
                        ++rep_.shells;
                        surfaces.push_back(w_.add("SHELL_BASED_SURFACE_MODEL", {V::string(std::string(m.name(s))), V::refs({open})}));
                    }
                }
                const auto origin = origin_axis();
                auto items = [origin](std::vector<std::uint64_t> v) {
                    v.push_back(origin);
                    return V::refs(v);
                };
                if (surfaces.empty())
                    return {w_.add("ADVANCED_BREP_SHAPE_REPRESENTATION", {V::string(""), items(solids), V::ref(ctx_)}), origin};
                if (solids.empty())
                    return {w_.add("MANIFOLD_SURFACE_SHAPE_REPRESENTATION", {V::string(""), items(surfaces), V::ref(ctx_)}), origin};
                // both: a shape representation linked to one representation of each kind
                const auto sr = w_.add("SHAPE_REPRESENTATION", {V::string(""), V::refs({origin}), V::ref(ctx_)});
                const auto a = w_.add("ADVANCED_BREP_SHAPE_REPRESENTATION", {V::string(""), items(solids), V::ref(ctx_)});
                const auto b = w_.add("MANIFOLD_SURFACE_SHAPE_REPRESENTATION", {V::string(""), items(surfaces), V::ref(ctx_)});
                w_.add("SHAPE_REPRESENTATION_RELATIONSHIP", {V::string(""), V::string(""), V::ref(sr), V::ref(a)});
                w_.add("SHAPE_REPRESENTATION_RELATIONSHIP", {V::string(""), V::string(""), V::ref(sr), V::ref(b)});
                return {sr, origin};
            }

            static bool leaf(const brep::ShapeId &s)
            {
                return std::holds_alternative<brep::SolidId>(s) || std::holds_alternative<brep::ShellId>(s) || std::holds_alternative<brep::FaceId>(s);
            }

            std::array<std::size_t, 3> signature(const Model &m, const brep::ShapeId &s) const
            {
                return {brep::explore<brep::FaceId>(m, s).size(), brep::explore<brep::EdgeId>(m, s).size(), brep::explore<brep::VertexId>(m, s).size()};
            }

            /// The shape in the coordinates of its part: itself if its placement is the identity,
            /// otherwise a copy brought back by the inverse placement, in the scratch model.
            std::pair<const Model *, brep::ShapeId> in_part_coordinates(const brep::ShapeId &s, const brep_place &p)
            {
                if (p.is_identity(1e-12))
                    return {&m0_, s};
                if (!scratch_)
                {
                    scratch_.emplace();
                    scratch_->setUnitScale(m0_.unitScale());
                    scratch_map_ = scratch_->append(m0_);
                }
                return {&*scratch_, brep::copy_shape(*scratch_, scratch_map_.map(s), p.inverse())};
            }

            /// A part (solid, shell or face): written once per part key, with a fresh product otherwise.
            Written part(const brep::ShapeId &s)
            {
                const auto inst = m0_.instance(s);
                const auto sig = signature(m0_, s);
                if (inst)
                    if (auto it = parts_by_key_.find(inst->part); it != parts_by_key_.end() && it->second.signature == sig)
                        return it->second;
                const auto [m, local] = in_part_coordinates(s, inst ? inst->place : brep_place{});
                const auto [rep, origin] = part_geometry(*m, {local});
                const auto pd = product(name_of(m0_, s, "part"));
                shape_definition(pd, rep);
                ++rep_.parts;
                Written w{pd, rep, origin, sig};
                if (inst)
                    parts_by_key_.emplace(inst->part, w);
                return w;
            }

            /// An assembly from a compound: its members are components placed relatively to it.
            Written assembly(brep::CompoundId cid, int depth = 0)
            {
                if (depth > 64)
                    throw BRepError("compound nested deeper than 64 levels");
                const brep::ShapeId s{cid};
                const auto inst = m0_.instance(s);
                const auto sig = signature(m0_, s);
                if (inst)
                    if (auto it = parts_by_key_.find(inst->part); it != parts_by_key_.end() && it->second.signature == sig)
                        return it->second;
                const brep_place own = inst ? inst->place : brep_place{};
                struct Member
                {
                    Written w;
                    std::uint64_t placement;
                };
                std::vector<Member> members;
                for (const auto &ms : m0_.compound(cid).shapes)
                {
                    std::optional<Written> w;
                    if (leaf(ms))
                        w = part(ms);
                    else if (auto sub = std::get_if<brep::CompoundId>(&ms))
                        w = assembly(*sub, depth + 1);
                    else
                    {
                        rep_.skipped.push_back(ms);
                        continue;
                    }
                    // placement of the member in the assembly: own^-1 o (member part -> model)
                    const auto mi = m0_.instance(ms);
                    const auto rel = own.inverse() * (mi ? mi->place : brep_place{});
                    members.push_back({*w, axis_of(rel)});
                }
                std::vector<std::uint64_t> items{origin_axis()};
                for (const auto &mb : members)
                    items.push_back(mb.placement);
                const auto rep = w_.add("SHAPE_REPRESENTATION", {V::string(""), V::refs(items), V::ref(ctx_)});
                const auto pd = product(name_of(m0_, s, "assembly"));
                shape_definition(pd, rep);
                std::size_t k = 0;
                for (const auto &mb : members)
                {
                    const auto label = "occurrence" + std::to_string(++k);
                    const auto nauo = w_.add("NEXT_ASSEMBLY_USAGE_OCCURRENCE", {V::string(label), V::string(label), V::string(""), V::ref(pd), V::ref(mb.w.pd), V::omitted()});
                    const auto pds = w_.add("PRODUCT_DEFINITION_SHAPE", {V::string(""), V::string(""), V::ref(nauo)});
                    const auto idt = w_.add("ITEM_DEFINED_TRANSFORMATION", {V::string(""), V::string(""), V::ref(mb.w.origin), V::ref(mb.placement)});
                    const auto rel = w_.add_complex({{"REPRESENTATION_RELATIONSHIP", {V::string(""), V::string(""), V::ref(mb.w.rep), V::ref(rep)}},
                                                     {"REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION", {V::ref(idt)}},
                                                     {"SHAPE_REPRESENTATION_RELATIONSHIP", {}}});
                    w_.add("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION", {V::ref(rel), V::ref(pds)});
                    ++rep_.occurrences;
                }
                ++rep_.assemblies;
                Written w{pd, rep, items.front(), sig};
                if (inst)
                    parts_by_key_.emplace(inst->part, w);
                return w;
            }

            void flatten_into(const brep::ShapeId &s, std::vector<brep::ShapeId> &leaves)
            {
                if (leaf(s))
                    leaves.push_back(s);
                else if (auto c = std::get_if<brep::CompoundId>(&s))
                    for (const auto &x : m0_.compound(*c).shapes)
                        flatten_into(x, leaves);
                else
                    rep_.skipped.push_back(s);
            }

        public:
            TopologyWriter(const Model &m, const StepWriteOptions &o) : m0_{m}, o_{o}, w_{}, g_{w_, o.approx_tol} {}

            std::pair<std::string, StepWriteReport> run(std::span<const brep::ShapeId> shapes)
            {
                for (const auto &s : shapes)
                    if (!m0_.alive(s))
                        throw BRepError("write_step: dead or invalid shape");
                const bool ap242 = o_.schema == "AP242";
                if (!ap242 && o_.schema != "AP214")
                    throw std::invalid_argument("write_step: schema must be AP214 or AP242");
                auto &h = w_.header();
                h.schemas = {std::string(ap242 ? schema::ap242 : schema::ap214)};
                h.author = {o_.author};
                h.organization = {o_.organization};

                ctx_ = g_.context(m0_.unitScale(), o_.uncertainty.value_or(gbs::brep_default_tolerance<T>));
                ctx2d_ = g_.context2d();
                app_ = w_.add("APPLICATION_CONTEXT", {V::string(ap242 ? "managed model based 3d engineering" : "core data for automotive mechanical design processes")});
                w_.add("APPLICATION_PROTOCOL_DEFINITION", {V::string("international standard"),
                                                           V::string(ap242 ? "ap242_managed_model_based_3d_engineering" : "automotive_design"),
                                                           V::integer(ap242 ? 2014 : 2000), V::ref(app_)});
                pc_ = w_.add("PRODUCT_CONTEXT", {V::string(""), V::ref(app_), V::string("mechanical")});
                pdc_ = w_.add("PRODUCT_DEFINITION_CONTEXT", {V::string("part definition"), V::ref(app_), V::string("design")});

                if (o_.flatten)
                {
                    std::vector<brep::ShapeId> leaves;
                    for (const auto &s : shapes)
                        flatten_into(s, leaves);
                    if (!leaves.empty())
                    {
                        const auto [rep, origin] = part_geometry(m0_, leaves);
                        const auto pd = product(shapes.size() == 1 ? name_of(m0_, shapes.front(), "part") : (o_.product_name.empty() ? "part" : o_.product_name));
                        shape_definition(pd, rep);
                        ++rep_.parts;
                    }
                }
                else
                    for (const auto &s : shapes)
                    {
                        if (leaf(s))
                            part(s);
                        else if (auto c = std::get_if<brep::CompoundId>(&s))
                            assembly(*c);
                        else
                            rep_.skipped.push_back(s);
                    }
                if (!products_.empty())
                    w_.add("PRODUCT_RELATED_PRODUCT_CATEGORY", {V::string("part"), V::omitted(), V::refs(products_)});

                rep_.entities = w_.size();
                rep_.approximated_curves = g_.stats().approximated_curves;
                rep_.approximated_surfaces = g_.stats().approximated_surfaces;
                rep_.max_deviation = g_.stats().max_deviation;
                return {w_.str(), std::move(rep_)};
            }
        };
    } // namespace detail

    /**
     * @brief Writes shapes of `m` as a STEP file content.
     *
     * Solids, shells and faces become parts, compounds assemblies (see the file comment);
     * wires, edges and vertices are skipped and reported. Fails (InvalidId, InvalidFile) on a
     * dead shape, an unknown schema or invalid geometry; nothing is written then.
     */
    inline auto write_step(const brep::Model<double> &m, std::span<const brep::ShapeId> shapes, const StepWriteOptions &opts = {})
        -> brep::BuildResult<std::pair<std::string, StepWriteReport>>
    {
        try
        {
            detail::TopologyWriter w{m, opts};
            return w.run(shapes);
        }
        catch (const BRepError &e)
        {
            return std::unexpected(brep::BuildError{brep::BuildErrc::InvalidId, e.what(), {}});
        }
        catch (const std::exception &e)
        {
            return std::unexpected(brep::BuildError{brep::BuildErrc::InvalidFile, e.what(), {}});
        }
    }

    inline auto write_step(const brep::Model<double> &m, const brep::ShapeId &shape, const StepWriteOptions &opts = {})
        -> brep::BuildResult<std::pair<std::string, StepWriteReport>>
    {
        return write_step(m, std::span<const brep::ShapeId>{&shape, 1}, opts);
    }

    /// Writes shapes of `m` to a STEP file; see write_step().
    inline auto write_step_file(const brep::Model<double> &m, std::span<const brep::ShapeId> shapes, const std::filesystem::path &path,
                                const StepWriteOptions &opts = {}) -> brep::BuildResult<StepWriteReport>
    {
        auto r = write_step(m, shapes, opts);
        if (!r)
            return std::unexpected(r.error());
        std::ofstream f(path, std::ios::binary);
        f.write(r->first.data(), static_cast<std::streamsize>(r->first.size()));
        if (!f)
            return std::unexpected(brep::BuildError{brep::BuildErrc::InvalidFile, "cannot write " + path.string(), {}});
        return std::move(r->second);
    }

    inline auto write_step_file(const brep::Model<double> &m, const brep::ShapeId &shape, const std::filesystem::path &path,
                                const StepWriteOptions &opts = {}) -> brep::BuildResult<StepWriteReport>
    {
        return write_step_file(m, std::span<const brep::ShapeId>{&shape, 1}, path, opts);
    }
}
