// Python bindings of the native BREP core (gbs-brep), sub-module gbs.brep.
// Design: docs/sources/design/brep_core.md § 8.2; architecture notes
// docs/sources/design/brep_pr10_python.md and step_pr10_python.md (STEP reader).
#include <pybind11/pybind11.h>
#include <cstdint>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <gbs-brep/brep>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>

#include <algorithm>
#include <sstream>
#include <string>

namespace py = pybind11;
using namespace gbs;
using namespace gbs::brep;

namespace
{
    using T = double;

    // Builders return std::expected in C++; in Python a failure raises BRepError.
    template <typename V>
    V value_or_raise(BuildResult<V> r)
    {
        return unwrap(std::move(r));
    }

    template <typename Id>
    void declare_id(py::module_ &m, const char *name)
    {
        py::class_<Id>(m, name, "Typed identifier of an entity of a Model (an index, owns nothing)")
            .def(py::init<>())
            .def(py::init([](std::uint32_t i) { return Id{i}; }), py::arg("index"))
            .def_readonly("index", &Id::index)
            .def("valid", &Id::valid)
            .def("__int__", [](const Id &h) { return h.index; })
            .def("__index__", [](const Id &h) { return h.index; })
            .def("__eq__", [](const Id &a, const Id &b) { return a == b; })
            .def("__lt__", [](const Id &a, const Id &b) { return a < b; })
            .def("__hash__", [](const Id &h) { return std::hash<Id>{}(h); })
            .def("__repr__", [name](const Id &h) {
                return std::string(name) + "(" + (h.valid() ? std::to_string(h.index) : std::string("invalid")) + ")";
            });
    }

    template <typename Id>
    void declare_explore(py::module_ &m, const char *name)
    {
        m.def(name, [](const Model<T> &mdl, const ShapeId &s) { return explore<Id>(mdl, s); },
              "Sub-entities of a shape by type, without duplicates", py::arg("model"), py::arg("shape"));
    }
}

void gbs_bind_brep(py::module &parent)
{
    auto m = parent.def_submodule("brep", "Native BREP core: model, builders, sewing, check, solids");

    py::register_exception<BRepError>(m, "BRepError", PyExc_RuntimeError);

    // ---- enums --------------------------------------------------------------
    py::enum_<ShapeType>(m, "ShapeType")
        .value("Vertex", ShapeType::Vertex)
        .value("Edge", ShapeType::Edge)
        .value("Wire", ShapeType::Wire)
        .value("Face", ShapeType::Face)
        .value("Shell", ShapeType::Shell)
        .value("Solid", ShapeType::Solid)
        .value("Compound", ShapeType::Compound);
    py::enum_<Orientation>(m, "Orientation")
        .value("Forward", Orientation::Forward)
        .value("Reversed", Orientation::Reversed);
    m.def("reverse", &gbs::brep::reverse, py::arg("orientation"));

    // ---- identifiers ----------------------------------------------------------
    declare_id<VertexId>(m, "VertexId");
    declare_id<EdgeId>(m, "EdgeId");
    declare_id<WireId>(m, "WireId");
    declare_id<FaceId>(m, "FaceId");
    declare_id<ShellId>(m, "ShellId");
    declare_id<SolidId>(m, "SolidId");
    declare_id<CompoundId>(m, "CompoundId");
    m.def("shape_type", &gbs::brep::shape_type, py::arg("shape"));

    // ---- entities (plain data) -----------------------------------------------
    py::class_<Vertex<T>>(m, "Vertex")
        .def(py::init([](const point<T, 3> &p, T tol) { return Vertex<T>{p, tol}; }), py::arg("pnt"), py::arg("tol") = brep_default_tolerance<T>)
        .def_readwrite("pnt", &Vertex<T>::pnt)
        .def_readwrite("tol", &Vertex<T>::tol);
    py::class_<Edge<T>>(m, "Edge")
        .def_readwrite("curve", &Edge<T>::curve)
        .def_readwrite("u1", &Edge<T>::u1)
        .def_readwrite("u2", &Edge<T>::u2)
        .def_readwrite("v1", &Edge<T>::v1)
        .def_readwrite("v2", &Edge<T>::v2)
        .def_readwrite("tol", &Edge<T>::tol)
        .def_readwrite("degenerate", &Edge<T>::degenerate)
        .def_readwrite("same_parameter", &Edge<T>::same_parameter);
    py::class_<CoEdge<T>>(m, "CoEdge")
        .def_readwrite("edge", &CoEdge<T>::edge)
        .def_readwrite("orient", &CoEdge<T>::orient)
        .def_readwrite("pcurve", &CoEdge<T>::pcurve);
    py::class_<Wire<T>>(m, "Wire")
        .def_readwrite("coedges", &Wire<T>::coedges)
        .def_readwrite("closed", &Wire<T>::closed);
    py::class_<Face<T>>(m, "Face")
        .def_readwrite("surface", &Face<T>::surface)
        .def_readwrite("wires", &Face<T>::wires)
        .def_readwrite("tol", &Face<T>::tol)
        .def_readwrite("natural_bounds", &Face<T>::natural_bounds)
        .def_property_readonly(
            "surface_id", [](const Face<T> &f) { return reinterpret_cast<std::uintptr_t>(f.surface.get()); },
            "Identity of the surface object: equal for faces sharing one surface (pieces of a split face)");
    py::class_<FaceUse>(m, "FaceUse")
        .def(py::init([](FaceId f, Orientation o) { return FaceUse{f, o}; }), py::arg("face"), py::arg("orient") = Orientation::Forward)
        .def_readwrite("face", &FaceUse::face)
        .def_readwrite("orient", &FaceUse::orient);
    py::class_<Shell>(m, "Shell")
        .def_readwrite("faces", &Shell::faces)
        .def_readwrite("closed", &Shell::closed);
    py::class_<Solid>(m, "Solid")
        .def_readwrite("outer", &Solid::outer)
        .def_readwrite("voids", &Solid::voids);
    py::class_<Compound>(m, "Compound")
        .def_readwrite("shapes", &Compound::shapes);

    // ---- model ---------------------------------------------------------------
    py::class_<IdRemap>(m, "IdRemap")
        .def("map", [](const IdRemap &r, const ShapeId &old) { return r.map(old); }, py::arg("old"));

    using M = Model<T>;
    // Getters return copies and set_* write back: a Python reference into the arena
    // would dangle as soon as an entity is added (the tables are std::vector).
    py::class_<M, std::shared_ptr<M>>(m, "Model", "Arena owning every topological entity of a BREP model")
        .def(py::init<>())
        .def("vertex", [](const M &x, VertexId i) { return x.vertex(i); }, py::arg("id"), "Copy of the vertex")
        .def("edge", [](const M &x, EdgeId i) { return x.edge(i); }, py::arg("id"), "Copy of the edge")
        .def("wire", [](const M &x, WireId i) { return x.wire(i); }, py::arg("id"), "Copy of the wire")
        .def("face", [](const M &x, FaceId i) { return x.face(i); }, py::arg("id"), "Copy of the face")
        .def("shell", [](const M &x, ShellId i) { return x.shell(i); }, py::arg("id"), "Copy of the shell")
        .def("solid", [](const M &x, SolidId i) { return x.solid(i); }, py::arg("id"), "Copy of the solid")
        .def("compound", [](const M &x, CompoundId i) { return x.compound(i); }, py::arg("id"), "Copy of the compound")
        .def("set_vertex", [](M &x, VertexId i, const Vertex<T> &v) { x.vertex(i) = v; }, py::arg("id"), py::arg("vertex"))
        .def("set_edge", [](M &x, EdgeId i, const Edge<T> &v) { x.edge(i) = v; }, py::arg("id"), py::arg("edge"))
        .def("set_wire", [](M &x, WireId i, const Wire<T> &v) { x.wire(i) = v; }, py::arg("id"), py::arg("wire"))
        .def("set_face", [](M &x, FaceId i, const Face<T> &v) { x.face(i) = v; }, py::arg("id"), py::arg("face"))
        .def("set_shell", [](M &x, ShellId i, const Shell &v) { x.shell(i) = v; }, py::arg("id"), py::arg("shell"))
        .def("set_solid", [](M &x, SolidId i, const Solid &v) { x.solid(i) = v; }, py::arg("id"), py::arg("solid"))
        .def("set_compound", [](M &x, CompoundId i, const Compound &v) { x.compound(i) = v; }, py::arg("id"), py::arg("compound"))
        .def("set_name", [](M &x, const ShapeId &s, std::string n) { x.setName(s, std::move(n)); }, py::arg("shape"), py::arg("name"),
             "Names an entity (an empty name removes it)")
        .def("name", [](const M &x, const ShapeId &s) { return std::string(x.name(s)); }, py::arg("shape"))
        .def("set_external_id", [](M &x, const ShapeId &s, std::optional<std::int64_t> id) { x.setExternalId(s, id); },
             py::arg("shape"), py::arg("id"), "Identifier in an external source, e.g. a STEP #id; None removes it")
        .def("external_id", [](const M &x, const ShapeId &s) { return x.externalId(s); }, py::arg("shape"))
        .def_property("unit_scale", &M::unitScale, &M::setUnitScale, "Length of one model unit in millimetres")
        .def("ids", [](const M &x, ShapeType t) {
            std::vector<ShapeId> r;
            auto add = [&](auto ids) { r.insert(r.end(), ids.begin(), ids.end()); };
            switch (t)
            {
            case ShapeType::Vertex: add(x.template ids<VertexId>()); break;
            case ShapeType::Edge: add(x.template ids<EdgeId>()); break;
            case ShapeType::Wire: add(x.template ids<WireId>()); break;
            case ShapeType::Face: add(x.template ids<FaceId>()); break;
            case ShapeType::Shell: add(x.template ids<ShellId>()); break;
            case ShapeType::Solid: add(x.template ids<SolidId>()); break;
            case ShapeType::Compound: add(x.template ids<CompoundId>()); break;
            }
            return r; }, py::arg("type"), "Identifiers of the live entities of a type")
        .def("alive", [](const M &x, const ShapeId &s) { return x.alive(s); }, py::arg("shape"))
        .def("count", [](const M &x, ShapeType t) { return x.count(t); }, py::arg("type"))
        .def("empty", &M::empty)
        .def("erase", [](M &x, const ShapeId &s) { x.erase(s); }, py::arg("shape"))
        .def("compact", &M::compact)
        .def("append", &M::append, py::arg("other"))
        .def("__repr__", [](const M &mdl) {
            std::ostringstream os;
            os << "Model(vertices=" << mdl.count(ShapeType::Vertex) << ", edges=" << mdl.count(ShapeType::Edge)
               << ", wires=" << mdl.count(ShapeType::Wire) << ", faces=" << mdl.count(ShapeType::Face)
               << ", shells=" << mdl.count(ShapeType::Shell) << ", solids=" << mdl.count(ShapeType::Solid)
               << ", compounds=" << mdl.count(ShapeType::Compound) << ")";
            return os.str();
        });

    // ---- explorer and queries ------------------------------------------------
    declare_explore<VertexId>(m, "explore_vertices");
    declare_explore<EdgeId>(m, "explore_edges");
    declare_explore<WireId>(m, "explore_wires");
    declare_explore<FaceId>(m, "explore_faces");
    declare_explore<ShellId>(m, "explore_shells");
    declare_explore<SolidId>(m, "explore_solids");
    declare_explore<CompoundId>(m, "explore_compounds");

    py::class_<TopologyIndex<T>>(m, "TopologyIndex", "Upward adjacency under a root shape; invalidated by any change of the model")
        .def(py::init<const M &, const ShapeId &>(), py::arg("model"), py::arg("root"))
        .def("faces_of", [](const TopologyIndex<T> &x, EdgeId e) { auto s = x.faces_of(e); return std::vector<FaceId>(s.begin(), s.end()); }, py::arg("edge"))
        .def("edges_of", [](const TopologyIndex<T> &x, VertexId v) { auto s = x.edges_of(v); return std::vector<EdgeId>(s.begin(), s.end()); }, py::arg("vertex"))
        .def("coedges_of", [](const TopologyIndex<T> &x, EdgeId e) {
            std::vector<std::pair<WireId, std::uint32_t>> r;
            for (const auto &c : x.coedges_of(e))
                r.emplace_back(c.wire, c.index);
            return r; }, py::arg("edge"), "(wire, index in wire.coedges) of every co-edge of the edge")
        .def("face_of", &TopologyIndex<T>::face_of, py::arg("wire"))
        .def("shells_of", [](const TopologyIndex<T> &x, FaceId f) { auto s = x.shells_of(f); return std::vector<ShellId>(s.begin(), s.end()); }, py::arg("face"));

    m.def("is_closed", py::overload_cast<const M &, WireId>(&is_closed<T>), py::arg("model"), py::arg("wire"));
    m.def("is_closed", py::overload_cast<const M &, ShellId>(&is_closed<T>), py::arg("model"), py::arg("shell"));
    m.def("is_chained", &is_chained<T>, py::arg("model"), py::arg("wire"));
    m.def("is_manifold", &is_manifold<T>, py::arg("model"), py::arg("shell"));
    m.def("is_orientable", &is_orientable<T>, py::arg("model"), py::arg("shell"));
    m.def("free_edges", &free_edges<T>, py::arg("model"), py::arg("shell"));
    m.def("bounding_box", [](const M &mdl, const ShapeId &s, std::size_t n) {
        auto b = bounding_box(mdl, s, n);
        return std::make_pair(b.min, b.max); }, py::arg("model"), py::arg("shape"), py::arg("n_samples") = 10,
          "(min, max) corners, tolerances included");
    m.def("signed_volume", &signed_volume<T>, py::arg("model"), py::arg("shell"), py::arg("n") = 64);
    m.def("uv_signed_area", py::overload_cast<const M &, WireId, std::size_t>(&uv_signed_area<T>),
          py::arg("model"), py::arg("wire"), py::arg("n_per_coedge") = 16);

    py::class_<SurfaceClosure>(m, "SurfaceClosure")
        .def_readonly("closed_u", &SurfaceClosure::closed_u)
        .def_readonly("closed_v", &SurfaceClosure::closed_v)
        .def_readonly("degenerate_u1", &SurfaceClosure::degenerate_u1)
        .def_readonly("degenerate_u2", &SurfaceClosure::degenerate_u2)
        .def_readonly("degenerate_v1", &SurfaceClosure::degenerate_v1)
        .def_readonly("degenerate_v2", &SurfaceClosure::degenerate_v2);
    m.def("surface_closure", &surface_closure<T>, py::arg("surface"), py::arg("tol"), py::arg("n_samples") = 9);

    // ---- builders ------------------------------------------------------------
    using Crv3 = std::shared_ptr<Curve<T, 3>>;
    using Crv2 = std::shared_ptr<Curve<T, 2>>;
    using Srf = std::shared_ptr<Surface<T, 3>>;
    const T dtol = brep_default_tolerance<T>;

    m.def("make_vertex", [](M &mdl, const point<T, 3> &p, T tol) { return value_or_raise(make_vertex(mdl, p, tol)); },
          py::arg("model"), py::arg("point"), py::arg("tol") = dtol);
    m.def("make_edge", [](M &mdl, Crv3 c, T tol) { return value_or_raise(make_edge(mdl, std::move(c), tol)); },
          py::arg("model"), py::arg("curve"), py::arg("tol") = dtol, "Edge on the whole range of a bounded curve");
    m.def("make_edge", [](M &mdl, Crv3 c, T u1, T u2, T tol) { return value_or_raise(make_edge(mdl, std::move(c), u1, u2, tol)); },
          py::arg("model"), py::arg("curve"), py::arg("u1"), py::arg("u2"), py::arg("tol") = dtol);
    m.def("make_edge", [](M &mdl, Crv3 c, T u1, T u2, VertexId v1, VertexId v2, T tol) { return value_or_raise(make_edge(mdl, std::move(c), u1, u2, v1, v2, tol)); },
          py::arg("model"), py::arg("curve"), py::arg("u1"), py::arg("u2"), py::arg("v1"), py::arg("v2"), py::arg("tol") = dtol);
    m.def("make_edge", [](M &mdl, const point<T, 3> &p1, const point<T, 3> &p2, T tol) { return value_or_raise(make_edge(mdl, p1, p2, tol)); },
          py::arg("model"), py::arg("p1"), py::arg("p2"), py::arg("tol") = dtol, "Straight edge between two points");
    m.def("make_edge", [](M &mdl, VertexId v1, VertexId v2, T tol) { return value_or_raise(make_edge(mdl, v1, v2, tol)); },
          py::arg("model"), py::arg("v1"), py::arg("v2"), py::arg("tol") = dtol, "Straight edge between two vertices");
    m.def("make_degenerate_edge", [](M &mdl, VertexId v, T u1, T u2) { return value_or_raise(make_degenerate_edge(mdl, v, u1, u2)); },
          py::arg("model"), py::arg("vertex"), py::arg("u1") = 0., py::arg("u2") = 1.);
    m.def("make_wire", [](M &mdl, const std::vector<EdgeId> &e, T tol) { return value_or_raise(make_wire(mdl, e, tol)); },
          py::arg("model"), py::arg("edges"), py::arg("tol") = dtol, "Wire chaining edges given in any order and sense");
    m.def("make_wire_ordered",
          [](M &mdl, const std::vector<std::pair<EdgeId, Orientation>> &ce, T tol) {
              std::vector<OrientedEdge> oe;
              for (const auto &[e, o] : ce)
                  oe.push_back({e, o});
              return value_or_raise(make_wire_ordered(mdl, oe, tol));
          },
          py::arg("model"), py::arg("coedges"), py::arg("tol") = dtol,
          "Wire of (edge, orientation) given in order; an edge may be used twice in opposite senses (seam)");

    m.def("split_edge",
          [](M &mdl, EdgeId e, T u) {
              auto r = value_or_raise(split_edge(mdl, e, u));
              return std::make_tuple(r.first, r.second, r.vertex);
          },
          py::arg("model"), py::arg("edge"), py::arg("u"),
          "Splits an edge at its parameter u: returns (first, second, vertex); every wire using it is updated");

    py::class_<MakeFaceOptions<T>>(m, "MakeFaceOptions")
        .def(py::init([](T tol, T pcurve_tol, std::size_t degree, std::size_t n_min, std::size_t n_max) {
                 return MakeFaceOptions<T>{tol, pcurve_tol, degree, n_min, n_max}; }),
             py::arg("tol") = dtol, py::arg("pcurve_tol") = brep_pcurve_approx_tol<T>, py::arg("pcurve_degree") = 3,
             py::arg("n_samples_min") = 9, py::arg("n_samples_max") = 1025)
        .def_readwrite("tol", &MakeFaceOptions<T>::tol)
        .def_readwrite("pcurve_tol", &MakeFaceOptions<T>::pcurve_tol)
        .def_readwrite("pcurve_degree", &MakeFaceOptions<T>::pcurve_degree)
        .def_readwrite("n_samples_min", &MakeFaceOptions<T>::n_samples_min)
        .def_readwrite("n_samples_max", &MakeFaceOptions<T>::n_samples_max);

    m.def("make_face", [](M &mdl, Srf s, T tol) { return value_or_raise(make_face(mdl, std::move(s), tol)); },
          py::arg("model"), py::arg("surface"), py::arg("tol") = dtol, "Face on the whole parametric rectangle of the surface");
    m.def("make_face", [](M &mdl, Srf s, WireId outer, const std::vector<WireId> &holes, const MakeFaceOptions<T> &o) {
              return value_or_raise(make_face(mdl, std::move(s), outer, holes, o)); },
          py::arg("model"), py::arg("surface"), py::arg("outer"), py::arg("holes") = std::vector<WireId>{},
          py::arg("options") = MakeFaceOptions<T>{}, "Face bounded by free wires of 3D edges lying on the surface");
    m.def("make_face_use", [](M &mdl, Srf s, WireId outer, const std::vector<WireId> &holes, const MakeFaceOptions<T> &o) {
              return value_or_raise(make_face_use(mdl, std::move(s), outer, holes, o)); },
          py::arg("model"), py::arg("surface"), py::arg("outer"), py::arg("holes") = std::vector<WireId>{},
          py::arg("options") = MakeFaceOptions<T>{},
          "Face bounded by wires given in the sense of the face, as the FaceUse keeping those senses");
    m.def("make_face_from_pcurves", [](M &mdl, Srf s, const std::vector<Crv2> &outer, const std::vector<std::vector<Crv2>> &holes, const MakeFaceOptions<T> &o) {
              return value_or_raise(make_face(mdl, std::move(s), outer, holes, o)); },
          py::arg("model"), py::arg("surface"), py::arg("outer"), py::arg("holes") = std::vector<std::vector<Crv2>>{},
          py::arg("options") = MakeFaceOptions<T>{}, "Face bounded by loops of 2D curves in the surface's parametric space");

    py::class_<SewReport>(m, "SewReport")
        .def_readonly("shells", &SewReport::shells)
        .def_readonly("merged", &SewReport::merged)
        .def_readonly("free_edges", &SewReport::free_edges)
        .def_readonly("rejected", &SewReport::rejected)
        .def_readonly("ambiguous", &SewReport::ambiguous)
        .def_readonly("orientable", &SewReport::orientable);
    m.def("sew", [](M &mdl, const std::vector<FaceId> &faces, T tol, std::size_t n_samples) {
              return value_or_raise(sew(mdl, faces, SewOptions<T>{tol, n_samples})); },
          py::arg("model"), py::arg("faces"), py::arg("tol") = dtol, py::arg("n_samples") = 10,
          "Sews faces into consistently oriented shells (no edge cutting)");
    m.def("make_solid", [](M &mdl, ShellId outer, const std::vector<ShellId> &voids) { return value_or_raise(make_solid(mdl, outer, voids)); },
          py::arg("model"), py::arg("outer"), py::arg("voids") = std::vector<ShellId>{},
          "Solid from closed shells, outer turned outwards and cavities inwards");
    m.def("make_compound", [](M &mdl, const std::vector<ShapeId> &shapes) { return value_or_raise(make_compound(mdl, shapes)); },
          py::arg("model"), py::arg("shapes"));

    // ---- check ---------------------------------------------------------------
    py::enum_<Issue> issue(m, "Issue");
    for (int i = 0; i <= static_cast<int>(Issue::VoidNotInward); ++i)
    {
        std::string n = to_string(static_cast<Issue>(i));
        std::ranges::replace(n, ' ', '_');
        std::ranges::replace(n, '-', '_');
        issue.value(n.c_str(), static_cast<Issue>(i));
    }
    py::class_<CheckEntry>(m, "CheckEntry")
        .def_readonly("shape", &CheckEntry::shape)
        .def_readonly("issue", &CheckEntry::issue)
        .def_readonly("detail", &CheckEntry::detail)
        .def("__repr__", [](const CheckEntry &e) { return std::string(to_string(e.issue)) + " on #" + std::to_string(shape_index(e.shape)) + (e.detail.empty() ? "" : " (" + e.detail + ")"); });
    py::class_<CheckReport>(m, "CheckReport")
        .def_readonly("entries", &CheckReport::entries)
        .def("ok", &CheckReport::ok)
        .def("count", &CheckReport::count, py::arg("issue"))
        .def("has", &CheckReport::has, py::arg("issue"), py::arg("shape"))
        .def("__bool__", &CheckReport::ok);
    m.def("check", [](const M &mdl, const ShapeId &s, bool geometry, std::size_t n_samples) {
              return check(mdl, s, CheckOptions{.n_samples = n_samples, .geometry = geometry}); },
          py::arg("model"), py::arg("shape"), py::arg("geometry") = true, py::arg("n_samples") = 9,
          "Every violated invariant of the shape; never modifies the model");

    // ---- placed copies --------------------------------------------------------
    using P = Placement<T>;
    py::class_<P>(m, "Placement", "Rigid placement p -> R p + t")
        .def(py::init<>(), "Identity")
        .def_static("frame", &P::frame, py::arg("origin"), py::arg("x"), py::arg("y"), py::arg("z"),
                    "Placement of the frame (origin, x, y, z): local coordinates -> global ones")
        .def_readwrite("r", &P::r, "Rows of the rotation")
        .def_readwrite("t", &P::t, "Translation")
        .def("apply", &P::apply, py::arg("point"))
        .def("inverse", &P::inverse)
        .def("is_identity", &P::is_identity, py::arg("eps") = 1e-14)
        .def("__mul__", [](const P &a, const P &b) { return a * b; }, "a * b applies b first");
    py::class_<InstanceOf<T>>(m, "InstanceOf", "Provenance of an instance of a part: part key and placement")
        .def(py::init([](std::uint64_t part, const P &place) { return InstanceOf<T>{part, place}; }), py::arg("part"), py::arg("place") = P{})
        .def_readwrite("part", &InstanceOf<T>::part)
        .def_readwrite("place", &InstanceOf<T>::place);
    m.def("new_part", [](M &mdl) { return mdl.newPart(); }, py::arg("model"), "A new part key, unique in the model");
    m.def("set_instance", [](M &mdl, const ShapeId &s, std::optional<InstanceOf<T>> i) { mdl.setInstance(s, i); },
          py::arg("model"), py::arg("shape"), py::arg("instance"), "Records (or removes with None) the part a shape is an instance of");
    m.def("instance", [](const M &mdl, const ShapeId &s) { return mdl.instance(s); }, py::arg("model"), py::arg("shape"),
          "InstanceOf of the shape, or None");
    m.def("copy_shape", &copy_shape<T>, py::arg("model"), py::arg("shape"), py::arg("placement"),
          "Copy of a shape and its sub-shapes moved by a rigid placement (shared topology and pcurves kept)");

    // ---- STEP reader ----------------------------------------------------------
    using step::StepIssue;
    using step::StepReadReport;
    py::class_<StepIssue>(m, "StepIssue", "An entity of the STEP file that was not read")
        .def_readonly("id", &StepIssue::id, "#id of the entity")
        .def_readonly("type", &StepIssue::type)
        .def_readonly("message", &StepIssue::message)
        .def("__repr__", [](const StepIssue &i) { return "#" + std::to_string(i.id) + " " + i.type + ": " + i.message; });
    py::class_<StepReadReport>(m, "StepReadReport", "Report of a STEP reading")
        .def_readonly("schemas", &StepReadReport::schemas)
        .def_readonly("length_unit", &StepReadReport::length_unit)
        .def_readonly("angle_unit", &StepReadReport::angle_unit)
        .def_readonly("length_factor", &StepReadReport::length_factor, "factor applied to the lengths of the file")
        .def_readonly("file_tolerance", &StepReadReport::file_tolerance, "uncertainty of the file, in model units")
        .def_readonly("solids", &StepReadReport::solids)
        .def_readonly("shells", &StepReadReport::shells)
        .def_readonly("faces", &StepReadReport::faces)
        .def_readonly("edges", &StepReadReport::edges)
        .def_readonly("vertices", &StepReadReport::vertices)
        .def_readonly("unsupported", &StepReadReport::unsupported)
        .def_readonly("failed", &StepReadReport::failed)
        .def_readonly("max_tolerance", &StepReadReport::max_tolerance)
        .def_readonly("raised_tolerances", &StepReadReport::raised_tolerances)
        .def_readonly("sense_mismatches", &StepReadReport::sense_mismatches)
        .def_readonly("check", &StepReadReport::check)
        .def("__repr__", [](const StepReadReport &r) {
            std::ostringstream os;
            os << "StepReadReport(solids=" << r.solids << ", shells=" << r.shells << ", faces=" << r.faces << ", edges=" << r.edges
               << ", unit='" << r.length_unit << "', unsupported=" << r.unsupported.size() << ", failed=" << r.failed.size()
               << ", valid=" << (r.check.ok() ? "True" : "False") << ")";
            return os.str();
        });

    auto options = [](std::optional<T> target_unit_mm, bool strict, T max_tolerance, const MakeFaceOptions<T> &face) {
        return step::StepReadOptions{.target_unit_mm = target_unit_mm, .strict = strict, .face = face, .max_tolerance = max_tolerance};
    };
    m.def("read_step",
          [options](M &mdl, const std::string &path, std::optional<T> target_unit_mm, bool strict, T max_tolerance, const MakeFaceOptions<T> &face) {
              auto r = value_or_raise(step::read_step_file(mdl, path, options(target_unit_mm, strict, max_tolerance, face)));
              return std::make_pair(r.root, std::move(r.report));
          },
          py::arg("model"), py::arg("path"), py::arg("target_unit_mm") = std::nullopt, py::arg("strict") = false,
          py::arg("max_tolerance") = 1e-2, py::arg("face_options") = MakeFaceOptions<T>{},
          "Reads the BREP shapes of a STEP file into the model; returns (root, report). Lengths are converted "
          "to the model unit (target_unit_mm sets it on an empty model). Raises BRepError if the file cannot be "
          "read, or in strict mode on the first face or entity that cannot be read (the model is then unchanged).");
    m.def("read_step_text",
          [options](M &mdl, const std::string &text, std::optional<T> target_unit_mm, bool strict, T max_tolerance, const MakeFaceOptions<T> &face) {
              auto r = value_or_raise(step::read_step(mdl, text, options(target_unit_mm, strict, max_tolerance, face)));
              return std::make_pair(r.root, std::move(r.report));
          },
          py::arg("model"), py::arg("text"), py::arg("target_unit_mm") = std::nullopt, py::arg("strict") = false,
          py::arg("max_tolerance") = 1e-2, py::arg("face_options") = MakeFaceOptions<T>{},
          "Same as read_step, from the content of a STEP file.");

    // ---- STEP writer ----------------------------------------------------------
    using step::StepWriteReport;
    py::class_<StepWriteReport>(m, "StepWriteReport", "Report of a STEP writing")
        .def_readonly("entities", &StepWriteReport::entities)
        .def_readonly("parts", &StepWriteReport::parts)
        .def_readonly("assemblies", &StepWriteReport::assemblies)
        .def_readonly("occurrences", &StepWriteReport::occurrences)
        .def_readonly("solids", &StepWriteReport::solids)
        .def_readonly("shells", &StepWriteReport::shells)
        .def_readonly("faces", &StepWriteReport::faces)
        .def_readonly("edges", &StepWriteReport::edges)
        .def_readonly("vertices", &StepWriteReport::vertices)
        .def_readonly("approximated_curves", &StepWriteReport::approximated_curves)
        .def_readonly("approximated_surfaces", &StepWriteReport::approximated_surfaces)
        .def_readonly("max_deviation", &StepWriteReport::max_deviation)
        .def_readonly("skipped", &StepWriteReport::skipped)
        .def("__repr__", [](const StepWriteReport &r) {
            std::ostringstream os;
            os << "StepWriteReport(parts=" << r.parts << ", assemblies=" << r.assemblies << ", solids=" << r.solids << ", faces=" << r.faces
               << ", approximated=" << r.approximated_curves + r.approximated_surfaces << ", skipped=" << r.skipped.size() << ")";
            return os.str();
        });

    auto write_options = [](const std::string &schema, T approx_tol, bool pcurves, bool flatten, std::optional<T> uncertainty,
                            const std::string &product_name, const std::string &author, const std::string &organization) {
        return step::StepWriteOptions{.schema = schema, .approx_tol = approx_tol, .pcurves = pcurves, .flatten = flatten,
                                      .uncertainty = uncertainty, .product_name = product_name, .author = author, .organization = organization};
    };
    // a shape or a list of shapes
    auto as_list = [](const py::object &shapes) {
        std::vector<ShapeId> v;
        if (py::isinstance<py::list>(shapes) || py::isinstance<py::tuple>(shapes))
            for (const auto &x : shapes)
                v.push_back(x.cast<ShapeId>());
        else
            v.push_back(shapes.cast<ShapeId>());
        return v;
    };
    m.def("write_step",
          [write_options, as_list](const M &mdl, const py::object &shapes, const std::string &path, const std::string &schema, T approx_tol, bool pcurves,
                                   bool flatten, std::optional<T> uncertainty, const std::string &product_name, const std::string &author,
                                   const std::string &organization) {
              const auto v = as_list(shapes);
              return value_or_raise(step::write_step_file(mdl, std::span<const ShapeId>{v}, path,
                                                          write_options(schema, approx_tol, pcurves, flatten, uncertainty, product_name, author, organization)));
          },
          py::arg("model"), py::arg("shapes"), py::arg("path"), py::arg("schema") = "AP214", py::arg("approx_tol") = 1e-5, py::arg("pcurves") = true,
          py::arg("flatten") = false, py::arg("uncertainty") = std::nullopt, py::arg("product_name") = "", py::arg("author") = "",
          py::arg("organization") = "",
          "Writes a shape or a list of shapes to a STEP file (solids, shells, faces as parts; compounds as assemblies sharing their parts); "
          "returns the report. Raises BRepError on a dead shape, an unknown schema or a file that cannot be written.");
    m.def("write_step_text",
          [write_options, as_list](const M &mdl, const py::object &shapes, const std::string &schema, T approx_tol, bool pcurves, bool flatten,
                                   std::optional<T> uncertainty, const std::string &product_name, const std::string &author, const std::string &organization) {
              const auto v = as_list(shapes);
              auto r = value_or_raise(step::write_step(mdl, std::span<const ShapeId>{v},
                                                       write_options(schema, approx_tol, pcurves, flatten, uncertainty, product_name, author, organization)));
              return std::make_pair(std::move(r.first), std::move(r.second));
          },
          py::arg("model"), py::arg("shapes"), py::arg("schema") = "AP214", py::arg("approx_tol") = 1e-5, py::arg("pcurves") = true,
          py::arg("flatten") = false, py::arg("uncertainty") = std::nullopt, py::arg("product_name") = "", py::arg("author") = "",
          py::arg("organization") = "", "Same as write_step, returning (text, report).");
}
