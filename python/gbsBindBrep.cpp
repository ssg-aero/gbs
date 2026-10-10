// Python bindings of the native BREP core (gbs-brep), sub-module gbs.brep.
// Design: docs/sources/design/brep_core.md § 8.2; architecture note
// docs/sources/design/brep_pr10_python.md.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <gbs-brep/brep>

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
        .def_readwrite("natural_bounds", &Face<T>::natural_bounds);
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
}
