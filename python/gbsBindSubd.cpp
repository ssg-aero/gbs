// Python bindings of the Catmull-Clark subdivision (gbs/subd), sub-module gbs.subd.
// Design: docs/sources/design/subd_catmull_clark.md; architecture note docs/sources/design/subd_pr01_cage.md.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gbs/subd/limit.h>
#include "np_points.h"

#include <exception>
#include <stdexcept>
#include <vector>

namespace py = pybind11;
using namespace gbs;
using namespace gbs::subd;

namespace
{
    using T = double;
    using release_gil = py::call_guard<py::gil_scoped_release>;

    // Evaluates f at the points (patches[i], us[i], vs[i]): an (n, 3) array.
    template <typename F>
    auto bulk(const std::vector<std::size_t> &patches, const std::vector<T> &us, const std::vector<T> &vs, F &&f) -> py::array_t<T>
    {
        if (us.size() != patches.size() || vs.size() != patches.size())
            throw std::invalid_argument{"subdivision: patches, us and vs differ in length"};
        points_vector<T, 3> out(patches.size());
        {
            py::gil_scoped_release release;
            for (std::size_t i{}; i < patches.size(); ++i)
                out[i] = f(patches[i], us[i], vs[i]);
        }
        return points_to_numpy(out);
    }

    // gbs.subd.CageError, for the translator (a function pointer, without capture)
    py::handle cage_error_type;

    // A CageError raises gbs.subd.CageError (a ValueError) with its issue and the index at fault.
    void translate_cage_error(std::exception_ptr p)
    {
        try
        {
            if (p)
                std::rethrow_exception(p);
        }
        catch (const CageError &e)
        {
            const auto type = py::reinterpret_borrow<py::object>(cage_error_type);
            py::object err = type(e.what());
            err.attr("issue") = py::cast(e.issue());
            err.attr("index") = e.index() == none ? py::object(py::none()) : py::cast(e.index());
            PyErr_SetObject(type.ptr(), err.ptr());
        }
    }
}

void gbs_bind_subd(py::module &parent)
{
    auto m = parent.def_submodule("subd", "Catmull-Clark subdivision: control cage, refinement, exact limit surface");

    auto &error = py::register_exception<CageError>(m, "CageError", PyExc_ValueError);
    error.attr("__doc__") = "Refusal of a subdivision cage: issue (a CageIssue) and index (the face, vertex, crease or corner at fault, or None)";
    error.attr("issue") = py::none();
    error.attr("index") = py::none();
    cage_error_type = error;
    py::register_exception_translator(&translate_cage_error);

    py::enum_<CageIssue>(m, "CageIssue", "What makes a cage unfit for the subdivision")
        .value("NoFace", CageIssue::NoFace)
        .value("FaceTooSmall", CageIssue::FaceTooSmall)
        .value("VertexOutOfRange", CageIssue::VertexOutOfRange)
        .value("RepeatedVertex", CageIssue::RepeatedVertex)
        .value("NonFiniteVertex", CageIssue::NonFiniteVertex)
        .value("NonManifoldEdge", CageIssue::NonManifoldEdge)
        .value("InconsistentOrientation", CageIssue::InconsistentOrientation)
        .value("NonManifoldVertex", CageIssue::NonManifoldVertex)
        .value("IsolatedVertex", CageIssue::IsolatedVertex)
        .value("InteriorValenceTwo", CageIssue::InteriorValenceTwo)
        .value("CreaseNotAnEdge", CageIssue::CreaseNotAnEdge)
        .value("CornerOutOfRange", CageIssue::CornerOutOfRange);
    py::enum_<VertexRule>(m, "VertexRule", "Subdivision rule of a vertex, from its sharp edges")
        .value("Smooth", VertexRule::Smooth)
        .value("Dart", VertexRule::Dart)
        .value("Crease", VertexRule::Crease)
        .value("Corner", VertexRule::Corner);

    py::class_<Cage<T>>(m, "Cage", "Control cage: vertices, faces (counterclockwise seen from outside), creases, corners")
        .def(py::init<>())
        .def(py::init([](points_vector<T, 3> vertices, std::vector<std::vector<std::size_t>> faces,
                         std::vector<std::array<std::size_t, 2>> creases, std::vector<std::size_t> corners) {
                 return Cage<T>{std::move(vertices), std::move(faces), std::move(creases), std::move(corners)};
             }),
             py::arg("vertices"), py::arg("faces"), py::arg("creases") = std::vector<std::array<std::size_t, 2>>{},
             py::arg("corners") = std::vector<std::size_t>{})
        .def_readwrite("vertices", &Cage<T>::vertices)
        .def_readwrite("faces", &Cage<T>::faces)
        .def_readwrite("creases", &Cage<T>::creases, "Infinitely sharp edges, as pairs of vertices")
        .def_readwrite("corners", &Cage<T>::corners, "Vertices interpolated by the limit surface")
        .def("__repr__", [](const Cage<T> &c) {
            return "Cage(" + std::to_string(c.vertices.size()) + " vertices, " + std::to_string(c.faces.size()) + " faces, " +
                   std::to_string(c.creases.size()) + " creases, " + std::to_string(c.corners.size()) + " corners)";
        });

    py::class_<Topology>(m, "Topology", "Half-edges of a cage, its edges and the rule of each vertex (checked: CageError)")
        .def(py::init([](const Cage<T> &c) { return Topology{c}; }), py::arg("cage"))
        .def("n_vertices", &Topology::n_vertices)
        .def("n_faces", &Topology::n_faces)
        .def("n_edges", &Topology::n_edges)
        .def("face_size", &Topology::face_size, py::arg("face"))
        .def("edge_vertices", &Topology::edge_vertices, py::arg("edge"))
        .def("is_boundary_edge", &Topology::is_boundary_edge, py::arg("edge"))
        .def("is_crease", &Topology::is_crease, py::arg("edge"))
        .def("is_sharp", &Topology::is_sharp, py::arg("edge"), "Crease or boundary")
        .def("valence", &Topology::valence, py::arg("vertex"), "Number of edges")
        .def("face_count", &Topology::face_count, py::arg("vertex"))
        .def("is_boundary", &Topology::is_boundary, py::arg("vertex"))
        .def("is_corner_tagged", &Topology::is_corner_tagged, py::arg("vertex"))
        .def("rule", &Topology::rule, py::arg("vertex"));

    m.def("refine", py::overload_cast<const Cage<T> &, std::size_t>(&refine<T>),
          "Cage refined `levels` times by Catmull-Clark (OpenSubdiv's Catmark rules); vertex v stays v",
          py::arg("cage"), py::arg("levels") = 1);

    py::class_<LimitSurface<T>>(m, "LimitSurface",
                                "Limit surface of a Catmull-Clark cage, evaluated exactly: one patch per quadrangle, n per face of n sides, "
                                "each on [0, 1]^2")
        .def(py::init<Cage<T>>(), py::arg("cage"))
        .def("cage", &LimitSurface<T>::cage, py::return_value_policy::reference_internal)
        .def("topology", &LimitSurface<T>::topology, py::return_value_policy::reference_internal)
        .def("patch_count", &LimitSurface<T>::patch_count)
        .def("patch", &LimitSurface<T>::patch, "Patch of a face: corner of an n-gon (0 for a quadrangle)", py::arg("face"), py::arg("corner") = 0)
        .def("patch_face", &LimitSurface<T>::patch_face, "[face, corner of an n-gon] of a patch", py::arg("patch"))
        .def("value", &LimitSurface<T>::value,
             "Point of a patch at (u, v), or its derivative of orders (du, dv); ValueError for a derivative at an irregular vertex",
             py::arg("patch"), py::arg("u"), py::arg("v"), py::arg("du") = 0, py::arg("dv") = 0, release_gil())
        .def("normal", &LimitSurface<T>::normal,
             "Unit normal of a patch at (u, v); at an irregular vertex, that of the tangent plane on the side of the patch, "
             "ValueError where none is given",
             py::arg("patch"), py::arg("u"), py::arg("v"), release_gil())
        .def(
            "values",
            [](const LimitSurface<T> &s, const std::vector<std::size_t> &patches, const std::vector<T> &us, const std::vector<T> &vs,
               std::size_t du, std::size_t dv) {
                return bulk(patches, us, vs, [&](std::size_t p, T u, T v) { return s.value(p, u, v, du, dv); });
            },
            "value at the points (patches[i], us[i], vs[i]), without the GIL: an (n, 3) array", py::arg("patches"), py::arg("us"),
            py::arg("vs"), py::arg("du") = 0, py::arg("dv") = 0)
        .def(
            "normals",
            [](const LimitSurface<T> &s, const std::vector<std::size_t> &patches, const std::vector<T> &us, const std::vector<T> &vs) {
                return bulk(patches, us, vs, [&](std::size_t p, T u, T v) { return s.normal(p, u, v); });
            },
            "normal at the points (patches[i], us[i], vs[i]), without the GIL: an (n, 3) array", py::arg("patches"), py::arg("us"),
            py::arg("vs"))
        .def("vertex_limit", &LimitSurface<T>::vertex_limit, "Limit position of a vertex of the cage", py::arg("vertex"), release_gil());
}
