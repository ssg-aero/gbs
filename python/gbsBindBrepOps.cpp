// Python bindings of stage 2 of the native BREP kernel: sub-module gbs.intersect (curve-curve,
// curve-surface, surface-surface) and the operations of gbs.brep (classification, section,
// split, booleans, sweeps, offset). Design: docs/sources/design/brep_palier2.md; architecture
// note docs/sources/design/brep2_pr10_occt_python.md.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <gbs-brep/brep>
#include <gbs/intersect/curves.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/intersect/surface_surface.h>

#include <array>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace gbs;
using namespace gbs::brep;

namespace
{
    using T = double;
    using Crv3 = std::shared_ptr<Curve<T, 3>>;
    using Crv2 = std::shared_ptr<Curve<T, 2>>;
    using Srf = std::shared_ptr<Surface<T, 3>>;
    using M = Model<T>;

    template <typename V>
    V value_or_raise(BuildResult<V> r)
    {
        return unwrap(std::move(r));
    }

    IntersectOptions<T> options(T tol, T angular_tol, std::size_t max_pairs, bool closed_forms)
    {
        IntersectOptions<T> o;
        o.tol = tol;
        o.angular_tol = angular_tol;
        o.max_pairs = max_pairs;
        o.closed_forms = closed_forms;
        return o;
    }

    template <std::size_t dim>
    void declare_curve_curve(py::module_ &m, const char *suffix)
    {
        using Hit = CurveCurveHit<T, dim>;
        using Res = CurveCurveResult<T, dim>;
        py::class_<Hit>(m, (std::string("CurveCurveHit") + suffix).c_str(), "A crossing of two curves: c1(t1) = c2(t2) = p")
            .def_readonly("t1", &Hit::t1)
            .def_readonly("t2", &Hit::t2)
            .def_readonly("p", &Hit::p)
            .def_readonly("tangent", &Hit::tangent)
            .def("__repr__", [](const Hit &h) { return "CurveCurveHit(t1=" + std::to_string(h.t1) + ", t2=" + std::to_string(h.t2) + ")"; });
        py::class_<Res>(m, (std::string("CurveCurveResult") + suffix).c_str())
            .def_readonly("hits", &Res::hits)
            .def_readonly("overlaps", &Res::overlaps)
            .def_readonly("complete", &Res::complete);
        m.def(
            "curves",
            [](const std::shared_ptr<Curve<T, dim>> &c1, const std::shared_ptr<Curve<T, dim>> &c2, std::optional<std::array<T, 2>> r1,
               std::optional<std::array<T, 2>> r2, T tol, T angular_tol, std::size_t max_pairs) {
                const auto b1 = r1.value_or(c1->bounds()), b2 = r2.value_or(c2->bounds());
                return intersect(*c1, b1[0], b1[1], *c2, b2[0], b2[1], options(tol, angular_tol, max_pairs, true));
            },
            py::arg("c1"), py::arg("c2"), py::arg("range1") = std::nullopt, py::arg("range2") = std::nullopt, py::arg("tol") = extrema_tol<T>,
            py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000,
            "Every crossing and common portion of two curves (2D or 3D) on their ranges (their bounds by default)");
    }
}

void gbs_bind_brep_ops(py::module &parent)
{
    // ---------------------------------------------------------------- gbs.intersect
    auto mi = parent.def_submodule("intersect", "Intersections of curves and surfaces: every solution on bounded ranges");
    py::class_<CurveOverlap<T>>(mi, "CurveOverlap", "c1 on [t1a, t1b] lies on c2 from t2a to t2b")
        .def_readonly("t1a", &CurveOverlap<T>::t1a)
        .def_readonly("t1b", &CurveOverlap<T>::t1b)
        .def_readonly("t2a", &CurveOverlap<T>::t2a)
        .def_readonly("t2b", &CurveOverlap<T>::t2b);
    declare_curve_curve<2>(mi, "2d");
    declare_curve_curve<3>(mi, "3d");

    py::class_<CurveSurfaceHit<T>>(mi, "CurveSurfaceHit", "A crossing of a curve and a surface: c(t) = s(u, v) = p")
        .def_readonly("t", &CurveSurfaceHit<T>::t)
        .def_readonly("u", &CurveSurfaceHit<T>::u)
        .def_readonly("v", &CurveSurfaceHit<T>::v)
        .def_readonly("p", &CurveSurfaceHit<T>::p)
        .def_readonly("tangent", &CurveSurfaceHit<T>::tangent);
    py::class_<CurveSurfacePortion<T>>(mi, "CurveSurfacePortion", "A portion of the curve lying on the surface")
        .def_readonly("ta", &CurveSurfacePortion<T>::ta)
        .def_readonly("tb", &CurveSurfacePortion<T>::tb)
        .def_readonly("ua", &CurveSurfacePortion<T>::ua)
        .def_readonly("va", &CurveSurfacePortion<T>::va)
        .def_readonly("ub", &CurveSurfacePortion<T>::ub)
        .def_readonly("vb", &CurveSurfacePortion<T>::vb);
    py::class_<CurveSurfaceResult<T>>(mi, "CurveSurfaceResult")
        .def_readonly("hits", &CurveSurfaceResult<T>::hits)
        .def_readonly("portions", &CurveSurfaceResult<T>::portions)
        .def_readonly("complete", &CurveSurfaceResult<T>::complete);
    mi.def(
        "curve_surface",
        [](const Crv3 &c, const Srf &s, std::optional<std::array<T, 2>> range, std::optional<std::array<T, 4>> box, T tol, T angular_tol,
           std::size_t max_pairs) {
            const auto r = range.value_or(c->bounds());
            return intersect(*c, r[0], r[1], *s, box.value_or(s->bounds()), options(tol, angular_tol, max_pairs, true));
        },
        py::arg("curve"), py::arg("surface"), py::arg("range") = std::nullopt, py::arg("box") = std::nullopt, py::arg("tol") = extrema_tol<T>,
        py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000,
        "Every crossing of a curve with a surface and every portion of the curve lying on it");

    py::class_<IntersectionCurve<T>>(mi, "IntersectionCurve", "A curve of the intersection of two surfaces, with its pcurves on both")
        .def_readonly("curve", &IntersectionCurve<T>::curve)
        .def_readonly("a", &IntersectionCurve<T>::a)
        .def_readonly("b", &IntersectionCurve<T>::b)
        .def_readonly("pc1", &IntersectionCurve<T>::pc1)
        .def_readonly("pc2", &IntersectionCurve<T>::pc2)
        .def_readonly("deviation", &IntersectionCurve<T>::deviation)
        .def_readonly("exact", &IntersectionCurve<T>::exact)
        .def_readonly("closed", &IntersectionCurve<T>::closed)
        .def_readonly("tangent", &IntersectionCurve<T>::tangent);
    py::class_<SurfaceSurfaceHit<T>>(mi, "SurfaceSurfaceHit", "An isolated contact point of two surfaces")
        .def_readonly("u1", &SurfaceSurfaceHit<T>::u1)
        .def_readonly("v1", &SurfaceSurfaceHit<T>::v1)
        .def_readonly("u2", &SurfaceSurfaceHit<T>::u2)
        .def_readonly("v2", &SurfaceSurfaceHit<T>::v2)
        .def_readonly("p", &SurfaceSurfaceHit<T>::p);
    py::class_<SurfaceSurfaceResult<T>>(mi, "SurfaceSurfaceResult")
        .def_readonly("curves", &SurfaceSurfaceResult<T>::curves)
        .def_readonly("points", &SurfaceSurfaceResult<T>::points)
        .def_readonly("coincident", &SurfaceSurfaceResult<T>::coincident)
        .def_readonly("singular", &SurfaceSurfaceResult<T>::singular)
        .def_readonly("complete", &SurfaceSurfaceResult<T>::complete);
    mi.def(
        "surfaces",
        [](const Srf &s1, const Srf &s2, std::optional<std::array<T, 4>> box1, std::optional<std::array<T, 4>> box2, T tol, T angular_tol,
           std::size_t max_pairs, bool closed_forms) {
            return intersect(*s1, box1.value_or(s1->bounds()), *s2, box2.value_or(s2->bounds()), options(tol, angular_tol, max_pairs, closed_forms));
        },
        py::arg("s1"), py::arg("s2"), py::arg("box1") = std::nullopt, py::arg("box2") = std::nullopt, py::arg("tol") = extrema_tol<T>,
        py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000, py::arg("closed_forms") = true,
        "Intersection of two surfaces on their parameter boxes: closed forms for pairs of primitives, marching otherwise");

    // ---------------------------------------------------------------- gbs.brep: stage 2
    py::module_ m = parent.attr("brep");

    py::enum_<Location>(m, "Location")
        .value("In", Location::In)
        .value("Out", Location::Out)
        .value("On", Location::On);
    m.def("classify", [](const M &mdl, FaceId f, const point<T, 2> &uv, T tol) { return classify(mdl, f, uv, tol); }, py::arg("model"),
          py::arg("face"), py::arg("uv"), py::arg("tol") = brep_default_tolerance<T>, "Position of the parameters uv against the face");
    m.def("classify", [](const M &mdl, SolidId s, const point<T, 3> &p, T tol) { return classify(mdl, s, p, tol); }, py::arg("model"),
          py::arg("solid"), py::arg("point"), py::arg("tol") = brep_default_tolerance<T>, "Position of the point against the solid");
    m.def("classify_face", [](const M &mdl, FaceId f, SolidId s, T tol) { return classify(mdl, f, s, tol); }, py::arg("model"),
          py::arg("face"), py::arg("solid"), py::arg("tol") = brep_default_tolerance<T>,
          "Position of the face against the solid, by an interior point (None for a face without inside)");
    m.def("interior_point", [](const M &mdl, FaceId f, T tol) { return interior_point(mdl, f, tol); }, py::arg("model"), py::arg("face"),
          py::arg("tol") = brep_default_tolerance<T>, "Parameters of a point inside the face, away from its boundary");

    py::class_<SectionEdge<T>>(m, "SectionEdge", "An edge of a section and the faces carrying it")
        .def_readonly("edge", &SectionEdge<T>::edge)
        .def_readonly("face_a", &SectionEdge<T>::face_a)
        .def_readonly("face_b", &SectionEdge<T>::face_b)
        .def_readonly("pcurve_a", &SectionEdge<T>::pcurve_a)
        .def_readonly("pcurve_b", &SectionEdge<T>::pcurve_b)
        .def_readonly("tangent", &SectionEdge<T>::tangent);
    py::class_<SectionResult<T>>(m, "SectionResult")
        .def_readonly("compound", &SectionResult<T>::compound)
        .def_readonly("wires", &SectionResult<T>::wires)
        .def_readonly("edges", &SectionResult<T>::edges)
        .def_readonly("boundary_edges", &SectionResult<T>::boundary_edges)
        .def_readonly("points", &SectionResult<T>::points)
        .def_readonly("coincident", &SectionResult<T>::coincident)
        .def_readonly("singular", &SectionResult<T>::singular)
        .def_readonly("complete", &SectionResult<T>::complete);
    m.def("section", [](M &mdl, const ShapeId &a, const ShapeId &b, T tol) { return value_or_raise(section(mdl, a, b, SectionOptions<T>{tol})); },
          py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>,
          "Section of two shapes: the edges where their faces meet, chained into wires");
    m.def(
        "section",
        [](M &mdl, const ShapeId &a, const Srf &tool, std::optional<std::array<T, 4>> box, T tol) {
            return value_or_raise(section(mdl, a, tool, box, SectionOptions<T>{tol}));
        },
        py::arg("model"), py::arg("shape"), py::arg("surface"), py::arg("box") = std::nullopt, py::arg("tol") = brep_default_tolerance<T>,
        "Section of a shape by a surface (on its parameter box)");

    py::class_<SplitFaceResult>(m, "SplitFaceResult")
        .def_readonly("faces", &SplitFaceResult::faces)
        .def_readonly("dangling", &SplitFaceResult::dangling);
    m.def(
        "split_face",
        [](M &mdl, FaceId f, const std::vector<std::pair<EdgeId, Crv2>> &cuts, T tol) {
            std::vector<FaceCut<T>> c;
            for (const auto &[e, pc] : cuts)
                c.push_back({e, pc});
            return value_or_raise(split_face(mdl, f, c, tol));
        },
        py::arg("model"), py::arg("face"), py::arg("cuts"), py::arg("tol") = brep_default_tolerance<T>,
        "Splits a face by edges crossing it, given as (edge, pcurve on the face)");
    py::class_<SplitResult<T>>(m, "SplitResult")
        .def_readonly("faces", &SplitResult<T>::faces)
        .def_readonly("dangling", &SplitResult<T>::dangling)
        .def_readonly("section", &SplitResult<T>::section);
    m.def("split", [](M &mdl, const ShapeId &s, const ShapeId &tool, T tol) { return value_or_raise(split(mdl, s, tool, SectionOptions<T>{tol})); },
          py::arg("model"), py::arg("shape"), py::arg("tool"), py::arg("tol") = brep_default_tolerance<T>,
          "Splits the faces of a shape along its section by another shape");
    m.def(
        "split",
        [](M &mdl, const ShapeId &s, const Srf &tool, std::optional<std::array<T, 4>> box, T tol) {
            return value_or_raise(split(mdl, s, tool, box, SectionOptions<T>{tol}));
        },
        py::arg("model"), py::arg("shape"), py::arg("surface"), py::arg("box") = std::nullopt, py::arg("tol") = brep_default_tolerance<T>,
        "Splits the faces of a shape along its section by a surface");

    py::enum_<BooleanOp>(m, "BooleanOp")
        .value("Fuse", BooleanOp::Fuse)
        .value("Cut", BooleanOp::Cut)
        .value("Common", BooleanOp::Common);
    py::enum_<FaceOrigin>(m, "FaceOrigin")
        .value("A", FaceOrigin::A)
        .value("B", FaceOrigin::B);
    py::class_<BooleanResult<T>>(m, "BooleanResult")
        .def_readonly("shape", &BooleanResult<T>::shape)
        .def_readonly("solids", &BooleanResult<T>::solids)
        .def_readonly("faces", &BooleanResult<T>::faces);
    auto bool_op = [](BooleanOp op) {
        return [op](M &mdl, const ShapeId &a, const ShapeId &b, T tol) {
            return value_or_raise(boolean(mdl, op, a, b, BooleanOptions<T>{tol}));
        };
    };
    m.def("boolean", [](M &mdl, BooleanOp op, const ShapeId &a, const ShapeId &b, T tol) { return value_or_raise(boolean(mdl, op, a, b, BooleanOptions<T>{tol})); },
          py::arg("model"), py::arg("op"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>,
          "Boolean operation on solids; raises BRepError (non-generic configuration) on tangencies and singular intersections");
    m.def("fuse", bool_op(BooleanOp::Fuse), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, "Union of solids");
    m.def("cut", bool_op(BooleanOp::Cut), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, "a minus b");
    m.def("common", bool_op(BooleanOp::Common), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>,
          "Intersection of solids");

    py::class_<SweepResult<T>>(m, "SweepResult")
        .def_readonly("shape", &SweepResult<T>::shape)
        .def_readonly("approximated", &SweepResult<T>::approximated);
    m.def("make_prism",
          [](M &mdl, const ShapeId &profile, const point<T, 3> &d, T tol) { return value_or_raise(make_prism(mdl, profile, d, SweepOptions<T>{tol})); },
          py::arg("model"), py::arg("profile"), py::arg("direction"), py::arg("tol") = brep_default_tolerance<T>,
          "Linear extrusion: an edge gives a face, a wire a shell, a face a solid");
    m.def(
        "make_revolution",
        [](M &mdl, const ShapeId &profile, const point<T, 3> &origin, const point<T, 3> &direction, T angle, T tol) {
            return value_or_raise(make_revolution(mdl, profile, ax1<T, 3>{{origin, direction}}, angle, SweepOptions<T>{tol}));
        },
        py::arg("model"), py::arg("profile"), py::arg("origin"), py::arg("direction"), py::arg("angle") = 2. * std::numbers::pi,
        py::arg("tol") = brep_default_tolerance<T>, "Revolution around an axis by an angle in (0, 2 pi]");
    m.def("offset_face", [](M &mdl, FaceId f, T d, T tol) { return value_or_raise(offset_face(mdl, f, d, tol)); }, py::arg("model"),
          py::arg("face"), py::arg("distance"), py::arg("tol") = brep_default_tolerance<T>,
          "The face at a distance along the normal of its surface (same parametrization, pcurves reused)");
}
