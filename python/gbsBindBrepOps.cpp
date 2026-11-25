// Python bindings of stage 2 of the native BREP kernel: sub-module gbs.intersect (curve-curve,
// curve-surface, surface-surface) and the operations of gbs.brep (classification, section,
// split, booleans, sweeps, offset). Design: docs/sources/design/brep_palier2.md; architecture
// note docs/sources/design/brep2_pr10_occt_python.md.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <pybind11/functional.h>

#include <gbs-brep/brep>
#include <gbs-io/mesh_writer.h>
#include <gbs-io/obj.h>
#include "np_points.h"
#include <gbs/intersect/curves.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/intersect/surface_surface.h>
#include "gbsBindBrep.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace gbs;
using namespace gbs::brep;
using gbs_bind::release_gil;
using gbs_bind::value_or_raise;

namespace
{
    using T = double;
    using Crv3 = std::shared_ptr<Curve<T, 3>>;
    using Crv2 = std::shared_ptr<Curve<T, 2>>;
    using Srf = std::shared_ptr<Surface<T, 3>>;
    using M = Model<T>;

    // The rule of the faces a fillet between faces swallows (stage 8, PR 2, § 5)
    constexpr const char *FILLET_FACES_DOC = R"(The solid with a fillet of a radius between two faces, or two sets of faces joined tangentially: a new solid, the
solid given unchanged.

The fillet is that of the sharp edge where both sides meet once continued over the faces between them (the virtual
edge). Which faces it swallows:

- the strip: the faces between the sides, connected, each face of both sides next to it;
- it ends at end faces: a face that touches both sides is an end face (kept, continued to meet the virtual edge, as
  the end of a box), except
- a face under the common edge: where the sides already meet along an edge e, a face G is of the strip when (1) G is
  of neither side, (2) G touches both sides, (3) its edges with them end on a vertex of e. A face tapering to a point
  on e, the triangle closing a chamfer stopped on e: the virtual edge continues e, both fused into one edge, one
  fillet in one piece. A common edge from the strip to the strip (the rest of the foot of a boss chamfered over
  its seam) lies under the virtual edge from end to end.

Faces that meet along edges with no strip running on from them: the fillet of those edges.

Refused (BRepError) when the fillet spills over the strip: its contacts off the sides, a face of the strip cut by it,
a point of the strip not taken off (or not filled in, for a concave fillet). A larger radius swallows it: in a U of a
wall and a floor at a right angle (a concave corner, the fillet filling it),
the arc of the fillet is r (sqrt 2 - 1) from the corner, so a face tapering to a point on the corner, w wide at its
widest point from it, is swallowed from r = w / (sqrt 2 - 1) on (w = 0.085: refused at r = 0.2, swallowed at 0.3).

Refused as well: strips on both sides of a common edge; sides that meet, once continued, on the far side of the
strip (the flanks of a rib in undercut, which meet below it, inside the plate), or a strip in part beyond the
sides continued (a bump past them). Where several strips are tried, the refusal is that of the one that went
furthest.

fillet of an edge does not cross a strip: where one runs on from the edge, fillet is refused at its end and names
fillet_faces.)";

    // A full round fillet (stage 8, PR 3)
    constexpr const char *FULL_ROUND_DOC = R"(The solid with its face `middle` replaced by the full round fillet tangent to it and to its neighbours
`faces1` and `faces2`: a new solid, the solid given unchanged.

No radius is given: a ball touching the three faces rolls along the middle face, its radius following them (half the
width of a rib with parallel flanks, (w / 2) cos(a) / (1 - sin(a)) with a draft a). The middle face is taken away, the
sides are cut along the contacts of the ball, the faces at the ends of the middle face along the curve where they meet
the fillet. The surface is a cylinder along a straight edge (a straight rib, with or without draft), a cone between
three planes meeting at a point (a rib whose top narrows), a rational B-spline within approx_tol otherwise (the blunt
trailing edge of a blade). Along a middle face concave with its sides (the bottom of a slot) the fillet adds material.

One face on each side and one in the middle; the same face on both sides (a blade of one patch) is split between its
two edges with the middle face and merged again. Refused (BRepError): a middle face that has no edge with a side; no
ball touching the three faces, a contact out of its face (a rib lower than the radius: fillet its two edges instead);
an end face that does not cut the fillet; a middle face closing on itself, with holes, narrowing to a point, with
several faces at an end or ending on a face tangent to a side (a rounded end), all not handled; a conic or G2 section.)";

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
            py::arg("c1").none(false), py::arg("c2").none(false), py::arg("range1") = std::nullopt, py::arg("range2") = std::nullopt,
            py::arg("tol") = extrema_tol<T>, py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000,
            "Every crossing and common portion of two curves (2D or 3D) on their ranges (their bounds by default)", release_gil());
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
        py::arg("curve").none(false), py::arg("surface").none(false), py::arg("range") = std::nullopt, py::arg("box") = std::nullopt,
        py::arg("tol") = extrema_tol<T>, py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000,
        "Every crossing of a curve with a surface and every portion of the curve lying on it", release_gil());

    py::class_<IntersectionCurve<T>>(mi, "IntersectionCurve", "A curve of the intersection of two surfaces, with its pcurves on both")
        .def_readonly("curve", &IntersectionCurve<T>::curve)
        .def_readonly("a", &IntersectionCurve<T>::a)
        .def_readonly("b", &IntersectionCurve<T>::b)
        .def_readonly("pc1", &IntersectionCurve<T>::pc1)
        .def_readonly("pc2", &IntersectionCurve<T>::pc2)
        .def_readonly("deviation", &IntersectionCurve<T>::deviation)
        .def_readonly("exact", &IntersectionCurve<T>::exact)
        .def_readonly("closed", &IntersectionCurve<T>::closed)
        .def_readonly("tangent", &IntersectionCurve<T>::tangent)
        .def_readonly("branch_spans", &IntersectionCurve<T>::branch_spans, "parameter intervals stepping over a branch point");
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
        py::arg("s1").none(false), py::arg("s2").none(false), py::arg("box1") = std::nullopt, py::arg("box2") = std::nullopt,
        py::arg("tol") = extrema_tol<T>, py::arg("angular_tol") = 1e-3, py::arg("max_pairs") = 200000, py::arg("closed_forms") = true,
        "Intersection of two surfaces on their parameter boxes: closed forms for pairs of primitives, marching otherwise", release_gil());

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
        .def_readonly("edges_on_faces", &SectionResult<T>::edges_on_faces, "(edge of one shape, face of the other it lies on)")
        .def_readonly("points", &SectionResult<T>::points)
        .def_readonly("coincident", &SectionResult<T>::coincident)
        .def_readonly("singular", &SectionResult<T>::singular)
        .def_readonly("complete", &SectionResult<T>::complete);
    m.def("section", [](M &mdl, const ShapeId &a, const ShapeId &b, T tol) { return value_or_raise(section(mdl, a, b, SectionOptions<T>{tol})); },
          py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>,
          "Section of two shapes: the edges where their faces meet, chained into wires", release_gil());
    m.def(
        "section",
        [](M &mdl, const ShapeId &a, const Srf &tool, std::optional<std::array<T, 4>> box, T tol) {
            return value_or_raise(section(mdl, a, tool, box, SectionOptions<T>{tol}));
        },
        py::arg("model"), py::arg("shape"), py::arg("surface"), py::arg("box") = std::nullopt, py::arg("tol") = brep_default_tolerance<T>,
        "Section of a shape by a surface (on its parameter box)", release_gil());

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
        "Splits a face by edges crossing it, given as (edge, pcurve on the face)", release_gil());
    py::class_<SplitResult<T>>(m, "SplitResult")
        .def_readonly("faces", &SplitResult<T>::faces)
        .def_readonly("dangling", &SplitResult<T>::dangling)
        .def_readonly("section", &SplitResult<T>::section);
    m.def("split", [](M &mdl, const ShapeId &s, const ShapeId &tool, T tol) { return value_or_raise(split(mdl, s, tool, SectionOptions<T>{tol})); },
          py::arg("model"), py::arg("shape"), py::arg("tool"), py::arg("tol") = brep_default_tolerance<T>,
          "Splits the faces of a shape along its section by another shape", release_gil());
    m.def(
        "split",
        [](M &mdl, const ShapeId &s, const Srf &tool, std::optional<std::array<T, 4>> box, T tol) {
            return value_or_raise(split(mdl, s, tool, box, SectionOptions<T>{tol}));
        },
        py::arg("model"), py::arg("shape"), py::arg("surface"), py::arg("box") = std::nullopt, py::arg("tol") = brep_default_tolerance<T>,
        "Splits the faces of a shape along its section by a surface", release_gil());

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
    auto options = [](T tol, bool unify) {
        BooleanOptions<T> o{tol};
        o.unify = unify;
        return o;
    };
    auto bool_op = [options](BooleanOp op) {
        return [op, options](M &mdl, const ShapeId &a, const ShapeId &b, T tol, bool unify) {
            return value_or_raise(boolean(mdl, op, a, b, options(tol, unify)));
        };
    };
    m.def("boolean", [options](M &mdl, BooleanOp op, const ShapeId &a, const ShapeId &b, T tol, bool unify) { return value_or_raise(boolean(mdl, op, a, b, options(tol, unify))); },
          py::arg("model"), py::arg("op"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, py::arg("unify") = true,
          "Boolean operation on solids; raises BRepError (code BuildErrc.non_generic_configuration) on tangencies and singular intersections. "
          "With unify, adjacent faces of the result on one surface are merged into one face",
          release_gil());
    m.def("fuse", bool_op(BooleanOp::Fuse), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, py::arg("unify") = true,
          "Union of solids", release_gil());
    m.def("cut", bool_op(BooleanOp::Cut), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, py::arg("unify") = true,
          "a minus b", release_gil());
    m.def("common", bool_op(BooleanOp::Common), py::arg("model"), py::arg("a"), py::arg("b"), py::arg("tol") = brep_default_tolerance<T>, py::arg("unify") = true,
          "Intersection of solids", release_gil());

    py::class_<UnifyReport>(m, "UnifyReport", "What unify_faces did")
        .def_readonly("faces_merged", &UnifyReport::faces_merged, "pairs of faces merged into one")
        .def_readonly("edges_joined", &UnifyReport::edges_joined, "pairs of edge pieces joined into one edge")
        .def_readonly("faces_planed", &UnifyReport::faces_planed, "faces on a plane in B-spline form put on the Plane primitive")
        .def_readonly("poles_reframed", &UnifyReport::poles_reframed, "faces of a sphere holding a pole put on a turned sphere")
        .def_readonly("surfaces_joined", &UnifyReport::surfaces_joined, "pairs of faces put on the join of their B-spline surfaces")
        .def_readonly("history", &UnifyReport::history,
                      "[new face, first piece, second piece], in order; the second piece is invalid for a face only put on another surface");
    m.def(
        "unify_faces",
        [](M &mdl, std::vector<FaceUse> faces, T tol) {
            auto r = unify_faces(mdl, faces, tol);
            return std::pair{std::move(faces), std::move(r)};
        },
        py::arg("model"), py::arg("faces"), py::arg("tol") = brep_default_tolerance<T>,
        "Adjacent faces of one domain among `faces` (every face using their edges: a whole shell or solid) merged into one face, "
        "in the shells too; returns (the faces, UnifyReport)",
        release_gil());

    py::class_<SweepResult<T>>(m, "SweepResult")
        .def_readonly("shape", &SweepResult<T>::shape)
        .def_readonly("approximated", &SweepResult<T>::approximated)
        .def_readonly("deviation", &SweepResult<T>::deviation);
    m.def("make_prism",
          [](M &mdl, const ShapeId &profile, const point<T, 3> &d, T tol) { return value_or_raise(make_prism(mdl, profile, d, SweepOptions<T>{tol})); },
          py::arg("model"), py::arg("profile"), py::arg("direction"), py::arg("tol") = brep_default_tolerance<T>,
          "Linear extrusion: an edge gives a face, a wire a shell, a face a solid", release_gil());
    m.def(
        "make_revolution",
        [](M &mdl, const ShapeId &profile, const point<T, 3> &origin, const point<T, 3> &direction, T angle, T tol) {
            return value_or_raise(make_revolution(mdl, profile, ax1<T, 3>{{origin, direction}}, angle, SweepOptions<T>{tol}));
        },
        py::arg("model"), py::arg("profile"), py::arg("origin"), py::arg("direction"), py::arg("angle") = 2. * std::numbers::pi,
        py::arg("tol") = brep_default_tolerance<T>, "Revolution around an axis by an angle in (0, 2 pi]", release_gil());
    py::enum_<PipeFrame>(m, "PipeFrame", "Moving frame carrying the profile of a pipe")
        .value("RotationMinimizing", PipeFrame::RotationMinimizing)
        .value("Frenet", PipeFrame::Frenet)
        .value("FixedBinormal", PipeFrame::FixedBinormal);
    py::class_<PipeOptions<T>>(m, "PipeOptions", "Options of make_pipe")
        .def(py::init([](T tol, T approx_tol, PipeFrame frame, const point<T, 3> &binormal) {
                 return PipeOptions<T>{tol, approx_tol, frame, binormal};
             }),
             py::arg("tol") = brep_default_tolerance<T>, py::arg("approx_tol") = T(1e-5), py::arg("frame") = PipeFrame::RotationMinimizing,
             py::arg("binormal") = point<T, 3>{})
        .def_readwrite("tol", &PipeOptions<T>::tol)
        .def_readwrite("approx_tol", &PipeOptions<T>::approx_tol)
        .def_readwrite("frame", &PipeOptions<T>::frame)
        .def_readwrite("binormal", &PipeOptions<T>::binormal);
    m.def(
        "make_pipe",
        [](M &mdl, const ShapeId &profile, const ShapeId &path, const PipeOptions<T> &o) { return value_or_raise(make_pipe(mdl, profile, path, o)); },
        py::arg("model"), py::arg("profile"), py::arg("path"), py::arg("options") = PipeOptions<T>{},
        "Sweep of a profile along a G1 path (an edge or a wire), carried by a moving frame: an edge gives a face or a shell, a wire a "
        "shell, a face a solid (no caps along a closed path)",
        release_gil());
    m.def("offset_face", [](M &mdl, FaceId f, T d, T tol) { return value_or_raise(offset_face(mdl, f, d, tol)); }, py::arg("model"),
          py::arg("face"), py::arg("distance"), py::arg("tol") = brep_default_tolerance<T>,
          "The face at a distance along the normal of its surface (same parametrization, pcurves reused)", release_gil());
    py::class_<OffsetOptions<T>>(m, "OffsetOptions", "Options of offset_solid and make_thick")
        .def(py::init([](T tol, T max_extension) {
                 OffsetOptions<T> o;
                 o.tol = tol;
                 o.max_extension = max_extension;
                 return o;
             }),
             py::arg("tol") = brep_default_tolerance<T>, py::arg("max_extension") = T(100))
        .def_readwrite("tol", &OffsetOptions<T>::tol)
        .def_readwrite("max_extension", &OffsetOptions<T>::max_extension);
    m.def(
        "offset_solid", [](M &mdl, SolidId s, T d, const OffsetOptions<T> &o) { return value_or_raise(offset_solid(mdl, s, d, o)); }, py::arg("model"),
        py::arg("solid"), py::arg("distance"), py::arg("options") = OffsetOptions<T>{},
        "The solid offset by a distance along its outward normal (> 0 grows it), with sharp joints", release_gil());
    m.def(
        "make_thick",
        [](M &mdl, SolidId s, const std::vector<FaceId> &removed, T t, const OffsetOptions<T> &o) { return value_or_raise(make_thick(mdl, s, removed, t, o)); },
        py::arg("model"), py::arg("solid"), py::arg("removed"), py::arg("thickness"), py::arg("options") = OffsetOptions<T>{},
        "The solid hollowed out to a thickness (> 0 inward), opened by the removed faces (none: a closed hollow solid)", release_gil());

    py::class_<Setback<T>>(m, "Setback", "A setback of a fillet at a vertex blend: the fillet of the edge stops at the distance from the vertex along the edge")
        .def(py::init([](VertexId v, EdgeId e, T d) { return Setback<T>{v, e, d}; }), py::arg("vertex"), py::arg("edge"), py::arg("distance"))
        .def_readwrite("vertex", &Setback<T>::vertex)
        .def_readwrite("edge", &Setback<T>::edge)
        .def_readwrite("distance", &Setback<T>::distance);
    py::enum_<SectionKind>(m, "SectionKind", "The shape of the section of a fillet")
        .value("Circular", SectionKind::Circular)
        .value("Conic", SectionKind::Conic)
        .value("G2", SectionKind::G2);
    py::class_<BlendSection<T>>(m, "BlendSection", "The section of a fillet: the arc of circle, the conic of parameter rho (0.5 a parabola) "
                                                   "tangent to the faces at the contacts of the ball of the radius, or the quartic of continuous "
                                                   "curvature (g2: contacts spread for the sag of the arc, 1.305 r at right angles)")
        .def(py::init<>())
        .def_static("circular", &BlendSection<T>::circular)
        .def_static("conic", py::overload_cast<T>(&BlendSection<T>::conic), py::arg("rho"))
        .def_static("conic", py::overload_cast<const BlendLaw<T> &>(&BlendSection<T>::conic), py::arg("law"), "rho following a law along the edges")
        .def_static("g2", &BlendSection<T>::g2, "the quartic of continuous curvature with the faces (stage 7, rule C)")
        .def_readwrite("kind", &BlendSection<T>::kind)
        .def_readwrite("rho", &BlendSection<T>::rho)
        .def_readwrite("rho_law", &BlendSection<T>::rho_law, "rho following a law (None: constant); a copy, set by assigning a law (changing the copy changes nothing)");
    py::class_<FilletOptions<T>>(m, "FilletOptions", "Options of fillet, chamfer, fillet_geometry and chamfer_geometry")
        .def(py::init([](T tol, T approx_tol, T angular_tol, T extension) {
                 FilletOptions<T> o;
                 o.tol = tol;
                 o.approx_tol = approx_tol;
                 o.angular_tol = angular_tol;
                 o.extension = extension;
                 return o;
             }),
             py::arg("tol") = brep_default_tolerance<T>, py::arg("approx_tol") = T(1e-5), py::arg("angular_tol") = T(1e-5), py::arg("extension") = T(1))
        .def_readwrite("tol", &FilletOptions<T>::tol, "tolerance of the result; the contacts of a chamfer within tol / 10")
        .def_readwrite("approx_tol", &FilletOptions<T>::approx_tol, "approximated surfaces: largest distance from the exact blend")
        .def_readwrite("angular_tol", &FilletOptions<T>::angular_tol, "approximated fillets: sine of the angle of their normal to the exact one")
        .def_readwrite("extension", &FilletOptions<T>::extension, "length of the blend past the ends of an open edge, in radius or largest distance")
        .def_readwrite("section", &FilletOptions<T>::section, "the section of a fillet (BlendSection: circular by default, or conic)")
        // (a tuple: appending to a list converted from the vector would change a copy, silently)
        .def_property(
            "setbacks",
            [](const FilletOptions<T> &o) {
                py::tuple t(o.setbacks.size());
                for (std::size_t i = 0; i < o.setbacks.size(); ++i)
                    t[i] = py::cast(o.setbacks[i]);
                return t;
            },
            [](FilletOptions<T> &o, std::vector<Setback<T>> s) { o.setbacks = std::move(s); },
            "setbacks of fillets at vertex blends (fillet only), each at least where the contacts cross: a tuple of Setback, set by assigning a "
            "list of them");
    py::class_<BlendContact<T>>(m, "BlendContact", "A contact curve of a blend on a face of the edge")
        .def_readonly("face", &BlendContact<T>::face)
        .def_readonly("support", &BlendContact<T>::support, "the surface of the face, continued past the ends of the edge")
        .def_readonly("pcurve", &BlendContact<T>::pcurve, "on the support")
        .def_readonly("curve", &BlendContact<T>::curve, "in space")
        .def_readonly("on_blend", &BlendContact<T>::on_blend, "on the surface of the blend")
        .def_readonly("deviation", &BlendContact<T>::deviation);
    py::class_<EdgeBlend<T>>(m, "EdgeBlend", "The geometry of a fillet or a chamfer along an edge, every curve on [a, b]")
        .def_readonly("edge", &EdgeBlend<T>::edge)
        .def_readonly("convex", &EdgeBlend<T>::convex)
        .def_readonly("closed", &EdgeBlend<T>::closed)
        .def_readonly("surface", &EdgeBlend<T>::surface)
        .def_readonly("exact", &EdgeBlend<T>::exact)
        .def_readonly("deviation", &EdgeBlend<T>::deviation)
        .def_readonly("spine", &EdgeBlend<T>::spine, "the centre line of a fillet, the edge of a chamfer")
        .def_readonly("a", &EdgeBlend<T>::a)
        .def_readonly("b", &EdgeBlend<T>::b)
        .def_readonly("t1", &EdgeBlend<T>::t1, "parameter at the vertex v1 of the edge (a closed edge: an end of its loop, a or b)")
        .def_readonly("t2", &EdgeBlend<T>::t2, "parameter at the vertex v2 of the edge (t1 > t2: the spine runs the other way; a closed edge: the other end of its loop)")
        .def_readonly("variable", &EdgeBlend<T>::variable, "a fillet whose radius follows a law")
        .def_readonly("contacts", &EdgeBlend<T>::contacts, "on the faces of the edge (a chamfer: the face of d1 first)")
        .def_readonly("pieces", &EdgeBlend<T>::pieces, "a fillet whose radius law kinks: its surface in pieces (surface None)")
        .def_readonly("lunes", &EdgeBlend<T>::lunes, "the lunes of the ball between pieces (a hollow of the law)")
        .def_readonly("creases", &EdgeBlend<T>::creases, "the creases between pieces (a peak of the law)")
        .def_readonly("circular", &EdgeBlend<T>::circular, "sections arcs of circle (false: a conic section)");
    py::class_<BlendPiece<T>>(m, "BlendPiece", "A piece of the surface of a fillet in pieces, on [a, b]")
        .def_readonly("surface", &BlendPiece<T>::surface)
        .def_readonly("a", &BlendPiece<T>::a)
        .def_readonly("b", &BlendPiece<T>::b)
        .def_readonly("on_blend", &BlendPiece<T>::on_blend, "the contacts on the surface")
        .def_readonly("deviation", &BlendPiece<T>::deviation);
    py::class_<BlendLune<T>>(m, "BlendLune", "The lune of the ball between two pieces of a fillet, at t")
        .def_readonly("t", &BlendLune<T>::t)
        .def_readonly("sphere", &BlendLune<T>::sphere)
        .def_readonly("arcs", &BlendLune<T>::arcs, "the characteristic circles before and after t, from contact 0 to contact 1")
        .def_readonly("on_sphere", &BlendLune<T>::on_sphere);
    py::class_<BlendCrease<T>>(m, "BlendCrease", "The crease where two pieces of a fillet meet, at t")
        .def_readonly("t", &BlendCrease<T>::t)
        .def_readonly("curve", &BlendCrease<T>::curve, "from contact 0 to contact 1, on [0, 1]")
        .def_readonly("on_pieces", &BlendCrease<T>::on_pieces)
        .def_readonly("deviation", &BlendCrease<T>::deviation);
    py::enum_<BlendVariation>(m, "BlendVariation", "How a value of a blend (a radius, a distance) varies between two control points")
        .value("Linear", BlendVariation::Linear)
        .value("Smooth", BlendVariation::Smooth);
    m.attr("RadiusVariation") = m.attr("BlendVariation");
    py::class_<BlendPoint<T>>(m, "BlendPoint", "A value of a blend (a radius, a distance) at a point of an edge: at, the fraction of its length from its first vertex; "
                                               "or at a point given by position (BlendPoint.at_point), put at its foot on the nearest edge blended")
        .def(py::init([](EdgeId e, T at, T v) { return BlendPoint<T>{e, at, v}; }), py::arg("edge"), py::arg("at"), py::arg("value"))
        // (the keyword of stage 5, RadiusPoint(edge, at, radius=...))
        .def(py::init([](EdgeId e, T at, T r) { return BlendPoint<T>{e, at, r}; }), py::arg("edge"), py::arg("at"), py::arg("radius"))
        .def_readwrite("edge", &BlendPoint<T>::edge)
        .def_readwrite("at", &BlendPoint<T>::at)
        .def_readwrite("value", &BlendPoint<T>::value)
        .def_readwrite("radius", &BlendPoint<T>::value, "the value, for a radius")
        .def_readwrite("where", &BlendPoint<T>::where, "the point given by position (None: on `edge` at `at`)")
        .def_static("at_point", &BlendPoint<T>::at_point, py::arg("point"), py::arg("value"),
                    "The value at the point of the edges blended nearest `point` (its foot): on the edges given for a fillet of edges, on the "
                    "virtual edge for a fillet of faces. Refused: a point farther than ten times its value and the length of the edge, two "
                    "points at the same foot with different values");
    m.attr("RadiusPoint") = m.attr("BlendPoint");
    py::class_<BlendLaw<T>>(m, "BlendLaw", "A value of a blend (the radius of a fillet, a distance of a chamfer) along its edges: control points, a linear or "
                                           "smooth variation between them along each chain of tangent edges")
        .def(py::init([](const std::vector<BlendPoint<T>> &pts, BlendVariation v) { return BlendLaw<T>{pts, v}; }), py::arg("points"),
             py::arg("variation") = BlendVariation::Smooth)
        .def_readwrite("points", &BlendLaw<T>::points)
        .def_readwrite("variation", &BlendLaw<T>::variation)
        .def_static("linear", &BlendLaw<T>::linear, py::arg("edge"), py::arg("v1"), py::arg("v2"),
                    "From v1 at the first vertex of the edge to v2 at its second, linearly")
        // (the keywords of stage 5, RadiusLaw.linear(edge, r1=..., r2=...))
        .def_static("linear", &BlendLaw<T>::linear, py::arg("edge"), py::arg("r1"), py::arg("r2"), "As linear(edge, v1, v2), with the names of a radius")
        .def_static("constant", &BlendLaw<T>::constant, py::arg("edge"), py::arg("value"), "The same value along the edge");
    m.attr("RadiusLaw") = m.attr("BlendLaw");
    m.def(
        "fillet_geometry",
        [](const M &mdl, SolidId s, EdgeId e, T r, const FilletOptions<T> &o) { return value_or_raise(fillet_geometry(mdl, s, e, r, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("radius"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the fillet of a radius along an edge of a solid: centre line, contacts, surface (the model is not modified)",
        release_gil());
    m.def(
        "fillet_geometry",
        [](const M &mdl, SolidId s, EdgeId e, const RadiusLaw<T> &law, const FilletOptions<T> &o) { return value_or_raise(fillet_geometry(mdl, s, e, law, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("law"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the fillet along an edge whose radius follows a law (its points on the edge): the envelope of the balls, "
        "in pieces with lunes or creases where the law kinks (the model is not modified)",
        release_gil());
    m.def(
        "chamfer_geometry",
        [](const M &mdl, SolidId s, EdgeId e, FaceId f, T d1, T d2, const FilletOptions<T> &o) {
            return value_or_raise(chamfer_geometry(mdl, s, e, f, d1, d2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("face"), py::arg("d1"), py::arg("d2"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the chamfer along an edge of a solid, d1 on the face `face` of the edge, d2 on the other (the model is not "
        "modified)",
        release_gil());
    m.def(
        "chamfer_geometry",
        [](const M &mdl, SolidId s, EdgeId e, FaceId f, const BlendLaw<T> &d1, const BlendLaw<T> &d2, const FilletOptions<T> &o) {
            return value_or_raise(chamfer_geometry(mdl, s, e, f, d1, d2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("face"), py::arg("d1"), py::arg("d2"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the chamfer along an edge at distances following laws, d1 on `face`, d2 on the other: ruled, in pieces with "
        "creases where a law kinks (the model is not modified)",
        release_gil());
    m.def(
        "fillet_chordal_geometry",
        [](const M &mdl, SolidId s, EdgeId e, T c, const FilletOptions<T> &o) { return value_or_raise(fillet_chordal_geometry(mdl, s, e, c, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("chord"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the fillet of constant chord along an edge: its contacts that far apart at every station, its radius following "
        "(the model is not modified)",
        release_gil());
    m.def(
        "fillet_chordal",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, T c, const FilletOptions<T> &o) { return value_or_raise(fillet_chordal(mdl, s, edges, c, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("chord"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges filleted at a constant chord (the distance between the contacts): a new solid, the solid given unchanged",
        release_gil());
    m.def(
        "chamfer_angle_geometry",
        [](const M &mdl, SolidId s, EdgeId e, FaceId f, T d1, T angle, const FilletOptions<T> &o) {
            return value_or_raise(chamfer_angle_geometry(mdl, s, e, f, d1, angle, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edge"), py::arg("face"), py::arg("d1"), py::arg("angle"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the chamfer along an edge at the distance d1 on `face` and the angle (radians) from it (the model is not modified)",
        release_gil());
    m.def(
        "fillet",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, T r, const FilletOptions<T> &o) { return value_or_raise(fillet(mdl, s, edges, r, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("radius"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges filleted at a radius: a new solid, the solid given unchanged (edges filleted together meet "
        "along chains of tangent edges, or three at a vertex of three edges: a sphere at the corner)",
        release_gil());
    m.def(
        "fillet_faces",
        [](M &mdl, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &faces2, T r, const FilletOptions<T> &o) {
            return value_or_raise(fillet_faces(mdl, s, faces1, faces2, r, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("faces1"), py::arg("faces2"), py::arg("radius"), py::arg("options") = FilletOptions<T>{},
        FILLET_FACES_DOC, release_gil());
    m.def(
        "fillet_faces",
        [](M &mdl, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &faces2, const RadiusLaw<T> &law, const FilletOptions<T> &o) {
            return value_or_raise(fillet_faces(mdl, s, faces1, faces2, law, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("faces1"), py::arg("faces2"), py::arg("law"), py::arg("options") = FilletOptions<T>{},
        "fillet_faces at a radius following a law along the virtual edge: its points on edges of the solid or by position "
        "(BlendPoint.at_point), each put at its foot on the virtual edge",
        release_gil());
    m.def(
        "fillet_faces_chordal",
        [](M &mdl, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &faces2, T c, const FilletOptions<T> &o) {
            return value_or_raise(fillet_faces_chordal(mdl, s, faces1, faces2, c, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("faces1"), py::arg("faces2"), py::arg("chord"), py::arg("options") = FilletOptions<T>{},
        "fillet_faces at a constant chord (as fillet_chordal) along the virtual edge", release_gil());
    m.def(
        "fillet_full_round",
        [](M &mdl, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &middle, const std::vector<FaceId> &faces2, const FilletOptions<T> &o) {
            return value_or_raise(fillet_full_round(mdl, s, faces1, middle, faces2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("faces1"), py::arg("middle"), py::arg("faces2"), py::arg("options") = FilletOptions<T>{}, FULL_ROUND_DOC,
        release_gil());
    m.def(
        "full_round_geometry",
        [](const M &mdl, SolidId s, FaceId face1, FaceId middle, FaceId face2, const FilletOptions<T> &o) {
            return value_or_raise(full_round_geometry(mdl, s, face1, middle, face2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("face1"), py::arg("middle"), py::arg("face2"), py::arg("options") = FilletOptions<T>{},
        "The geometry of the full round fillet of the face `middle` between `face1` and `face2`: the surface, the centres of the balls (spine; the "
        "radius at t is its distance to a contact) and the contacts on both faces, on the parameter of the edge between `middle` and `face1` (the "
        "model is not modified)",
        release_gil());
    m.def(
        "fillet",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, const RadiusLaw<T> &law, const FilletOptions<T> &o) {
            return value_or_raise(fillet(mdl, s, edges, law, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("law"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges filleted at a radius following a law along each chain of tangent edges given: a new solid, the "
        "solid given unchanged",
        release_gil());
    m.def(
        "chamfer",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, const BlendLaw<T> &d, const FilletOptions<T> &o) {
            return value_or_raise(chamfer(mdl, s, edges, d, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("law"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered at the same distance on both faces, following a law along each chain of tangent edges",
        release_gil());
    m.def(
        "chamfer",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, FaceId f, const BlendLaw<T> &d1, const BlendLaw<T> &d2, const FilletOptions<T> &o) {
            return value_or_raise(chamfer(mdl, s, edges, f, d1, d2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("face"), py::arg("d1"), py::arg("d2"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered at distances following laws, d1 on `face`, d2 on the other face of each edge",
        release_gil());
    m.def(
        "chamfer_angle",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, FaceId f, T d1, T angle, const FilletOptions<T> &o) {
            return value_or_raise(chamfer_angle(mdl, s, edges, f, d1, angle, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("face"), py::arg("d1"), py::arg("angle"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered at the distance d1 on `face` and the angle (radians, within ]0, pi/2[) from it",
        release_gil());
    m.def(
        "chamfer_angle",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, FaceId f, const BlendLaw<T> &d1, T angle, const FilletOptions<T> &o) {
            return value_or_raise(chamfer_angle(mdl, s, edges, f, d1, angle, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("face"), py::arg("d1"), py::arg("angle"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered at a distance following a law on `face` and the angle (radians) from it",
        release_gil());
    m.def(
        "chamfer",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, T d, const FilletOptions<T> &o) { return value_or_raise(chamfer(mdl, s, edges, d, o)); },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("distance"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered at a distance on both faces: a new solid, the solid given unchanged", release_gil());
    m.def(
        "chamfer",
        [](M &mdl, SolidId s, const std::vector<EdgeId> &edges, FaceId f, T d1, T d2, const FilletOptions<T> &o) {
            return value_or_raise(chamfer(mdl, s, edges, f, d1, d2, o));
        },
        py::arg("model"), py::arg("solid"), py::arg("edges"), py::arg("face"), py::arg("d1"), py::arg("d2"), py::arg("options") = FilletOptions<T>{},
        "The solid with its edges chamfered, d1 on the face `face` (a face of every edge), d2 on the other: a new solid", release_gil());

    // ---------------------------------------------------------------- gbs.brep: stage 3 (meshing)
    py::class_<MeshOptions<T>>(m, "MeshOptions", "Precision of a mesh")
        .def(py::init([](T deflection, T angle, T max_edge, std::size_t max_points, T min_angle) {
                 return MeshOptions<T>{deflection, angle, max_edge, max_points, min_angle};
             }),
             py::arg("deflection") = 0., py::arg("angle") = 0.5, py::arg("max_edge") = 0., py::arg("max_points") = 200000,
             py::arg("min_angle") = MeshOptions<T>{}.min_angle)
        .def_readwrite("deflection", &MeshOptions<T>::deflection, "largest distance between the mesh and the geometry (0: 1e-3 of the diagonal of the exact box of the shape, exact_bounding_box)")
        .def_readwrite("angle", &MeshOptions<T>::angle, "largest turn of the tangent along a segment, of the normal across a triangle (radians)")
        .def_readwrite("max_edge", &MeshOptions<T>::max_edge, "longest segment or side of a triangle (0: none)")
        .def_readwrite("max_points", &MeshOptions<T>::max_points, "most points of a face")
        .def_readwrite("min_angle", &MeshOptions<T>::min_angle, "triangles under it refined, in the scaled parameter space (radians; 0: none)");
    using FaceRange = typename Mesh<T>::FaceRange;
    py::class_<FaceRange>(m, "MeshFaceRange", "The triangles of a face: triangles[first:first + count]")
        .def_readonly("face", &FaceRange::face)
        .def_readonly("first", &FaceRange::first)
        .def_readonly("count", &FaceRange::count);
    py::class_<Mesh<T>>(m, "Mesh", "A triangle mesh of a shape")
        .def_property_readonly("vertices", [](const Mesh<T> &x) { return points_to_numpy(x.vertices); }, "(n, 3) array of the vertices")
        .def_property_readonly("triangles", [](const Mesh<T> &x) { return points_to_numpy(x.triangles); },
                               "(m, 3) array of the vertices of the triangles, counter-clockwise seen from outside")
        .def_property_readonly("uv", [](const Mesh<T> &x) { return points_to_numpy(x.uv); }, "(3 m, 2) array: (u, v) of each corner on its face")
        .def_property_readonly("normals", [](const Mesh<T> &x) { return points_to_numpy(x.normals); }, "(3 m, 3) array: unit normal at each corner, oriented")
        .def_readonly("faces", &Mesh<T>::faces)
        .def_readonly("incomplete", &Mesh<T>::incomplete, "faces whose refinement stopped at max_points")
        .def_readonly("folded", &Mesh<T>::folded, "faces with a triangle against the normal of the surface at a corner (a defect); not those incomplete")
        .def_property_readonly(
            "edges",
            [](const Mesh<T> &x) {
                py::list out;
                for (const auto &e : x.edges)
                    out.append(py::make_tuple(e.edge, py::array_t<std::uint32_t>(py::ssize_t(e.vertices.size()), e.vertices.data())));
                return out;
            },
            "the non-degenerate edges of the shape, to draw them: (EdgeId, array of the vertices along it, from its v1 to its v2)")
        .def("area", [](const Mesh<T> &x) { return area(x); })
        .def("volume", [](const Mesh<T> &x) { return volume(x); }, "volume enclosed by a closed mesh (positive facing outwards)")
        .def("open_sides", [](const Mesh<T> &x) { return open_sides(x); }, "sides not used once in each sense (0: closed)");
    m.def("mesh", [](const M &mdl, const ShapeId &s, const MeshOptions<T> &o) { return value_or_raise(mesh(mdl, s, o)); }, py::arg("model"),
          py::arg("shape"), py::arg_v("options", MeshOptions<T>{}, "MeshOptions()"),
          "Triangle mesh of a shape, watertight along the shared edges, within the deflection and the angle of the options", release_gil());
    m.def(
        "write_stl",
        [](const Mesh<T> &x, const std::filesystem::path &path, bool binary) {
            if (auto r = io::write_stl_file(x, path, binary); !r)
                throw BuildFailure(std::move(r.error()));
        },
        py::arg("mesh"), py::arg("path"), py::arg("binary") = true, "Writes the STL file of a mesh (binary or ASCII)", release_gil());
    m.def(
        "write_obj",
        [](const Mesh<T> &x, const std::filesystem::path &path) {
            if (auto r = io::write_obj_file(x, path); !r)
                throw BuildFailure(std::move(r.error()));
        },
        py::arg("mesh"), py::arg("path"), "Writes the Wavefront OBJ file of a mesh (a normal per corner, a group per face)", release_gil());

    // polygon meshes in OBJ (gbs-io/obj.h): faces of any number of sides, lines and points
    using PM = io::PolygonMesh<T>;
    using Lists = std::vector<std::vector<std::uint32_t>>;
    const auto lists = [](const io::Polygons &p) {
        Lists r;
        r.reserve(p.size());
        for (const auto f : p)
            r.emplace_back(f.begin(), f.end());
        return r;
    };
    py::class_<PM>(m, "PolygonMesh",
                   "A polygon mesh: vertices, faces of 3 sides or more, polylines (OBJ l) and points (OBJ p), by the index of their vertices "
                   "from 0. Its attributes are copies: assign a whole new value to change one (mesh.faces = faces), not an element of it "
                   "(mesh.faces.append(f) changes a copy only)")
        .def(py::init([](const points_vector<T, 3> &vertices, const Lists &faces, const Lists &lines, std::vector<std::uint32_t> points) {
                 return PM{vertices, io::Polygons{faces}, io::Polygons{lines}, std::move(points)};
             }),
             py::arg("vertices") = points_vector<T, 3>{}, py::arg("faces") = Lists{}, py::arg("lines") = Lists{},
             py::arg("points") = std::vector<std::uint32_t>{})
        .def_property(
            "vertices", [](const PM &x) { return points_to_numpy(x.vertices); }, [](PM &x, const points_vector<T, 3> &v) { x.vertices = v; },
            "the vertices, an (n, 3) array (a copy)")
        .def_property(
            "faces", [lists](const PM &x) { return lists(x.faces); }, [](PM &x, const Lists &f) { x.faces = io::Polygons{f}; },
            "the faces, lists of vertex indices (a copy)")
        .def_property(
            "lines", [lists](const PM &x) { return lists(x.lines); }, [](PM &x, const Lists &l) { x.lines = io::Polygons{l}; },
            "the polylines, lists of vertex indices (a copy)")
        .def_property(
            "points", [](const PM &x) { return x.points; }, [](PM &x, std::vector<std::uint32_t> p) { x.points = std::move(p); },
            "the points, vertex indices (a copy)")
        .def("__eq__", [](const PM &a, const PM &b) { return a == b; }, py::is_operator());
    // ObjError: a ValueError that carries the line (from 1) and the reason
    const py::object obj_error = py::exception<io::ObjError>(m, "ObjError", PyExc_ValueError);
    const auto raise = [obj_error](const io::ObjError &e) {
        py::object x = obj_error(e.what());
        x.attr("line") = e.line;
        x.attr("reason") = std::string{e.reason()};
        PyErr_SetObject(obj_error.ptr(), x.ptr());
        throw py::error_already_set();
    };
    const auto os_error = [](const std::filesystem::filesystem_error &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        throw py::error_already_set();
    };
    // read without the GIL, the mesh moved into Python (not copied)
    const auto read = [raise, os_error](auto &&f) {
        std::optional<io::ObjRead<T>> r;
        try
        {
            py::gil_scoped_release nogil;
            r = f();
        }
        catch (const io::ObjError &e)
        {
            raise(e);
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            os_error(e);
        }
        return py::make_tuple(std::move(r->mesh), std::move(r->skipped));
    };
    m.def(
        "read_obj", [read](const std::filesystem::path &path) { return read([&] { return io::read_obj_file<T>(path); }); }, py::arg("path"),
        "The polygon mesh of a Wavefront OBJ file and the statements skipped, counted by keyword: (PolygonMesh, [(keyword, count)]); "
        "ObjError (a ValueError, with line and reason) for a fault, OSError if the file cannot be read");
    m.def(
        "read_obj_text", [read](const std::string &text) { return read([&] { return io::read_obj<T>(text); }); }, py::arg("text"),
        "Same as read_obj, from the text of an OBJ file");
    m.def(
        "write_obj",
        [os_error](const PM &x, const std::filesystem::path &path, const std::string &comment) {
            try
            {
                py::gil_scoped_release nogil;
                io::write_obj_file(x, path, comment);
            }
            catch (const std::filesystem::filesystem_error &e)
            {
                os_error(e);
            }
        },
        py::arg("mesh"), py::arg("path"), py::arg("comment") = "",
        "Writes the OBJ file of a polygon mesh, coordinates exact, whole or not at all; ValueError for what read_obj would refuse, "
        "OSError if the file cannot be written");
    m.def(
        "write_obj_text", [](const PM &x, const std::string &comment) { return io::write_obj(x, comment); }, py::arg("mesh"),
        py::arg("comment") = "", "The OBJ text of a polygon mesh", release_gil());
}
