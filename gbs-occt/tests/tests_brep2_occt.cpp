// Stage 2 of the native BREP kernel checked by OpenCascade: booleans, sections and sweeps of
// gbs on shapes built by OCCT (passed through STEP) are compared with BRepAlgoAPI and
// BRepPrimAPI; the results of gbs are written to STEP by gbs and read back by OCCT, which must
// find them valid, with the same volumes, and as many faces as the results of OCCT with their
// faces of one domain merged (ShapeUpgrade_UnifySameDomain).
#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include <gbs-occt/curvesbuild.h>

#include <blade_fixture.h>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
#include <NCollection_List.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <Geom_BSplineCurve.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Pnt.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <gp_Pln.hxx>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>

using namespace gbs;
using namespace gbs::brep;

namespace
{
    constexpr double pi = std::numbers::pi;

    std::string tmp_file(const std::string &name) { return (std::filesystem::temp_directory_path() / name).string(); }

    /// An OCCT shape read by gbs, through a STEP file written by OCCT.
    ShapeId to_gbs(Model<double> &m, const TopoDS_Shape &s, const std::string &name)
    {
        const auto file = tmp_file("gbs_brep2_" + name + ".stp");
        STEPControl_Writer w;
        Interface_Static::SetCVal("write.step.schema", "AP214IS");
        REQUIRE(w.Transfer(s, STEPControl_AsIs) == IFSelect_RetDone);
        REQUIRE(w.Write(file.c_str()) == IFSelect_RetDone);
        auto r = gbs::step::read_step_file(m, file);
        REQUIRE(r.has_value());
        std::string issues;
        for (const auto &e : r->report.check.entries)
            issues += std::string(to_string(e.issue)) + " " + e.detail + "; ";
        INFO(issues);
        REQUIRE(r->report.check.ok());
        return r->root;
    }

    /// A gbs shape read by OCCT, through a STEP text written by gbs.
    TopoDS_Shape to_occt(const Model<double> &m, const ShapeId &s, const std::string &name)
    {
        auto w = gbs::step::write_step(m, s, {.flatten = true});
        REQUIRE(w.has_value());
        const auto file = tmp_file("gbs_brep2_written_" + name + ".stp");
        {
            std::ofstream f(file, std::ios::binary);
            f << w->first;
        }
        STEPControl_Reader rd;
        REQUIRE(rd.ReadFile(file.c_str()) == IFSelect_RetDone);
        rd.TransferRoots();
        return rd.OneShape();
    }

    double volume(const TopoDS_Shape &s)
    {
        GProp_GProps p;
        BRepGProp::VolumeProperties(s, p, 1e-9, true);
        return p.Mass();
    }

    double length(const TopoDS_Shape &s)
    {
        GProp_GProps p;
        BRepGProp::LinearProperties(s, p);
        return p.Mass();
    }

    int count(const TopoDS_Shape &s, TopAbs_ShapeEnum t)
    {
        NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> mp;
        TopExp::MapShapes(s, t, mp);
        return mp.Extent();
    }

    /// The number of faces of a gbs result against the OCCT reference with its faces of one domain merged.
    enum class Faces
    {
        same,  ///< as many
        fewer, ///< as many or fewer: gbs merges pieces of a face that OCCT keeps apart
    };

    /// The faces of `s` on one domain merged by OCCT: the reference of the unified booleans of gbs.
    TopoDS_Shape unified(const TopoDS_Shape &s)
    {
        ShapeUpgrade_UnifySameDomain u(s, true, true, false);
        u.Build();
        return u.Shape();
    }

    /// The gbs result, read by OCCT, against the OCCT result: valid, same volume, and as many
    /// faces as the OCCT result with its faces of one domain merged.
    void expect_same(const Model<double> &m, const ShapeId &g, const TopoDS_Shape &ref, const std::string &name, double rel = 2e-3, Faces faces = Faces::same,
                     bool same_solids = true)
    {
        INFO(name);
        const auto back = to_occt(m, g, name);
        BRepCheck_Analyzer ana(back);
        CHECK(ana.IsValid());
        const double vr = volume(ref), vb = volume(back);
        INFO("volume occt " << vr << ", gbs read by occt " << vb);
        CHECK(std::abs(vb - vr) <= rel * std::max(1., std::abs(vr)));
        if (same_solids)
            CHECK(count(back, TopAbs_SOLID) == count(ref, TopAbs_SOLID));
        else if (count(back, TopAbs_SOLID) != count(ref, TopAbs_SOLID))
            return; // (other solids: their faces are not comparable)
        const int nb = count(back, TopAbs_FACE), nr = count(unified(ref), TopAbs_FACE);
        INFO("faces occt " << count(ref, TopAbs_FACE) << ", unified " << nr << ", gbs " << nb);
        switch (faces)
        {
        case Faces::same: CHECK(nb == nr); break;
        case Faces::fewer: CHECK(nb <= nr); break;
        }
    }

    TopoDS_Shape nurbs(const TopoDS_Shape &s) { return BRepBuilderAPI_NurbsConvert(s, true).Shape(); }

    /// The block -1..1 x -1..1 under z = -bend max(x, 0)^3 (a cubic B-spline with a knot at 0),
    /// down to z = -2: its profile in y = -1 extruded along y, converted to NURBS.
    TopoDS_Shape bspline_block(double bend)
    {
        NCollection_Array1<gp_Pnt> poles(1, 5);
        const double gx[]{-1., -2. / 3., 0., 2. / 3., 1.}, zx[]{0., 0., 0., 0., -bend};
        for (int i = 0; i < 5; ++i)
            poles.SetValue(i + 1, gp_Pnt(gx[i], -1., zx[i]));
        NCollection_Array1<double> knots(1, 3);
        NCollection_Array1<int> mults(1, 3);
        knots.SetValue(1, -1.), knots.SetValue(2, 0.), knots.SetValue(3, 1.);
        mults.SetValue(1, 4), mults.SetValue(2, 1), mults.SetValue(3, 4);
        Handle(Geom_BSplineCurve) top = new Geom_BSplineCurve(poles, knots, mults, 3);
        BRepBuilderAPI_MakeWire w;
        w.Add(BRepBuilderAPI_MakeEdge(top).Edge());
        w.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1., -1., -bend), gp_Pnt(1., -1., -2.)).Edge());
        w.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1., -1., -2.), gp_Pnt(-1., -1., -2.)).Edge());
        w.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(-1., -1., -2.), gp_Pnt(-1., -1., 0.)).Edge());
        const auto face = BRepBuilderAPI_MakeFace(w.Wire(), true).Face();
        return nurbs(BRepPrimAPI_MakePrism(face, gp_Vec(0., 2., 0.)).Shape());
    }

    TopoDS_Shape moved(const TopoDS_Shape &s, const gp_Vec &v)
    {
        gp_Trsf t;
        t.SetTranslation(v);
        return BRepBuilderAPI_Transform(s, t, true).Shape();
    }
}

TEST(tests_brep2_occt, booleans)
{
    const TopoDS_Shape box = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
    struct Pair
    {
        std::string name;
        TopoDS_Shape a, b;
        std::array<Faces, 3> faces{}; ///< number of faces against OCCT, its faces of one domain merged: fuse, cut, common
        TopoDS_Shape ra{}, rb{}; ///< operands of the OCCT reference, when not a and b
        bool same_solids{true};  ///< same number of solids as OCCT
    };
    const std::vector<Pair> pairs{
        {"boxes", box, moved(box, gp_Vec(1., 1., 1.))},
        {"box_cylinder", box, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 1.5, -1.), gp_Dir(0., 0., 1.)), 0.8, 6.).Shape()},
        {"box_sphere", box, BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(2., 3., 4.), gp_Dir(1., 1., 1.)), 1.5).Shape()},
        {"crossed_cylinders", BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., -3.), gp_Dir(0., 0., 1.)), 1., 6.).Shape(),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-3., 0., 0.), gp_Dir(1., 0., 0.), gp_Dir(0., 0., 1.)), 0.6, 6.).Shape()},
        {"coplanar_boxes", box, moved(box, gp_Vec(1., 1., 0.))},
        {"stacked_boxes", box, moved(box, gp_Vec(1., 0., 4.))},
        {"coaxial_cylinders", BRepPrimAPI_MakeCylinder(1., 2.).Shape(), moved(BRepPrimAPI_MakeCylinder(1., 2.).Shape(), gp_Vec(0., 0., 1.))},
        // contacts: at a point, along a line, and two slivers touching along a line
        {"sphere_touching", box, BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(1., 1.5, 5.5), gp_Dir(1., 1., 1.)), 1.5).Shape()},
        // (the face under a line of contact: one face in gbs, two in OCCT)
        {"cylinder_lying", box, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -1., 5.), gp_Dir(0., 1., 0.), gp_Dir(0., 0., 1.)), 1., 5.).Shape(), {Faces::fewer, Faces::same, Faces::same}},
        {"tangent_slivers", moved(box, gp_Vec(3., -1.5, 3.)),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(std::cos(1.), std::sin(1.), 0.)), 5., 10.).Shape()},
        // edges lying on a face of the other operand
        {"box_on_edge", box, BRepPrimAPI_MakeBox(gp_Ax2(gp_Pnt(1., 3.5, 4.), gp_Dir(0., -1., 0.), gp_Dir(1., 0., 1.)), 2., 3., 4.).Shape(), {Faces::fewer, Faces::same, Faces::same}},
        {"boxes_on_edge", box, moved(box, gp_Vec(2., 3., 1.))},
        {"seam_on_face", box, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -1., 5.), gp_Dir(0., 1., 0.), gp_Dir(0., 0., -1.)), 1., 5.).Shape(),
         {Faces::fewer, Faces::same, Faces::same}},
        // intersection curves crossing at branch points: equal crossed cylinders, Viviani
        {"steinmetz", BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., -3.), gp_Dir(0., 0., 1.)), 1., 6.).Shape(),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-3., 0., 0.), gp_Dir(1., 0., 0.), gp_Dir(0., 0., 1.)), 1., 6.).Shape()},
        // (the axis of the sphere along y: its poles off the curve, which goes through (0, 0, +-2);
        // the side of the cylinder inside the sphere pinched at the double point of the curve: one
        // face in gbs, two in OCCT)
        {"viviani", BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.), gp_Dir(0., 0., 1.)), 2.).Shape(),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 0., -3.), gp_Dir(0., 0., 1.), gp_Dir(0., 1., 0.)), 1., 6.).Shape(), {Faces::same, Faces::fewer, Faces::same}},
        // curves through poles and apexes: Viviani with the axis of the sphere along the cylinder
        // (both poles on the curve), a cone cut through its apex
        {"viviani_poles", BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(0., 1., 0.)), 2.).Shape(),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 0., -3.), gp_Dir(0., 0., 1.), gp_Dir(0., 1., 0.)), 1., 6.).Shape()},
        {"half_cylinder", BRepPrimAPI_MakeCylinder(1., 2.).Shape(), BRepPrimAPI_MakeBox(gp_Pnt(-2., -2., -1.), 2., 4., 4.).Shape()},
        {"cone_apex", BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(std::cos(1.), std::sin(1.), 0.)), 2., 0., 3.).Shape(),
         BRepPrimAPI_MakeBox(gp_Pnt(0., -3., -1.), 3., 6., 5.).Shape()},
        // coincident faces of B-spline surfaces (the operands converted to NURBS by OCCT), and of
        // a primitive and a B-spline
        {"coplanar_nurbs", nurbs(box), nurbs(moved(box, gp_Vec(1., 1., 0.)))},
        {"stacked_nurbs", nurbs(box), nurbs(moved(box, gp_Vec(1., 0., 4.)))},
        {"plane_and_nurbs", box, nurbs(moved(box, gp_Vec(1., 1., 0.)))},
        // tangencies along a line between the B-spline faces of a box (converted to NURBS by
        // OCCT) and a cylinder, traced by marching (OCCT's NURBS cylinder is read with pcurves
        // 5e-6 off its edges, so the cylinder is kept)
        {"cylinder_lying_nurbs", nurbs(box), BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -1., 5.), gp_Dir(0., 1., 0.), gp_Dir(0., 0., 1.)), 1., 5.).Shape(),
         {Faces::fewer, Faces::same, Faces::same}},
        {"tangent_slivers_nurbs", nurbs(moved(box, gp_Vec(3., -1.5, 3.))),
         BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(std::cos(1.), std::sin(1.), 0.)), 5., 10.).Shape(), {},
         moved(box, gp_Vec(3., -1.5, 3.)), BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(std::cos(1.), std::sin(1.), 0.)), 5., 10.).Shape()},
        // a contact at a point between B-spline surfaces (the sphere, or both operands, converted
        // to NURBS by OCCT), certified by the second fundamental forms
        {"sphere_touching_nurbs", nurbs(box), nurbs(BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(1., 1.5, 5.5), gp_Dir(1., 1., 1.)), 1.5).Shape())},
        {"plane_touching_nurbs_sphere", box, nurbs(BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(1., 1.5, 5.5), gp_Dir(1., 1., 1.)), 1.5).Shape())},
        {"seam_on_face_nurbs", nurbs(box), BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., -1., 5.), gp_Dir(0., 1., 0.), gp_Dir(0., 0., -1.)), 1., 5.).Shape(),
         {Faces::fewer, Faces::same, Faces::same}},
        // faces coincident on a part of their overlap: the bottom of a box on the B-spline top of a
        // block over its patch x <= 0, parting from it beyond, or rising into the box
        // (OCCT does not see the contact of the parting case: its fusion is the two operands as
        // they are, two solids, where gbs merges them along the common region; the planar part of
        // the B-spline top merges with the coplanar faces next to it)
        {"partial_coincidence", bspline_block(1.), BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., 0.), 1.5, 2., 1.).Shape(), {}, {}, {}, false},
        {"partial_coincidence_rising", bspline_block(-1.), BRepPrimAPI_MakeBox(gp_Pnt(-1., -1., 0.), 1.5, 2., 1.).Shape()},
        {"cone_apex_seam", BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.), gp_Dir(0., 1., 0.)), 2., 0., 3.).Shape(),
         BRepPrimAPI_MakeBox(gp_Pnt(0., -3., -1.), 3., 6., 5.).Shape()},
    };
    for (const auto &[name, a, b, faces, ra_, rb_, same_solids] : pairs)
    {
        Model<double> m;
        const auto ga = to_gbs(m, a, name + "_a"), gb = to_gbs(m, b, name + "_b");
        const auto &ra = ra_.IsNull() ? a : ra_, &rb = rb_.IsNull() ? b : rb_;
        struct Op
        {
            const char *name;
            BooleanOp op;
            TopoDS_Shape ref;
        };
        for (const auto &[oname, op, ref] : {Op{"fuse", BooleanOp::Fuse, BRepAlgoAPI_Fuse(ra, rb).Shape()},
                                             Op{"cut", BooleanOp::Cut, BRepAlgoAPI_Cut(ra, rb).Shape()},
                                             Op{"common", BooleanOp::Common, BRepAlgoAPI_Common(ra, rb).Shape()}})
        {
            const auto label = name + "_" + oname;
            INFO(label);
            const auto t0 = std::chrono::steady_clock::now();
            auto r = boolean(m, op, ga, gb);
            const auto t1 = std::chrono::steady_clock::now();
            if (!r)
                FAIL(r.error().message);
            MESSAGE(label << ": " << std::chrono::duration<double>(t1 - t0).count() << " s, " << r->solids.size() << " solids");
            REQUIRE(check(m, r->shape).ok());
            if (r->solids.empty())
            {
                CHECK(volume(ref) < 1e-9);
                continue;
            }
            expect_same(m, r->shape, ref, label, 2e-3, faces[op == BooleanOp::Fuse ? 0 : op == BooleanOp::Cut ? 1 : 2], same_solids);
        }
    }
}

TEST(tests_brep2_occt, sections)
{
    const TopoDS_Shape box = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
    const TopoDS_Shape cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 2., 5.).Shape();
    const TopoDS_Shape sph = BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0.2, 0.3, 1.)), 2.).Shape();
    const gp_Pln oblique(gp_Pnt(1., 1.5, 2.), gp_Dir(0.3, 0.2, 1.));
    struct Case
    {
        std::string name;
        TopoDS_Shape s;
        gp_Pln pl;
    };
    for (const auto &[name, s, pl] : {Case{"box", box, oblique}, Case{"cylinder", cyl, gp_Pln(gp_Pnt(0., 0., 2.5), gp_Dir(0., 0.4, 1.))},
                                      Case{"sphere", sph, gp_Pln(gp_Pnt(0., 0., 0.5), gp_Dir(0.1, 0., 1.))}})
    {
        INFO(name);
        const double ref = length(BRepAlgoAPI_Section(s, pl).Shape());
        Model<double> m;
        const auto g = to_gbs(m, s, "section_" + name);
        const auto &o = pl.Location();
        const auto &n = pl.Axis().Direction();
        const auto &x = pl.XAxis().Direction();
        auto plane = std::make_shared<Plane<double>>(ax2<double, 3>{{{o.X(), o.Y(), o.Z()}, {n.X(), n.Y(), n.Z()}, {x.X(), x.Y(), x.Z()}}},
                                                     std::array<double, 4>{-50., 50., -50., 50.});
        auto r = section(m, g, plane);
        REQUIRE(r.has_value());
        double l{};
        for (const auto &se : r->edges)
        {
            const auto &e = m.edge(se.edge);
            constexpr int k = 2000;
            for (int i = 0; i < k; ++i)
                l += norm(e.curve->value(e.u1 + (e.u2 - e.u1) * (i + 1.) / k) - e.curve->value(e.u1 + (e.u2 - e.u1) * double(i) / k));
        }
        INFO("length occt " << ref << " gbs " << l);
        CHECK(std::abs(l - ref) <= 1e-4 * ref);
        for (auto w : r->wires)
            CHECK(m.wire(w).closed);
    }
}

TEST(tests_brep2_occt, sweeps)
{
    // a plate 4 x 2 with a hole of radius 0.5, built by OCCT
    BRepBuilderAPI_MakePolygon poly(gp_Pnt(0., 0., 0.), gp_Pnt(4., 0., 0.), gp_Pnt(4., 2., 0.), gp_Pnt(0., 2., 0.), true);
    const auto hole = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(1., 1., 0.), gp_Dir(0., 0., -1.)), 0.5))).Wire();
    BRepBuilderAPI_MakeFace mf(poly.Wire());
    mf.Add(hole);
    const TopoDS_Face plate = mf.Face();
    Model<double> m;
    const auto gp = to_gbs(m, plate, "plate");
    const auto f = explore<FaceId>(m, gp).front();

    // extrusion
    const auto prism = unwrap(make_prism(m, f, point<double, 3>{0.5, 0.2, 3.}));
    expect_same(m, prism.shape, BRepPrimAPI_MakePrism(plate, gp_Vec(0.5, 0.2, 3.)).Shape(), "prism");
    // revolution around an axis of the plane of the plate, away from it: full and partial
    const gp_Ax1 axis(gp_Pnt(0., -1., 0.), gp_Dir(1., 0., 0.));
    const ax1<double, 3> gaxis{{{0., -1., 0.}, {1., 0., 0.}}};
    const auto full = unwrap(make_revolution(m, f, gaxis));
    expect_same(m, full.shape, BRepPrimAPI_MakeRevol(plate, axis).Shape(), "revolution");
    const auto part = unwrap(make_revolution(m, f, gaxis, 2.));
    expect_same(m, part.shape, BRepPrimAPI_MakeRevol(plate, axis, 2.).Shape(), "partial_revolution");
}

namespace
{
    /// The OCCT pipe of a profile wire along a spine, corrected Frenet (OCCT's frame closest to a
    /// rotation-minimizing one), closed into a solid.
    TopoDS_Shape occt_pipe(const TopoDS_Wire &spine, const TopoDS_Wire &profile)
    {
        BRepOffsetAPI_MakePipeShell p(spine);
        p.SetMode(false);
        p.Add(profile);
        p.Build();
        REQUIRE(p.IsDone());
        REQUIRE(p.MakeSolid());
        return p.Shape();
    }

    /// The disk of radius r at c normal to n, in gbs.
    FaceId gbs_disk(Model<double> &m, double r, const point<double, 3> &c, const point<double, 3> &n, const point<double, 3> &x)
    {
        const ax2<double, 3> ax{{c, n, x}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(r, r, ax)));
        const auto w = unwrap(make_wire(m, std::vector<EdgeId>{e}));
        return unwrap(make_face(m, std::make_shared<Plane<double>>(ax, std::array<double, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), w));
    }
}

// make_pipe against BRepOffsetAPI_MakePipeShell: the volume does not depend on the frame for a
// profile normal to the path whose centroid is on it (Pappus), nor the number of faces.
TEST(tests_brep2_occt, pipes)
{
    // a disk r 0.5 along an elbow: segment, quarter of circle R 2, segment
    {
        Model<double> m;
        const auto d = gbs_disk(m, 0.5, {0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
        const auto arc = std::make_shared<EllipseCurve<double>>(2., 2., ax2<double, 3>{{{2., 0., 3.}, {0., 1., 0.}, {-1., 0., 0.}}}, 0., pi / 2);
        const auto end = arc->value(pi / 2);
        const auto e1 = unwrap(make_edge(m, point<double, 3>{0., 0., 0.}, point<double, 3>{0., 0., 3.}));
        const auto e2 = unwrap(make_edge(m, arc));
        const auto e3 = unwrap(make_edge(m, end, point<double, 3>{6., 0., 5.}));
        const auto g = unwrap(make_pipe(m, d, unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3}))));

        BRepBuilderAPI_MakeWire spine;
        spine.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(0., 0., 0.), gp_Pnt(0., 0., 3.)).Edge());
        spine.Add(BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(gp_Pnt(0., 0., 3.), gp_Vec(0., 0., 1.), gp_Pnt(end[0], end[1], end[2])).Value()).Edge());
        spine.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(end[0], end[1], end[2]), gp_Pnt(6., 0., 5.)).Edge());
        const auto profile = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 0.5))).Wire();
        expect_same(m, g.shape, occt_pipe(spine.Wire(), profile), "pipe_elbow");
    }
    // a square 0.4 x 0.4 along a B-spline helix (radius 2, pitch 1.5, two turns)
    {
        Model<double> m;
        std::vector<double> u(97);
        std::vector<constrType<double, 3, 1>> q(97);
        for (std::size_t i = 0; i <= 96; ++i)
        {
            u[i] = double(i) / 96.;
            const double t = 4 * pi * u[i];
            q[i] = {point<double, 3>{2 * std::cos(t), 2 * std::sin(t), 1.5 * t / (2 * pi)}};
        }
        const auto h = std::make_shared<BSCurve<double, 3>>(interpolate<double, 3>(q, u, std::size_t{3}));
        const auto t0 = h->value(0., 1) / norm(h->value(0., 1));
        const point<double, 3> c{2., 0., 0.};
        auto x = point<double, 3>{1., 0., 0.} - t0[0] * t0; // normal to the path (the spline's end tangent is a little off the helix's)
        x = x / norm(x);
        const auto y = cross(t0, x);
        std::array<point<double, 3>, 4> sq{c - 0.2 * x - 0.2 * y, c + 0.2 * x - 0.2 * y, c + 0.2 * x + 0.2 * y, c - 0.2 * x + 0.2 * y};
        std::vector<EdgeId> es;
        BRepBuilderAPI_MakePolygon poly;
        for (int i = 0; i < 4; ++i)
        {
            es.push_back(unwrap(make_edge(m, sq[i], sq[(i + 1) % 4])));
            poly.Add(gp_Pnt(sq[i][0], sq[i][1], sq[i][2]));
        }
        poly.Close();
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{c, t0, x}}, std::array<double, 4>{-1., 1., -1., 1.}),
                                        unwrap(make_wire(m, es))));
        const auto g = unwrap(make_pipe(m, f, unwrap(make_edge(m, h, 0., 1.))));
        const auto spine = BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(occt_utils::BSplineCurve(*h)).Edge()).Wire();
        expect_same(m, g.shape, occt_pipe(spine, poly.Wire()), "pipe_helix");
    }
}

namespace
{
    /// BRepOffsetAPI_MakeOffsetShape with sharp joints (GeomAbs_Intersection).
    TopoDS_Shape occt_offset(const TopoDS_Shape &s, double d)
    {
        BRepOffsetAPI_MakeOffsetShape mk;
        mk.PerformByJoin(s, d, 1e-7, BRepOffset_Skin, false, false, GeomAbs_Intersection);
        REQUIRE(mk.IsDone());
        return mk.Shape();
    }

    /// The face of `s` whose box has its centre at c.
    TopoDS_Face occt_face_at(const TopoDS_Shape &s, const point<double, 3> &c)
    {
        for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next())
        {
            Bnd_Box b;
            BRepBndLib::Add(ex.Current(), b);
            double x1, y1, z1, x2, y2, z2;
            b.Get(x1, y1, z1, x2, y2, z2);
            if (std::hypot((x1 + x2) / 2 - c[0], (y1 + y2) / 2 - c[1], (z1 + z2) / 2 - c[2]) < 1e-5)
                return TopoDS::Face(ex.Current());
        }
        FAIL("no face there");
        return {};
    }

    /// BRepOffsetAPI_MakeThickSolid with sharp joints, gbs's thickness (positive inward).
    TopoDS_Shape occt_thick(const TopoDS_Shape &s, const TopoDS_Face &removed, double t)
    {
        NCollection_List<TopoDS_Shape> faces;
        faces.Append(removed);
        BRepOffsetAPI_MakeThickSolid mk;
        mk.MakeThickSolidByJoin(s, faces, -t, 1e-7, BRepOffset_Skin, false, false, GeomAbs_Intersection);
        REQUIRE(mk.IsDone());
        return mk.Shape();
    }

    FaceId gbs_polygon(Model<double> &m, const std::vector<point<double, 3>> &p, const point<double, 3> &n, const point<double, 3> &x)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{p[0], n, x}}, std::array<double, 4>{-10., 10., -10., 10.}),
                                unwrap(make_wire(m, e))));
    }

    SolidId gbs_solid(Model<double> &m, const BuildResult<SweepResult<double>> &r) { return explore<SolidId>(m, unwrap(r).shape).front(); }
}

// offset_solid and make_thick against BRepOffsetAPI_MakeOffsetShape and BRepOffsetAPI_MakeThickSolid
// with sharp joints (GeomAbs_Intersection), on the same solids (gbs's, read by OCCT)
TEST(tests_brep2_occt, offsets)
{
    const point<double, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
    // a box
    {
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1.5, 0.}, {0., 1.5, 0.}}, Z, X), point<double, 3>{0., 0., 1.}));
        const auto ref = to_occt(m, s, "offset_box_source");
        for (double d : {0.25, -0.2})
            expect_same(m, unwrap(offset_solid(m, s, d)), occt_offset(ref, d), "offset_box " + std::to_string(d));
        FaceId top;
        for (auto g : explore<FaceId>(m, ShapeId{s}))
            if (const auto *p = dynamic_cast<const Plane<double> *>(m.face(g).surface.get()); p && p->frame().o[2] > 0.5)
                top = g;
        expect_same(m, unwrap(make_thick(m, s, std::vector<FaceId>{top}, 0.1)),
                    occt_thick(ref, occt_face_at(ref, {1., 0.75, 1.}), 0.1), "thick_box");
    }
    // a cylinder, and a cup
    {
        Model<double> m;
        const ax2<double, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(1., 1., ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax, std::array<double, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        const auto s = gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., 2.}));
        const auto ref = to_occt(m, s, "offset_cylinder_source");
        for (double d : {0.3, -0.25})
            expect_same(m, unwrap(offset_solid(m, s, d)), occt_offset(ref, d), "offset_cylinder " + std::to_string(d));
        FaceId top;
        for (auto g : explore<FaceId>(m, ShapeId{s}))
            if (const auto *p = dynamic_cast<const Plane<double> *>(m.face(g).surface.get()); p && p->frame().o[2] > 1.)
                top = g;
        expect_same(m, unwrap(make_thick(m, s, std::vector<FaceId>{top}, 0.1)), occt_thick(ref, occt_face_at(ref, {0., 0., 2.}), 0.1), "thick_cup");
    }
    // a slab on a stadium: planes and half cylinders tangent along the vertical edges
    {
        Model<double> m;
        const double L = 3., R = 1.;
        const auto v = std::array{unwrap(make_vertex(m, point<double, 3>{0., -R, 0.})), unwrap(make_vertex(m, point<double, 3>{L, -R, 0.})),
                                  unwrap(make_vertex(m, point<double, 3>{L, R, 0.})), unwrap(make_vertex(m, point<double, 3>{0., R, 0.}))};
        const auto right = std::make_shared<EllipseCurve<double>>(R, R, ax2<double, 3>{{{L, 0., 0.}, Z, X}}, -pi / 2, pi / 2);
        const auto left = std::make_shared<EllipseCurve<double>>(R, R, ax2<double, 3>{{{0., 0., 0.}, Z, X}}, pi / 2, 3 * pi / 2);
        const std::vector<EdgeId> e{unwrap(make_edge(m, v[0], v[1])), unwrap(make_edge(m, right, -pi / 2, pi / 2, v[1], v[2])),
                                    unwrap(make_edge(m, v[2], v[3])), unwrap(make_edge(m, left, pi / 2, 3 * pi / 2, v[3], v[0]))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{{0., 0., 0.}, Z, X}}, std::array<double, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, e))));
        const auto s = gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., 0.8}));
        const auto ref = to_occt(m, s, "offset_slab_source");
        for (double d : {0.2, -0.3})
            expect_same(m, unwrap(offset_solid(m, s, d)), occt_offset(ref, d), "offset_slab " + std::to_string(d));
    }
    // a square along a helix: B-spline faces, sharp edges between them
    {
        Model<double> m;
        std::vector<double> u(97);
        std::vector<constrType<double, 3, 1>> q(97);
        for (std::size_t i = 0; i <= 96; ++i)
        {
            u[i] = double(i) / 96.;
            const double t = 2 * pi * u[i];
            q[i] = {point<double, 3>{2 * std::cos(t), 2 * std::sin(t), 1.5 * t / (2 * pi)}};
        }
        const auto h = std::make_shared<BSCurve<double, 3>>(interpolate<double, 3>(q, u, std::size_t{3}));
        const auto t0 = h->value(0., 1) / norm(h->value(0., 1));
        const point<double, 3> c{2., 0., 0.};
        auto x = X - t0[0] * t0;
        x = x / norm(x);
        const auto y = cross(t0, x);
        const auto sq = gbs_polygon(m, {c - 0.2 * x - 0.2 * y, c + 0.2 * x - 0.2 * y, c + 0.2 * x + 0.2 * y, c - 0.2 * x + 0.2 * y}, t0, x);
        const auto s = gbs_solid(m, make_pipe(m, sq, unwrap(make_edge(m, h, 0., 1.))));
        const auto ref = to_occt(m, s, "offset_helix_source");
        for (double d : {0.05, -0.05})
            expect_same(m, unwrap(offset_solid(m, s, d)), occt_offset(ref, d), "offset_helix " + std::to_string(d));
    }
    // OCCT's shapes read by gbs, where OCCT's own offset fails (MakeOffsetShape gives no solid on the
    // wedges of a ball, is not done on the cylinder of NURBS): gbs's offset read back by OCCT, valid,
    // with the volume in closed form (a half ball and a cylinder) or integrated (the wedge: the ball of
    // radius 2.7 in either half-space shifted by 0.3, over z the area of the disk in them, 45.423926)
    auto expect_read_back = [](const Model<double> &m, const ShapeId &g, double v, const std::string &name) {
        INFO(name);
        const auto back = to_occt(m, g, name);
        CHECK(BRepCheck_Analyzer(back).IsValid());
        CHECK(count(back, TopAbs_SOLID) == 1);
        INFO("volume of gbs read by occt " << volume(back) << ", expected " << v);
        CHECK(std::abs(volume(back) - v) <= 1e-5 * v);
    };
    {
        Model<double> m;
        const auto s = explore<SolidId>(m, to_gbs(m, BRepPrimAPI_MakeSphere(3., pi).Shape(), "offset_half_ball")).front();
        expect_read_back(m, unwrap(offset_solid(m, s, -0.3)), pi * 2.4 * 2.4 * (2 * 2.7 + 0.3) / 3., "offset_half_ball");
    }
    {
        Model<double> m;
        const auto s = explore<SolidId>(m, to_gbs(m, BRepPrimAPI_MakeSphere(3., 4.).Shape(), "offset_sphere_wedge")).front();
        expect_read_back(m, unwrap(offset_solid(m, s, -0.3)), 45.423926, "offset_sphere_wedge");
    }
    {
        Model<double> m;
        const auto s = explore<SolidId>(m, to_gbs(m, nurbs(BRepPrimAPI_MakeCylinder(1., 2.).Shape()), "offset_nurbs_cylinder_source")).front();
        for (double d : {0.1, -0.1})
            expect_read_back(m, unwrap(offset_solid(m, s, d)), pi * (1. + d) * (1. + d) * (2. + 2 * d), "offset_nurbs_cylinder " + std::to_string(d));
    }
}

namespace
{
    /// The edge of `s` nearest the point c.
    TopoDS_Edge occt_edge_near(const TopoDS_Shape &s, const point<double, 3> &c)
    {
        const auto v = BRepBuilderAPI_MakeVertex(gp_Pnt(c[0], c[1], c[2])).Vertex();
        TopoDS_Edge best;
        double gap = 1e300;
        for (TopExp_Explorer ex(s, TopAbs_EDGE); ex.More(); ex.Next())
        {
            BRepExtrema_DistShapeShape d(v, ex.Current());
            if (d.IsDone() && d.Value() < gap)
                gap = d.Value(), best = TopoDS::Edge(ex.Current());
        }
        REQUIRE(gap < 1e-6);
        return best;
    }

    /// The face of `s` nearest the point c (on it).
    TopoDS_Face occt_face_near(const TopoDS_Shape &s, const point<double, 3> &c)
    {
        const auto v = BRepBuilderAPI_MakeVertex(gp_Pnt(c[0], c[1], c[2])).Vertex();
        TopoDS_Face best;
        double gap = 1e300;
        for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next())
        {
            BRepExtrema_DistShapeShape d(v, ex.Current());
            if (d.IsDone() && d.Value() < gap)
                gap = d.Value(), best = TopoDS::Face(ex.Current());
        }
        REQUIRE(gap < 1e-6);
        return best;
    }

    /// The edge of the gbs solid whose middle is nearest p.
    EdgeId gbs_edge_near(const Model<double> &m, SolidId s, const point<double, 3> &p)
    {
        EdgeId best{};
        double gap = 1e300;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            for (int i = 0; i <= 16; ++i)
            {
                const double d = norm(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * i / 16.) - p);
                if (d < gap)
                    gap = d, best = e;
            }
        }
        return best;
    }

    /// The face of the edge e of the gbs solid whose outward normal at the middle of the edge is nearest the direction `dir`.
    FaceId gbs_face_toward(const Model<double> &m, SolidId s, EdgeId e, const point<double, 3> &dir)
    {
        const auto &ed = m.edge(e);
        FaceId best{};
        double dot = -2.;
        for (auto sh : explore<ShellId>(m, ShapeId{s}))
        {
            std::map<std::uint32_t, Orientation> orient;
            for (const auto &fu : m.shell(sh).faces)
                orient[fu.face.index] = fu.orient;
            const auto uses = shell_edge_uses(m, sh);
            if (const auto it = uses.find(e); it != uses.end())
                for (const auto &u : it->second)
                {
                    const auto &ce = m.wire(u.coedge.wire).coedges[u.coedge.index];
                    const auto uv = ce.pcurve->value((ed.u1 + ed.u2) / 2.);
                    auto n = normal_direction(*m.face(u.face).surface, uv[0], uv[1]);
                    if (orient.at(u.face.index) == Orientation::Reversed)
                        n = -1. * n;
                    if (const double d = n * dir / norm(dir); d > dot)
                        dot = d, best = u.face;
                }
        }
        return best;
    }

    /**
     * The surface of the blend of gbs over the edge (5 % short of its ends) lies on the faces
     * that OCCT made from the edge: the largest distance of 17 x 9 points between the contacts.
     */
    double gap_to_occt(const EdgeBlend<double> &b, const NCollection_List<TopoDS_Shape> &generated)
    {
        REQUIRE(generated.Extent() > 0);
        double worst = 0.;
        const double t1 = std::min(b.t1, b.t2), t2 = std::max(b.t1, b.t2);
        for (int i = 0; i <= 16; ++i)
        {
            const double t = t1 + (t2 - t1) * (0.05 + 0.9 * i / 16.);
            const auto q1 = b.contacts[0].on_blend->value(t), q2 = b.contacts[1].on_blend->value(t);
            for (int j = 0; j <= 8; ++j)
            {
                const auto q = q1 + (j / 8.) * (q2 - q1);
                const auto p = b.surface->value(q[0], q[1]);
                const auto v = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])).Vertex();
                double gap = 1e300;
                for (const auto &f : generated)
                {
                    BRepExtrema_DistShapeShape d(v, f);
                    if (d.IsDone())
                        gap = std::min(gap, d.Value());
                }
                worst = std::max(worst, gap);
            }
        }
        return worst;
    }

    /// The fillet of gbs along the edge near p against BRepFilletAPI_MakeFillet, within `tol`.
    void expect_fillet(const Model<double> &m, SolidId s, const TopoDS_Shape &ref, const point<double, 3> &p, double r, double tol, const std::string &name)
    {
        INFO(name);
        const auto b = unwrap(fillet_geometry(m, s, gbs_edge_near(m, s, p), r));
        const auto edge = occt_edge_near(ref, p);
        BRepFilletAPI_MakeFillet mk(ref);
        mk.Add(r, edge);
        mk.Build();
        REQUIRE(mk.IsDone());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        const double gap = gap_to_occt(b, mk.Generated(edge));
        INFO("gap to OCCT " << gap);
        CHECK(gap <= tol);
    }

    /// The chamfer of gbs, d1 on the face of the edge facing `dir`, against BRepFilletAPI_MakeChamfer on the same face.
    void expect_chamfer(const Model<double> &m, SolidId s, const TopoDS_Shape &ref, const point<double, 3> &p, const point<double, 3> &dir, double d1,
                        double d2, double tol, const std::string &name)
    {
        INFO(name);
        const auto e = gbs_edge_near(m, s, p);
        const auto b = unwrap(chamfer_geometry(m, s, e, gbs_face_toward(m, s, e, dir), d1, d2));
        const auto edge = occt_edge_near(ref, p);
        // the face of d1 in OCCT: where the first contact lies, inside the edge
        const auto c = b.contacts[0].curve->value((b.t1 + b.t2) / 2.);
        BRepFilletAPI_MakeChamfer mk(ref);
        mk.Add(d1, d2, edge, occt_face_near(ref, c));
        mk.Build();
        REQUIRE(mk.IsDone());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        const double gap = gap_to_occt(b, mk.Generated(edge));
        INFO("gap to OCCT " << gap);
        CHECK(gap <= tol);
    }
}

// fillet_geometry and chamfer_geometry against BRepFilletAPI_MakeFillet and BRepFilletAPI_MakeChamfer:
// the surface of gbs lies on the faces OCCT made from the edge, on the same solids (gbs's, read by OCCT)
TEST(tests_brep2_occt, fillet_geometry)
{
    const point<double, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
    // the edges of a box, convex and concave (an L-shaped prism): cylinders, planes
    {
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, Z, X),
                                               point<double, 3>{0., 0., 1.}));
        const auto ref = to_occt(m, s, "fillet_l_source");
        expect_fillet(m, s, ref, {2., 0.5, 1.}, 0.3, 1e-7, "fillet_convex");
        expect_fillet(m, s, ref, {1., 1., 0.5}, 0.25, 1e-7, "fillet_concave");
        expect_chamfer(m, s, ref, {2., 0.5, 1.}, Z, 0.2, 0.3, 1e-7, "chamfer_convex");
        expect_chamfer(m, s, ref, {1., 1., 0.5}, Y, 0.3, 0.15, 1e-7, "chamfer_concave");
    }
    // the top edge of a cylinder: a torus, a cone
    {
        Model<double> m;
        const ax2<double, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(1., 1., ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax, std::array<double, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        const auto s = gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., 2.}));
        const auto ref = to_occt(m, s, "fillet_cylinder_source");
        const point<double, 3> p{std::cos(1.), std::sin(1.), 2.};
        expect_fillet(m, s, ref, p, 0.2, 1e-7, "fillet_cylinder");
        expect_chamfer(m, s, ref, p, Z, 0.15, 0.1, 1e-7, "chamfer_cylinder");
    }
    // two crossed cylinders: a saddle, concave, a marched centre line. OCCT's fillet is the
    // further from the rolling ball: 1.1e-4 from gbs's surface, its points 8e-5 off the distance r
    // to gbs's centre line (gbs's contacts are within 1e-6 of it, along the normals of the faces);
    // its approximation tolerance is 1e-4 by default
    {
        Model<double> m;
        const ax2<double, 3> ax{{{0., 0., 0.}, Z, X}};
        auto cyl = [&](double r, double h) {
            const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(r, r, ax)));
            const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax, std::array<double, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
            return gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., h}));
        };
        const auto lying = std::get<SolidId>(copy_shape(m, ShapeId{cyl(1., 4.)}, Placement<double>::frame({-2., 0., 0.}, Y, Z, X)));
        const auto fu = unwrap(fuse(m, ShapeId{lying}, ShapeId{cyl(0.4, 1.5)}));
        REQUIRE(fu.solids.size() == 1u);
        const auto s = fu.solids.front();
        const auto ref = to_occt(m, ShapeId{s}, "fillet_crossed_source");
        const point<double, 3> p{0.4 * std::cos(1.), 0.4 * std::sin(1.), std::sqrt(1. - 0.16 * std::sin(1.) * std::sin(1.))};
        expect_fillet(m, s, ref, p, 0.15, 2e-4, "fillet_crossed");
        expect_chamfer(m, s, ref, p, Z, 0.1, 0.1, 1e-5, "chamfer_crossed");
    }
    // a B-spline bump: plane and B-spline faces, marched centre line, approximated surfaces
    {
        Model<double> m;
        const points_vector<double, 3> q{{0., 0., 0.}, {0.5, -0.3, 0.}, {1., 0.2, 0.}, {1.5, -0.1, 0.}, {2., 0., 0.}};
        const auto bump = std::make_shared<BSCurve<double, 3>>(interpolate<double, 3>(q, std::vector<double>{0., 0.25, 0.5, 0.75, 1.}, std::size_t{3}));
        const auto v0 = unwrap(make_vertex(m, point<double, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<double, 3>{2., 0., 0.}));
        const auto v2 = unwrap(make_vertex(m, point<double, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<double, 3>{0., 1., 0.}));
        const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                    unwrap(make_edge(m, v3, v0))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{{0., 0., 0.}, Z, X}}, std::array<double, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, w))));
        const auto s = gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., 0.8}));
        const auto ref = to_occt(m, s, "fillet_bump_source");
        const auto p = bump->value(0.5) + point<double, 3>{0., 0., 0.8};
        expect_fillet(m, s, ref, p, 0.1, 1e-5, "fillet_bump");
        expect_chamfer(m, s, ref, p, Z, 0.08, 0.05, 1e-5, "chamfer_bump");
    }
}

namespace
{
    /// The fillet (r) or chamfer (d, r < 0) of gbs of the edges through the points `at` against
    /// BRepFilletAPI_MakeFillet / MakeChamfer on the same solid: valid once read by OCCT, the
    /// same volume within `rel`, as many faces.
    void expect_blended(Model<double> &m, SolidId s, const std::vector<point<double, 3>> &at, double r, double d, double rel, const std::string &name,
                        Faces faces = Faces::same)
    {
        INFO(name);
        const auto ref = to_occt(m, s, name + "_source");
        std::vector<EdgeId> edges;
        for (const auto &p : at)
            edges.push_back(gbs_edge_near(m, s, p));
        const auto g = r > 0. ? fillet(m, s, edges, r) : chamfer(m, s, edges, d);
        INFO((g ? std::string{} : g.error().message));
        REQUIRE(g.has_value());
        TopoDS_Shape occt;
        if (r > 0.)
        {
            BRepFilletAPI_MakeFillet mk(ref);
            for (const auto &p : at)
                mk.Add(r, occt_edge_near(ref, p));
            mk.Build();
            REQUIRE(mk.IsDone());
            occt = mk.Shape();
        }
        else
        {
            BRepFilletAPI_MakeChamfer mk(ref);
            for (const auto &p : at)
                mk.Add(d, occt_edge_near(ref, p));
            mk.Build();
            REQUIRE(mk.IsDone());
            occt = mk.Shape();
        }
        CHECK(BRepCheck_Analyzer(occt).IsValid());
        expect_same(m, ShapeId{*g}, occt, name, rel, faces);
    }
}

// fillet and chamfer against BRepFilletAPI_MakeFillet and BRepFilletAPI_MakeChamfer: the solids
TEST(tests_brep2_occt, fillet_solids)
{
    const point<double, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
    auto cylinder = [&](Model<double> &m, double r, double h) {
        const ax2<double, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax, std::array<double, 4>{-2. * r, 2. * r, -2. * r, 2. * r}),
                                        unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., h}));
    };
    // the edges of an L-shaped prism: two parallel convex ones, the concave one
    {
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, Z, X),
                                               point<double, 3>{0., 0., 1.}));
        expect_blended(m, s, {{2., 0.5, 1.}, {0., 1., 1.}}, 0.3, 0., 1e-6, "fillet_l_convex");
        expect_blended(m, s, {{1., 1., 0.5}}, 0.25, 0., 1e-6, "fillet_l_concave");
        expect_blended(m, s, {{2., 0.5, 1.}, {1., 1., 0.5}}, -1., 0.2, 1e-6, "chamfer_l");
    }
    // the rims of a cylinder: tori, cones
    {
        Model<double> m;
        const auto s = cylinder(m, 1., 2.);
        const std::vector<point<double, 3>> rims{{std::cos(1.), std::sin(1.), 2.}, {std::cos(1.), std::sin(1.), 0.}};
        expect_blended(m, s, rims, 0.2, 0., 1e-6, "fillet_cylinder_rims");
        expect_blended(m, s, rims, -1., 0.15, 1e-6, "chamfer_cylinder_rims");
    }
    // a cylinder standing on a plate: its foot concave, its rim convex
    {
        Model<double> m;
        const auto plate = gbs_solid(m, make_prism(m, gbs_polygon(m, {{-2., -2., -1.}, {2., -2., -1.}, {2., 2., -1.}, {-2., 2., -1.}}, Z, X),
                                                   point<double, 3>{0., 0., 1.}));
        const auto fu = unwrap(fuse(m, ShapeId{plate}, ShapeId{cylinder(m, 0.5, 1.)}));
        REQUIRE(fu.solids.size() == 1u);
        expect_blended(m, fu.solids.front(), {{0.5 * std::cos(1.), 0.5 * std::sin(1.), 0.}, {0.5 * std::cos(1.), 0.5 * std::sin(1.), 1.}}, 0.15, 0., 1e-6,
                       "fillet_boss");
    }
    // the top rim of an oblong prism: a chain of four edges round it
    {
        Model<double> m;
        const double L = 2., R = 0.5;
        std::vector<EdgeId> e{unwrap(make_edge(m, point<double, 3>{0., -R, 0.}, point<double, 3>{L, -R, 0.})),
                              unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(R, R, ax2<double, 3>{{{L, 0., 0.}, Z, X}}, -std::numbers::pi / 2.,
                                                                                         std::numbers::pi / 2.))),
                              unwrap(make_edge(m, point<double, 3>{L, R, 0.}, point<double, 3>{0., R, 0.})),
                              unwrap(make_edge(m, std::make_shared<EllipseCurve<double>>(R, R, ax2<double, 3>{{{0., 0., 0.}, Z, X}}, std::numbers::pi / 2.,
                                                                                         3. * std::numbers::pi / 2.)))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{{0., 0., 0.}, Z, X}}, std::array<double, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, e))));
        const auto s = gbs_solid(m, make_prism(m, f, point<double, 3>{0., 0., 1.}));
        expect_blended(m, s, {{1., -R, 1.}, {L + R, 0., 1.}, {1., R, 1.}, {-R, 0., 1.}}, 0.1, 0., 1e-6, "fillet_oblong_chain");
    }
    // two crossed cylinders: a closed saddle edge, approximated fillets
    {
        Model<double> m;
        const auto lying = std::get<SolidId>(copy_shape(m, ShapeId{cylinder(m, 1., 4.)}, Placement<double>::frame({-2., 0., 0.}, Y, Z, X)));
        const auto fu = unwrap(fuse(m, ShapeId{lying}, ShapeId{cylinder(m, 0.4, 1.5)}));
        REQUIRE(fu.solids.size() == 1u);
        const point<double, 3> p{0.4 * std::cos(1.), 0.4 * std::sin(1.), std::sqrt(1. - 0.16 * std::sin(1.) * std::sin(1.))};
        expect_blended(m, fu.solids.front(), {p}, 0.15, 0., 1e-5, "fillet_crossed");
        expect_blended(m, fu.solids.front(), {p}, -1., 0.1, 1e-5, "chamfer_crossed");
    }
    // a tee: a cylinder of radius 0.99 standing on one of radius 1, the faces nearly tangent at
    // two points (OCCT cuts its fillet in more faces: 8 against 6). The tee of equal radii, where
    // the fillet vanishes, is not compared: OCCT's result is not valid there (16.98 against 16.19)
    {
        Model<double> m;
        const auto lying = std::get<SolidId>(copy_shape(m, ShapeId{cylinder(m, 1., 4.)}, Placement<double>::frame({-2., 0., 0.}, Y, Z, X)));
        const auto fu = unwrap(fuse(m, ShapeId{lying}, ShapeId{cylinder(m, 0.99, 2.)}));
        REQUIRE(fu.solids.size() == 1u);
        const point<double, 3> p{0.99 * std::cos(1.), 0.99 * std::sin(1.), std::sqrt(1. - 0.99 * 0.99 * std::sin(1.) * std::sin(1.))};
        expect_blended(m, fu.solids.front(), {p}, 0.05, 0., 1e-5, "fillet_tee_099", Faces::fewer);
    }
    // the top edge of the prism of a B-spline bump, square at its ends: free ends on planes
    {
        Model<double> m;
        const auto bump = std::make_shared<BSCurve<double, 3>>(
            points_vector<double, 3>{{0., 0., 0.}, {0.4, 0., 0.}, {0.8, -0.3, 0.}, {1.2, 0.2, 0.}, {1.6, 0., 0.}, {2., 0., 0.}},
            std::vector<double>{0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.}, std::size_t{3});
        const auto v0 = unwrap(make_vertex(m, point<double, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<double, 3>{2., 0., 0.}));
        const auto v2 = unwrap(make_vertex(m, point<double, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<double, 3>{0., 1., 0.}));
        const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                    unwrap(make_edge(m, v3, v0))};
        const auto pf = unwrap(make_face(m, std::make_shared<Plane<double>>(ax2<double, 3>{{{0., 0., 0.}, Z, X}}, std::array<double, 4>{-10., 10., -10., 10.}),
                                         unwrap(make_wire(m, w))));
        const auto s = gbs_solid(m, make_prism(m, pf, point<double, 3>{0., 0., 0.8}));
        const auto p = bump->value(0.5) + point<double, 3>{0., 0., 0.8};
        expect_blended(m, s, {p}, 0.1, 0., 1e-5, "fillet_bump_solid");
        expect_blended(m, s, {p}, -1., 0.05, 1e-5, "chamfer_bump_solid");
    }
}


// the corners of fillets (spheres) and the fillet at the root of a blade against
// BRepFilletAPI_MakeFillet
TEST(tests_brep2_occt, fillet_corners_and_blade)
{
    const point<double, 3> X{1., 0., 0.}, Z{0., 0., 1.};
    auto box = [&](Model<double> &m, double a, double b, double c, const point<double, 3> &o) {
        return gbs_solid(m, make_prism(m, gbs_polygon(m, {o, o + point<double, 3>{a, 0., 0.}, o + point<double, 3>{a, b, 0.}, o + point<double, 3>{0., b, 0.}}, Z, X),
                                       point<double, 3>{0., 0., c}));
    };
    // the three edges at a corner of a box, then all twelve
    {
        Model<double> m;
        const auto s = box(m, 1., 1., 1., {0., 0., 0.});
        expect_blended(m, s, {{0.5, 1., 1.}, {1., 0.5, 1.}, {1., 1., 0.5}}, 0.2, 0., 1e-6, "fillet_box_corner");
        std::vector<point<double, 3>> all;
        for (double a : {0., 1.})
            for (double b : {0., 1.})
                all.insert(all.end(), {{0.5, a, b}, {a, 0.5, b}, {a, b, 0.5}});
        expect_blended(m, s, all, 0.2, 0., 1e-6, "fillet_box_all");
    }
    // the inner corner of a pocket: three concave edges
    {
        Model<double> m;
        const auto block = box(m, 2., 2., 2., {0., 0., 0.});
        const auto pocket = box(m, 1., 1., 1., {1., 1., 1.});
        const auto cu = unwrap(cut(m, ShapeId{block}, ShapeId{pocket}));
        REQUIRE(cu.solids.size() == 1u);
        expect_blended(m, cu.solids.front(), {{1.5, 1., 1.}, {1., 1.5, 1.}, {1., 1., 1.5}}, 0.2, 0., 1e-6, "fillet_pocket_corner");
    }
    // the top corner of a triangular prism: faces at 60 degrees
    {
        Model<double> m;
        const double h = std::sqrt(3.) / 2.;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {1., 0., 0.}, {0.5, h, 0.}}, Z, X), point<double, 3>{0., 0., 1.}));
        expect_blended(m, s, {{0.75, h / 2., 1.}, {0.25, h / 2., 1.}, {0.5, h, 0.5}}, 0.1, 0., 1e-6, "fillet_prism_corner");
    }
    // the root of a blade on its hub (a closed B-spline face, not periodic, its sections joined
    // with C2 knots) filleted to the cylinder: OCCT's fillet is not valid (BRepCheck_Analyzer;
    // 7 faces, the volume 0.014 above the solid filleted where the fillet adds 3.8e-4), it is
    // only reported; gbs's read by OCCT valid, of the same volume and faces as in gbs
    {
        Model<double> m;
        const auto [s, root] = blade_fixture::blade_on_hub(m);
        REQUIRE(root.size() == 1u);
        const auto g = fillet(m, s, root, 0.03);
        REQUIRE(g.has_value());
        const auto back = to_occt(m, ShapeId{*g}, "fillet_blade_root");
        CHECK(BRepCheck_Analyzer(back).IsValid());
        CHECK(count(back, TopAbs_FACE) == 6);
        const double v0 = signed_volume(m, m.solid(s).outer, 256), v1 = signed_volume(m, m.solid(*g).outer, 256);
        CHECK(std::abs(volume(back) - v1) <= 1e-6);
        CHECK(v1 - v0 > 3e-4);
        CHECK(v1 - v0 < 4.2e-4);
        const auto ref = to_occt(m, ShapeId{s}, "fillet_blade_root_source");
        CHECK(BRepCheck_Analyzer(ref).IsValid());
        BRepFilletAPI_MakeFillet mk(ref);
        const auto &ed = m.edge(root.front());
        mk.Add(0.03, occt_edge_near(ref, ed.curve->value(ed.u1 + 0.3 * (ed.u2 - ed.u1))));
        mk.Build();
        std::string occt = "not done";
        if (mk.IsDone())
            occt = std::string{BRepCheck_Analyzer(mk.Shape()).IsValid() ? "valid" : "not valid"} + ", volume " + std::to_string(volume(mk.Shape()));
        MESSAGE("OCCT fillet of the blade root: " << occt << "; gbs: volume " << v1);
    }
}

// fillets of a radius varying linearly along an edge beside BRepFilletAPI_MakeFillet::Add(R1,
// R2, E): the edges of a box (straight: the parameter along the length), the radius from R1 at
// the first vertex of the edge to R2 at its second
TEST(tests_brep2_occt, variable_fillet_solids)
{
    const point<double, 3> X{1., 0., 0.}, Z{0., 0., 1.};
    struct Case
    {
        point<double, 3> at;
        double r1, r2;
        std::string name;
    };
    for (const auto &c : {Case{{1., 0., 1.}, 0.1, 0.3, "variable_box_top"}, Case{{2., 0.5, 0.}, 0.2, 0.05, "variable_box_bottom"}})
    {
        INFO(c.name);
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {0., 1., 0.}}, Z, X), point<double, 3>{0., 0., 1.}));
        const auto e = gbs_edge_near(m, s, c.at);
        const auto ref = to_occt(m, ShapeId{s}, c.name + "_source");
        const auto g = fillet(m, s, std::vector{e}, RadiusLaw<double>::linear(e, c.r1, c.r2));
        INFO((g ? std::string{} : g.error().message));
        REQUIRE(g.has_value());
        // the first vertex of gbs's edge on OCCT's edge
        const auto oe = occt_edge_near(ref, c.at);
        TopoDS_Vertex first, last;
        TopExp::Vertices(oe, first, last, true);
        const auto p1 = BRep_Tool::Pnt(first);
        const auto v1 = m.vertex(m.edge(e).v1).pnt;
        const bool same = norm(point<double, 3>{p1.X(), p1.Y(), p1.Z()} - v1) < 1e-9;
        BRepFilletAPI_MakeFillet mk(ref);
        mk.Add(same ? c.r1 : c.r2, same ? c.r2 : c.r1, oe);
        mk.Build();
        REQUIRE(mk.IsDone());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        // the result of gbs read by OCCT valid, as many faces; the volumes only reported: OCCT's
        // linear law is not the rolling ball of gbs (its volume checked against the union of the
        // balls in tests_brep_variable_fillet_solid) nor arcs of radius r(x) in the planes normal
        // to the edge (on the top edge, 0.1 to 0.3: OCCT takes off 0.01895, the arcs
        // (1 - pi / 4) int r^2 = 0.01860, the rolling ball 0.01854)
        const auto back = to_occt(m, ShapeId{*g}, c.name);
        CHECK(BRepCheck_Analyzer(back).IsValid());
        CHECK(count(back, TopAbs_FACE) == count(unified(mk.Shape()), TopAbs_FACE));
        CHECK(std::abs(volume(back) - signed_volume(m, m.solid(*g).outer, 256)) <= 1e-6);
        MESSAGE(c.name << ": volume taken off by OCCT " << 2. - volume(mk.Shape()) << ", by gbs " << 2. - volume(back));
    }
}

// blends meeting at a vertex (stage 6, PR 1) against BRepFilletAPI_MakeFillet and
// BRepFilletAPI_MakeChamfer: mitres of two fillets or chamfers at a corner (convex, concave),
// three chamfers at a corner, the eight edges of a pyramid (a sphere of four arcs at its top), one
// side edge of it (a free end at a vertex of four edges)
TEST(tests_brep2_occt, vertex_blends)
{
    const point<double, 3> X{1., 0., 0.}, Z{0., 0., 1.};
    auto box = [&](Model<double> &m, double a, double b, double c, const point<double, 3> &o) {
        return gbs_solid(m, make_prism(m, gbs_polygon(m, {o, o + point<double, 3>{a, 0., 0.}, o + point<double, 3>{a, b, 0.}, o + point<double, 3>{0., b, 0.}}, Z, X),
                                       point<double, 3>{0., 0., c}));
    };
    {
        Model<double> m;
        const auto s = box(m, 2., 1., 1.5, {0., 0., 0.});
        expect_blended(m, s, {{1., 0., 1.5}, {2., 0., 0.75}}, 0.2, 0., 1e-6, "fillet_mitre");
        expect_blended(m, s, {{1., 0., 1.5}, {2., 0., 0.75}}, -1., 0.2, 1e-6, "chamfer_mitre");
        // three chamfers at a corner: gbs mitres them (they meet at one point), OCCT cuts the
        // corner off by the plane through the three points where the contacts cross, a triangle
        // more: the tetrahedron d^3 / 12 below gbs's volume
        {
            const double d = 0.2;
            const std::vector<point<double, 3>> at{{1., 1., 1.5}, {2., 0.5, 1.5}, {2., 1., 0.75}};
            const auto ref = to_occt(m, s, "chamfer_corner_source");
            std::vector<EdgeId> edges;
            for (const auto &p : at)
                edges.push_back(gbs_edge_near(m, s, p));
            const auto g = unwrap(chamfer(m, s, edges, d));
            BRepFilletAPI_MakeChamfer mk(ref);
            for (const auto &p : at)
                mk.Add(d, occt_edge_near(ref, p));
            mk.Build();
            REQUIRE(mk.IsDone());
            const auto back = to_occt(m, ShapeId{g}, "chamfer_corner");
            CHECK(BRepCheck_Analyzer(back).IsValid());
            CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
            CHECK(std::abs(volume(back) - volume(mk.Shape()) - d * d * d / 12.) <= 1e-9);
            CHECK(count(back, TopAbs_FACE) + 1 == count(mk.Shape(), TopAbs_FACE));
        }
    }
    {
        Model<double> m;
        const auto block = box(m, 2., 2., 2., {0., 0., 0.});
        const auto pocket = box(m, 1., 1., 1., {1., 1., 1.});
        const auto cu = unwrap(cut(m, ShapeId{block}, ShapeId{pocket}));
        REQUIRE(cu.solids.size() == 1u);
        expect_blended(m, cu.solids.front(), {{1.5, 1., 1.}, {1., 1.5, 1.}}, 0.2, 0., 1e-6, "fillet_pocket_mitre");
    }
    {
        Model<double> m;
        const double h = 1.5;
        const std::array<point<double, 3>, 4> b{{{-1., -1., 0.}, {1., -1., 0.}, {1., 1., 0.}, {-1., 1., 0.}}};
        const point<double, 3> top{0., 0., h};
        auto face = [&](const std::vector<point<double, 3>> &p) {
            auto n = cross(p[1] - p[0], p[2] - p[0]);
            auto x = p[1] - p[0];
            return gbs_polygon(m, p, n / norm(n), x / norm(x));
        };
        std::vector<FaceId> faces{face({b[3], b[2], b[1], b[0]})};
        for (std::size_t i = 0; i < 4; ++i)
            faces.push_back(face({b[i], b[(i + 1) % 4], top}));
        const auto sw = unwrap(sew(m, faces, {.tol = 1e-9}));
        REQUIRE(sw.shells.size() == 1u);
        const auto s = unwrap(make_solid(m, sw.shells.front()));
        REQUIRE(signed_volume(m, m.solid(s).outer, 64) > 0.);
        std::vector<point<double, 3>> all;
        for (std::size_t i = 0; i < 4; ++i)
            all.insert(all.end(), {(b[i] + b[(i + 1) % 4]) / 2., (b[i] + top) / 2.});
        expect_blended(m, s, all, 0.1, 0., 1e-6, "fillet_pyramid");
        expect_blended(m, s, {(b[2] + top) / 2.}, 0.05, 0., 1e-6, "fillet_pyramid_side");
    }
    // a vertex blend (PR 3): the corner of a box, its three edges filleted at three radii. Both
    // valid, the same faces; OCCT's corner surface is not gbs's: gbs takes off 0.014 more (26.0370
    // against 26.0513, about 4 % of the volume of the corner)
    {
        Model<double> m;
        const auto s = box(m, 3., 3., 3., {0., 0., 0.});
        const std::vector<std::pair<point<double, 3>, double>> at{{{1.5, 3., 3.}, 0.5}, {{3., 1.5, 3.}, 0.7}, {{3., 3., 1.5}, 0.9}};
        const auto ref = to_occt(m, s, "vertex_blend_source");
        RadiusLaw<double> law;
        std::vector<EdgeId> edges;
        BRepFilletAPI_MakeFillet mk(ref);
        for (const auto &[p, r] : at)
        {
            edges.push_back(gbs_edge_near(m, s, p));
            law.points.push_back({edges.back(), 0.5, r});
            mk.Add(r, occt_edge_near(ref, p));
        }
        const auto g = unwrap(fillet(m, s, edges, law));
        mk.Build();
        REQUIRE(mk.IsDone());
        const auto back = to_occt(m, ShapeId{g}, "vertex_blend");
        CHECK(BRepCheck_Analyzer(back).IsValid());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        CHECK(count(back, TopAbs_FACE) == count(mk.Shape(), TopAbs_FACE));
        CHECK(volume(back) < volume(mk.Shape()));
        CHECK(volume(mk.Shape()) - volume(back) <= 0.02);
    }
}

TEST(tests_brep2_occt, variable_chamfers)
{
    // stage 7, PR 1: chamfers of two distances (AddDD) and of a distance and an angle (AddDA), at
    // right angles and at 120 degrees; gbs builds them through its laws (constant here, as OCCT
    // has none): the same volumes, the same faces, both valid
    const point<double, 3> X{1., 0., 0.}, Z{0., 0., 1.};
    const double c = std::cos(2. * std::numbers::pi / 3.), s3 = std::sin(2. * std::numbers::pi / 3.);
    struct Case
    {
        std::vector<point<double, 3>> base;
        point<double, 3> edge, face; // a point of the edge, a point of the face of d1
        bool angle;
        double d1, d2; // d2: the angle in radians
        const char *name;
    };
    const std::vector<Case> cases{
        {{{0., 0., 0.}, {3., 0., 0.}, {3., 1., 0.}, {0., 1., 0.}}, {1.5, 1., 1.}, {1.5, 0.5, 1.}, false, 0.1, 0.3, "chamfer_dd"},
        {{{0., 0., 0.}, {3., 0., 0.}, {3., 1., 0.}, {0., 1., 0.}}, {1.5, 1., 1.}, {1.5, 0.5, 1.}, true, 0.3, std::numbers::pi / 6., "chamfer_da"},
        {{{0., 0., 0.}, {3., 0., 0.}, {3., 2. * s3, 0.}, {2. * c, 2. * s3, 0.}}, {0., 0., 0.5}, {1., 0., 0.5}, true, 0.5, 40. * std::numbers::pi / 180., "chamfer_da_120"},
    };
    for (const auto &k : cases)
    {
        INFO(k.name);
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, k.base, Z, X), point<double, 3>{0., 0., 1.}));
        const auto ref = to_occt(m, s, std::string(k.name) + "_source");
        const auto e = gbs_edge_near(m, s, k.edge);
        auto bounds_edge = [&](FaceId x) {
            for (auto w : m.face(x).wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (ce.edge == e)
                        return true;
            return false;
        };
        FaceId f{};
        for (auto x : explore<FaceId>(m, ShapeId{s}))
            if (const auto *pl = dynamic_cast<const Plane<double> *>(m.face(x).surface.get()))
            {
                const auto uv = pl->parameters(k.face);
                if (norm(pl->value(uv[0], uv[1]) - k.face) < 1e-9 && bounds_edge(x))
                    f = x;
            }
        REQUIRE(f.valid());
        const auto g = k.angle ? unwrap(chamfer_angle(m, s, std::vector{e}, f, k.d1, k.d2))
                               : unwrap(chamfer(m, s, std::vector{e}, f, BlendLaw<double>::constant(e, k.d1), BlendLaw<double>::constant(e, k.d2)));
        BRepFilletAPI_MakeChamfer mk(ref);
        if (k.angle)
            mk.AddDA(k.d1, k.d2, occt_edge_near(ref, k.edge), occt_face_near(ref, k.face));
        else
            mk.Add(k.d1, k.d2, occt_edge_near(ref, k.edge), occt_face_near(ref, k.face));
        mk.Build();
        REQUIRE(mk.IsDone());
        const auto back = to_occt(m, ShapeId{g}, k.name);
        CHECK(BRepCheck_Analyzer(back).IsValid());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        CHECK(count(back, TopAbs_FACE) == count(mk.Shape(), TopAbs_FACE));
        CHECK(std::abs(volume(back) - volume(mk.Shape())) <= 1e-9);
    }
}

// fillets between two faces over a chamfer (stage 8, PR 1): OCCT has none; the fillet of the
// edge before the chamfer by BRepFilletAPI_MakeFillet is the same solid
TEST(tests_brep2_occt, fillet_faces_over_chamfers)
{
    const point<double, 3> X{1., 0., 0.}, Z{0., 0., 1.};
    auto face_where = [](const Model<double> &m, SolidId s, auto &&keep) {
        for (auto f : explore<FaceId>(m, ShapeId{s}))
        {
            bool all = true;
            for (auto v : explore<VertexId>(m, ShapeId{f}))
                all = all && keep(m.vertex(v).pnt);
            if (all)
                return f;
        }
        return FaceId{};
    };
    // the box 2 x 1 x 1, its edge y = 1, z = 1 chamfered at 0.1, r = 0.3
    {
        Model<double> m;
        const auto s = gbs_solid(m, make_prism(m, gbs_polygon(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {0., 1., 0.}}, Z, X), point<double, 3>{0., 0., 1.}));
        const auto ref = to_occt(m, s, "fillet_faces_box_source");
        const auto ch = unwrap(chamfer(m, s, std::vector<EdgeId>{gbs_edge_near(m, s, {1., 1., 1.})}, 0.1));
        const auto top = face_where(m, ch, [](const point<double, 3> &p) { return std::abs(p[2] - 1.) < 1e-12; });
        const auto side = face_where(m, ch, [](const point<double, 3> &p) { return std::abs(p[1] - 1.) < 1e-12; });
        const auto g = fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, 0.3);
        INFO((g ? std::string{} : g.error().message));
        REQUIRE(g.has_value());
        BRepFilletAPI_MakeFillet mk(ref);
        mk.Add(0.3, occt_edge_near(ref, {1., 1., 1.}));
        mk.Build();
        REQUIRE(mk.IsDone());
        CHECK(BRepCheck_Analyzer(mk.Shape()).IsValid());
        expect_same(m, ShapeId{*g}, mk.Shape(), "fillet_faces_box", 1e-9, Faces::same);
    }
}
