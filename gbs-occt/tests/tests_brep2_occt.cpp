// Stage 2 of the native BREP kernel checked by OpenCascade: booleans, sections and sweeps of
// gbs on shapes built by OCCT (passed through STEP) are compared with BRepAlgoAPI and
// BRepPrimAPI; the results of gbs are written to STEP by gbs and read back by OCCT, which must
// find them valid, with the same volumes, and as many faces as the results of OCCT with their
// faces of one domain merged (ShapeUpgrade_UnifySameDomain).
#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>

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
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
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
        BRepGProp::VolumeProperties(s, p, 1e-9, Standard_True);
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
        TopTools_IndexedMapOfShape mp;
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
        TColgp_Array1OfPnt poles(1, 5);
        const double gx[]{-1., -2. / 3., 0., 2. / 3., 1.}, zx[]{0., 0., 0., 0., -bend};
        for (int i = 0; i < 5; ++i)
            poles.SetValue(i + 1, gp_Pnt(gx[i], -1., zx[i]));
        TColStd_Array1OfReal knots(1, 3);
        TColStd_Array1OfInteger mults(1, 3);
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
