// Comparison of the gbs STEP reader (gbs-io/step/reader.h) with OpenCascade's:
// shapes built by OCCT (primitives, booleans, fillets, sweeps) are written by
// STEPControl_Writer, read back by gbs and by STEPControl_Reader, and compared
// by their numbers of faces, edges and vertices and by their volumes; the gbs
// result must be complete and valid (check()).
#include <doctest_gtest.hpp>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include <gbs/primitives.h>
#include <BRepCheck_Analyzer.hxx>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <TopoDS_Compound.hxx>
#include <TopLoc_Location.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <gp_Circ.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Trsf.hxx>

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
    std::string tmp_file(const std::string &name)
    {
        return (std::filesystem::temp_directory_path() / name).string();
    }

    std::string dump(const gbs::step::StepReadReport &r)
    {
        std::ostringstream os;
        for (const auto &i : r.unsupported)
            os << "unsupported #" << i.id << " " << i.type << ": " << i.message << "\n";
        for (const auto &i : r.failed)
            os << "failed #" << i.id << " " << i.type << ": " << i.message << "\n";
        for (const auto &e : r.check.entries)
            os << "check: " << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        os << "sense mismatches: " << r.sense_mismatches << "\n";
        return os.str();
    }

    struct Counts
    {
        std::size_t faces{}, edges{}, vertices{};
    };

    // faces, non-degenerate edges and vertices of an OCCT shape (STEP has no degenerate edges)
    Counts occt_counts(const TopoDS_Shape &s)
    {
        TopTools_IndexedMapOfShape f, e, v;
        TopExp::MapShapes(s, TopAbs_FACE, f);
        TopExp::MapShapes(s, TopAbs_EDGE, e);
        TopExp::MapShapes(s, TopAbs_VERTEX, v);
        std::size_t edges = 0;
        for (int i = 1; i <= e.Extent(); ++i)
            edges += !BRep_Tool::Degenerated(TopoDS::Edge(e(i)));
        return {static_cast<std::size_t>(f.Extent()), edges, static_cast<std::size_t>(v.Extent())};
    }

    double occt_volume(const TopoDS_Shape &s)
    {
        GProp_GProps p;
        BRepGProp::VolumeProperties(s, p, 1e-9, Standard_True);
        return p.Mass();
    }

    // non-degenerate edges of a gbs shape (degenerate ones are inserted at poles by the reader)
    std::size_t gbs_edges(const Model<double> &m, const ShapeId &s)
    {
        std::size_t n = 0;
        for (auto e : explore<EdgeId>(m, s))
            n += !m.edge(e).degenerate;
        return n;
    }

    double gbs_volume(const Model<double> &m, const ShapeId &s)
    {
        double v = 0.;
        for (auto so : explore<SolidId>(m, s))
        {
            v += signed_volume(m, m.solid(so).outer, 256);
            for (auto vo : m.solid(so).voids)
                v += signed_volume(m, vo, 256);
        }
        return v;
    }


    /// OCCT reads a STEP text written by gbs: the shape must be valid.
    TopoDS_Shape occt_read_text(const std::string &text, const std::string &name)
    {
        const auto file = tmp_file("gbs_written_" + name + ".stp");
        {
            std::ofstream f(file, std::ios::binary);
            f << text;
        }
        STEPControl_Reader rd;
        REQUIRE(rd.ReadFile(file.c_str()) == IFSelect_RetDone);
        rd.TransferRoots();
        return rd.OneShape();
    }

    /// Double round trip: the shape read by gbs is written by gbs and read by OCCT again.
    void check_written_by_gbs(const Model<double> &m, const ShapeId &root, const TopoDS_Shape &reference, const std::string &name,
                              double rel_volume_tol)
    {
        auto w = gbs::step::write_step(m, root);
        if (!w)
            FAIL(w.error().message);
        CHECK(w->second.approximated_curves == 0u);
        CHECK(w->second.approximated_surfaces == 0u);
        const auto back = occt_read_text(w->first, name);
        BRepCheck_Analyzer ana(back);
        CHECK(ana.IsValid());
        const auto a = occt_counts(reference), b = occt_counts(back);
        CHECK(b.faces == a.faces);
        CHECK(b.edges == a.edges);
        CHECK(b.vertices == a.vertices);
        const double vr = occt_volume(reference), vb = occt_volume(back);
        INFO("volume written by gbs, read by occt " << vb << " (" << vr << ")");
        CHECK(std::abs(vb - vr) <= rel_volume_tol * std::abs(vr));
    }

    /**
     * Writes `shape` to STEP with OCCT, reads it with gbs and OCCT, compares.
     */
    void compare(const TopoDS_Shape &shape, const std::string &name, double rel_volume_tol = 2e-3)
    {
        INFO(name);
        const auto file = tmp_file("gbs_occt_" + name + ".stp");
        STEPControl_Writer w;
        Interface_Static::SetCVal("write.step.schema", "AP214IS");
        REQUIRE(w.Transfer(shape, STEPControl_AsIs) == IFSelect_RetDone);
        REQUIRE(w.Write(file.c_str()) == IFSelect_RetDone);

        STEPControl_Reader rd;
        REQUIRE(rd.ReadFile(file.c_str()) == IFSelect_RetDone);
        rd.TransferRoots();
        const TopoDS_Shape occt = rd.OneShape();
        const auto oc = occt_counts(occt);

        Model<double> m;
        auto r = gbs::step::read_step_file(m, file);
        if (!r)
            FAIL(r.error().message);
        INFO(dump(r->report));
        CHECK(r->report.failed.empty());
        CHECK(r->report.unsupported.empty());
        CHECK(r->report.check.ok());
        CHECK(r->report.faces == oc.faces);
        CHECK(gbs_edges(m, r->root) == oc.edges);
        CHECK(r->report.vertices == oc.vertices);

        const double vo = occt_volume(occt), vg = gbs_volume(m, r->root);
        INFO("volume occt " << vo << " gbs " << vg);
        CHECK(std::abs(vg - vo) <= rel_volume_tol * std::abs(vo));
        TopTools_IndexedMapOfShape solids;
        TopExp::MapShapes(occt, TopAbs_SOLID, solids);
        CHECK(r->report.solids == static_cast<std::size_t>(solids.Extent()));
        // and back: written by gbs, read by OCCT
        check_written_by_gbs(m, r->root, occt, name, rel_volume_tol);
    }

    // closed planar profile in the xz plane: a B-spline arc closed by two segments and the axis
    TopoDS_Face profile_face()
    {
        TColgp_Array1OfPnt pts(1, 5);
        pts(1) = gp_Pnt(2., 0., 0.);
        pts(2) = gp_Pnt(3., 0., 1.);
        pts(3) = gp_Pnt(2.5, 0., 2.);
        pts(4) = gp_Pnt(3.2, 0., 3.);
        pts(5) = gp_Pnt(2., 0., 4.);
        Handle(Geom_BSplineCurve) bs = GeomAPI_PointsToBSpline(pts).Curve();
        BRepBuilderAPI_MakeWire mw;
        mw.Add(BRepBuilderAPI_MakeEdge(bs).Edge());
        mw.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(2., 0., 4.), gp_Pnt(1., 0., 4.)).Edge());
        mw.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1., 0., 4.), gp_Pnt(1., 0., 0.)).Edge());
        mw.Add(BRepBuilderAPI_MakeEdge(gp_Pnt(1., 0., 0.), gp_Pnt(2., 0., 0.)).Edge());
        return BRepBuilderAPI_MakeFace(mw.Wire(), Standard_True).Face();
    }
}

TEST(tests_step_occt, primitives)
{
    compare(BRepPrimAPI_MakeBox(gp_Pnt(-1., 2., 3.), 10., 20., 30.).Shape(), "box");
    compare(BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(1., 2., 3.), gp_Dir(1., 1., 0.)), 5., 10.).Shape(), "cylinder");
    compare(BRepPrimAPI_MakeSphere(gp_Pnt(1., -2., 0.5), 4.).Shape(), "sphere");
    compare(BRepPrimAPI_MakeCone(3., 0., 6.).Shape(), "cone");
    compare(BRepPrimAPI_MakeCone(3., 1., 6.).Shape(), "frustum");
    compare(BRepPrimAPI_MakeTorus(6., 2.).Shape(), "torus");
    compare(BRepPrimAPI_MakeCylinder(5., 10., std::numbers::pi / 2).Shape(), "cylinder_quarter");
    compare(BRepPrimAPI_MakeSphere(4., -0.5, 1.).Shape(), "sphere_zone");
}

TEST(tests_step_occt, booleans)
{
    const auto box = BRepPrimAPI_MakeBox(gp_Pnt(-5., -5., 0.), 10., 10., 10.).Shape();
    const auto cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., -1.), gp_Dir(0., 0., 1.)), 2., 12.).Shape();
    compare(BRepAlgoAPI_Cut(box, cyl).Shape(), "box_minus_cylinder");
    const auto sph = BRepPrimAPI_MakeSphere(gp_Pnt(5., 5., 10.), 3.).Shape();
    compare(BRepAlgoAPI_Cut(box, sph).Shape(), "box_minus_sphere");
    compare(BRepAlgoAPI_Fuse(box, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 10.), gp_Dir(0., 0., 1.)), 2., 4.).Shape()).Shape(), "box_plus_boss");
    // a cylinder crossing the seam of another one
    const auto c1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 5., 10.).Shape();
    const auto c2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-10., 0., 5.), gp_Dir(1., 0., 0.)), 1.5, 20.).Shape();
    compare(BRepAlgoAPI_Cut(c1, c2).Shape(), "cylinder_cross_hole");
}

TEST(tests_step_occt, fillets)
{
    const auto box = BRepPrimAPI_MakeBox(10., 8., 6.).Shape();
    BRepFilletAPI_MakeFillet fillet(box);
    for (TopExp_Explorer ex(box, TopAbs_EDGE); ex.More(); ex.Next())
        fillet.Add(1., TopoDS::Edge(ex.Current()));
    compare(fillet.Shape(), "box_fillets", 5e-3);
}

TEST(tests_step_occt, sweeps)
{
    const auto face = profile_face();
    compare(BRepPrimAPI_MakeRevol(face, gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.))).Shape(), "revolution");
    compare(BRepPrimAPI_MakeRevol(face, gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 2.).Shape(), "revolution_partial");
    compare(BRepPrimAPI_MakePrism(face, gp_Vec(0., 5., 0.)).Shape(), "prism");
}

TEST(tests_step_occt, harder_cases)
{
    // loft between a circle and a square: B-spline faces
    BRepOffsetAPI_ThruSections loft(Standard_True);
    loft.AddWire(BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 3.)).Edge()).Wire());
    loft.AddWire(BRepBuilderAPI_MakePolygon(gp_Pnt(-2., -2., 5.), gp_Pnt(2., -2., 5.), gp_Pnt(2., 2., 5.), gp_Pnt(-2., 2., 5.), Standard_True).Wire());
    loft.Build();
    compare(loft.Shape(), "loft", 5e-3);
    // tilted torus
    compare(BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(1., 2., 3.), gp_Dir(1., 2., 2.)), 6., 2.).Shape(), "torus_tilted");
    // sphere drilled through its poles: the hole meets the seam and the poles disappear
    const auto sph = BRepPrimAPI_MakeSphere(5.).Shape();
    compare(BRepAlgoAPI_Cut(sph, BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., -6.), gp_Dir(0., 0., 1.)), 1., 12.).Shape()).Shape(), "sphere_drilled");
    // revolution of a circle off the axis (surface of revolution of a conic)
    BRepBuilderAPI_MakeWire mw(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(gp_Pnt(4., 0., 0.), gp_Dir(0., 1., 0.)), 1.)).Edge());
    compare(BRepPrimAPI_MakeRevol(BRepBuilderAPI_MakeFace(mw.Wire()).Face(), gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.)), 3.).Shape(), "revolved_circle");
}

TEST(tests_step_occt, assembly)
{
    // a shared part instanced three times, one of them in a placed sub-assembly
    const auto part = BRepPrimAPI_MakeBox(4., 3., 2.).Shape();
    const auto pin = BRepPrimAPI_MakeCylinder(0.5, 6.).Shape();
    auto located = [](const TopoDS_Shape &s, gp_Vec t, double angle, gp_Dir axis) {
        gp_Trsf r, tr;
        r.SetRotation(gp_Ax1(gp_Pnt(0., 0., 0.), axis), angle);
        tr.SetTranslation(t);
        return s.Moved(TopLoc_Location(tr * r));
    };
    BRep_Builder b;
    TopoDS_Compound sub, top;
    b.MakeCompound(sub);
    b.Add(sub, part);
    b.Add(sub, located(pin, gp_Vec(2., 1.5, 2.), 0., gp_Dir(0., 0., 1.)));
    b.MakeCompound(top);
    b.Add(top, located(part, gp_Vec(10., 0., 0.), std::numbers::pi / 2, gp_Dir(0., 0., 1.)));
    b.Add(top, located(part, gp_Vec(0., 10., 0.), std::numbers::pi / 3, gp_Dir(1., 1., 0.)));
    b.Add(top, located(sub, gp_Vec(-10., -5., 3.), -std::numbers::pi / 4, gp_Dir(0., 1., 0.)));

    const auto file = tmp_file("gbs_occt_assembly.stp");
    STEPControl_Writer w;
    Interface_Static::SetCVal("write.step.schema", "AP214IS");
    Interface_Static::SetIVal("write.step.assembly", 1);
    REQUIRE(w.Transfer(top, STEPControl_AsIs) == IFSelect_RetDone);
    REQUIRE(w.Write(file.c_str()) == IFSelect_RetDone);
    Interface_Static::SetIVal("write.step.assembly", 0);

    STEPControl_Reader rd;
    REQUIRE(rd.ReadFile(file.c_str()) == IFSelect_RetDone);
    rd.TransferRoots();
    const TopoDS_Shape occt = rd.OneShape();

    Model<double> m;
    auto r = gbs::step::read_step_file(m, file);
    if (!r)
        FAIL(r.error().message);
    INFO(dump(r->report));
    CHECK(r->report.failed.empty());
    CHECK(r->report.check.ok());
    TopTools_IndexedMapOfShape solids;
    TopExp::MapShapes(occt, TopAbs_SOLID, solids);
    // MapShapes merges instances of a shared solid: count them by exploration
    std::size_t n_occt = 0;
    for (TopExp_Explorer ex(occt, TopAbs_SOLID); ex.More(); ex.Next())
        ++n_occt;
    CHECK(r->report.solids == n_occt);
    CHECK(n_occt == 4u);
    const double vo = occt_volume(occt), vg = gbs_volume(m, r->root);
    INFO("volume occt " << vo << " gbs " << vg);
    CHECK(std::abs(vg - vo) <= 2e-3 * vo);
    Bnd_Box bb;
    BRepBndLib::AddOptimal(occt, bb, Standard_False, Standard_False);
    double x0, y0, z0, x1, y1, z1;
    bb.Get(x0, y0, z0, x1, y1, z1);
    const auto gb = bounding_box(m, r->root, 37);
    INFO("occt box " << x0 << " " << y0 << " " << z0 << " / " << x1 << " " << y1 << " " << z1);
    INFO("gbs box " << gb.min[0] << " " << gb.min[1] << " " << gb.min[2] << " / " << gb.max[0] << " " << gb.max[1] << " " << gb.max[2]);
    for (auto [a, g] : {std::pair{x0, gb.min[0]}, {y0, gb.min[1]}, {z0, gb.min[2]}, {x1, gb.max[0]}, {y1, gb.max[1]}, {z1, gb.max[2]}})
        CHECK(std::abs(a - g) < 1e-3);

    // written back by gbs as an assembly sharing its parts, read by OCCT: same solids, volume, box
    auto gw = gbs::step::write_step(m, r->root);
    REQUIRE(gw.has_value());
    CHECK(gw->second.solids == 2u); // the box and the pin, each written once
    const auto back = occt_read_text(gw->first, "assembly");
    CHECK(BRepCheck_Analyzer(back).IsValid());
    std::size_t n_back = 0;
    for (TopExp_Explorer ex(back, TopAbs_SOLID); ex.More(); ex.Next())
        ++n_back;
    CHECK(n_back == n_occt);
    CHECK(std::abs(occt_volume(back) - vo) <= 2e-3 * vo);
    Bnd_Box bb2;
    BRepBndLib::AddOptimal(back, bb2, Standard_False, Standard_False);
    double u0, v0, w0, u1, v1, w1;
    bb2.Get(u0, v0, w0, u1, v1, w1);
    for (auto [a, g] : {std::pair{x0, u0}, {y0, v0}, {z0, w0}, {x1, u1}, {y1, v1}, {z1, w1}})
        CHECK(std::abs(a - g) < 1e-3);
}

TEST(tests_step_occt, shapes_built_in_gbs_read_by_occt)
{
    // a sewn box of bilinear B-spline faces, and solids of natural faces on primitives
    Model<double> m;
    auto quad = [](point<double, 3> a, point<double, 3> b, point<double, 3> c, point<double, 3> d) {
        std::vector<double> k{0., 0., 1., 1.};
        return std::make_shared<BSSurface<double, 3>>(points_vector<double, 3>{a, b, c, d}, k, k, 1, 1);
    };
    std::vector<FaceId> faces;
    auto corner = [](unsigned i) { return point<double, 3>{double(i & 1), double((i >> 1) & 1), double((i >> 2) & 1)}; };
    for (unsigned axis = 0; axis < 3; ++axis)
        for (unsigned side = 0; side < 2; ++side)
        {
            const unsigned a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
            auto idx = [&](unsigned u, unsigned v) { return (side << axis) | (u << a1) | (v << a2); };
            faces.push_back(unwrap(make_face(m, quad(corner(idx(0, 0)), corner(idx(1, 0)), corner(idx(0, 1)), corner(idx(1, 1))))));
        }
    const auto box = unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()));
    const double pi = std::numbers::pi;
    const ax2<double, 3> z{{{0., 0., 5.}, {0., 0., 1.}, {1., 0., 0.}}};
    auto natural_solid = [&](std::shared_ptr<Surface<double, 3>> s) {
        return unwrap(make_solid(m, m.add(Shell{{FaceUse{unwrap(make_face(m, s))}}})));
    };
    const auto ball = natural_solid(std::make_shared<SphericalSurface<double>>(2., z, std::array<double, 4>{0., 2 * pi, -pi / 2, pi / 2}));
    const auto torus = natural_solid(std::make_shared<ToroidalSurface<double>>(4., 1., z, std::array<double, 4>{0., 2 * pi, 0., 2 * pi}));
    // a closed cylinder: lateral natural face sewn to two disks is not built here; its lateral face is an open shell
    const auto tube = m.add(Shell{{FaceUse{unwrap(make_face(m, std::make_shared<CylindricalSurface<double>>(1., z, std::array<double, 4>{0., 2 * pi, 0., 3.})))}}});
    struct Case
    {
        ShapeId s;
        double volume;
        int faces;
    };
    for (auto [s, volume, nfaces] : {Case{box, 1., 6}, Case{ball, 4. / 3. * pi * 8., 1}, Case{torus, 2 * pi * pi * 4., 1}, Case{tube, 0., 1}})
    {
        auto w = gbs::step::write_step(m, s);
        REQUIRE(w.has_value());
        const auto back = occt_read_text(w->first, "gbs_built");
        CHECK(BRepCheck_Analyzer(back).IsValid());
        TopTools_IndexedMapOfShape f;
        TopExp::MapShapes(back, TopAbs_FACE, f);
        CHECK(f.Extent() == nfaces);
        if (volume > 0.)
            CHECK(std::abs(occt_volume(back) - volume) <= 1e-3 * volume);
    }
}
