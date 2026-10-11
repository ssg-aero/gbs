// Comparison of the gbs STEP reader (gbs-io/step/reader.h) with OpenCascade's:
// shapes built by OCCT (primitives, booleans, fillets, sweeps) are written by
// STEPControl_Writer, read back by gbs and by STEPControl_Reader, and compared
// by their numbers of faces, edges and vertices and by their volumes; the gbs
// result must be complete and valid (check()).
#include <doctest_gtest.hpp>
#include <gbs-io/step/reader.h>

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
