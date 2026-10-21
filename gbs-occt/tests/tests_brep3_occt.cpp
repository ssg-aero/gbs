// Stage 3 of the native BREP kernel checked by OpenCascade: the meshes of gbs (gbs::brep::mesh)
// of shapes built by OCCT (passed through STEP), and of booleans of gbs on them, against the
// meshes of BRepMesh_IncrementalMesh at the same deflection and angle: both closed, their
// volumes close to the exact volume (BRepGProp), their triangles within the deflection of the
// faces of OCCT (sampled by BRepExtrema), and their numbers of triangles of the same order.
#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBndLib.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax2.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

using namespace gbs;
using namespace gbs::brep;

namespace
{
    std::string tmp_file(const std::string &name) { return (std::filesystem::temp_directory_path() / name).string(); }

    /// An OCCT shape read by gbs, through a STEP file written by OCCT.
    ShapeId to_gbs(Model<double> &m, const TopoDS_Shape &s, const std::string &name)
    {
        const auto file = tmp_file("gbs_brep3_" + name + ".stp");
        STEPControl_Writer w;
        Interface_Static::SetCVal("write.step.schema", "AP214IS");
        REQUIRE(w.Transfer(s, STEPControl_AsIs) == IFSelect_RetDone);
        REQUIRE(w.Write(file.c_str()) == IFSelect_RetDone);
        auto r = gbs::step::read_step_file(m, file);
        REQUIRE(r.has_value());
        REQUIRE(r->report.check.ok());
        return r->root;
    }

    double volume(const TopoDS_Shape &s)
    {
        GProp_GProps p;
        BRepGProp::VolumeProperties(s, p, 1e-9, Standard_True);
        return p.Mass();
    }

    double diagonal(const TopoDS_Shape &s)
    {
        Bnd_Box b;
        BRepBndLib::AddOptimal(s, b, Standard_False);
        return std::sqrt(b.SquareExtent());
    }

    /// A triangle mesh as points and oriented triangles (indices), for both meshers.
    struct Triangles
    {
        std::vector<point<double, 3>> points;
        std::vector<std::array<std::uint32_t, 3>> triangles;
        double seconds{};

        double volume() const
        {
            double v{};
            for (const auto &t : triangles)
                v += points[t[0]] * cross(points[t[1]], points[t[2]]);
            return v / 6.;
        }
        double area() const
        {
            double a{};
            for (const auto &t : triangles)
                a += norm(cross(points[t[1]] - points[t[0]], points[t[2]] - points[t[0]]));
            return a / 2.;
        }
        /// Closed: the vertices of OCCT are merged by position first (it has a node per face).
        std::size_t open_sides() const
        {
            std::map<std::array<double, 3>, std::uint32_t> merged;
            std::vector<std::uint32_t> id(points.size());
            for (std::size_t k = 0; k < points.size(); ++k)
                id[k] = merged.emplace(points[k], static_cast<std::uint32_t>(merged.size())).first->second;
            std::vector<std::array<std::uint32_t, 3>> t;
            for (const auto &x : triangles)
                if (id[x[0]] != id[x[1]] && id[x[1]] != id[x[2]] && id[x[2]] != id[x[0]])
                    t.push_back({id[x[0]], id[x[1]], id[x[2]]});
            return brep::open_sides(t);
        }
        /// The largest distance of the centroids of up to n triangles (spread) to the faces of s.
        double deviation(const TopoDS_Shape &s, std::size_t n = 300) const
        {
            TopoDS_Compound faces;
            BRep_Builder b;
            b.MakeCompound(faces);
            for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
                b.Add(faces, e.Current());
            double d{};
            const std::size_t step = std::max<std::size_t>(1, triangles.size() / n);
            for (std::size_t k = 0; k < triangles.size(); k += step)
            {
                const auto &t = triangles[k];
                const auto g = (points[t[0]] + points[t[1]] + points[t[2]]) / 3.;
                BRepExtrema_DistShapeShape dist(BRepBuilderAPI_MakeVertex(gp_Pnt(g[0], g[1], g[2])).Vertex(), faces);
                REQUIRE(dist.IsDone());
                d = std::max(d, dist.Value());
            }
            return d;
        }
    };

    Triangles of_gbs(const Model<double> &m, const ShapeId &s, double deflection, double angle)
    {
        const auto t0 = std::chrono::steady_clock::now();
        const auto mesh = unwrap(brep::mesh(m, s, MeshOptions<double>{.deflection = deflection, .angle = angle}));
        Triangles r{mesh.vertices, mesh.triangles};
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        REQUIRE(mesh.incomplete.empty());
        REQUIRE(mesh.folded.empty());
        return r;
    }

    Triangles of_occt(const TopoDS_Shape &shape, double deflection, double angle)
    {
        TopoDS_Shape s = shape; // (BRepMesh stores the triangulations in the shape)
        const auto t0 = std::chrono::steady_clock::now();
        BRepMesh_IncrementalMesh mesher(s, deflection, Standard_False, angle, Standard_False);
        REQUIRE(mesher.IsDone());
        Triangles r;
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next())
        {
            const TopoDS_Face &f = TopoDS::Face(e.Current());
            TopLoc_Location loc;
            const auto tri = BRep_Tool::Triangulation(f, loc);
            REQUIRE(!tri.IsNull());
            const auto first = static_cast<std::uint32_t>(r.points.size());
            for (int i = 1; i <= tri->NbNodes(); ++i)
            {
                const gp_Pnt p = tri->Node(i).Transformed(loc.Transformation());
                r.points.push_back({p.X(), p.Y(), p.Z()});
            }
            for (int i = 1; i <= tri->NbTriangles(); ++i)
            {
                int a, b, c;
                tri->Triangle(i).Get(a, b, c);
                if (f.Orientation() == TopAbs_REVERSED)
                    std::swap(b, c);
                r.triangles.push_back({first + a - 1, first + b - 1, first + c - 1});
            }
        }
        return r;
    }
}

TEST(tests_brep3_occt, meshes)
{
    const TopoDS_Shape box = BRepPrimAPI_MakeBox(2., 3., 4.).Shape();
    const TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(gp_Ax2(gp_Pnt(2., 3., 4.), gp_Dir(1., 1., 1.)), 1.5).Shape();
    const TopoDS_Shape cyl1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0., 0., -3.), gp_Dir(0., 0., 1.)), 1., 6.).Shape();
    const TopoDS_Shape cyl2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-3., 0., 0.), gp_Dir(1., 0., 0.), gp_Dir(0., 0., 1.)), 0.6, 6.).Shape();
    struct Case
    {
        std::string name;
        TopoDS_Shape a, b{};      ///< b: the second operand of a boolean of gbs
        BooleanOp op{BooleanOp::Fuse};
        TopoDS_Shape reference{}; ///< the boolean of OCCT (a when none)
    };
    const std::vector<Case> cases{
        {"box", box},
        {"cylinder", BRepPrimAPI_MakeCylinder(5., 10.).Shape()},
        {"sphere", BRepPrimAPI_MakeSphere(3.).Shape()},
        {"cone", BRepPrimAPI_MakeCone(2., 0., 3.).Shape()},
        {"torus", BRepPrimAPI_MakeTorus(3., 1.).Shape()},
        {"box_minus_sphere", box, sphere, BooleanOp::Cut, BRepAlgoAPI_Cut(box, sphere).Shape()},
        {"crossed_cylinders", cyl1, cyl2, BooleanOp::Fuse, BRepAlgoAPI_Fuse(cyl1, cyl2).Shape()},
    };
    for (const auto &c : cases)
    {
        INFO(c.name);
        Model<double> m;
        ShapeId s = to_gbs(m, c.a, c.name + "_a");
        if (!c.b.IsNull())
            s = unwrap(boolean(m, c.op, s, to_gbs(m, c.b, c.name + "_b"))).shape;
        const TopoDS_Shape &ref = c.reference.IsNull() ? c.a : c.reference;
        const double exact = volume(ref);
        for (double relative : {1e-2, 1e-3})
        {
            const double deflection = relative * diagonal(ref), angle = 0.5;
            INFO("deflection " << deflection);
            const auto g = of_gbs(m, s, deflection, angle);
            const auto o = of_occt(ref, deflection, angle);
            CHECK(g.open_sides() == 0u);
            CHECK(o.open_sides() == 0u);
            const double dg = g.deviation(ref), d_occt = o.deviation(ref);
            MESSAGE(c.name << " " << relative << ": gbs " << g.triangles.size() << " triangles, " << g.seconds << " s, volume " << (g.volume() - exact) / exact
                           << ", deviation " << dg / deflection << " | occt " << o.triangles.size() << " triangles, " << o.seconds << " s, volume "
                           << (o.volume() - exact) / exact << ", deviation " << d_occt / deflection);
            // the volume: within the area times the deflection, inscribed in the convex parts
            CHECK(std::abs(g.volume() - exact) <= g.area() * deflection);
            // the centroids within the deflection (the vertices are on the faces)
            CHECK(dg <= 1.05 * deflection);
            // the same order of triangles
            CHECK(g.triangles.size() <= 6 * o.triangles.size() + 12);
            CHECK(o.triangles.size() <= 6 * g.triangles.size() + 12);
        }
    }
}
