// Round trip of IgesWriter::add_geometry through OpenCascade: rational and 2D
// NURBS written by gbs must be read back with the right geometry.
#include <doctest_gtest.hpp>
#include <gbs-io/iges.h>
#include <gbs/bscbuild.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <IGESControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <cmath>
#include <filesystem>

using namespace gbs;

namespace
{
    std::string tmp_file(const std::string &name)
    {
        return (std::filesystem::temp_directory_path() / name).string();
    }

    TopoDS_Shape read_back(const std::string &file)
    {
        IGESControl_Reader reader;
        REQUIRE(reader.ReadFile(file.c_str()) == IFSelect_RetDone);
        reader.TransferRoots();
        return reader.OneShape();
    }

    // every sampled point of every edge of the shape lies at `radius` from (cx, cy) in the plane z = 0
    void check_circle(const TopoDS_Shape &shape, double cx, double cy, double radius)
    {
        int edges = 0;
        for (TopExp_Explorer ex(shape, TopAbs_EDGE); ex.More(); ex.Next(), ++edges)
        {
            BRepAdaptor_Curve c(TopoDS::Edge(ex.Current()));
            for (int i = 0; i <= 20; ++i)
            {
                const auto p = c.Value(c.FirstParameter() + (c.LastParameter() - c.FirstParameter()) * i / 20.);
                ASSERT_NEAR(std::hypot(p.X() - cx, p.Y() - cy), radius, 1e-7);
                ASSERT_NEAR(p.Z(), 0., 1e-12);
            }
        }
        ASSERT_GE(edges, 1);
    }
}

TEST(tests_iges_writer_occt, rational_circle_3d)
{
    IgesWriter<double> w;
    w.model().SetUnitsFlag(UNIT_MILLIMETER);
    w.add_geometry(build_circle<double, 3>(2., {1., 0., 0.}));
    const auto f = tmp_file("gbs_writer_circle3d_occt.igs");
    w.write(f);
    check_circle(read_back(f), 1., 0., 2.);
}

TEST(tests_iges_writer_occt, rational_circle_2d)
{
    IgesWriter<double> w;
    w.model().SetUnitsFlag(UNIT_MILLIMETER);
    w.add_geometry(build_circle<double, 2>(1.5, {0.5, -1.}));
    const auto f = tmp_file("gbs_writer_circle2d_occt.igs");
    w.write(f);
    check_circle(read_back(f), 0.5, -1., 1.5);
}

TEST(tests_iges_writer_occt, rational_surface)
{
    const auto circle = build_circle<double, 3>(1.);
    points_vector<double, 4> poles;
    for (double z : {0., 2.})
        for (auto p : circle.poles())
        {
            p[2] += p[3] * z;
            poles.push_back(p);
        }
    IgesWriter<double> w;
    w.model().SetUnitsFlag(UNIT_MILLIMETER);
    w.add_geometry(BSSurfaceRational<double, 3>{poles, circle.knotsFlats(), std::vector<double>{0., 0., 2., 2.}, circle.degree(), 1});
    const auto f = tmp_file("gbs_writer_cylinder_occt.igs");
    w.write(f);

    int faces = 0;
    for (TopExp_Explorer ex(read_back(f), TopAbs_FACE); ex.More(); ex.Next(), ++faces)
    {
        BRepAdaptor_Surface s(TopoDS::Face(ex.Current()));
        for (int i = 0; i <= 10; ++i)
            for (int j = 0; j <= 4; ++j)
            {
                const double u = s.FirstUParameter() + (s.LastUParameter() - s.FirstUParameter()) * i / 10.;
                const double v = s.FirstVParameter() + (s.LastVParameter() - s.FirstVParameter()) * j / 4.;
                const auto p = s.Value(u, v);
                ASSERT_NEAR(std::hypot(p.X(), p.Y()), 1., 1e-7); // a true cylinder of radius 1
            }
    }
    ASSERT_EQ(faces, 1);
}
