#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-io/iges_brep.h>
#include "tests_brep_geom.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <string>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using namespace brep_tests;

    // Entity type -> count, read from the Directory Entry section of an IGES file
    // (fixed 80-column lines, 'D' in column 73, two lines per entity, type in columns 1-8).
    std::map<int, int> entity_counts(const std::string &file)
    {
        std::map<int, int> counts;
        std::ifstream in(file);
        std::string line;
        bool first = true;
        while (std::getline(in, line))
        {
            if (line.size() < 73 || line[72] != 'D')
                continue;
            if (first)
                ++counts[std::stoi(line.substr(0, 8))];
            first = !first;
        }
        return counts;
    }

    std::string tmp_file(const std::string &name)
    {
        return (std::filesystem::temp_directory_path() / name).string();
    }

    std::shared_ptr<Surface<T, 3>> plane_z(T z)
    {
        return quad({-2., -2., z}, {2., -2., z}, {-2., 2., z}, {2., 2., z});
    }

    SolidId box_solid(Model<T> &m)
    {
        auto r = unwrap(sew(m, box_faces(m)));
        return unwrap(make_solid(m, r.shells.front()));
    }

    // Cylinder closed by two disks, all exact rational NURBS: the rational cylinder and rational circles.
    SolidId nurbs_cylinder_solid(Model<T> &m, T h = 2.)
    {
        auto lateral = unwrap(make_face(m, nurbs_cylinder(h)));
        auto disk = [&](T z) {
            auto rc = build_circle<T, 3>(1.);
            auto poles = rc.poles();
            for (auto &p : poles)
                p[2] += p[3] * z;
            auto c = std::make_shared<BSCurveRational<T, 3>>(poles, rc.knotsFlats(), rc.degree());
            return unwrap(make_face(m, plane_z(z), unwrap(make_wire(m, std::vector{unwrap(make_edge(m, c))}))));
        };
        auto r = unwrap(sew(m, std::vector{lateral, disk(0.), disk(h)}, {.tol = 1e-5}));
        return unwrap(make_solid(m, r.shells.front()));
    }
}

TEST(tests_brep_iges, box_solid)
{
    Model<T> m;
    auto so = box_solid(m);
    IgesWriter<T> w;
    auto rep = add_brep(w, m, so, "box");
    ASSERT_EQ(rep.faces, 6);
    ASSERT_EQ(rep.wires, 6);
    ASSERT_EQ(rep.free_edges, 0);
    ASSERT_EQ(rep.approximated_surfaces, 0); // bilinear B-spline faces
    // the pcurves rebuilt by the sewing are B-splines; the edge curves are CurveOnSurface: approximated
    ASSERT_LE(rep.max_deviation, IgesExportOptions<T>{}.approx_tol);

    const auto file = tmp_file("gbs_brep_box.igs");
    w.write(file);
    auto n = entity_counts(file);
    ASSERT_EQ(n[144], 6);
    ASSERT_EQ(n[128], 6);
    ASSERT_EQ(n[142], 6);
    ASSERT_EQ(n[102], 12);     // parameter space and model space of each wire
    ASSERT_EQ(n[126], 6 * 8);  // 4 co-edges per face, a 2D and a 3D curve each

    DLL_IGES back;
    ASSERT_TRUE(back.Read(file.c_str()));
}

TEST(tests_brep_iges, curved_and_trimmed_faces)
{
    Model<T> m;
    // natural faces with seams and poles: an analytic sphere (approximated) and an exact NURBS cylinder
    auto sphere_face = unwrap(make_face(m, sphere(1.)));
    auto cyl_face = unwrap(make_face(m, nurbs_cylinder()));
    // a trimmed face with a hole, from 2D loops
    auto seg = [](point<T, 2> a, point<T, 2> b) -> std::shared_ptr<Curve<T, 2>> { return std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>(a, b)); };
    auto holed = unwrap(make_face(m, plane_z(0.),
                                  std::vector{seg({0.1, 0.1}, {0.9, 0.1}), seg({0.9, 0.1}, {0.9, 0.9}), seg({0.9, 0.9}, {0.1, 0.9}), seg({0.1, 0.9}, {0.1, 0.1})},
                                  {std::vector{seg({0.4, 0.4}, {0.4, 0.6}), seg({0.4, 0.6}, {0.6, 0.4}), seg({0.6, 0.4}, {0.4, 0.4})}}));
    // a free edge
    auto free_e = unwrap(make_edge(m, point<T, 3>{5., 0., 0.}, point<T, 3>{6., 1., 0.}));
    auto c = unwrap(make_compound(m, std::vector<ShapeId>{sphere_face, cyl_face, holed, free_e}));

    IgesWriter<T> w;
    auto rep = add_brep(w, m, c, "", IgesExportOptions<T>{.approx_tol = 1e-5});
    ASSERT_EQ(rep.faces, 3);
    ASSERT_EQ(rep.wires, 4);
    ASSERT_EQ(rep.free_edges, 1);
    ASSERT_EQ(rep.approximated_surfaces, 1); // the analytic sphere only
    ASSERT_GT(rep.approximated_curves, 0);
    ASSERT_LE(rep.max_deviation, 1e-5);

    const auto file = tmp_file("gbs_brep_mixed.igs");
    w.write(file);
    auto n = entity_counts(file);
    ASSERT_EQ(n[144], 3);
    ASSERT_EQ(n[142], 4);
    ASSERT_EQ(n[102], 8);
    // co-edges: sphere 4, cylinder 4, holed 4 + 3 -> 15, two 126 each, plus the free edge
    ASSERT_EQ(n[126], 15 * 2 + 1);
    DLL_IGES back;
    ASSERT_TRUE(back.Read(file.c_str()));
}

TEST(tests_brep_iges, rational_solid_and_scale)
{
    Model<T> m;
    auto so = nurbs_cylinder_solid(m);
    IgesWriter<T> w;
    auto rep = add_brep(w, m, so, "cyl", IgesExportOptions<T>{.scale = 1000.});
    ASSERT_EQ(rep.faces, 3);
    ASSERT_EQ(rep.approximated_surfaces, 0); // rational cylinder and planes written exactly
    const auto file = tmp_file("gbs_brep_cylinder.igs");
    w.write(file);
    auto n = entity_counts(file);
    ASSERT_EQ(n[144], 3);
    ASSERT_EQ(n[128], 3);
    DLL_IGES back;
    ASSERT_TRUE(back.Read(file.c_str()));
}
