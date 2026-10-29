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
#include <vector>

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

    // Parameter data of the entities of a type, by directory entry (columns 1-64 of the P lines, 'P' in column 73)
    std::vector<std::vector<std::string>> parameters(const std::string &file, int type)
    {
        std::map<int, std::string> text;
        std::ifstream in(file);
        std::string line;
        while (std::getline(in, line))
            if (line.size() >= 73 && line[72] == 'P')
                text[std::stoi(line.substr(64, 8))] += line.substr(0, 64);
        std::vector<std::vector<std::string>> r;
        for (const auto &[de, t] : text)
        {
            std::vector<std::string> fields{""};
            for (char c : t)
                if (c == ',' || c == ';')
                    fields.emplace_back();
                else if (c != ' ')
                    fields.back() += c;
            if (std::stoi(fields[0]) == type)
                r.push_back(fields);
        }
        return r;
    }

    // Global section field k (1-based), Hollerith strings kept as they are
    std::string global_field(const std::string &file, std::size_t k)
    {
        std::string g;
        std::ifstream in(file);
        std::string line;
        while (std::getline(in, line))
            if (line.size() >= 73 && line[72] == 'G')
                g += line.substr(0, 72);
        std::vector<std::string> fields{""};
        for (std::size_t i = 0; i < g.size(); ++i)
        {
            const char c = g[i];
            if (c == 'H' && !fields.back().empty() && std::ranges::all_of(fields.back(), [](char d) { return d >= '0' && d <= '9'; }))
            {
                const auto n = std::stoul(fields.back());
                fields.back() += g.substr(i, n + 1);
                i += n;
            }
            else if (c == ',' || c == ';')
                fields.emplace_back();
            else
                fields.back() += c;
        }
        return k - 1 < fields.size() ? fields[k - 1] : "";
    }

    // Largest coordinate of the control points of the 128 entities
    double largest_pole_coordinate(const std::string &file)
    {
        double r = 0.;
        for (const auto &f : parameters(file, 128))
        {
            const int k1 = std::stoi(f[1]), k2 = std::stoi(f[2]), m1 = std::stoi(f[3]), m2 = std::stoi(f[4]);
            const int knots = (k1 + m1 + 2) + (k2 + m2 + 2), poles = (k1 + 1) * (k2 + 1);
            for (int i = 0; i < 3 * poles; ++i)
            {
                auto x = f[std::size_t(10 + knots + poles + i)];
                std::ranges::replace(x, 'D', 'E');
                r = std::max(r, std::abs(std::stod(x)));
            }
        }
        return r;
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
    m.setName(explore<FaceId>(m, so)[2], "LID"); // a face with its own name keeps it as label
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

    // labels are in columns 57-64 of the second line of each directory entry, right-justified
    std::ifstream in(file);
    std::string line;
    int lid = 0, box = 0;
    while (std::getline(in, line))
        if (line.size() >= 73 && line[72] == 'D')
        {
            auto label = std::string_view{line}.substr(56, 8); // right-justified in its 8 columns
            label.remove_prefix(std::min(label.find_first_not_of(' '), label.size()));
            lid += label == "LID";
            box += label.starts_with("box");
        }
    ASSERT_EQ(lid, 1);
    ASSERT_EQ(box, 5);

    DLL_IGES back;
    ASSERT_TRUE(back.Read(file.c_str()));
}

TEST(tests_brep_iges, long_names)
{
    // a name longer than the 8 columns of a label (libIGES truncates them: "pressure_side_0001" and
    // "pressure_side_0002" both became "pressure") is also written whole, as a Name property (406, form 15)
    Model<T> m;
    auto so = box_solid(m);
    const auto faces = explore<FaceId>(m, so);
    m.setName(faces[0], "pressure_side_0001");
    m.setName(faces[1], "pressure_side_0002");
    m.setName(faces[2], "LID");
    IgesWriter<T> w;
    add_brep(w, m, so, "box");
    const auto file = tmp_file("gbs_brep_long_names.igs");
    w.write(file);
    ASSERT_EQ(entity_counts(file)[406], 2);
    // the parameter data: columns 1-64 of the 'P' lines
    std::ifstream in(file);
    std::string line, data;
    while (std::getline(in, line))
        if (line.size() >= 73 && line[72] == 'P')
            data += line.substr(0, 64);
    ASSERT_NE(data.find("18Hpressure_side_0001"), std::string::npos);
    ASSERT_NE(data.find("18Hpressure_side_0002"), std::string::npos);

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

TEST(tests_brep_iges, units_and_resolution)
{
    // a unit box in a model in inches: written in millimetres (the default unit of the file), or in inches when
    // the file is in inches; the coordinates were written in the unit of the model under the millimetre flag
    Model<T> m;
    auto so = box_solid(m);
    m.setUnitScale(25.4);
    for (const bool inches : {false, true})
    {
        INFO("file in inches " << inches);
        IgesWriter<T> w;
        if (inches)
            w.model().SetUnitsFlag(UNIT_INCH);
        add_brep(w, m, so, "box");
        const auto file = tmp_file("gbs_brep_units.igs");
        w.write(file);
        ASSERT_EQ(global_field(file, 14), inches ? "1" : "2");
        ASSERT_NEAR(largest_pole_coordinate(file), inches ? 1. : 25.4, 1e-12);
        // the minimum resolution covers the tolerances of the shape (1e-6 in), not the 1e-8 of libIGES
        auto res = global_field(file, 19);
        std::ranges::replace(res, 'D', 'E');
        ASSERT_GE(std::stod(res), (inches ? 1. : 25.4) * 1e-6 * (1. - 1e-9));
    }
}

TEST(tests_brep_iges, primitives_written_exactly)
{
    // a cylinder, its circles and generatrices as primitives: the cylinder and the circles are written as their
    // exact rational NURBS (they were interpolated), the pcurves carried to the parameters of the NURBS
    Model<T> m;
    const auto cyl = std::make_shared<CylindricalSurface<T>>(2., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., 2. * pi, 0., 3.});
    const auto f = unwrap(make_face(m, cyl));
    IgesWriter<T> w;
    auto rep = add_brep(w, m, f, "cyl");
    ASSERT_EQ(rep.approximated_surfaces, 0u);
    ASSERT_EQ(rep.over_tolerance, 0u);
    ASSERT_LE(rep.max_deviation, 1e-5);
    const auto file = tmp_file("gbs_brep_primitives.igs");
    w.write(file);
    // PROP3 (0: rational): the 128 entity, and the 126 of the circles
    const auto srf = parameters(file, 128);
    ASSERT_EQ(srf.size(), 1u);
    ASSERT_EQ(srf[0][7], "0");
    const auto crv = parameters(file, 126);
    ASSERT_EQ(std::ranges::count_if(crv, [](const auto &c) { return c[5] == "0"; }), 2);
}
