#include <doctest_gtest.hpp>
#include <gbs-io/step/reader.h>

#include <cmath>
#include <numbers>
#include <sstream>
#include <string>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;
using namespace gbs::step;

namespace
{
    constexpr double pi = std::numbers::pi;

    std::string file(const std::string &name) { return std::string(GBS_TESTS_IN_DIR) + "/step/" + name; }

    std::string dump(const StepReadReport &r)
    {
        std::ostringstream os;
        for (const auto &i : r.unsupported)
            os << "unsupported #" << i.id << " " << i.type << ": " << i.message << "\n";
        for (const auto &i : r.failed)
            os << "failed #" << i.id << " " << i.type << ": " << i.message << "\n";
        for (const auto &e : r.check.entries)
            os << "check: " << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        return os.str();
    }

    StepReadResult read_ok(Model<double> &m, const std::string &name, StepReadOptions o = {})
    {
        auto r = read_step_file(m, file(name), o);
        if (!r)
            FAIL(r.error().message);
        INFO(dump(r->report));
        REQUIRE(r->report.failed.empty());
        REQUIRE(r->report.unsupported.empty());
        REQUIRE(r->report.check.ok());
        return std::move(*r);
    }

    double volume(const Model<double> &m, const ShapeId &s)
    {
        return signed_volume(m, m.solid(std::get<SolidId>(s)).outer, 256);
    }
}

TEST(tests_step_read, box)
{
    Model<double> m;
    auto r = read_ok(m, "box.stp");
    ASSERT_TRUE(std::holds_alternative<SolidId>(r.root));
    ASSERT_EQ(r.report.solids, 1u);
    ASSERT_EQ(r.report.faces, 6u);
    ASSERT_EQ(r.report.edges, 12u);
    ASSERT_EQ(r.report.vertices, 8u);
    ASSERT_EQ(r.report.sense_mismatches, 0u); // three planes against the outward normal, same_sense .F.
    ASSERT_NEAR(volume(m, r.root), 24., 1e-9);
    // names and STEP ids
    ASSERT_EQ(m.name(r.root), "BOX");
    ASSERT_TRUE(m.externalId(r.root).has_value());
    std::size_t named = 0;
    for (auto f : explore<FaceId>(m, r.root))
    {
        ASSERT_TRUE(m.externalId(f).has_value());
        named += m.name(f) == "TOP";
    }
    ASSERT_EQ(named, 1u);
    ASSERT_EQ(r.report.length_unit, "MILLI METRE");
    ASSERT_NEAR(*r.report.file_tolerance, 1e-7, 1e-20);
    ASSERT_EQ(m.unitScale(), 1.);
    ASSERT_TRUE(r.report.schemas.front().starts_with("AUTOMOTIVE_DESIGN"));
}

TEST(tests_step_read, curved_solids)
{
    // cylinder: seam used twice; sphere: seam only, poles inserted; cone: apex inserted, angles in degrees
    struct Case
    {
        const char *name;
        double volume, tol;
        std::size_t faces, degenerate;
    };
    const double H = 2. / std::tan(pi / 6);
    for (const auto &c : {Case{"cylinder.stp", pi * 25. * 10., 0.2, 3, 0}, Case{"sphere.stp", 4. / 3. * pi * 27., 0.2, 1, 2},
                          Case{"cone.stp", pi * 4. * H / 3., 0.05, 2, 1}})
    {
        INFO(c.name);
        Model<double> m;
        auto r = read_ok(m, c.name);
        ASSERT_EQ(r.report.faces, c.faces);
        ASSERT_NEAR(volume(m, r.root), c.volume, c.tol);
        std::size_t deg = 0;
        for (auto e : explore<EdgeId>(m, r.root))
            deg += m.edge(e).degenerate;
        ASSERT_EQ(deg, c.degenerate);
    }
}

TEST(tests_step_read, inches_hole_and_open_shell)
{
    Model<double> m;
    auto r = read_ok(m, "plate_inch.stp");
    ASSERT_TRUE(std::holds_alternative<ShellId>(r.root)); // shell based surface model: an open shell
    ASSERT_EQ(r.report.length_unit, "INCH");
    ASSERT_NEAR(r.report.length_factor, 25.4, 1e-12);
    ASSERT_NEAR(*r.report.file_tolerance, 25.4e-7, 1e-18);
    auto faces = explore<FaceId>(m, r.root);
    ASSERT_EQ(faces.size(), 1u);
    const auto &f = m.face(faces[0]);
    ASSERT_EQ(f.wires.size(), 2u); // outer and hole
    const double in2 = 25.4 * 25.4;
    ASSERT_NEAR(uv_signed_area(m, f.wires[0]), 8. * in2, 1e-6);           // counter-clockwise, in mm
    ASSERT_NEAR(uv_signed_area(m, f.wires[1], 2048), -pi * 0.25 * in2, 1e-3); // clockwise (bound orientation .F.)
    ASSERT_EQ(m.name(faces[0]), "PLATE");
    ASSERT_EQ(m.name(r.root), "PLATE");

    // the same file read in inches: coordinates kept
    Model<double> mi;
    auto ri = read_ok(mi, "plate_inch.stp", {.target_unit_mm = 25.4});
    ASSERT_EQ(mi.unitScale(), 25.4);
    ASSERT_NEAR(ri.report.length_factor, 1., 1e-15);
    ASSERT_NEAR(uv_signed_area(mi, mi.face(explore<FaceId>(mi, ri.root)[0]).wires[0]), 8., 1e-9);
}

TEST(tests_step_read, rational_bspline_face)
{
    Model<double> m;
    auto r = read_ok(m, "bspline_patch.stp");
    auto faces = explore<FaceId>(m, r.root);
    ASSERT_EQ(faces.size(), 1u);
    // the surface is the quarter cylinder of radius 1
    const auto &s = *m.face(faces[0]).surface;
    for (double u : {0., 0.3, 0.7, 1.})
        for (double v : {0., 0.5, 1.})
        {
            const auto p = s(u, v);
            ASSERT_NEAR(std::hypot(p[0], p[1]), 1., 1e-12);
        }
    ASSERT_EQ(r.report.edges, 4u);
}

TEST(tests_step_read, partial_and_strict)
{
    // the top face of the box lies on an unsupported surface
    Model<double> m;
    auto r = read_step_file(m, file("box_unsupported_face.stp"));
    ASSERT_TRUE(r.has_value());
    INFO(dump(r->report));
    ASSERT_TRUE(std::holds_alternative<ShellId>(r->root)); // no solid: its shell is returned
    ASSERT_EQ(r->report.faces, 5u);
    ASSERT_EQ(r->report.unsupported.size(), 1u);
    ASSERT_EQ(r->report.unsupported[0].type, "SURFACE_REPLICA");
    ASSERT_EQ(r->report.failed.size(), 2u); // the face, then the solid
    ASSERT_EQ(r->report.failed[1].type, "MANIFOLD_SOLID_BREP");
    ASSERT_FALSE(is_closed(m, std::get<ShellId>(r->root)));

    // strict: the first failure fails the reading, the model is unchanged
    Model<double> ms;
    auto v = unwrap(make_vertex(ms, {0., 0., 0.}));
    auto rs = read_step_file(ms, file("box_unsupported_face.stp"), {.strict = true});
    ASSERT_FALSE(rs.has_value());
    ASSERT_EQ(rs.error().code, BuildErrc::UnsupportedEntity);
    ASSERT_EQ(ms.ids<VertexId>().size(), 1u);
    ASSERT_TRUE(ms.alive(v));
    ASSERT_TRUE(ms.ids<EdgeId>().empty());
    ASSERT_TRUE(ms.ids<FaceId>().empty());
    ASSERT_TRUE(ms.ids<ShellId>().empty());
}

TEST(tests_step_read, errors)
{
    Model<double> m;
    // syntax error located
    auto r = read_step(m, "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.))\nENDSEC;\nEND-ISO-10303-21;\n");
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::InvalidFile);
    ASSERT_NE(r.error().message.find("line"), std::string::npos);
    // no shape
    auto r2 = read_step(m, "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\nENDSEC;\nEND-ISO-10303-21;\n");
    ASSERT_FALSE(r2.has_value());
    ASSERT_EQ(r2.error().code, BuildErrc::InvalidFile);
    // missing file
    ASSERT_FALSE(read_step_file(m, file("no_such_file.stp")).has_value());
    // a model in another unit
    Model<double> mi;
    mi.setUnitScale(25.4);
    unwrap(make_vertex(mi, {0., 0., 0.}));
    auto r3 = read_step_file(mi, file("box.stp"), {.target_unit_mm = 1.});
    ASSERT_FALSE(r3.has_value());
    // reading twice into the same model: two independent solids
    Model<double> m2;
    auto a = read_ok(m2, "box.stp");
    auto b = read_ok(m2, "cylinder.stp");
    ASSERT_NE(shape_index(a.root), shape_index(b.root));
    ASSERT_EQ(m2.ids<SolidId>().size(), 2u);
}
