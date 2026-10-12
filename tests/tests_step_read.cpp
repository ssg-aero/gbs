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

TEST(tests_step_read, seam_insertion)
{
    // faces turning around their surface without seam edge: a seam is inserted
    struct Case
    {
        const char *name;
        double volume, tol;
        std::size_t edges;
    };
    const double H = 2. / std::tan(pi / 6), hc = 1.5;
    for (const auto &c : {Case{"band_cylinder.stp", pi * 25. * 10., 0.2, 4},   // two circles (one split at the seam angle), seam
                          Case{"sphere_cap.stp", pi * hc * hc * (9. - hc) / 3., 0.05, 3}, // circle, seam, pole
                          Case{"cone_noseam.stp", pi * 4. * H / 3., 0.05, 3}})
    {
        INFO(c.name);
        Model<double> m;
        auto r = read_ok(m, c.name);
        ASSERT_TRUE(std::holds_alternative<SolidId>(r.root));
        ASSERT_EQ(r.report.edges, c.edges);
        ASSERT_NEAR(volume(m, r.root), c.volume, c.tol);
    }
    // the seam avoids a hole across the angle 0: both circles are split at the seam angle
    Model<double> m;
    auto r = read_ok(m, "band_with_hole.stp");
    ASSERT_EQ(r.report.faces, 1u);
    ASSERT_EQ(r.report.edges, 2u * 2u + 1u + 4u);
    const auto &f = m.face(explore<FaceId>(m, r.root)[0]);
    ASSERT_EQ(f.wires.size(), 2u); // the band (with its seam) and the hole
}

TEST(tests_step_read, files_written_by_occt)
{
    // partial torus (revolution of a circle by 3 rad) written by OpenCascade: torus axis -z, so the
    // seam edge lies at v = 3 pi / 2, and its loop runs against same_sense
    {
        Model<double> m;
        auto r = read_ok(m, "occt_partial_torus.stp");
        ASSERT_EQ(r.report.solids, 1u);
        ASSERT_EQ(r.report.faces, 3u);
        ASSERT_EQ(r.report.sense_mismatches, 1u);
        ASSERT_NEAR(volume(m, r.root), 3. / (2 * pi) * 2 * pi * pi * 4., 1e-2);
    }
    // loft between a circle and a square: B-spline faces whose edges lie on their borders
    {
        Model<double> m;
        auto r = read_ok(m, "occt_loft.stp");
        ASSERT_EQ(r.report.solids, 1u);
        ASSERT_EQ(r.report.faces, 7u);
    }
}

TEST(tests_step_read, assemblies)
{
    auto box_of = [](const Model<double> &m, const ShapeId &s) { return bounding_box(m, s, 37); };
    // the bounding boxes include the edge tolerances
    auto near = [](const point<double, 3> &a, const point<double, 3> &b) { return norm(a - b) < 1e-4; };
    // products and placed occurrences, flattened: ASM { BOX, BOX, SUB { CYL } }
    {
        Model<double> m;
        auto r = read_ok(m, "assembly.stp");
        ASSERT_TRUE(std::holds_alternative<CompoundId>(r.root));
        ASSERT_EQ(m.name(r.root), "ASM");
        ASSERT_EQ(r.report.solids, 3u);
        const auto &c = m.compound(std::get<CompoundId>(r.root));
        ASSERT_EQ(c.shapes.size(), 3u);
        // first box moved to (10, 0, 0)
        auto b1 = box_of(m, c.shapes[0]);
        ASSERT_TRUE(near(b1.min, {10., 0., 0.}));
        ASSERT_TRUE(near(b1.max, {12., 3., 4.}));
        ASSERT_EQ(m.name(c.shapes[0]), "BOX");
        // second box turned a quarter about z at (0, 20, 0): (a, b, c) -> (-b, 20 + a, c)
        auto b2 = box_of(m, c.shapes[1]);
        ASSERT_TRUE(near(b2.min, {-3., 20., 0.}));
        ASSERT_TRUE(near(b2.max, {0., 22., 4.}));
        // sub-assembly at (0, 0, 10) holding the cylinder at (5, 5, 0)
        ASSERT_EQ(m.name(c.shapes[2]), "SUB");
        auto b3 = box_of(m, c.shapes[2]);
        ASSERT_TRUE(near(b3.min, {4., 4., 10.}));
        ASSERT_TRUE(near(b3.max, {6., 6., 15.}));
        double v = 0.;
        for (auto so : explore<SolidId>(m, r.root))
            v += signed_volume(m, m.solid(so).outer, 256);
        ASSERT_NEAR(v, 2 * 24. + pi * 5., 1e-2);
        // the part read once is not left in the model besides its instances
        ASSERT_EQ(m.ids<SolidId>().size(), 3u);
    }
    // mapped items
    {
        Model<double> m;
        auto r = read_ok(m, "mapped_items.stp");
        ASSERT_EQ(r.report.solids, 2u);
        const auto &c = m.compound(std::get<CompoundId>(r.root));
        auto b1 = box_of(m, c.shapes[0]);
        ASSERT_TRUE(near(b1.min, {-5., 0., 0.}));
        ASSERT_TRUE(near(b1.max, {-3., 3., 4.}));
        auto b2 = box_of(m, c.shapes[1]); // (a, b, c) -> (5 + c, a, b)
        ASSERT_TRUE(near(b2.min, {5., 0., 0.}));
        ASSERT_TRUE(near(b2.max, {9., 2., 3.}));
        ASSERT_EQ(m.ids<SolidId>().size(), 2u);
    }
}

TEST(tests_step_read, elementary_geometry_stays_analytic)
{
    Model<double> m;
    auto r = read_ok(m, "cylinder.stp");
    std::size_t cylinders = 0, planes = 0, circles = 0, lines = 0;
    for (auto f : explore<FaceId>(m, r.root))
    {
        cylinders += std::dynamic_pointer_cast<CylindricalSurface<double>>(m.face(f).surface) != nullptr;
        planes += std::dynamic_pointer_cast<Plane<double>>(m.face(f).surface) != nullptr;
    }
    for (auto e : explore<EdgeId>(m, r.root))
    {
        circles += std::dynamic_pointer_cast<EllipseCurve<double>>(m.edge(e).curve) != nullptr;
        lines += std::dynamic_pointer_cast<LineCurve<double>>(m.edge(e).curve) != nullptr;
    }
    ASSERT_EQ(cylinders, 1u);
    ASSERT_EQ(planes, 2u);
    ASSERT_EQ(circles, 2u);
    ASSERT_EQ(lines, 1u);
    // the edge parameters are those of the file: a full circle runs over 2 pi
    for (auto e : explore<EdgeId>(m, r.root))
        if (std::dynamic_pointer_cast<EllipseCurve<double>>(m.edge(e).curve))
            ASSERT_NEAR(m.edge(e).u2 - m.edge(e).u1, 2 * pi, 1e-12);
    // B-spline geometry stays B-spline
    Model<double> mb;
    auto rb = read_ok(mb, "bspline_patch.stp");
    const auto patch = mb.face(explore<FaceId>(mb, rb.root)[0]).surface;
    using Rational = BSSurfaceRational<double, 3>;
    ASSERT_TRUE(std::dynamic_pointer_cast<Rational>(patch) != nullptr);
}

TEST(tests_step_read, assembly_provenance)
{
    Model<double> m;
    auto r = read_ok(m, "assembly.stp");
    const auto &c = m.compound(std::get<CompoundId>(r.root));
    const auto i1 = m.instance(c.shapes[0]), i2 = m.instance(c.shapes[1]), isub = m.instance(c.shapes[2]);
    ASSERT_TRUE(i1 && i2 && isub);
    ASSERT_EQ(i1->part, i2->part); // two instances of BOX
    ASSERT_NE(isub->part, i1->part);
    ASSERT_TRUE(m.instance(r.root)->place.is_identity()); // the root assembly
    // placements: BOX at (10, 0, 0); BOX turned a quarter about z at (0, 20, 0)
    ASSERT_LT(norm(i1->place.apply({0., 0., 0.}) - point<double, 3>{10., 0., 0.}), 1e-12);
    ASSERT_LT(norm(i2->place.apply({1., 0., 0.}) - point<double, 3>{0., 21., 0.}), 1e-12);
    // both instances are the same part: their vertices brought back to the part coincide
    auto part_vertices = [&](const ShapeId &s, const InstanceOf<double> &i) {
        std::vector<point<double, 3>> v;
        for (auto x : explore<VertexId>(m, s))
            v.push_back(i.place.inverse().apply(m.vertex(x).pnt));
        std::ranges::sort(v);
        return v;
    };
    const auto v1 = part_vertices(c.shapes[0], *i1), v2 = part_vertices(c.shapes[1], *i2);
    ASSERT_EQ(v1.size(), 8u);
    for (std::size_t k = 0; k < v1.size(); ++k)
        ASSERT_LT(norm(v1[k] - v2[k]), 1e-9);
    // the cylinder inside SUB: an instance of CYL, placed by both occurrences
    const auto cyl = m.compound(std::get<CompoundId>(c.shapes[2])).shapes.front();
    ASSERT_LT(norm(m.instance(cyl)->place.apply({0., 0., 0.}) - point<double, 3>{5., 5., 10.}), 1e-12);
    // a second reading into the same model gets its own part keys
    auto r2 = read_ok(m, "assembly.stp");
    const auto first = m.compound(std::get<CompoundId>(r2.root)).shapes[0];
    ASSERT_NE(m.instance(first)->part, i1->part);
}
