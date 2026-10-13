#include <doctest_gtest.hpp>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include "tests_brep_geom.h"

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
    using brep_tests::pi;

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

    struct Summary
    {
        std::size_t solids, faces, edges, vertices;
        double volume;
    };

    Summary summary(const Model<double> &m, const ShapeId &s)
    {
        Summary r{explore<SolidId>(m, s).size(), explore<FaceId>(m, s).size(), 0, explore<VertexId>(m, s).size(), 0.};
        for (auto e : explore<EdgeId>(m, s))
            r.edges += !m.edge(e).degenerate;
        for (auto so : explore<SolidId>(m, s))
            r.volume += signed_volume(m, m.solid(so).outer, 256);
        return r;
    }

    // reads back a text written by gbs: complete and valid
    StepReadResult read_back(Model<double> &m, const std::string &text, StepReadOptions o = {})
    {
        auto r = read_step(m, text, o);
        if (!r)
            FAIL(r.error().message);
        INFO(dump(r->report));
        REQUIRE(r->report.failed.empty());
        REQUIRE(r->report.unsupported.empty());
        REQUIRE(r->report.check.ok());
        return std::move(*r);
    }
}

TEST(tests_step_write, round_trip_of_the_test_files)
{
    // every reading test file: read, write, read again, same shape
    for (const char *name : {"box.stp", "cylinder.stp", "sphere.stp", "cone.stp", "plate_inch.stp", "bspline_patch.stp", "band_cylinder.stp",
                             "sphere_cap.stp", "cone_noseam.stp", "band_with_hole.stp", "occt_partial_torus.stp", "occt_loft.stp", "assembly.stp",
                             "mapped_items.stp"})
    {
        INFO(std::string(name));
        Model<double> m1;
        auto r1 = read_step_file(m1, file(name));
        REQUIRE(r1.has_value());
        auto w = write_step(m1, r1->root);
        if (!w)
            FAIL(w.error().message);
        ASSERT_TRUE(w->second.skipped.empty());
        ASSERT_EQ(w->second.approximated_surfaces, 0u); // primitives and B-splines: everything exact
        ASSERT_EQ(w->second.approximated_curves, 0u);
        Model<double> m2;
        auto r2 = read_back(m2, w->first);
        const auto a = summary(m1, r1->root), b = summary(m2, r2.root);
        ASSERT_EQ(b.solids, a.solids);
        ASSERT_EQ(b.faces, a.faces);
        ASSERT_EQ(b.edges, a.edges);
        ASSERT_EQ(b.vertices, a.vertices);
        ASSERT_NEAR(b.volume, a.volume, 1e-6 * (1. + std::abs(a.volume)));
        if (!m1.name(r1->root).empty()) // an unnamed root gets a generated product name
            ASSERT_EQ(m2.name(r2.root), m1.name(r1->root));
        ASSERT_EQ(m2.unitScale(), m1.unitScale());
    }
}

TEST(tests_step_write, analytic_entities_and_pcurves)
{
    Model<double> m;
    auto r = read_step_file(m, file("cylinder.stp"));
    REQUIRE(r.has_value());
    auto w = write_step(m, r->root);
    REQUIRE(w.has_value());
    const auto f = parse_p21(w->first);
    REQUIRE(f.has_value());
    ASSERT_EQ(f->instances_of("CYLINDRICAL_SURFACE").size(), 1u);
    ASSERT_EQ(f->instances_of("PLANE").size(), 2u);
    ASSERT_EQ(f->instances_of("CIRCLE").size(), 2u);
    ASSERT_EQ(f->instances_of("SEAM_CURVE").size(), 1u);     // the seam, used twice by the lateral face
    ASSERT_EQ(f->instances_of("SURFACE_CURVE").size(), 2u);  // the circles, shared by two faces
    ASSERT_EQ(f->instances_of("PCURVE").size(), 6u);
    ASSERT_EQ(f->instances_of("EDGE_CURVE").size(), 3u);
    ASSERT_EQ(f->instances_of("MANIFOLD_SOLID_BREP").size(), 1u);
    ASSERT_TRUE(f->instances_of("B_SPLINE_SURFACE_WITH_KNOTS").empty());
    ASSERT_TRUE(f->schemas().front().starts_with("AUTOMOTIVE_DESIGN"));
    // without pcurves
    auto np = write_step(m, r->root, {.pcurves = false});
    REQUIRE(np.has_value());
    ASSERT_TRUE(parse_p21(np->first)->instances_of("PCURVE").empty());
    Model<double> m2;
    read_back(m2, np->first);
    // AP242 header
    auto a242 = write_step(m, r->root, {.schema = "AP242"});
    REQUIRE(a242.has_value());
    ASSERT_TRUE(parse_p21(a242->first)->schemas().front().starts_with("AP242"));
    ASSERT_FALSE(write_step(m, r->root, {.schema = "AP203"}).has_value());
}

TEST(tests_step_write, assemblies_share_their_parts)
{
    Model<double> m;
    auto r = read_step_file(m, file("assembly.stp"));
    REQUIRE(r.has_value());
    auto w = write_step(m, r->root);
    REQUIRE(w.has_value());
    // ASM { BOX, BOX, SUB { CYL } }: BOX written once, 4 products, 4 occurrences
    ASSERT_EQ(w->second.parts, 2u);
    ASSERT_EQ(w->second.assemblies, 2u);
    ASSERT_EQ(w->second.occurrences, 4u);
    ASSERT_EQ(w->second.solids, 2u);
    const auto f = parse_p21(w->first);
    ASSERT_EQ(f->instances_of("PRODUCT").size(), 4u);
    ASSERT_EQ(f->instances_of("NEXT_ASSEMBLY_USAGE_OCCURRENCE").size(), 4u);
    ASSERT_EQ(f->instances_of("MANIFOLD_SOLID_BREP").size(), 2u);
    // read back: same hierarchy, names and positions
    Model<double> m2;
    auto r2 = read_back(m2, w->first);
    ASSERT_EQ(m2.name(r2.root), "ASM");
    const auto &c1 = m.compound(std::get<CompoundId>(r->root)).shapes;
    const auto &c2 = m2.compound(std::get<CompoundId>(r2.root)).shapes;
    ASSERT_EQ(c2.size(), 3u);
    for (std::size_t k = 0; k < 3; ++k)
    {
        ASSERT_EQ(m2.name(c2[k]), m.name(c1[k]));
        const auto b1 = bounding_box(m, c1[k], 37), b2 = bounding_box(m2, c2[k], 37);
        ASSERT_LT(norm(b1.min - b2.min), 1e-6);
        ASSERT_LT(norm(b1.max - b2.max), 1e-6);
    }
    ASSERT_EQ(m2.instance(c2[0])->part, m2.instance(c2[1])->part); // shared again
    // flattened: one part, three solids
    auto fl = write_step(m, r->root, {.flatten = true});
    REQUIRE(fl.has_value());
    ASSERT_EQ(fl->second.parts, 1u);
    ASSERT_EQ(fl->second.assemblies, 0u);
    Model<double> m3;
    auto r3 = read_back(m3, fl->first);
    ASSERT_EQ(r3.report.solids, 3u);
}

TEST(tests_step_write, shapes_built_in_gbs)
{
    using namespace brep_tests;
    Model<double> m;
    // a box sewn from bilinear B-spline faces, a sphere and a torus of natural faces on primitives
    auto sewn = unwrap(sew(m, box_faces(m)));
    const auto box = unwrap(make_solid(m, sewn.shells.front()));
    m.setName(box, "BOX");
    const ax2<double, 3> z{{{0., 0., 5.}, {0., 0., 1.}, {1., 0., 0.}}};
    auto sphere_face = unwrap(make_face(m, std::make_shared<SphericalSurface<double>>(1., z, std::array<double, 4>{0., 2 * pi, -pi / 2, pi / 2})));
    const auto ball = unwrap(make_solid(m, m.add(Shell{{FaceUse{sphere_face}}})));
    auto torus_face = unwrap(make_face(m, std::make_shared<ToroidalSurface<double>>(3., 1., z, std::array<double, 4>{0., 2 * pi, 0., 2 * pi})));
    const auto torus = unwrap(make_solid(m, m.add(Shell{{FaceUse{torus_face}}})));
    // an analytic surface given by a formula: approximated
    auto free_face = unwrap(make_face(m, brep_tests::sphere(2.)));
    const std::vector<ShapeId> shapes{box, ball, torus, free_face, unwrap(make_vertex(m, {0., 0., 0.}))};
    auto w = write_step(m, std::span<const ShapeId>{shapes}, {.approx_tol = 1e-6});
    REQUIRE(w.has_value());
    ASSERT_EQ(w->second.skipped.size(), 1u); // the vertex
    ASSERT_EQ(w->second.approximated_surfaces, 1u);
    ASSERT_EQ(w->second.parts, 4u);
    const auto f = parse_p21(w->first);
    ASSERT_EQ(f->instances_of("SPHERICAL_SURFACE").size(), 1u);
    ASSERT_EQ(f->instances_of("TOROIDAL_SURFACE").size(), 1u);
    ASSERT_EQ(f->instances_of("SEAM_CURVE").size(), 4u); // sphere, torus (two), the sphere given by a formula
    Model<double> m2;
    auto r2 = read_back(m2, w->first);
    ASSERT_EQ(r2.report.solids, 3u);
    double v = 0.;
    for (auto so : explore<SolidId>(m2, r2.root))
        v += signed_volume(m2, m2.solid(so).outer, 256);
    ASSERT_NEAR(v, 1. + 4. / 3. * pi + 2 * pi * pi * 3., 2e-2);
}

TEST(tests_step_write, units_files_and_errors)
{
    // a model in inches stays in inches
    Model<double> mi;
    auto ri = read_step_file(mi, file("plate_inch.stp"), {.target_unit_mm = 25.4});
    REQUIRE(ri.has_value());
    auto w = write_step(mi, ri->root);
    REQUIRE(w.has_value());
    Model<double> mm;
    auto rm = read_back(mm, w->first); // read into millimetres
    ASSERT_EQ(rm.report.length_unit, "INCH");
    ASSERT_NEAR(rm.report.length_factor, 25.4, 1e-12);
    // file output
    const auto path = (std::filesystem::temp_directory_path() / "gbs_step_write.stp").string();
    auto wf = write_step_file(mi, ri->root, path);
    REQUIRE(wf.has_value());
    Model<double> m2;
    ASSERT_TRUE(read_step_file(m2, path).has_value());
    // errors: dead shape, unwritable path
    ASSERT_FALSE(write_step(mi, ShapeId{SolidId{99}}).has_value());
    ASSERT_FALSE(write_step_file(mi, ri->root, "/no/such/dir/x.stp").has_value());
}
