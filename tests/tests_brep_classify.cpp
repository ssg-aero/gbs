#include <doctest_gtest.hpp>
#include <gbs-brep/classify.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <string>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using brep_tests::pi;
    using T = double;

    std::string file(const std::string &name) { return std::string(GBS_TESTS_IN_DIR) + "/step/" + name; }

    SolidId read_solid(Model<T> &m, const std::string &name)
    {
        auto r = step::read_step_file(m, file(name));
        REQUIRE(r.has_value());
        const auto solids = explore<SolidId>(m, r->root);
        REQUIRE(solids.size() == 1u);
        return solids.front();
    }

    // the face of a solid whose surface contains p
    FaceId face_at(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        for (auto f : explore<FaceId>(m, ShapeId{s}))
        {
            const auto &srf = *m.face(f).surface;
            const auto b = srf.bounds();
            const auto [u, v, d] = intersect_detail::project_on_surface(srf, b, p);
            if (d < 1e-9)
                return f;
        }
        FAIL("no face at the point");
        return {};
    }

    // a hash of the parameters in [-1, 1): a roughness at every scale
    T noise(T u, T v = 0.)
    {
        std::uint64_t h = std::bit_cast<std::uint64_t>(u) * 0x9E3779B97F4A7C15ull ^ std::bit_cast<std::uint64_t>(v) * 0xC2B2AE3D27D4EB4Full;
        h ^= h >> 31;
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 29;
        return T(h >> 11) * 0x1p-52 - T(1);
    }

    // the plane z = h on [0, a]² (parameters in [0, 1]²), rough by 0.005: never flat, its
    // subdivision against a ray reaches IntersectOptions::max_pairs
    class RoughPlane : public Surface<T, 3>
    {
        T a_, h_;

    public:
        RoughPlane(T a, T h) : a_{a}, h_{h} {}
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            if (du + dv == 0)
                return {a_ * u, a_ * v, h_ + T(0.005) * noise(u, v)};
            if (du + dv == 1)
                return {du == 1 ? a_ : T(0), dv == 1 ? a_ : T(0), T(0)};
            return {};
        }
        auto bounds() const -> std::array<T, 4> override { return {0., 1., 0., 1.}; }
        auto has_second_derivatives() const -> bool override { return false; }
    };

    // a pcurve made rough by 0.5 (in the parameter space of the face): never flat either
    class RoughCurve : public Curve<T, 2>
    {
        std::shared_ptr<Curve<T, 2>> c_;

    public:
        explicit RoughCurve(std::shared_ptr<Curve<T, 2>> c) : c_{std::move(c)} {}
        auto value(T t, size_t d = 0) const -> std::array<T, 2> override
        {
            auto p = c_->value(t, d);
            if (d == 0)
                p[0] += T(0.5) * noise(t);
            return p;
        }
        auto bounds() const -> std::array<T, 2> override { return c_->bounds(); }
        auto has_second_derivatives() const -> bool override { return false; }
    };
}

TEST(tests_brep_classify, points_in_faces)
{
    // the plate (4 x 2 in, hole of radius 0.5 in at (1, 1)), read in millimetres
    Model<T> m;
    auto r = step::read_step_file(m, file("plate_inch.stp"));
    REQUIRE(r.has_value());
    const auto faces = explore<FaceId>(m, r->root);
    REQUIRE(faces.size() == 1u);
    const auto f = faces.front();
    const auto &plane = dynamic_cast<const Plane<T> &>(*m.face(f).surface);
    auto uv = [&](T x, T y) {
        const auto q = plane.parameters(point<T, 3>{x * 25.4, y * 25.4, 0.});
        return point<T, 2>{q[0], q[1]};
    };
    ASSERT_EQ(classify(m, f, uv(3., 1.)), Location::In);
    ASSERT_EQ(classify(m, f, uv(1., 1.)), Location::Out);   // in the hole
    ASSERT_EQ(classify(m, f, uv(1.25, 1.)), Location::Out); // in the hole, near its edge
    ASSERT_EQ(classify(m, f, uv(1.6, 1.)), Location::In);
    ASSERT_EQ(classify(m, f, uv(1.5, 1.)), Location::On);   // on the hole
    ASSERT_EQ(classify(m, f, uv(2., 0.)), Location::On);    // on the outer boundary
    ASSERT_EQ(classify(m, f, uv(4., 2.)), Location::On);    // at a vertex
    ASSERT_EQ(classify(m, f, uv(5., 1.)), Location::Out);
    ASSERT_EQ(classify(m, f, uv(2., 1e-3)), Location::In);  // 25 µm inside, tolerance 1e-6
    // ray directions: points on the horizontal and vertical lines of the vertices and of the hole
    for (T y : {0.5, 1., 1.5})
        ASSERT_EQ(classify(m, f, uv(0.25, y)), Location::In);
    const auto c = interior_point(m, f);
    REQUIRE(c.has_value());
    ASSERT_EQ(classify(m, f, *c), Location::In);
}

TEST(tests_brep_classify, faces_with_seams_and_poles)
{
    // lateral face of the cylinder (seam used twice), sphere bounded by its seam (poles)
    Model<T> m;
    const auto cyl = read_solid(m, "cylinder.stp");
    const auto lat = face_at(m, cyl, {5., 0., 5.});
    const auto box = m.face(lat).surface->bounds();
    for (T u : {0.1, 1., 3., 6.2})
        ASSERT_EQ(classify(m, lat, point<T, 2>{u, 5.}), Location::In);
    ASSERT_EQ(classify(m, lat, point<T, 2>{0., 5.}), Location::On); // on the seam
    ASSERT_EQ(classify(m, lat, point<T, 2>{1., 10.}), Location::On);
    ASSERT_EQ(classify(m, lat, point<T, 2>{1., box[3] + 1.}), Location::Out);
    const auto sph = read_solid(m, "sphere.stp");
    const auto sf = explore<FaceId>(m, ShapeId{sph}).front();
    for (T v : {-1.5, -0.5, 0., 1.2, 1.5})
        ASSERT_EQ(classify(m, sf, point<T, 2>{2., v}), Location::In);
    ASSERT_EQ(classify(m, sf, point<T, 2>{2., pi / 2}), Location::On); // the north pole
    // interior points of every face of the test solids
    for (const char *name : {"box.stp", "cylinder.stp", "sphere.stp", "cone.stp", "band_with_hole.stp", "occt_partial_torus.stp", "occt_loft.stp"})
    {
        INFO(std::string(name));
        Model<T> mm;
        auto rr = step::read_step_file(mm, file(name));
        REQUIRE(rr.has_value());
        for (auto f : explore<FaceId>(mm, rr->root))
        {
            const auto p = interior_point(mm, f);
            REQUIRE(p.has_value());
            ASSERT_EQ(classify(mm, f, *p), Location::In);
        }
    }
}

TEST(tests_brep_classify, points_in_solids)
{
    Model<T> m;
    // box 2 x 3 x 4
    const auto box = read_solid(m, "box.stp");
    ASSERT_EQ(classify(m, box, {1., 1.5, 2.}), Location::In);
    ASSERT_EQ(classify(m, box, {1e-3, 1e-3, 1e-3}), Location::In);
    ASSERT_EQ(classify(m, box, {-1e-3, 1.5, 2.}), Location::Out);
    ASSERT_EQ(classify(m, box, {5., 5., 5.}), Location::Out);
    ASSERT_EQ(classify(m, box, {2.5, 1.5, 2.}), Location::Out);
    ASSERT_EQ(classify(m, box, {0., 1.5, 2.}), Location::On);   // on a face
    ASSERT_EQ(classify(m, box, {0., 0., 2.}), Location::On);    // on an edge
    ASSERT_EQ(classify(m, box, {2., 3., 4.}), Location::On);    // at a vertex
    // points aligned with edges and faces in every axis direction
    for (T x : {0.5, 1., 1.5})
        for (T z : {1., 2., 3.})
            ASSERT_EQ(classify(m, box, {x, 1.5, z}), Location::In);
    // cylinder R 5 h 10 (seam), sphere R 3 (poles), cone
    const auto cyl = read_solid(m, "cylinder.stp");
    ASSERT_EQ(classify(m, cyl, {0., 0., 5.}), Location::In);
    ASSERT_EQ(classify(m, cyl, {4.9, 0., 0.1}), Location::In);
    ASSERT_EQ(classify(m, cyl, {4., 4., 5.}), Location::Out); // inside the box, outside the cylinder
    ASSERT_EQ(classify(m, cyl, {5., 0., 5.}), Location::On);  // on the seam
    ASSERT_EQ(classify(m, cyl, {3., 4., 10.}), Location::On); // on the top circle
    const auto sph = read_solid(m, "sphere.stp");
    ASSERT_EQ(classify(m, sph, {0., 0., 0.}), Location::In);
    ASSERT_EQ(classify(m, sph, {0., 0., 2.99}), Location::In);
    ASSERT_EQ(classify(m, sph, {2.5, 2.5, 0.}), Location::Out);
    ASSERT_EQ(classify(m, sph, {0., 0., 3.}), Location::On); // pole
    const auto cone = read_solid(m, "cone.stp");
    const T H = 2. / std::tan(pi / 6);
    ASSERT_EQ(classify(m, cone, {0., 0., -H / 2}), Location::In);
    ASSERT_EQ(classify(m, cone, {1.5, 0., -H / 2}), Location::Out);
    ASSERT_EQ(classify(m, cone, {0., 0., -H}), Location::On); // apex
}

TEST(tests_brep_classify, cavities_and_faces)
{
    using namespace brep_tests;
    // a box 0..4 with a cavity 1..2: the inner shell is reversed by make_solid
    Model<T> m;
    auto outer = unwrap(sew(m, box_faces(m, {0., 0., 0.}, 4.)));
    auto inner = unwrap(sew(m, box_faces(m, {1., 1., 1.}, 1.)));
    const auto solid = unwrap(make_solid(m, outer.shells.front(), std::vector<ShellId>{inner.shells.front()}));
    ASSERT_EQ(classify(m, solid, {3., 3., 3.}), Location::In);
    ASSERT_EQ(classify(m, solid, {1.5, 1.5, 1.5}), Location::Out); // in the cavity
    ASSERT_EQ(classify(m, solid, {1., 1.5, 1.5}), Location::On);   // on the cavity
    ASSERT_EQ(classify(m, solid, {0.5, 1.5, 1.5}), Location::In);

    // the faces of the solid against a box 0.5..2.5 around the cavity: the six faces of the
    // cavity are inside, the six outer faces outside
    auto other = unwrap(sew(m, box_faces(m, {0.5, 0.5, 0.5}, 2.)));
    const auto probe = unwrap(make_solid(m, other.shells.front()));
    std::size_t in{}, out{};
    for (auto f : explore<FaceId>(m, ShapeId{solid}))
    {
        const auto where = classify(m, f, probe);
        REQUIRE(where.has_value());
        (*where == Location::In ? in : out) += 1;
    }
    ASSERT_EQ(in, 6u);
    ASSERT_EQ(out, 6u);
}

TEST(tests_brep_classify, capped_intersections)
{
    // an intersection capped by max_pairs may miss crossings: classify used its partial count of
    // the first ray (here none on the rough face, the parity wrong) instead of another direction
    using namespace brep_tests;
    Model<T> m;
    auto shell = unwrap(sew(m, box_faces(m, {0., 0., 0.}, 4.)));
    const auto solid = unwrap(make_solid(m, shell.shells.front()));
    for (auto f : explore<FaceId>(m, ShapeId{solid}))
        if (std::abs(m.face(f).surface->value(0.5, 0.5)[2] - 4.) < 1e-12)
            m.face(f).surface = std::make_shared<RoughPlane>(4., 4.); // the top face, crossed by the first rays
    ASSERT_EQ(classify(m, solid, {2., 2., 0.5}), Location::In);
    ASSERT_EQ(classify(m, solid, {2., 2., 5.}), Location::Out);
    const LineCurve<T> up{{2., 2., 0.5}, classify_detail::directions_3d<T>(0), 0., 10.};
    ASSERT_FALSE(intersect(up, 0., 10., RoughPlane{4., 4.}, {0., 1., 0., 1.}).complete);

    // the right side of the plate (x = 4 in) made rough: the first ray from the right half crosses it
    Model<T> mp;
    auto r = step::read_step_file(mp, file("plate_inch.stp"));
    REQUIRE(r.has_value());
    const auto f = explore<FaceId>(mp, r->root).front();
    const auto &plane = dynamic_cast<const Plane<T> &>(*mp.face(f).surface);
    auto uv = [&](T x, T y) {
        const auto q = plane.parameters(point<T, 3>{x * 25.4, y * 25.4, 0.});
        return point<T, 2>{q[0], q[1]};
    };
    const auto right = uv(4., 1.);
    for (auto w : mp.face(f).wires)
        for (auto &ce : mp.wire(w).coedges)
        {
            const auto &e = mp.edge(ce.edge);
            if (const auto mid = ce.pcurve->value((e.u1 + e.u2) / 2.); norm(mid - right) < 1e-6)
                ce.pcurve = std::make_shared<RoughCurve>(ce.pcurve);
        }
    ASSERT_EQ(classify(mp, f, uv(3., 1.)), Location::In);
    ASSERT_EQ(classify(mp, f, uv(1.6, 1.)), Location::In);
}

TEST(tests_brep_classify, near_the_extreme_of_a_curved_boundary)
{
    // a disk of radius 1 whose circle starts at 10 degrees: 17 samples of the circle missed its
    // points at 0 degree (box x <= cos 10° = 0.985), so points of the face beyond were outside it
    Model<T> m;
    points_vector<T, 3> poles{{-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}};
    const std::vector<T> k{-2., -2., 2., 2.};
    const auto plane = std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1); // parameters x and y
    const T a = 10. * pi / 180.;
    const auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(1., a, a + 2. * pi, ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}));
    const auto disk = unwrap(make_face(m, plane, unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle))}))));
    ASSERT_EQ(classify(m, disk, {0.99, 0.}), Location::In);
    ASSERT_EQ(classify(m, disk, {0.999, 0.}), Location::In);
    ASSERT_EQ(classify(m, disk, {1.001, 0.}), Location::Out);
    ASSERT_EQ(classify(m, disk, {0., -0.995}), Location::In);
    const auto box = classify_detail::face_uv_box(m, disk); // enclosing the circle, by 1e-3 of it at most
    ASSERT_GE(box[1], 1.);
    ASSERT_LE(box[1], 1.003);

    // its prism: the first ray from a point near the rim and the top crosses the top there
    const auto cyl = explore<SolidId>(m, unwrap(make_prism(m, ShapeId{disk}, point<T, 3>{0., 0., 1.})).shape).front();
    ASSERT_EQ(classify(m, cyl, {0.988, 0., 0.99}), Location::In);
    ASSERT_EQ(classify(m, cyl, {0.992, 0.001, 0.98}), Location::In);
    ASSERT_EQ(classify(m, cyl, {0.988, 0., 1.01}), Location::Out);
}
