#include <doctest_gtest.hpp>
#include <gbs-brep/sweep.h>
#include <gbs-brep/offset.h>
#include <gbs-brep/check.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

#include <array>
#include <cmath>
#include <typeinfo>
#include <numbers>
#include <sstream>
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

    ShapeId read(Model<T> &m, const std::string &name)
    {
        auto r = step::read_step_file(m, file(name));
        REQUIRE(r.has_value());
        return r->root;
    }

    void expect_valid(const Model<T> &m, const ShapeId &s)
    {
        const auto rep = check(m, s);
        std::ostringstream os;
        for (const auto &e : rep.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        INFO(os.str());
        REQUIRE(rep.ok());
    }

    // volume by the sampling of the faces (signed_volume): a few 1e-4 relative on trimmed faces
    T volume(const Model<T> &m, const ShapeId &s)
    {
        T v{};
        for (auto so : explore<SolidId>(m, s))
            v += signed_volume(m, m.solid(so).outer, 160);
        return v;
    }

    template <typename S>
    std::size_t count(const Model<T> &m, const ShapeId &s)
    {
        std::size_t n{};
        for (auto f : explore<FaceId>(m, s))
            n += dynamic_cast<const S *>(m.face(f).surface.get()) != nullptr;
        return n;
    }

    FaceId face_on(const Model<T> &m, const ShapeId &s, const point<T, 3> &p)
    {
        for (auto f : explore<FaceId>(m, s))
        {
            const auto &srf = *m.face(f).surface;
            if (intersect_detail::project_on_surface(srf, srf.bounds(), p)[2] < 1e-9)
                return f;
        }
        FAIL("no face");
        return {};
    }
}

TEST(tests_brep_sweep, prisms)
{
    // the plate 4 x 2 in with a hole, extruded by 10 mm: 4 planes, a cylinder, 2 caps
    Model<T> m;
    const auto plate = read(m, "plate_inch.stp");
    const auto f = explore<FaceId>(m, plate).front();
    auto r = unwrap(make_prism(m, f, point<T, 3>{0., 0., 10.}));
    ASSERT_EQ(r.approximated, 0u);
    expect_valid(m, r.shape);
    ASSERT_EQ(explore<FaceId>(m, r.shape).size(), 7u);
    ASSERT_EQ(count<Plane<T>>(m, r.shape), 6u);
    ASSERT_EQ(count<CylindricalSurface<T>>(m, r.shape), 1u);
    const T area = (4. * 2. - pi * 0.25) * 25.4 * 25.4;
    ASSERT_NEAR(volume(m, r.shape), area * 10., 2e-3 * area * 10.);
    // obliquely: the sides are exact NURBS extrusions, same volume (base x height)
    auto o = unwrap(make_prism(m, f, point<T, 3>{3., 2., 10.}));
    expect_valid(m, o.shape);
    ASSERT_NEAR(volume(m, o.shape), area * 10., 2e-3 * area * 10.);

    // the cap of the cylinder (a disc R 5) extruded: a cylinder of volume 25 pi h
    const auto cyl = read(m, "cylinder.stp");
    const auto disc = face_on(m, cyl, {0., 0., 0.});
    auto c = unwrap(make_prism(m, disc, point<T, 3>{0., 0., -4.}));
    expect_valid(m, c.shape);
    ASSERT_NEAR(volume(m, c.shape), 100. * pi, 2e-3 * 100. * pi);
    // an edge gives a face, a wire a shell
    const auto circle = explore<EdgeId>(m, ShapeId{disc}).front();
    auto e = unwrap(make_prism(m, circle, point<T, 3>{0., 0., 1.}));
    ASSERT_TRUE(std::holds_alternative<FaceId>(e.shape));
    ASSERT_TRUE(dynamic_cast<const CylindricalSurface<T> *>(m.face(std::get<FaceId>(e.shape)).surface.get()) != nullptr);
    const auto wire = m.face(f).wires.front();
    auto w = unwrap(make_prism(m, wire, point<T, 3>{0., 0., 1.}));
    ASSERT_TRUE(std::holds_alternative<ShellId>(w.shape));
    ASSERT_EQ(m.shell(std::get<ShellId>(w.shape)).faces.size(), 4u);
    // errors
    ASSERT_FALSE(make_prism(m, f, point<T, 3>{0., 0., 0.}).has_value());
    ASSERT_FALSE(make_prism(m, ShapeId{SolidId{999}}, point<T, 3>{0., 0., 1.}).has_value());
}

TEST(tests_brep_sweep, revolutions)
{
    using namespace brep_tests;
    Model<T> m;
    const ax1<T, 3> z_axis{{{0., 0., 0.}, {0., 0., 1.}}};
    // a rectangle x in [2, 3], z in [0, 1] of the plane y = 0 around z: a ring of volume 5 pi
    const auto rect = unwrap(make_face(m, quad({2., 0., 0.}, {3., 0., 0.}, {2., 0., 1.}, {3., 0., 1.})));
    auto ring = unwrap(make_revolution(m, rect, z_axis));
    expect_valid(m, ring.shape);
    ASSERT_EQ(count<CylindricalSurface<T>>(m, ring.shape), 2u);
    ASSERT_NEAR(volume(m, ring.shape), 5. * pi, 2e-3 * 5. * pi);
    // half a turn: closed by the profile and its turned copy
    auto half = unwrap(make_revolution(m, rect, z_axis, pi));
    expect_valid(m, half.shape);
    ASSERT_NEAR(volume(m, half.shape), 2.5 * pi, 2e-3 * 2.5 * pi);

    // the disc R 5 of the cylinder cap around an axis of its plane at 10 from its centre: a torus
    const auto cyl = read(m, "cylinder.stp");
    const auto disc = face_on(m, cyl, {0., 0., 0.});
    auto tor = unwrap(make_revolution(m, disc, ax1<T, 3>{{{10., 0., 0.}, {0., 1., 0.}}}));
    expect_valid(m, tor.shape);
    ASSERT_EQ(count<ToroidalSurface<T>>(m, tor.shape), 1u);
    ASSERT_NEAR(volume(m, tor.shape), 2. * pi * pi * 10. * 25., 2e-3 * 2. * pi * pi * 250.);

    // a right triangle with a side on the axis: a cone of volume pi r² h / 3 (apex on the axis)
    const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{2., 0., 0.}));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{2., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto e3 = unwrap(make_edge(m, point<T, 3>{0., 0., 3.}, point<T, 3>{0., 0., 0.}));
    const auto tw = unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3}));
    const auto xz = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, {1., 0., 0.}}}, std::array<T, 4>{-5., 5., -5., 5.});
    const auto tri = unwrap(make_face(m, xz, tw));
    auto cone = unwrap(make_revolution(m, tri, z_axis));
    expect_valid(m, cone.shape);
    ASSERT_EQ(count<ConicalSurface<T>>(m, cone.shape), 1u);
    ASSERT_NEAR(volume(m, cone.shape), 4. * pi, 2e-3 * 4. * pi);
    // a half disc R 2 of the plane y = 0 around its diameter on the axis: a sphere (two poles)
    const auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, {1., 0., 0.}}}, -pi / 2, pi / 2);
    const auto ea = unwrap(make_edge(m, arc));
    const auto ed = unwrap(make_edge(m, point<T, 3>{0., 0., 2.}, point<T, 3>{0., 0., -2.}));
    const auto hd = unwrap(make_face(m, xz, unwrap(make_wire(m, std::vector<EdgeId>{ea, ed}))));
    auto ball = unwrap(make_revolution(m, hd, z_axis));
    expect_valid(m, ball.shape);
    ASSERT_EQ(count<SphericalSurface<T>>(m, ball.shape), 1u);
    ASSERT_NEAR(volume(m, ball.shape), 32. / 3. * pi, 2e-3 * 32. / 3. * pi);
    // errors
    ASSERT_FALSE(make_revolution(m, rect, z_axis, 0.).has_value());
    ASSERT_FALSE(make_revolution(m, rect, ax1<T, 3>{{{0., 0., 0.}, {0., 0., 0.}}}).has_value());
}

TEST(tests_brep_sweep, offset_faces)
{
    Model<T> m;
    const auto cyl = read(m, "cylinder.stp"); // R 5, h 10
    for (auto f : explore<FaceId>(m, cyl))
    {
        INFO("cylinder face " << typeid(*m.face(f).surface).name());
        auto g = unwrap(offset_face(m, f, 1.));
        expect_valid(m, ShapeId{g});
        // every point of the offset face is at 1 from the face, along its normal
        const auto &s = *m.face(f).surface, &so = *m.face(g).surface;
        const auto uv = interior_point(m, f);
        REQUIRE(uv.has_value());
        const auto n = cross(s.value((*uv)[0], (*uv)[1], 1, 0), s.value((*uv)[0], (*uv)[1], 0, 1));
        const auto q = s.value((*uv)[0], (*uv)[1]) + (1. / norm(n)) * n;
        ASSERT_LT(norm(so.value((*uv)[0], (*uv)[1]) - q), 1e-9);
        ASSERT_TRUE(dynamic_cast<const ElementarySurface<T> *>(&so) != nullptr); // exact primitive
    }
    // the sphere, outwards and inwards
    const auto sph = read(m, "sphere.stp");
    const auto sf = explore<FaceId>(m, sph).front();
    for (T d : {0.5, -1.})
    {
        auto g = unwrap(offset_face(m, sf, d));
        expect_valid(m, ShapeId{g});
        ASSERT_NEAR(dynamic_cast<const SphericalSurface<T> &>(*m.face(g).surface).radius(), 3. + d, 1e-12);
    }
    // a cone (pcurves shifted in v)
    const auto cone = read(m, "cone.stp");
    for (auto f : explore<FaceId>(m, cone))
    {
        INFO("cone face " << typeid(*m.face(f).surface).name());
        expect_valid(m, ShapeId{unwrap(offset_face(m, f, 0.3))});
    }
    // a B-spline face: the generic offset surface
    using namespace brep_tests;
    const auto bump = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.5}, {0., 2., 0.3}, {2., 2., 1.})));
    auto g = unwrap(offset_face(m, bump, 0.2));
    expect_valid(m, ShapeId{g});
    const auto &s = *m.face(bump).surface, &so = *m.face(g).surface;
    for (T u : {0.2, 0.7})
        for (T v : {0.3, 0.9})
        {
            const auto n = cross(s.value(u, v, 1, 0), s.value(u, v, 0, 1));
            ASSERT_LT(norm(so.value(u, v) - (s.value(u, v) + (0.2 / norm(n)) * n)), 1e-12);
        }
    ASSERT_FALSE(offset_face(m, FaceId{999}, 1.).has_value());
}

namespace
{
    /// A closed curve r = 1 + a sin(k t) in z = 0, neither a primitive nor a B-spline.
    std::shared_ptr<Curve<T, 3>> flower(T a, int k)
    {
        return std::make_shared<brep_tests::AnalyticCurve>(
            [a, k](T t) {
                const T r = 1. + a * std::sin(k * t);
                return point<T, 3>{r * std::cos(t), r * std::sin(t), 0.};
            },
            std::array<T, 2>{0., 2. * pi});
    }

    std::array<std::size_t, 7> counts(const Model<T> &m)
    {
        std::array<std::size_t, 7> n{};
        for (std::size_t t = 0; t < n.size(); ++t)
            n[t] = m.count(ShapeType(t));
        return n;
    }
}

TEST(tests_brep_sweep, curves_approximated_within_approx_tol)
{
    // a curve neither primitive nor B-spline was interpolated through 65 points whatever
    // approx_tol: r = 1 + 0.2 sin 7t swept 1.1e-2 away from it, approx_tol 1e-5
    Model<T> m;
    const auto c = flower(0.2, 7);
    const auto w = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, c))}));
    const auto r = unwrap(make_prism(m, ShapeId{w}, point<T, 3>{0., 0., 1.}));
    ASSERT_EQ(r.approximated, 1u);
    ASSERT_GT(r.deviation, 0.);
    ASSERT_LE(r.deviation, 1e-5);
    // the bottom line of the swept surface against the curve, densely, both ways
    const auto &srf = *m.face(explore<FaceId>(m, r.shape).front()).surface;
    const auto b = srf.bounds();
    std::vector<point<T, 3>> bottom, curve;
    for (int i = 0; i <= 8000; ++i)
    {
        bottom.push_back(srf.value(b[0] + (b[1] - b[0]) * i / 8000., b[2]));
        curve.push_back(c->value(2. * pi * i / 8000.));
    }
    auto far = [](const std::vector<point<T, 3>> &x, const std::vector<point<T, 3>> &y) {
        T worst{};
        for (std::size_t i = 0; i < x.size(); i += 8)
        {
            T d = 1e300;
            for (const auto &q : y)
                d = std::min(d, norm(x[i] - q));
            worst = std::max(worst, d);
        }
        return worst;
    };
    ASSERT_LE(std::max(far(curve, bottom), far(bottom, curve)), 1e-5);
    // not reachable within the cap of the refinement: refused, nothing left
    SweepOptions<T> fine;
    fine.approx_tol = 1e-14;
    const auto before = counts(m);
    ASSERT_FALSE(make_prism(m, ShapeId{w}, point<T, 3>{0., 0., 1.}, fine).has_value());
    ASSERT_EQ(counts(m), before);
}

TEST(tests_brep_sweep, face_profile_of_any_curve)
{
    // a disk bounded by a circle given as such a curve: its caps are copies of the profile, which
    // copy_shape could not move (an exception out of make_prism, the side faces left behind);
    // the sides sew to the caps within the deviation of their approximation
    Model<T> m;
    points_vector<T, 3> poles{{-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}};
    const std::vector<T> k{-2., -2., 2., 2.};
    const auto xy = std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);
    const auto disk = unwrap(make_face(m, xy, unwrap(make_wire(m, std::vector{unwrap(make_edge(m, flower(0., 1)))}))));
    const auto r = unwrap(make_prism(m, ShapeId{disk}, point<T, 3>{0., 0., 2.}));
    expect_valid(m, r.shape);
    ASSERT_NEAR(volume(m, r.shape), 2. * pi, 2e-3 * 2. * pi);
    ASSERT_EQ(r.approximated, 1u);
}

TEST(tests_brep_sweep, failure_leaves_the_model)
{
    // a prism of a square along its own plane encloses no volume: refused, and the faces swept
    // so far no longer stay in the model
    Model<T> m;
    const auto sq = unwrap(make_face(m, brep_tests::quad({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.})));
    const auto before = counts(m);
    const auto r = make_prism(m, ShapeId{sq}, point<T, 3>{1., 0., 0.});
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::ZeroVolume);
    ASSERT_EQ(counts(m), before);
}

TEST(tests_brep_sweep, offset_of_an_apex_and_beyond_the_curvature)
{
    // the apex of a half cone becomes an arc: its two ends, and the generatrices ending there,
    // shared one vertex (vertex off curve twice)
    Model<T> m;
    const auto cone = std::make_shared<ConicalSurface<T>>(0., 0.4, ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., pi, 0., 2.});
    const auto f = unwrap(make_face(m, cone));
    const auto o = unwrap(offset_face(m, f, 0.1));
    expect_valid(m, ShapeId{o});
    std::size_t arcs{};
    for (auto e : explore<EdgeId>(m, ShapeId{o}))
        arcs += !m.edge(e).degenerate && m.edge(e).v1 != m.edge(e).v2 && norm(m.vertex(m.edge(e).v1).pnt - m.vertex(m.edge(e).v2).pnt) < 0.2 &&
                std::abs(m.vertex(m.edge(e).v1).pnt[2] - m.vertex(m.edge(e).v2).pnt[2]) < 1e-9;
    ASSERT_GE(arcs, 1u); // the arc of the apex, between two vertices

    // a sphere of radius 1 offset by -1 (a point) or -1.5 (turned inside out): refused, where it
    // gave a face of an offset surface, valid for check
    const auto sph = std::make_shared<SphericalSurface<T>>(1., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., 2. * pi, -pi / 2., pi / 2.});
    const auto s = unwrap(make_face(m, sph));
    const auto before = counts(m);
    for (T d : {-1., -1.5})
    {
        INFO("d = " << d);
        const auto r = offset_face(m, s, d);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::InvalidBounds);
        ASSERT_EQ(counts(m), before);
    }
    ASSERT_TRUE(offset_face(m, s, -0.5).has_value());
}

TEST(tests_brep_sweep, dead_reference_keeps_its_reason)
{
    // a profile referring to an erased vertex: InvalidId (it was PCurveApproximation)
    Model<T> m;
    const auto f = brep_tests::box_faces(m).front();
    const auto v = explore<VertexId>(m, f).front();
    m.erase(v);
    auto r = make_prism(m, f, point<T, 3>{0., 0., 1.});
    ASSERT_FALSE(r.has_value());
    ASSERT_TRUE(r.error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(r.error().shapes.front() == ShapeId{v});
}
