#include <doctest_gtest.hpp>
#include <gbs-brep/thick.h>
#include <gbs-brep/sweep.h>
#include <gbs-brep/pipe.h>
#include <gbs/bscinterp.h>
#include <gbs-brep/check.h>
#include <gbs-brep/mesh.h>

#include <array>
#include <cmath>
#include <numbers>
#include <sstream>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    void expect_valid(const Model<T> &m, const ShapeId &s)
    {
        const auto rep = check(m, s);
        std::ostringstream os;
        for (const auto &e : rep.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        INFO(os.str());
        REQUIRE(rep.ok());
    }

    /// The volume of the mesh of the shape, at a deflection of 5e-5 (signed_volume is not exact
    /// on faces with holes).
    T volume(const Model<T> &m, const ShapeId &s)
    {
        MeshOptions<T> o;
        o.deflection = 5e-5;
        o.max_points = 2'000'000;
        return gbs::brep::volume(unwrap(mesh(m, s, o)));
    }

    SolidId solid_of(Model<T> &m, const BuildResult<SweepResult<T>> &r)
    {
        return explore<SolidId>(m, unwrap(r).shape).front();
    }

    /// The volume within 3e-4 (the mesh is inside the curved faces by up to its deflection).
    void expect_volume(const Model<T> &m, const ShapeId &s, T v)
    {
        ASSERT_NEAR(volume(m, s), v, 3e-4 * std::abs(v));
    }

    /// The planar face bounded by the closed chain of edges `e`, in the plane of normal n through o.
    FaceId planar(Model<T> &m, const std::vector<EdgeId> &e, const point<T, 3> &o, const point<T, 3> &n, const point<T, 3> &x)
    {
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{o, n, x}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    }

    /// A circle (or arc) of radius r centred at c in the plane (x, y), from t1 to t2.
    std::shared_ptr<EllipseCurve<T>> arc(T r, const point<T, 3> &c, const point<T, 3> &x, const point<T, 3> &y, T t1 = 0., T t2 = 2 * pi)
    {
        return std::make_shared<EllipseCurve<T>>(r, r, ax2<T, 3>{{c, cross(x, y), x}}, t1, t2);
    }

    const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};

    SolidId revolved(Model<T> &m, const std::vector<EdgeId> &e)
    {
        return solid_of(m, make_revolution(m, planar(m, e, {}, {0., -1., 0.}, X), ax1<T, 3>{{{0., 0., 0.}, Z}}));
    }

    SolidId sphere(Model<T> &m, T r)
    {
        const auto e1 = unwrap(make_edge(m, arc(r, {}, X, Z, -pi / 2, pi / 2)));
        const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., r}, point<T, 3>{0., 0., -r}));
        return revolved(m, {e1, e2});
    }

    SolidId torus(Model<T> &m, T R, T r)
    {
        return revolved(m, {unwrap(make_edge(m, arc(r, {R, 0., 0.}, X, Z)))});
    }

    /// The frustum of radii r0 at z = 0 and r1 at z = h.
    SolidId frustum(Model<T> &m, T r0, T r1, T h)
    {
        const point<T, 3> a{0., 0., 0.}, b{r0, 0., 0.}, c{r1, 0., h}, d{0., 0., h};
        return revolved(m, {unwrap(make_edge(m, a, b)), unwrap(make_edge(m, b, c)), unwrap(make_edge(m, c, d)), unwrap(make_edge(m, d, a))});
    }

    /// A slab of height H on the stadium of straight sides L and half circles of radius R.
    SolidId slab(Model<T> &m, T L, T R, T H)
    {
        const auto v = std::array{unwrap(make_vertex(m, point<T, 3>{0., -R, 0.})), unwrap(make_vertex(m, point<T, 3>{L, -R, 0.})),
                                  unwrap(make_vertex(m, point<T, 3>{L, R, 0.})), unwrap(make_vertex(m, point<T, 3>{0., R, 0.}))};
        const auto right = arc(R, {L, 0., 0.}, X, Y, -pi / 2, pi / 2), left = arc(R, {0., 0., 0.}, X, Y, pi / 2, 3 * pi / 2);
        const std::vector<EdgeId> e{unwrap(make_edge(m, v[0], v[1])), unwrap(make_edge(m, right, -pi / 2, pi / 2, v[1], v[2])),
                                    unwrap(make_edge(m, v[2], v[3])), unwrap(make_edge(m, left, pi / 2, 3 * pi / 2, v[3], v[0]))};
        return solid_of(m, make_prism(m, planar(m, e, {}, Z, X), point<T, 3>{0., 0., H}));
    }

    std::array<std::size_t, 7> counts(const Model<T> &m)
    {
        std::array<std::size_t, 7> n{};
        for (std::size_t t = 0; t < n.size(); ++t)
            n[t] = m.count(ShapeType(t));
        return n;
    }

    template <typename S>
    std::size_t count(const Model<T> &m, const ShapeId &s)
    {
        std::size_t n{};
        for (auto f : explore<FaceId>(m, s))
            n += dynamic_cast<const S *>(m.face(f).surface.get()) != nullptr;
        return n;
    }

    /// The planar face bounded by the polygon `p` (counter-clockwise around z).
    FaceId polygon(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., p[0][2]}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                unwrap(make_wire(m, e))));
    }

    /// The planar face bounded by the polygon `p`, in the plane of normal n.
    FaceId polygon_in(Model<T> &m, const std::vector<point<T, 3>> &p, const point<T, 3> &n, const point<T, 3> &x)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return planar(m, e, p[0], n, x);
    }

    SolidId box(Model<T> &m, T a, T b, T c)
    {
        return solid_of(m, make_prism(m, polygon(m, {{0., 0., 0.}, {a, 0., 0.}, {a, b, 0.}, {0., b, 0.}}), point<T, 3>{0., 0., c}));
    }

    SolidId cylinder(Model<T> &m, T r, T h)
    {
        const ax2<T, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }
}

TEST(tests_brep_thick, box)
{
    for (T d : {0.25, -0.2})
    {
        INFO("d = " << d);
        Model<T> m;
        const auto s = box(m, 2., 1.5, 1.);
        const auto before = counts(m);
        const auto o = unwrap(offset_solid(m, s, d));
        expect_valid(m, ShapeId{o});
        ASSERT_EQ(explore<FaceId>(m, ShapeId{o}).size(), 6u);
        ASSERT_EQ(count<Plane<T>>(m, ShapeId{o}), 6u);
        ASSERT_NEAR(volume(m, ShapeId{o}), (2. + 2 * d) * (1.5 + 2 * d) * (1. + 2 * d), 1e-9);
        (void)before;
    }
    // shrunk beyond its half thickness: a face vanishes
    Model<T> m;
    const auto s = box(m, 2., 1.5, 1.);
    const auto before = counts(m);
    const auto r = offset_solid(m, s, -0.6);
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    ASSERT_EQ(counts(m), before);
}

TEST(tests_brep_thick, cylinder)
{
    for (T d : {0.3, -0.25})
    {
        INFO("d = " << d);
        Model<T> m;
        const auto s = cylinder(m, 1., 2.);
        const auto o = unwrap(offset_solid(m, s, d));
        expect_valid(m, ShapeId{o});
        ASSERT_EQ(count<CylindricalSurface<T>>(m, ShapeId{o}), 1u);
        expect_volume(m, ShapeId{o}, pi * (1. + d) * (1. + d) * (2. + 2 * d));
    }
}

TEST(tests_brep_thick, sphere)
{
    for (T d : {0.3, -0.2})
    {
        INFO("d = " << d);
        Model<T> m;
        const auto s = sphere(m, 1.);
        const auto o = unwrap(offset_solid(m, s, d));
        expect_valid(m, ShapeId{o});
        expect_volume(m, ShapeId{o}, 4. / 3. * pi * std::pow(1. + d, 3));
    }
    // a sphere shrunk to nothing: the offset folds
    Model<T> m;
    const auto s = sphere(m, 1.);
    ASSERT_FALSE(offset_solid(m, s, -1.2).has_value());
}

TEST(tests_brep_thick, torus)
{
    for (T d : {0.3, -0.2})
    {
        INFO("d = " << d);
        Model<T> m;
        const auto t = torus(m, 2., 0.5);
        const auto ot = unwrap(offset_solid(m, t, d));
        expect_valid(m, ShapeId{ot});
        expect_volume(m, ShapeId{ot}, 2. * pi * pi * 2. * (0.5 + d) * (0.5 + d));
    }
}

TEST(tests_brep_thick, frustum)
{
    // the cone offset by d along its normal, its radius by d / cos a; the disks moved by d
    for (T d : {0.3, -0.2})
    {
        INFO("d = " << d);
        Model<T> m;
        const T r0 = 2., r1 = 1., h = 2., k = (r0 - r1) / h;
        const auto f = frustum(m, r0, r1, h);
        const auto of = unwrap(offset_solid(m, f, d));
        expect_valid(m, ShapeId{of});
        ASSERT_EQ(count<ConicalSurface<T>>(m, ShapeId{of}), 1u);
        const T a = r0 + d * std::sqrt(1. + k * k);
        auto F = [&](T z) { return -std::pow(a - k * z, 3) / (3. * k); };
        expect_volume(m, ShapeId{of}, pi * (F(h + d) - F(-d)));
    }
}

TEST(tests_brep_thick, tangent_edges)
{
    // a slab on a stadium: planes and half cylinders tangent along the vertical edges
    for (T d : {0.2, -0.3})
    {
        INFO("d = " << d);
        Model<T> m;
        const T L = 3., R = 1., H = 0.8;
        const auto s = slab(m, L, R, H);
        const auto o = unwrap(offset_solid(m, s, d));
        expect_valid(m, ShapeId{o});
        ASSERT_EQ(count<CylindricalSurface<T>>(m, ShapeId{o}), 2u);
        INFO("volume " << volume(m, ShapeId{o}) << " expected " << (2. * (R + d) * L + pi * (R + d) * (R + d)) * (H + 2 * d));
        expect_volume(m, ShapeId{o}, (2. * (R + d) * L + pi * (R + d) * (R + d)) * (H + 2 * d));
    }
}

TEST(tests_brep_thick, thick_solids)
{
    const T a = 2., b = 1.5, c = 1., t = 0.1;
    {
        // a box opened at its top
        Model<T> m;
        const auto s = box(m, a, b, c);
        FaceId top;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
        {
            const auto &face = m.face(f);
            const auto bx = classify_detail::face_uv_box(m, f);
            const auto p = face.surface->value((bx[0] + bx[1]) / 2, (bx[2] + bx[3]) / 2);
            if (std::abs(p[2] - c) < 1e-9)
                top = f;
        }
        const auto before = counts(m);
        const auto r = unwrap(make_thick(m, s, std::vector<FaceId>{top}, t));
        expect_valid(m, ShapeId{r});
        INFO("open box " << volume(m, ShapeId{r}) << " expected " << a * b * c - (a - 2 * t) * (b - 2 * t) * (c - t));
        expect_volume(m, ShapeId{r}, a * b * c - (a - 2 * t) * (b - 2 * t) * (c - t));
        ASSERT_EQ(explore<FaceId>(m, ShapeId{r}).size(), 11u);
        // outward
        const auto g = unwrap(make_thick(m, s, std::vector<FaceId>{top}, -t));
        expect_valid(m, ShapeId{g});
        expect_volume(m, ShapeId{g}, (a + 2 * t) * (b + 2 * t) * (c + t) - a * b * c);
        // closed: a void
        const auto h = unwrap(make_thick(m, s, std::vector<FaceId>{}, t));
        expect_valid(m, ShapeId{h});
        ASSERT_EQ(m.solid(h).voids.size(), 1u);
        expect_volume(m, ShapeId{h}, a * b * c - (a - 2 * t) * (b - 2 * t) * (c - 2 * t));
        (void)before;
    }
    {
        // a cup: a cylinder opened at its top
        Model<T> m;
        const T r = 1., h = 2.;
        const auto s = cylinder(m, r, h);
        FaceId top;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (const auto *p = dynamic_cast<const Plane<T> *>(m.face(f).surface.get()); p && p->frame().o[2] > h / 2)
                top = f;
        const auto cup = unwrap(make_thick(m, s, std::vector<FaceId>{top}, t));
        expect_valid(m, ShapeId{cup});
        expect_volume(m, ShapeId{cup}, pi * r * r * h - pi * (r - t) * (r - t) * (h - t));
    }
}

TEST(tests_brep_thick, bspline_faces)
{
    // a square swept along a helix: four B-spline faces with sharp edges between them, two caps;
    // offset, the square grows by 2d and the path goes on by d at each end, about
    const T a = 0.4, R = 2., pitch = 1.5;
    for (T d : {0.05, -0.05})
    {
        INFO("d = " << d);
        Model<T> m;
        std::vector<T> u(97);
        std::vector<constrType<T, 3, 1>> q(97);
        for (std::size_t i = 0; i <= 96; ++i)
        {
            u[i] = T(i) / 96.;
            const T t = 2 * pi * u[i];
            q[i] = {point<T, 3>{R * std::cos(t), R * std::sin(t), pitch * t / (2 * pi)}};
        }
        const auto path = unwrap(make_edge(m, std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(q, u, std::size_t{3}))));
        const auto c = m.edge(path).curve;
        const auto t0 = c->value(0., 1) / norm(c->value(0., 1));
        auto x = X - t0[0] * t0;
        x = x / norm(x);
        const auto y = cross(t0, x);
        const point<T, 3> o{R, 0., 0.};
        const auto sq = polygon_in(m, {o - a / 2 * x - a / 2 * y, o + a / 2 * x - a / 2 * y, o + a / 2 * x + a / 2 * y, o - a / 2 * x + a / 2 * y}, t0, x);
        const auto s = solid_of(m, make_pipe(m, ShapeId{sq}, ShapeId{path}));
        const T L = std::hypot(2 * pi * R, pitch);
        const auto before = counts(m);
        const auto r = offset_solid(m, s, d);
        if (!r)
            FAIL(r.error().message);
        expect_valid(m, ShapeId{*r});
        ASSERT_EQ(explore<FaceId>(m, ShapeId{*r}).size(), 6u);
        // the section (a + 2d)² along the path, gone on by d at each end
        ASSERT_NEAR(volume(m, ShapeId{*r}), (a + 2 * d) * (a + 2 * d) * (L + 2 * d), 2e-3 * (a + 2 * d) * (a + 2 * d) * L);
        (void)before;
    }
}

TEST(tests_brep_thick, refusals)
{
    Model<T> m;
    // the apex of a cone: its offset is not a point
    const auto cone = revolved(m, {unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.})),
                                   unwrap(make_edge(m, point<T, 3>{1., 0., 0.}, point<T, 3>{0., 0., 2.})),
                                   unwrap(make_edge(m, point<T, 3>{0., 0., 2.}, point<T, 3>{0., 0., 0.}))});
    const auto before = counts(m);
    for (T d : {0.1, -0.1})
    {
        const auto r = offset_solid(m, cone, d);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
        ASSERT_EQ(counts(m), before);
    }
    // arguments
    const auto b = box(m, 1., 1., 1.);
    ASSERT_EQ(offset_solid(m, b, 1e-9).error().code, BuildErrc::InvalidBounds);
    ASSERT_EQ(offset_solid(m, b, std::numeric_limits<T>::quiet_NaN()).error().code, BuildErrc::InvalidBounds);
    ASSERT_EQ(offset_solid(m, SolidId{999}, 0.1).error().code, BuildErrc::InvalidId);
    const auto other = explore<FaceId>(m, ShapeId{cone}).front();
    ASSERT_EQ(make_thick(m, b, std::vector<FaceId>{other}, 0.1).error().code, BuildErrc::InvalidId);
    ASSERT_EQ(make_thick(m, b, std::vector<FaceId>{}, 0.).error().code, BuildErrc::InvalidBounds);
    // too thick: the inner solid vanishes
    const auto r = make_thick(m, b, std::vector<FaceId>{}, 0.6);
    ASSERT_FALSE(r.has_value());
}

TEST(tests_brep_thick, self_intersection)
{
    // a square ring cut open at its top, the cut 0.4 wide: offset outward by 0.3, its two ends
    // pass each other (no face vanishes: the offset intersects itself)
    Model<T> m;
    const auto f = polygon(m, {{0., 0., 0.}, {4., 0., 0.}, {4., 4., 0.}, {2.2, 4., 0.}, {2.2, 3., 0.}, {3., 3., 0.}, {3., 1., 0.}, {1., 1., 0.},
                               {1., 3., 0.}, {1.8, 3., 0.}, {1.8, 4., 0.}, {0., 4., 0.}});
    const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., 1.}));
    const auto before = counts(m);
    const auto r = offset_solid(m, s, 0.3);
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    ASSERT_NE(r.error().message.find("intersects itself"), std::string::npos);
    ASSERT_EQ(counts(m), before);
    // by 0.1 the cut stays open
    const auto g = unwrap(offset_solid(m, s, 0.1));
    expect_valid(m, ShapeId{g});
}
