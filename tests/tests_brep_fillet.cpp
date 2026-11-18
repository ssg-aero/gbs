#include <doctest_gtest.hpp>
#include <gbs-brep/fillet.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/sweep.h>
#include <gbs/bscinterp.h>

#include <array>
#include <cmath>
#include <numbers>
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
    const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};

    SolidId solid_of(Model<T> &m, const BuildResult<SweepResult<T>> &r)
    {
        return explore<SolidId>(m, unwrap(r).shape).front();
    }

    /// The planar face bounded by the polygon `p`, in the plane of normal n through p[0].
    FaceId polygon(Model<T> &m, const std::vector<point<T, 3>> &p, const point<T, 3> &n = Z, const point<T, 3> &x = X)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n, x}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    }

    SolidId prism(Model<T> &m, const std::vector<point<T, 3>> &p, T h)
    {
        return solid_of(m, make_prism(m, polygon(m, p), point<T, 3>{0., 0., h}));
    }

    SolidId box(Model<T> &m, T a, T b, T c, const point<T, 3> &o = {})
    {
        return prism(m, {o, o + point<T, 3>{a, 0., 0.}, o + point<T, 3>{a, b, 0.}, o + point<T, 3>{0., b, 0.}}, c);
    }

    SolidId cylinder(Model<T> &m, T r, T h, const point<T, 3> &o = {})
    {
        const ax2<T, 3> ax{{o, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    /// The edge of the solid whose middle is nearest p.
    EdgeId edge_near(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        EdgeId best{};
        T gap = 1e300;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            for (int i = 0; i <= 16; ++i)
            {
                const T d = norm(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / 16.) - p);
                if (d < gap)
                    gap = d, best = e;
            }
        }
        return best;
    }

    point<T, 3> unit_of(const point<T, 3> &v) { return v / norm(v); }

    /// The outward normal of the face f (in the solid s) at uv on its surface continued, `support`.
    point<T, 3> outward(const Model<T> &m, SolidId s, FaceId f, const Surface<T, 3> &support, const point<T, 2> &uv)
    {
        auto n = normal_direction(support, uv[0], uv[1]);
        for (auto sh : explore<ShellId>(m, ShapeId{s}))
            for (const auto &fu : m.shell(sh).faces)
                if (fu.face == f && fu.orient == Orientation::Reversed)
                    n = T(-1) * n;
        return n;
    }

    /**
     * The blend touches the faces: along the contacts, the surface of the blend runs through the
     * contact (within `gap`), and its normal is that of the face (sine under 1e-6); for a fillet,
     * the contacts are at r from the centre line along the normals, toward the inside of a convex
     * edge, and the surface at r from the centre line.
     */
    void expect_blend(const Model<T> &m, SolidId s, const EdgeBlend<T> &b, std::optional<T> r, T gap, T angle = 1e-9)
    {
        REQUIRE(b.surface);
        REQUIRE(b.b > b.a);
        for (int i = 0; i <= 32; ++i)
        {
            const T t = b.a + (b.b - b.a) * T(i) / 32.;
            const auto c = b.spine->value(t);
            for (const auto &ct : b.contacts)
            {
                const auto p = ct.curve->value(t);
                const auto q = ct.on_blend->value(t);
                ASSERT_LE(norm(b.surface->value(q[0], q[1]) - p), gap);
                const auto uv = ct.pcurve->value(t);
                ASSERT_LE(norm(ct.support->value(uv[0], uv[1]) - p), 1e-12 + 1e-12 * norm(p));
                const auto nf = outward(m, s, ct.face, *ct.support, uv);
                const auto nb = normal_direction(*b.surface, q[0], q[1]);
                ASSERT_LE(norm(cross(nf, nb)), r ? angle : 1.); // a chamfer is not tangent
                if (r)
                {
                    ASSERT_NEAR(norm(p - c), *r, std::max(1e-7, gap / 10.));
                    // c = p - r n toward the inside of a convex edge, outside a concave one
                    ASSERT_LE(norm(c - (p + (b.convex ? -*r : *r) * nf)), std::max(1e-7, gap / 10.));
                }
            }
            if (r)
                for (int j = 0; j <= 8; ++j)
                {
                    // across the blend: between the pcurves of the contacts on it
                    const auto q1 = b.contacts[0].on_blend->value(t), q2 = b.contacts[1].on_blend->value(t);
                    const T w = T(j) / 8.;
                    const auto q = q1 + w * (q2 - q1);
                    ASSERT_NEAR(norm(b.surface->value(q[0], q[1]) - c), *r, gap + 1e-9);
                }
        }
    }
}

TEST(tests_brep_fillet, box_edge)
{
    // the edge along y at x = 2, z = 1 of a box, convex: a cylinder of radius r
    Model<T> m;
    const auto s = box(m, 2., 1.5, 1.);
    const auto e = edge_near(m, s, {2., 0.75, 1.});
    const T r = 0.3;
    const auto b = unwrap(fillet_geometry(m, s, e, r));
    ASSERT_TRUE(b.convex);
    ASSERT_TRUE(b.exact);
    const auto *cyl = dynamic_cast<const CylindricalSurface<T> *>(b.surface.get());
    ASSERT_TRUE(cyl);
    ASSERT_NEAR(cyl->radius(), r, 1e-12);
    expect_blend(m, s, b, r, 1e-9);
    // the centre line at (2 - r, ., 1 - r), from y = 0 to 1.5 and past by r
    for (T t : {b.t1, b.t2})
    {
        const auto c = b.spine->value(t);
        ASSERT_NEAR(c[0], 2. - r, 1e-9);
        ASSERT_NEAR(c[2], 1. - r, 1e-9);
    }
    ASSERT_NEAR(std::abs(b.spine->value(b.t2)[1] - b.spine->value(b.t1)[1]), 1.5, 1e-9);
    ASSERT_NEAR(norm(b.spine->value(b.b) - b.spine->value(b.a)), 1.5 + 2 * r, 1e-9);
    // a quarter of a turn
    const auto bx = cyl->bounds();
    ASSERT_NEAR(bx[1] - bx[0], pi / 2, 1e-9);
}

namespace
{
    /// The chamfer contacts at d1 and d2 from the edge, in the plane across it.
    void expect_chamfer(const Model<T> &m, SolidId s, const EdgeBlend<T> &b, T d1, T d2, T gap)
    {
        expect_blend(m, s, b, std::nullopt, gap);
        const std::array<T, 2> d{d1, d2};
        for (int i = 0; i <= 32; ++i)
        {
            const T t = b.a + (b.b - b.a) * T(i) / 32.;
            const auto e = b.spine->value(t), te = unit_of(b.spine->value(t, 1));
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto p = b.contacts[k].curve->value(t);
                ASSERT_NEAR(norm(p - e), d[k], gap + 1e-9);
                ASSERT_NEAR((p - e) * te, 0., gap + 1e-9);
            }
        }
    }

    /// A cylinder of radius r and length h along the axis `z` from o.
    SolidId cylinder_along(Model<T> &m, T r, T h, const point<T, 3> &o, const point<T, 3> &x, const point<T, 3> &y, const point<T, 3> &z)
    {
        const auto c = cylinder(m, r, h);
        return std::get<SolidId>(copy_shape(m, ShapeId{c}, Placement<T>::frame(o, x, y, z)));
    }

    SolidId fused(Model<T> &m, SolidId a, SolidId b)
    {
        const auto r = unwrap(fuse(m, ShapeId{a}, ShapeId{b}));
        REQUIRE(r.solids.size() == 1u);
        return r.solids.front();
    }
}

TEST(tests_brep_fillet, concave_edge)
{
    // the inner vertical edge of an L-shaped prism: the ball outside, a cylinder
    Model<T> m;
    const auto s = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, 1.);
    const auto e = edge_near(m, s, {1., 1., 0.5});
    const T r = 0.25;
    const auto b = unwrap(fillet_geometry(m, s, e, r));
    ASSERT_FALSE(b.convex);
    ASSERT_TRUE(b.exact);
    ASSERT_TRUE(dynamic_cast<const CylindricalSurface<T> *>(b.surface.get()));
    expect_blend(m, s, b, r, 1e-9);
    const auto c = b.spine->value(b.t1);
    ASSERT_NEAR(c[0], 1. + r, 1e-9);
    ASSERT_NEAR(c[1], 1. + r, 1e-9);

    const auto ch = unwrap(chamfer_geometry(m, s, e, 0.2, 0.3));
    ASSERT_TRUE(ch.exact);
    ASSERT_TRUE(dynamic_cast<const Plane<T> *>(ch.surface.get()));
    expect_chamfer(m, s, ch, 0.2, 0.3, 1e-9);
}

TEST(tests_brep_fillet, cylinder_edges)
{
    // the top edge of a cylinder: a torus of major radius R - r (convex), and the edge at the
    // foot of a cylinder standing on a box: R + r (concave)
    const T R = 0.6, r = 0.2;
    {
        Model<T> m;
        const auto s = cylinder(m, R, 1.);
        const auto e = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), 1.});
        const auto b = unwrap(fillet_geometry(m, s, e, r));
        ASSERT_TRUE(b.convex);
        ASSERT_TRUE(b.closed);
        const auto *t = dynamic_cast<const ToroidalSurface<T> *>(b.surface.get());
        ASSERT_TRUE(t);
        ASSERT_NEAR(t->major_radius(), R - r, 1e-9);
        ASSERT_NEAR(t->minor_radius(), r, 1e-12);
        ASSERT_NEAR(b.b - b.a, 2 * pi, 1e-9);
        expect_blend(m, s, b, r, 1e-9);

        // the chamfer: a cone
        const auto ch = unwrap(chamfer_geometry(m, s, e, 0.1, 0.15));
        ASSERT_TRUE(ch.exact);
        ASSERT_TRUE(dynamic_cast<const ConicalSurface<T> *>(ch.surface.get()));
        expect_chamfer(m, s, ch, 0.1, 0.15, 1e-9);
    }
    {
        Model<T> m;
        const auto s = fused(m, box(m, 2., 2., 1.), cylinder(m, R, 1., {1., 1., 1.}));
        const auto e = edge_near(m, s, {1. + R * std::cos(1.), 1. + R * std::sin(1.), 1.});
        const auto b = unwrap(fillet_geometry(m, s, e, r));
        ASSERT_FALSE(b.convex);
        ASSERT_TRUE(b.closed);
        const auto *t = dynamic_cast<const ToroidalSurface<T> *>(b.surface.get());
        ASSERT_TRUE(t);
        ASSERT_NEAR(t->major_radius(), R + r, 1e-9);
        expect_blend(m, s, b, r, 1e-9);
    }
}

TEST(tests_brep_fillet, sphere_edge)
{
    // a ball cut by the plane z = h: the circle of radius a = sqrt(R^2 - h^2), a torus
    Model<T> m;
    const T R = 1., h = 0.4, r = 0.15;
    const T a = std::sqrt(R * R - h * h);
    const auto e1 = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(R, R, ax2<T, 3>{{{0., 0., 0.}, cross(X, Z), X}}, -pi / 2, pi / 2)));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., R}, point<T, 3>{0., 0., -R}));
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                    unwrap(make_wire(m, std::vector<EdgeId>{e1, e2}))));
    const auto ball = solid_of(m, make_revolution(m, f, ax1<T, 3>{{{0., 0., 0.}, Z}}));
    const auto cr = unwrap(cut(m, ShapeId{ball}, ShapeId{box(m, 4., 4., 2., {-2., -2., h})}));
    REQUIRE(cr.solids.size() == 1u);
    const auto s = cr.solids.front();
    const auto e = edge_near(m, s, {a * std::cos(1.), a * std::sin(1.), h});
    const auto b = unwrap(fillet_geometry(m, s, e, r));
    ASSERT_TRUE(b.convex);
    const auto *t = dynamic_cast<const ToroidalSurface<T> *>(b.surface.get());
    ASSERT_TRUE(t);
    // the centre circle: at h - r, on the sphere of radius R - r
    ASSERT_NEAR(t->major_radius(), std::sqrt((R - r) * (R - r) - (h - r) * (h - r)), 1e-9);
    expect_blend(m, s, b, r, 1e-9);
}

TEST(tests_brep_fillet, crossed_cylinders)
{
    // a cylinder of radius 0.4 standing on one of radius 1 lying along x: the edge where they
    // meet is a saddle, the centre line a marched curve, the surface a rational B-spline
    Model<T> m;
    const auto lying = cylinder_along(m, 1., 4., {-2., 0., 0.}, Y, Z, X);
    const auto standing = cylinder(m, 0.4, 1.5, {0., 0., 0.});
    const auto s = fused(m, lying, standing);
    const auto e = edge_near(m, s, {0.4, 0., std::sqrt(1. - 0.16)});
    const T r = 0.15;
    FilletOptions<T> o;
    const auto b = unwrap(fillet_geometry(m, s, e, r, o));
    ASSERT_FALSE(b.convex);
    ASSERT_TRUE(b.closed);
    ASSERT_FALSE(b.exact);
    ASSERT_TRUE((dynamic_cast<const BSSurfaceRational<T, 3> *>(b.surface.get())));
    ASSERT_LE(b.deviation, o.approx_tol);
    expect_blend(m, s, b, r, o.approx_tol, o.angular_tol);

    const auto ch = unwrap(chamfer_geometry(m, s, e, 0.1, 0.1, o));
    ASSERT_FALSE(ch.exact);
    expect_chamfer(m, s, ch, 0.1, 0.1, o.approx_tol);
}

TEST(tests_brep_fillet, bspline_face)
{
    // the top edge of the prism of a B-spline bump: a plane and a B-spline face
    Model<T> m;
    const points_vector<T, 3> q{{0., 0., 0.}, {0.5, -0.3, 0.}, {1., 0.2, 0.}, {1.5, -0.1, 0.}, {2., 0., 0.}};
    const auto bump = std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(q, std::vector<T>{0., 0.25, 0.5, 0.75, 1.}, std::size_t{3}));
    const auto v0 = unwrap(make_vertex(m, point<T, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<T, 3>{2., 0., 0.}));
    const auto v2 = unwrap(make_vertex(m, point<T, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<T, 3>{0., 1., 0.}));
    const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                unwrap(make_edge(m, v3, v0))};
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                    unwrap(make_wire(m, w))));
    const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., 0.8}));
    const auto e = edge_near(m, s, {1., bump->value(0.5)[1], 0.8});
    const T r = 0.1;
    FilletOptions<T> o;
    const auto b = unwrap(fillet_geometry(m, s, e, r, o));
    ASSERT_TRUE(b.convex);
    ASSERT_FALSE(b.closed);
    ASSERT_FALSE(b.exact);
    expect_blend(m, s, b, r, o.approx_tol, o.angular_tol);

    const auto ch = unwrap(chamfer_geometry(m, s, e, 0.08, 0.05, o));
    expect_chamfer(m, s, ch, 0.08, 0.05, o.approx_tol);
}

TEST(tests_brep_fillet, refusals)
{
    Model<T> m;
    const auto s = box(m, 2., 1.5, 1.);
    const auto e = edge_near(m, s, {2., 0.75, 1.});
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    ASSERT_EQ(code(fillet_geometry(m, s, e, 0.)), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet_geometry(m, s, e, -0.1)), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, 0.1, 0.)), BuildErrc::InvalidBounds);
    // wider than the side face (height 1)
    ASSERT_EQ(code(fillet_geometry(m, s, e, 1.2)), BuildErrc::NonGeneric);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, 1.2, 1.2)), BuildErrc::NonGeneric);
    // an edge of another solid
    const auto other = box(m, 1., 1., 1., {5., 0., 0.});
    ASSERT_EQ(code(fillet_geometry(m, s, edge_near(m, other, {6., 0.5, 1.}), 0.1)), BuildErrc::InvalidId);

    // beyond the radius of the top edge of a cylinder
    const auto c = cylinder(m, 0.5, 1., {0., 5., 0.});
    ASSERT_EQ(code(fillet_geometry(m, c, edge_near(m, c, {0.5 * std::cos(1.), 5. + 0.5 * std::sin(1.), 1.}), 0.6)), BuildErrc::NonGeneric);
    // its seam: used twice by the side face
    ASSERT_EQ(code(fillet_geometry(m, c, edge_near(m, c, {0.5, 5., 0.5}), 0.1)), BuildErrc::InvalidId);
}

namespace
{
    /// A cylinder of radius R2 standing on one of radius 1 lying along x, and the edge where they meet.
    std::pair<SolidId, EdgeId> tee(Model<T> &m, T R2)
    {
        const auto s = fused(m, cylinder_along(m, 1., 4., {-2., 0., 0.}, Y, Z, X), cylinder(m, R2, 2.));
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            const auto p = ed.curve->value((ed.u1 + ed.u2) / 2.);
            if (std::abs(std::hypot(p[0], p[1]) - R2) < 1e-6 && std::abs(std::hypot(p[1], p[2]) - 1.) < 1e-6)
                return {s, e};
        }
        FAIL("no edge where the cylinders meet");
        return {};
    }
}

TEST(tests_brep_fillet, near_tangent_faces)
{
    // a tee whose cylinders have nearly the same radius: on the sides the faces are nearly
    // tangent and the offsets meet at a grazing angle, where the marched centre line loses its
    // accuracy; the stations of the ball are put back on the rolling ball by Gauss-Newton
    FilletOptions<T> o;
    for (T R2 : {0.99, 0.999})
        for (T r : {0.05, 0.15})
        {
            INFO("R2 = " << R2 << ", r = " << r);
            Model<T> m;
            const auto [s, e] = tee(m, R2);
            const auto b = unwrap(fillet_geometry(m, s, e, r, o));
            ASSERT_FALSE(b.convex);
            ASSERT_TRUE(b.closed);
            ASSERT_LE(b.deviation, o.approx_tol);
            expect_blend(m, s, b, r, o.approx_tol, o.angular_tol);
        }
    // the same radius: the faces tangent at the vertices where the two arcs of the edge meet;
    // the fillet vanishes there, its contacts meeting
    {
        Model<T> m;
        const auto [s, e] = tee(m, 1.);
        const T r = 0.15;
        const auto b = unwrap(fillet_geometry(m, s, e, r, o));
        ASSERT_FALSE(b.closed);
        // at the end of the edge where the faces are tangent (y = +-1; the other may be a seam)
        int tangent_ends = 0;
        for (const auto &[t, v] : {std::pair{b.t1, m.edge(e).v1}, std::pair{b.t2, m.edge(e).v2}})
            if (std::abs(std::abs(m.vertex(v).pnt[1]) - 1.) < 1e-9)
            {
                ++tangent_ends;
                ASSERT_LE(norm(b.contacts[0].curve->value(t) - b.contacts[1].curve->value(t)), 1e-6);
            }
        ASSERT_GE(tangent_ends, 1);
        const T lo = std::min(b.t1, b.t2), hi = std::max(b.t1, b.t2);
        for (int i = 0; i <= 64; ++i)
        {
            const T t = lo + (hi - lo) * T(i) / 64.;
            const auto c = b.spine->value(t);
            for (const auto &ct : b.contacts)
                ASSERT_NEAR(norm(ct.curve->value(t) - c), r, o.approx_tol);
        }
    }
}

TEST(tests_brep_fillet, sharp_wedges)
{
    // a prism on a triangle of angle a at the origin: the contacts at r / tan(a / 2) from the
    // edge, within the faces (2 long) or not
    for (const auto &[deg, r, fits] : {std::tuple{20., 0.2, true}, std::tuple{5., 0.05, true}, std::tuple{5., 0.2, false}, std::tuple{1., 0.01, true},
                                       std::tuple{1., 0.05, false}})
    {
        INFO(deg << " degrees, r = " << r);
        Model<T> m;
        const T a = deg * pi / 180.;
        const auto s = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {2. * std::cos(a), 2. * std::sin(a), 0.}}, 1.);
        const auto e = edge_near(m, s, {0., 0., 0.5});
        const auto b = fillet_geometry(m, s, e, r);
        ASSERT_EQ(b.has_value(), fits);
        if (fits)
        {
            ASSERT_TRUE((dynamic_cast<const CylindricalSurface<T> *>(b->surface.get())));
            expect_blend(m, s, *b, r, 1e-9);
            ASSERT_NEAR(norm(b->contacts[0].curve->value(b->t1) - b->spine->value(b->t1)), r, 1e-9);
            ASSERT_NEAR(norm(b->contacts[0].curve->value(b->t1) - point<T, 3>{0., 0., b->contacts[0].curve->value(b->t1)[2]}), r / std::tan(a / 2), 1e-9);
        }
        else
            ASSERT_EQ(b.error().code, BuildErrc::NonGeneric);
    }
}
