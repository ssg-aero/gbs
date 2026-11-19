#include <doctest_gtest.hpp>
#include <gbs-brep/fillet.h>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/sweep.h>
#include <gbs/bscbuild.h>
#include <gbs/bscinterp.h>
#include <gbs/bselementary.h>

#include <blade_fixture.h>

#include <array>
#include <iomanip>
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

    /// The face of the edge e whose outward normal at the middle of the edge is nearest the direction `dir`.
    FaceId face_toward(const Model<T> &m, SolidId s, EdgeId e, const point<T, 3> &dir)
    {
        const auto &ed = m.edge(e);
        FaceId best{};
        T dot = -2.;
        for (auto sh : explore<ShellId>(m, ShapeId{s}))
        {
            std::map<std::uint32_t, Orientation> orient;
            for (const auto &fu : m.shell(sh).faces)
                orient[fu.face.index] = fu.orient;
            const auto uses = shell_edge_uses(m, sh);
            if (const auto it = uses.find(e); it != uses.end())
                for (const auto &u : it->second)
                {
                    const auto &ce = m.wire(u.coedge.wire).coedges[u.coedge.index];
                    const auto uv = ce.pcurve->value((ed.u1 + ed.u2) / 2.);
                    auto n = normal_direction(*m.face(u.face).surface, uv[0], uv[1]);
                    if (orient.at(u.face.index) == Orientation::Reversed)
                        n = T(-1) * n;
                    if (const T d = n * dir / norm(dir); d > dot)
                        dot = d, best = u.face;
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

    /// 37 points of [a, b]: away from the stations of the contacts and of the approximated surfaces (powers of 2) but at the ends.
    constexpr int samples = 37;

    /**
     * The blend touches the faces: along the contacts, the surface of the blend runs through the
     * contact (within `gap`), and its normal is that of the face (sine under 1e-6); for a fillet,
     * the contacts are at r from the centre line along the normals, toward the inside of a convex
     * edge, and the surface at r from the centre line. The spine runs along the edge from t1 to t2.
     */
    void expect_blend(const Model<T> &m, SolidId s, const EdgeBlend<T> &b, std::optional<T> r, T gap, T angle = 1e-9)
    {
        REQUIRE(b.surface);
        REQUIRE(b.b > b.a);
        for (const auto &ct : b.contacts)
            ASSERT_LE(ct.deviation, gap);
        const auto &ed = m.edge(b.edge);
        ASSERT_GT((b.t2 > b.t1 ? 1. : -1.) * (b.spine->value(b.t1, 1) * ed.curve->value(ed.u1, 1)), 0.);
        for (int i = 0; i <= samples; ++i)
        {
            const T t = b.a + (b.b - b.a) * T(i) / T(samples);
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
    /// The chamfer contacts at d1 (on the face f) and d2 from the edge, in the plane across it.
    void expect_chamfer(const Model<T> &m, SolidId s, const EdgeBlend<T> &b, FaceId f, T d1, T d2, T gap)
    {
        ASSERT_EQ(b.contacts[0].face, f);
        expect_blend(m, s, b, std::nullopt, gap);
        const std::array<T, 2> d{d1, d2};
        for (int i = 0; i <= samples; ++i)
        {
            const T t = b.a + (b.b - b.a) * T(i) / T(samples);
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

    /// A circle (or arc) of radius r centred at c in the plane (x, y), from t1 to t2.
    std::shared_ptr<EllipseCurve<T>> arc(T r, const point<T, 3> &c, const point<T, 3> &x, const point<T, 3> &y, T t1 = 0., T t2 = 2 * pi)
    {
        return std::make_shared<EllipseCurve<T>>(r, r, ax2<T, 3>{{c, cross(x, y), x}}, t1, t2);
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

    // 0.2 on the face y = 1 (outward normal +y), 0.3 on the face x = 1, either way
    for (const auto &[dir, d1, d2] : {std::tuple{Y, 0.2, 0.3}, std::tuple{X, 0.2, 0.3}})
    {
        const auto f = face_toward(m, s, e, dir);
        const auto ch = unwrap(chamfer_geometry(m, s, e, f, d1, d2));
        ASSERT_TRUE(ch.exact);
        ASSERT_TRUE(dynamic_cast<const Plane<T> *>(ch.surface.get()));
        expect_chamfer(m, s, ch, f, d1, d2, 1e-9);
        // the contact on f at d1 from the edge, on f
        const auto p = ch.contacts[0].curve->value((ch.t1 + ch.t2) / 2.);
        ASSERT_NEAR(dir == Y ? p[1] : p[0], 1., 1e-12);
        ASSERT_NEAR(dir == Y ? p[0] : p[1], 1. + d1, 1e-9);
    }
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
        expect_blend(m, s, b, r, 1e-7);

        // the chamfer: a cone, its contacts on it within tol
        const auto top = face_toward(m, s, e, Z);
        const auto ch = unwrap(chamfer_geometry(m, s, e, top, 0.1, 0.15));
        ASSERT_TRUE(ch.exact);
        ASSERT_TRUE(dynamic_cast<const ConicalSurface<T> *>(ch.surface.get()));
        expect_chamfer(m, s, ch, top, 0.1, 0.15, 1e-6);
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
        expect_blend(m, s, b, r, 1e-7);
    }
}

TEST(tests_brep_fillet, open_arc_edges)
{
    // the top edges of half cylinders: arcs on a circular centre line, from a vertex across the
    // start of the circle (x > 0: from -pi/2 to pi/2, in two edges) or not (y > 0: from 0 to pi)
    const T R = 1., r = 0.2;
    for (const auto &[cutter, p] : {std::pair{point<T, 3>{-4., -2., -1.}, point<T, 3>{std::cos(0.3), std::sin(0.3), 1.}},
                                    std::pair{point<T, 3>{-2., -4., -1.}, point<T, 3>{std::cos(1.), std::sin(1.), 1.}}})
    {
        Model<T> m;
        const auto cr = unwrap(cut(m, ShapeId{cylinder(m, R, 1.)}, ShapeId{box(m, 4., 4., 4., cutter)}));
        REQUIRE(cr.solids.size() == 1u);
        const auto s = cr.solids.front();
        const auto e = edge_near(m, s, p);
        const auto &ed = m.edge(e);
        const auto b = unwrap(fillet_geometry(m, s, e, r));
        ASSERT_TRUE(b.convex);
        ASSERT_FALSE(b.closed);
        const auto *t = dynamic_cast<const ToroidalSurface<T> *>(b.surface.get());
        ASSERT_TRUE(t);
        ASSERT_NEAR(t->major_radius(), R - r, 1e-9);
        expect_blend(m, s, b, r, 1e-7);
        // the centre line at the vertices, continued past them by r (a parameter of r / (R - r))
        for (const auto &[tv, v] : {std::pair{b.t1, ed.v1}, std::pair{b.t2, ed.v2}})
        {
            const auto q = m.vertex(v).pnt;
            const auto c = b.spine->value(tv);
            const T k = (R - r) / std::hypot(q[0], q[1]);
            ASSERT_LE(norm(c - point<T, 3>{k * q[0], k * q[1], 1. - r}), 1e-9);
        }
        ASSERT_NEAR(std::min(b.t1, b.t2) - b.a, r / (R - r), 1e-9);
        ASSERT_NEAR(b.b - std::max(b.t1, b.t2), r / (R - r), 1e-9);

        // the chamfer of an edge of a boolean (a curve on a surface), continued past its vertices too
        const auto top = face_toward(m, s, e, Z);
        const auto ch = unwrap(chamfer_geometry(m, s, e, top, 0.1, 0.15));
        ASSERT_TRUE(ch.exact);
        ASSERT_TRUE(dynamic_cast<const ConicalSurface<T> *>(ch.surface.get()));
        ASSERT_LT(ch.a, ed.u1);
        ASSERT_GT(ch.b, ed.u2);
        ASSERT_NEAR(norm(ch.spine->value(ch.a) - ch.spine->value(ed.u1)), 0.15, 1e-3);
        expect_chamfer(m, s, ch, top, 0.1, 0.15, 1e-6);
    }
}

TEST(tests_brep_fillet, acute_corners)
{
    // the top edges of an equilateral prism end at corners of 60 degrees: past them the contact on
    // the top leaves the face across the edges at its vertices, left to the faces at the ends
    Model<T> m;
    const auto s = prism(m, {{0., 0., 0.}, {1., 0., 0.}, {0.5, std::sqrt(0.75), 0.}}, 1.);
    const auto e = edge_near(m, s, {0.5, 0., 1.});
    const auto top = face_toward(m, s, e, Z);
    for (T r : {0.15, 0.3, 0.5})
    {
        const auto b = unwrap(fillet_geometry(m, s, e, r));
        expect_blend(m, s, b, r, 1e-9);
        const auto ch = unwrap(chamfer_geometry(m, s, e, top, r, r));
        expect_chamfer(m, s, ch, top, r, r, 1e-9);
    }
    // the contact on the side, of height 1, out of it at the middle
    ASSERT_EQ(fillet_geometry(m, s, e, 1.2).error().code, BuildErrc::NonGeneric);
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
    expect_blend(m, s, b, r, 1e-7);
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

    const auto f = face_toward(m, s, e, Z);
    const auto ch = unwrap(chamfer_geometry(m, s, e, f, 0.1, 0.1, o));
    ASSERT_FALSE(ch.exact);
    expect_chamfer(m, s, ch, f, 0.1, 0.1, o.approx_tol);
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

    const auto top = face_toward(m, s, e, Z);
    const auto ch = unwrap(chamfer_geometry(m, s, e, top, 0.08, 0.05, o));
    expect_chamfer(m, s, ch, top, 0.08, 0.05, o.approx_tol);
}

TEST(tests_brep_fillet, refusals)
{
    Model<T> m;
    const auto s = box(m, 2., 1.5, 1.);
    const auto e = edge_near(m, s, {2., 0.75, 1.});
    const auto top = face_toward(m, s, e, Z);
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    ASSERT_EQ(code(fillet_geometry(m, s, e, 0.)), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet_geometry(m, s, e, -0.1)), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, top, 0.1, 0.)), BuildErrc::InvalidBounds);
    // options
    FilletOptions<T> o;
    o.extension = -0.5;
    ASSERT_EQ(code(fillet_geometry(m, s, e, 0.1, o)), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, top, 0.1, 0.1, o)), BuildErrc::InvalidBounds);
    o = {};
    o.approx_tol = 0.;
    ASSERT_EQ(code(fillet_geometry(m, s, e, 0.1, o)), BuildErrc::InvalidBounds);
    // wider than the side face (height 1)
    ASSERT_EQ(code(fillet_geometry(m, s, e, 1.2)), BuildErrc::NonGeneric);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, top, 1.2, 1.2)), BuildErrc::NonGeneric);
    // an edge of another solid, no edge, a face not of the edge
    const auto other = box(m, 1., 1., 1., {5., 0., 0.});
    ASSERT_EQ(code(fillet_geometry(m, s, edge_near(m, other, {6., 0.5, 1.}), 0.1)), BuildErrc::InvalidId);
    ASSERT_EQ(code(fillet_geometry(m, s, EdgeId{}, 0.1)), BuildErrc::InvalidId);
    ASSERT_EQ(code(chamfer_geometry(m, s, e, face_toward(m, s, edge_near(m, s, {0., 0.75, 0.}), -1. * Z), 0.1, 0.1)), BuildErrc::InvalidId);

    // beyond the radius of the top edge of a cylinder
    const auto c = cylinder(m, 0.5, 1., {0., 5., 0.});
    ASSERT_EQ(code(fillet_geometry(m, c, edge_near(m, c, {0.5 * std::cos(1.), 5. + 0.5 * std::sin(1.), 1.}), 0.6)), BuildErrc::NonGeneric);
    // its seam: used twice by the side face
    ASSERT_EQ(code(fillet_geometry(m, c, edge_near(m, c, {0.5, 5., 0.5}), 0.1)), BuildErrc::InvalidId);

    // the edges of a slab on a stadium between its planes and its half cylinders: tangent faces
    {
        const T L = 2., R = 0.5;
        const auto v = std::array{unwrap(make_vertex(m, point<T, 3>{0., -R, 0.})), unwrap(make_vertex(m, point<T, 3>{L, -R, 0.})),
                                  unwrap(make_vertex(m, point<T, 3>{L, R, 0.})), unwrap(make_vertex(m, point<T, 3>{0., R, 0.}))};
        const std::vector<EdgeId> w{unwrap(make_edge(m, v[0], v[1])), unwrap(make_edge(m, arc(R, {L, 0., 0.}, X, Y, -pi / 2, pi / 2), -pi / 2, pi / 2, v[1], v[2])),
                                    unwrap(make_edge(m, v[2], v[3])), unwrap(make_edge(m, arc(R, {0., 0., 0.}, X, Y, pi / 2, 3 * pi / 2), pi / 2, 3 * pi / 2, v[3], v[0]))};
        const auto slab = solid_of(m, make_prism(m, unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                                                         unwrap(make_wire(m, w)))),
                                                 point<T, 3>{0., 0., 0.4}));
        const auto r = fillet_geometry(m, slab, edge_near(m, slab, {L, -R, 0.2}), 0.1);
        ASSERT_EQ(code(r), BuildErrc::NonGeneric);
        ASSERT_NE(r.error().message.find("tangent"), std::string::npos);
    }
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
                // it ends there, not continued past the vertex; the surface reduced to a point
                ASSERT_TRUE(t == b.a || t == b.b);
                const auto q1 = b.contacts[0].on_blend->value(t), q2 = b.contacts[1].on_blend->value(t);
                for (int j = 0; j <= 8; ++j)
                {
                    const auto q = q1 + (j / 8.) * (q2 - q1);
                    ASSERT_LE(norm(b.surface->value(q[0], q[1]) - b.contacts[0].curve->value(t)), 1e-6);
                }
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

// ---- topology (PR 6): fillet and chamfer of edges of a solid

namespace
{
    T volume(const Model<T> &m, SolidId s)
    {
        return signed_volume(m, m.solid(s).outer, 512);
    }

    std::size_t n_faces(const Model<T> &m, SolidId s)
    {
        return explore<FaceId>(m, ShapeId{s}).size();
    }

    /// A valid solid whose faces are those given, of the volume given.
    void expect_solid(const Model<T> &m, SolidId s, std::size_t faces, T vol, T vol_tol)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
        ASSERT_EQ(n_faces(m, s), faces);
        INFO(std::setprecision(15) << "volume " << volume(m, s) << " expected " << vol);
        ASSERT_NEAR(volume(m, s), vol, vol_tol);
    }

    /// The area of the corner left by a fillet of radius r in a right angle, and its centroid from the corner (along each side).
    std::pair<T, T> corner(T r)
    {
        const T a = (1. - pi / 4.) * r * r;
        const T x = (r * r * r / 2. - (pi * r * r / 4.) * (r - 4. * r / (3. * pi))) / a;
        return {a, x};
    }

    /// A prism of height h on the oblong of two half circles of radius R, centres (0, 0) and (L, 0).
    SolidId oblong(Model<T> &m, T L, T R, T h)
    {
        std::vector<EdgeId> e{unwrap(make_edge(m, point<T, 3>{0., -R, 0.}, point<T, 3>{L, -R, 0.})),
                              unwrap(make_edge(m, arc(R, {L, 0., 0.}, X, Y, -pi / 2., pi / 2.))),
                              unwrap(make_edge(m, point<T, 3>{L, R, 0.}, point<T, 3>{0., R, 0.})),
                              unwrap(make_edge(m, arc(R, {0., 0., 0.}, X, Y, pi / 2., 3. * pi / 2.)))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }
}

TEST(tests_brep_fillet, fillet_box_edges)
{
    // an edge of a box: a quarter cylinder, the corner (1 - pi / 4) r^2 off along it; the third
    // faces at its ends cut by quarter circles
    for (const T r : {0.1, 0.3})
    {
        INFO("r = " << r);
        Model<T> m;
        const auto s = box(m, 2., 1., 1.);
        const auto e = edge_near(m, s, {1., 0., 1.});
        const auto n0 = n_faces(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{e}, r));
        ASSERT_NE(f, s);
        expect_solid(m, f, 7, 2. - corner(r).first * 2., 1e-9);
        // the solid given unchanged
        expect_solid(m, s, n0, 2., 1e-12);
    }
    // two parallel edges of the top face, and the four vertical edges (independent)
    {
        Model<T> m;
        const auto s = box(m, 2., 1., 1.);
        const T r = 0.2;
        const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {1., 0., 1.}), edge_near(m, s, {1., 1., 1.})}, r));
        expect_solid(m, f, 8, 2. - corner(r).first * 4., 1e-9);
        std::vector<EdgeId> vertical;
        for (const auto &p : {point<T, 3>{0., 0., 0.5}, point<T, 3>{2., 0., 0.5}, point<T, 3>{2., 1., 0.5}, point<T, 3>{0., 1., 0.5}})
            vertical.push_back(edge_near(m, s, p));
        const auto g = unwrap(fillet(m, s, vertical, r));
        expect_solid(m, g, 10, 2. - corner(r).first * 4., 1e-9);
    }
}

TEST(tests_brep_fillet, fillet_concave_edge)
{
    // the inner edge of an L: material added, the end faces grow by the corner
    Model<T> m;
    const auto s = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, 1.);
    const T r = 0.25;
    const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {1., 1., 0.5})}, r));
    expect_solid(m, f, 9, 3. + corner(r).first, 1e-9);
}

TEST(tests_brep_fillet, fillet_closed_edges)
{
    // the rims of a cylinder: tori, the volume by Pappus (the corner turning at R - x)
    const T R = 1., h = 2., r = 0.2;
    const auto [a, x] = corner(r);
    Model<T> m;
    const auto s = cylinder(m, R, h);
    // (points away from the seam, at x = R)
    const auto top = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), h});
    const auto f = unwrap(fillet(m, s, std::vector{top}, r));
    // (the volume of the cylinder itself is integrated within 4e-7)
    expect_solid(m, f, 4, pi * R * R * h - 2. * pi * (R - x) * a, 1e-6);
    const auto g = unwrap(fillet(m, s, std::vector{top, edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), 0.})}, r));
    expect_solid(m, g, 5, pi * R * R * h - 4. * pi * (R - x) * a, 1e-6);
}

TEST(tests_brep_fillet, fillet_chain)
{
    // the top rim of an oblong prism: four edges tangent at their vertices, one chain round it;
    // the contacts on the top face a loop, those on the sides cut the edges between them
    const T L = 2., R = 0.5, h = 1., r = 0.1;
    const auto [a, x] = corner(r);
    Model<T> m;
    const auto s = oblong(m, L, R, h);
    std::vector<EdgeId> rim;
    for (const auto &p : {point<T, 3>{1., -R, h}, point<T, 3>{L + R, 0., h}, point<T, 3>{1., R, h}, point<T, 3>{-R, 0., h}})
        rim.push_back(edge_near(m, s, p));
    const T v0 = (pi * R * R + 2. * R * L) * h;
    const auto f = unwrap(fillet(m, s, rim, r));
    expect_solid(m, f, 6 + 4, v0 - a * (2. * L + 2. * pi * (R - x)), 1e-6);
    // part of the chain: refused, as no face ends it (the next edge goes on tangentially)
    const auto g = fillet(m, s, std::vector{rim[0], rim[1]}, r);
    REQUIRE_FALSE(g.has_value());
    ASSERT_EQ(g.error().code, BuildErrc::NonGeneric);
    ASSERT_NE(g.error().message.find("chain"), std::string::npos);
}

TEST(tests_brep_fillet, chamfer_edges)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto f = unwrap(chamfer(m, s, std::vector{e}, 0.2));
    expect_solid(m, f, 7, 2. - 0.2 * 0.2 / 2. * 2., 1e-9);
    const auto top = face_toward(m, s, e, Z);
    const auto g = unwrap(chamfer(m, s, std::vector{e}, top, 0.1, 0.3));
    expect_solid(m, g, 7, 2. - 0.1 * 0.3 / 2. * 2., 1e-9);
}

TEST(tests_brep_fillet, blend_refusals)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    auto alive = [&] { return std::array{m.ids<VertexId>().size(), m.ids<EdgeId>().size(), m.ids<WireId>().size(), m.ids<FaceId>().size(),
                                         m.ids<ShellId>().size(), m.ids<SolidId>().size()}; };
    const auto size = alive();
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    // two edges meeting at a corner, not tangentially: not handled
    const auto r = fillet(m, s, std::vector{edge_near(m, s, {1., 0., 1.}), edge_near(m, s, {2., 0.5, 1.})}, 0.2);
    ASSERT_EQ(code(r), BuildErrc::NonGeneric);
    ASSERT_EQ(alive(), size);
    ASSERT_EQ(code(fillet(m, s, std::vector<EdgeId>{}, 0.2)), BuildErrc::EmptyInput);
    const auto e = edge_near(m, s, {1., 0., 1.});
    ASSERT_EQ(code(fillet(m, s, std::vector{e, e}, 0.2)), BuildErrc::InvalidId);
    // wider than the faces next to it
    ASSERT_EQ(code(fillet(m, s, std::vector{e}, 1.2)), BuildErrc::NonGeneric);
    // the model as it was (but tombstones)
    ASSERT_EQ(alive(), size);
}

TEST(tests_brep_fillet, fillet_boss_and_corners)
{
    // a cylinder standing on a plate: its foot concave (a torus adding the corner, outside the
    // cylinder), its top rim convex, both closed
    {
        const T R = 0.5, r = 0.15;
        const auto [a, x] = corner(r);
        Model<T> m;
        const auto s = fused(m, box(m, 4., 4., 1., {-2., -2., -1.}), cylinder(m, R, 1.));
        const auto foot = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), 0.});
        const auto top = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), 1.});
        const T v0 = volume(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{foot, top}, r));
        expect_solid(m, f, n_faces(m, s) + 2, v0 + 2. * pi * (R + x) * a - 2. * pi * (R - x) * a, 1e-6);
    }
    // the top edges of a triangular prism: ends cut by the other sides at 60 degrees; the
    // corner swept between oblique planes, its volume along its centroid line
    {
        const T side = 2., h = 1., r = 0.2;
        const auto [a, x] = corner(r);
        Model<T> m;
        const auto s = prism(m, {{0., 0., 0.}, {side, 0., 0.}, {side / 2., side * std::sqrt(3.) / 2., 0.}}, h);
        const T v0 = volume(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {side / 2., 0., h})}, r));
        expect_solid(m, f, 6, v0 - a * (side - 2. * x / std::sqrt(3.)), 1e-9);
    }
}

TEST(tests_brep_fillet, chamfer_closed_and_concave)
{
    // the rim of a cylinder: a cone, the triangle d^2 / 2 turning at R - d / 3
    {
        const T R = 1., h = 2., d = 0.2;
        Model<T> m;
        const auto s = cylinder(m, R, h);
        const T v0 = volume(m, s);
        const auto top = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), h}), foot = edge_near(m, s, {R * std::cos(1.), R * std::sin(1.), 0.});
        const auto f = unwrap(chamfer(m, s, std::vector{top}, d));
        expect_solid(m, f, 4, v0 - 2. * pi * (R - d / 3.) * d * d / 2., 1e-6);
        const auto g = unwrap(chamfer(m, s, std::vector{top, foot}, d));
        expect_solid(m, g, 5, v0 - 4. * pi * (R - d / 3.) * d * d / 2., 1e-6);
    }
    // the inner edge of an L: material added
    {
        Model<T> m;
        const auto s = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, 1.);
        const auto e = edge_near(m, s, {1., 1., 0.5});
        const auto f = unwrap(chamfer(m, s, std::vector{e}, face_toward(m, s, e, Y), 0.2, 0.3));
        expect_solid(m, f, 9, 3. + 0.2 * 0.3 / 2., 1e-9);
    }
}

TEST(tests_brep_fillet, fillet_approximated)
{
    FilletOptions<T> o;
    // the saddle where a cylinder stands on another: a closed edge, a B-spline fillet round its
    // loop; the chamfer likewise
    {
        Model<T> m;
        const auto s = fused(m, cylinder_along(m, 1., 4., {-2., 0., 0.}, Y, Z, X), cylinder(m, 0.4, 1.5, {0., 0., 0.}));
        const auto e = edge_near(m, s, {0.4, 0., std::sqrt(1. - 0.16)});
        const T v0 = volume(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{e}, 0.15, o));
        expect_solid(m, f, n_faces(m, s) + 1, v0, 0.05);
        ASSERT_GT(volume(m, f), v0); // concave: material added
        const auto g = unwrap(chamfer(m, s, std::vector{e}, 0.1, o));
        expect_solid(m, g, n_faces(m, s) + 1, v0, 0.05);
        ASSERT_GT(volume(m, g), v0);
    }
    // a tee whose cylinders have nearly the same radius
    {
        Model<T> m;
        const auto [s, e] = tee(m, 0.99);
        const T v0 = volume(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{e}, 0.05, o));
        expect_solid(m, f, n_faces(m, s) + 1, v0, 0.01);
        ASSERT_GT(volume(m, f), v0);
    }
    // the top edge of a B-spline bump prism: free ends on the planes at x = 0 and x = 2 (the
    // bump leaving them square: an edge at a grazing angle to the face at its end needs the
    // blend continued far past it)
    {
        Model<T> m;
        const auto bump = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {0.4, 0., 0.}, {0.8, -0.3, 0.}, {1.2, 0.2, 0.}, {1.6, 0., 0.}, {2., 0., 0.}},
                                                          std::vector<T>{0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.}, std::size_t{3});
        const auto v0 = unwrap(make_vertex(m, point<T, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<T, 3>{2., 0., 0.}));
        const auto v2 = unwrap(make_vertex(m, point<T, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<T, 3>{0., 1., 0.}));
        const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                    unwrap(make_edge(m, v3, v0))};
        const auto pf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                         unwrap(make_wire(m, w))));
        const auto s = solid_of(m, make_prism(m, pf, point<T, 3>{0., 0., 0.8}));
        const auto e = edge_near(m, s, {bump->value(0.5)[0], bump->value(0.5)[1], 0.8});
        const T vol = volume(m, s);
        const auto f = unwrap(fillet(m, s, std::vector{e}, 0.1, o));
        // the corner, about (1 - pi / 4) r^2 along the bump (a little over 2 long)
        expect_solid(m, f, 7, vol - corner(0.1).first * 2.1, 2e-3);
        const auto top = face_toward(m, s, e, Z);
        const auto g = unwrap(chamfer(m, s, std::vector{e}, top, 0.08, 0.05, o));
        expect_solid(m, g, 7, vol - 0.08 * 0.05 / 2. * 2.1, 2e-3);
    }
}

TEST(tests_brep_fillet, fillet_vanishing)
{
    // a tee of two cylinders of the same radius: they meet along two half ellipses crossing at
    // the points where the cylinders are tangent; the fillet of one half ellipse (two edges, cut
    // at the seam of the standing cylinder) vanishes at both ends
    Model<T> m;
    const auto [s, e0] = tee(m, 1.);
    std::vector<EdgeId> half;
    for (auto e : explore<EdgeId>(m, ShapeId{s}))
    {
        const auto &ed = m.edge(e);
        if (ed.degenerate)
            continue;
        const auto p = ed.curve->value((ed.u1 + ed.u2) / 2.);
        if (std::abs(std::hypot(p[0], p[1]) - 1.) < 1e-6 && std::abs(std::hypot(p[1], p[2]) - 1.) < 1e-6 && p[0] > 0.)
            half.push_back(e);
    }
    REQUIRE(!half.empty());
    const T v0 = volume(m, s);
    const auto f = fillet(m, s, half, 0.15);
    INFO((f ? std::string{} : f.error().message));
    REQUIRE(f.has_value());
    expect_solid(m, *f, n_faces(m, s) + half.size(), v0, 0.05);
    ASSERT_GT(volume(m, *f), v0);
}


TEST(tests_brep_fillet, blend_refusals_more)
{
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    // a free end at a vertex of four edges: the top of a box cut along a diagonal
    {
        Model<T> m;
        const auto s = box(m, 2., 1., 1.);
        const auto top = face_toward(m, s, edge_near(m, s, {1., 0., 1.}), Z);
        VertexId a{}, b{};
        for (auto v : explore<VertexId>(m, ShapeId{s}))
        {
            if (norm(m.vertex(v).pnt - point<T, 3>{2., 0., 1.}) < 1e-12)
                a = v;
            if (norm(m.vertex(v).pnt - point<T, 3>{0., 1., 1.}) < 1e-12)
                b = v;
        }
        const auto d = unwrap(make_edge(m, a, b));
        const auto pc = plane_pcurve(m.edge(d), *m.face(top).surface, 1e-9);
        REQUIRE(pc.has_value());
        unwrap(split_face(m, top, std::vector<FaceCut<T>>{{d, pc->pcurve}}));
        REQUIRE(check(m, ShapeId{s}).ok());
        const auto r = fillet(m, s, std::vector{edge_near(m, s, {2., 0., 0.5})}, 0.1);
        ASSERT_EQ(code(r), BuildErrc::NonGeneric);
        ASSERT_NE(r.error().message.find("more than three edges"), std::string::npos);
    }
    // two fillets whose contacts meet on the face between them (a full round): refused as such
    {
        Model<T> m;
        const auto s = box(m, 2., 0.3, 1.);
        const auto r = fillet(m, s, std::vector{edge_near(m, s, {1., 0., 1.}), edge_near(m, s, {1., 0.3, 1.})}, 0.15);
        ASSERT_EQ(code(r), BuildErrc::NonGeneric);
        ASSERT_NE(r.error().message.find("full round"), std::string::npos);
    }
    // an edge of another solid, a dead face, a face that is not one of every edge
    {
        Model<T> m;
        const auto s = box(m, 2., 1., 1.), other = box(m, 1., 1., 1., {5., 0., 0.});
        const auto e = edge_near(m, s, {1., 0., 1.});
        ASSERT_EQ(code(fillet(m, s, std::vector{edge_near(m, other, {5.5, 0., 1.})}, 0.1)), BuildErrc::InvalidId);
        ASSERT_EQ(code(chamfer(m, s, std::vector{e}, FaceId{}, 0.1, 0.2)), BuildErrc::InvalidId);
        ASSERT_EQ(code(chamfer(m, s, std::vector{e, edge_near(m, s, {1., 1., 0.})}, face_toward(m, s, e, Z), 0.1, 0.2)), BuildErrc::InvalidId);
    }
}

TEST(tests_brep_fillet, fillet_solid_with_a_void)
{
    // a hollow box (make_thick, no face removed): an edge of its outer shell (convex, material
    // taken off) and one of its void (concave seen from the solid, material added)
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto h = unwrap(make_thick(m, s, std::vector<FaceId>{}, 0.1));
    REQUIRE(m.solid(h).voids.size() == 1u);
    auto total = [&](SolidId x) {
        T v = signed_volume(m, m.solid(x).outer, 512);
        for (auto sh : m.solid(x).voids)
            v += signed_volume(m, sh, 512);
        return v;
    };
    const T v0 = 2. - 1.8 * 0.8 * 0.8, r = 0.05, a = corner(r).first;
    ASSERT_NEAR(total(h), v0, 1e-9);
    const auto f = unwrap(fillet(m, h, std::vector{edge_near(m, h, {1., 0., 1.})}, r));
    expect_solid(m, f, 13, 2. - a * 2., 1e-9);
    ASSERT_NEAR(total(f), v0 - a * 2., 1e-9);
    const auto g = unwrap(fillet(m, h, std::vector{edge_near(m, h, {1., 0.1, 0.9})}, r));
    ASSERT_TRUE(check(m, ShapeId{g}).ok());
    ASSERT_NEAR(total(g), v0 + a * 1.8, 1e-9);
}

TEST(tests_brep_fillet, fillet_box_corner)
{
    // the three edges at a corner of a box: three quarter cylinders and the eighth of a sphere
    // at the corner; in the cube [1 - r, 1]^3 of the corner only the eighth of the ball is kept,
    // along each edge beyond it (1 - pi / 4) r^2
    for (const T r : {0.1, 0.3})
    {
        INFO("r = " << r);
        Model<T> m;
        const auto s = box(m, 1., 1., 1.);
        const std::vector<EdgeId> edges{edge_near(m, s, {0.5, 1., 1.}), edge_near(m, s, {1., 0.5, 1.}), edge_near(m, s, {1., 1., 0.5})};
        const auto f = unwrap(fillet(m, s, edges, r));
        const T removed = 3. * corner(r).first * (1. - r) + r * r * r * (1. - pi / 6.);
        expect_solid(m, f, 10, 1. - removed, 1e-9);
        // the sphere: of radius r, centred at (1 - r, 1 - r, 1 - r)
        const point<T, 3> c{1. - r, 1. - r, 1. - r};
        int spheres = 0;
        for (auto fc : explore<FaceId>(m, ShapeId{f}))
            if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(m.face(fc).surface.get()))
            {
                ++spheres;
                ASSERT_NEAR(sp->radius(), r, 1e-12);
                ASSERT_LE(norm(sp->frame().o - c), 1e-12);
            }
        ASSERT_EQ(spheres, 1);
    }
}

TEST(tests_brep_fillet, fillet_box_all_edges)
{
    // the twelve edges of a box: the rounded box, the core, slabs, quarter cylinders along the
    // edges and the eight eighths of a ball
    const T a = 2., b = 1.5, c = 1., r = 0.2;
    Model<T> m;
    const auto s = box(m, a, b, c);
    const auto all = explore<EdgeId>(m, ShapeId{s});
    ASSERT_EQ(all.size(), 12u);
    const auto f = unwrap(fillet(m, s, all, r));
    const T A = a - 2. * r, B = b - 2. * r, C = c - 2. * r;
    const T v = A * B * C + 2. * r * (A * B + B * C + A * C) + pi * r * r * (A + B + C) + 4. / 3. * pi * r * r * r;
    // (the volume of the spheres integrated within about 1e-9)
    expect_solid(m, f, 26, v, 1e-8);
}

TEST(tests_brep_fillet, fillet_concave_corner)
{
    // a pocket in a block: the three edges at its inner corner are concave; the ball rolled in
    // the corner adds the eighth of a cube less the eighth of the ball
    const T r = 0.2;
    Model<T> m;
    const auto block = box(m, 2., 2., 2.);
    const auto pocket = box(m, 1., 1., 1., {1., 1., 1.});
    const auto s = unwrap(cut(m, ShapeId{block}, ShapeId{pocket})).solids.front();
    const std::vector<EdgeId> edges{edge_near(m, s, {1.5, 1., 1.}), edge_near(m, s, {1., 1.5, 1.}), edge_near(m, s, {1., 1., 1.5})};
    const auto n0 = n_faces(m, s);
    const auto f = unwrap(fillet(m, s, edges, r));
    const T added = 3. * corner(r).first * (1. - r) + r * r * r * (1. - pi / 6.);
    expect_solid(m, f, n0 + 4, 7. + added, 1e-9);
}

TEST(tests_brep_fillet, corner_refusals)
{
    Model<T> m;
    const auto s = box(m, 1., 1., 1.);
    const std::vector<EdgeId> edges{edge_near(m, s, {0.5, 1., 1.}), edge_near(m, s, {1., 0.5, 1.}), edge_near(m, s, {1., 1., 0.5})};
    const auto n = m.ids<FaceId>().size();
    // chamfers at a corner: not handled
    const auto c = chamfer(m, s, edges, 0.2);
    ASSERT_FALSE(c.has_value());
    ASSERT_NE(c.error().message.find("three that are not fillets"), std::string::npos);
    ASSERT_EQ(m.ids<FaceId>().size(), n);
    // the top of the inner corner of an L: a concave edge and two convex ones
    const auto l = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}, 1.);
    const auto f = fillet(m, l, std::vector{edge_near(m, l, {1., 1., 0.5}), edge_near(m, l, {1.5, 1., 1.}), edge_near(m, l, {1., 1.5, 1.})}, 0.2);
    ASSERT_FALSE(f.has_value());
    ASSERT_NE(f.error().message.find("convex and concave"), std::string::npos);
}

TEST(tests_brep_fillet, fillet_blade_root)
{
    // a blade fused to its hub, its root filleted: the skinned blade closes across a seam without
    // being periodic, and its sections join pieces with C2 knots (multiplicity 7, degree 7)
    // where the fillet is no smoother; the loop of the fillet starts on the seam
    Model<T> m;
    const auto [s, root] = blade_fixture::blade_on_hub(m);
    ASSERT_EQ(root.size(), 1u);
    ASSERT_EQ(m.edge(root[0]).v1, m.edge(root[0]).v2);
    const T r = 0.03;
    FilletOptions<T> o;
    const auto b = unwrap(fillet_geometry(m, s, root[0], r, o));
    ASSERT_TRUE(b.closed);
    ASSERT_FALSE(b.convex);
    ASSERT_LE(b.deviation, o.approx_tol);
    expect_blend(m, s, b, r, o.approx_tol, o.angular_tol);
    const auto n0 = n_faces(m, s);
    const auto f = unwrap(fillet(m, s, root, r, o));
    // the corner (1 - pi / 4) r^2 added along the root, about (the hub turns away under it: 5 %
    // less); the volume against OCCT in tests_brep2_occt
    T length{};
    const auto &ed = m.edge(root[0]);
    for (int i = 0; i < 512; ++i)
        length += norm(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * (i + 1) / 512.) - ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * i / 512.));
    const T v0 = volume(m, s);
    expect_solid(m, f, n0 + 1, v0 + corner(r).first * length, 0.1 * corner(r).first * length);
    // the solid given unchanged
    expect_solid(m, s, n0, v0, 1e-12);
}

namespace
{
    /// The solid of revolution of radius r0 at z = 0 to r1 at z = 1, its side ruled between two
    /// rational B-spline circles (a closed B-spline, its seam at x > 0 below), its top circle
    /// turned by `twist` (the seam then oblique), closed by two planes.
    SolidId frustum(Model<T> &m, T r0, T r1, T twist = 0.)
    {
        const auto c0 = build_circle<T, 3>(r0, {0., 0., 0.});
        auto c1 = build_circle<T, 3>(r1, {0., 0., 1.});
        auto turned = c1.poles();
        for (auto &q : turned)
            q = {std::cos(twist) * q[0] - std::sin(twist) * q[1], std::sin(twist) * q[0] + std::cos(twist) * q[1], q[2], q[3]};
        c1 = BSCurveRational<T, 3>(turned, c1.knotsFlats(), c1.degree());
        const std::size_t nu = c0.poles().size();
        std::vector<std::array<T, 4>> poles(2 * nu);
        for (std::size_t j = 0; j < nu; ++j)
            poles[j] = c0.poles()[j], poles[j + nu] = c1.poles()[j];
        const auto side = std::make_shared<BSSurfaceRational<T, 3>>(poles, c0.knotsFlats(), std::vector<T>{0., 0., 1., 1.}, std::size_t{2}, std::size_t{1});
        std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{side}))};
        for (const auto &[c, z] : {std::pair{c0, 0.}, std::pair{c1, 1.}})
        {
            const auto e = unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<BSCurveRational<T, 3>>(c)}));
            const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., z}, Z, X}}, std::array<T, 4>{-3., 3., -3., 3.});
            faces.push_back(unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{plane}, unwrap(make_wire(m, std::vector<EdgeId>{e})))));
        }
        return unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()));
    }

    /// A cylinder of radius r on z from z0 to z0 + h whose side is a closed rational B-spline, its seam at the angle `turn`.
    SolidId nurbs_cylinder(Model<T> &m, T r, T h, T z0, T turn)
    {
        const ax2<T, 3> ax{{{0., 0., z0}, Z, X}};
        const auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(r, turn, turn + 2. * pi, ax));
        const auto e = unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{circle}));
        const auto plane = std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2. * r, 2. * r, -2. * r, 2. * r});
        const auto disc = unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{plane}, unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, disc, point<T, 3>{0., 0., h}));
    }
}

TEST(tests_brep_fillet, fillet_offset_of_cylinder)
{
    // a cylinder whose side is the offset of a cylinder (radius 0.8 + 0.2, its bounds a whole
    // turn: it closes, but repeats past its seam): its top edge filleted on a torus, exact, as
    // the offset of a primitive gives (a seamless view of it gave a rational B-spline)
    const T r = 0.1;
    Model<T> m;
    const auto base = std::make_shared<CylindricalSurface<T>>(0.8, ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{0., 2. * pi, 0., 1.});
    std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{std::make_shared<OffsetSurface<T>>(std::shared_ptr<Surface<T, 3>>{base}, 0.2)}))};
    for (const T z : {0., 1.})
    {
        const ax2<T, 3> ax{{{0., 0., z}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., ax)));
        faces.push_back(unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e})))));
    }
    const auto s = unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()));
    const auto top = edge_near(m, s, {-1., 0., 1.}); // (the seam of the side at x = 1)
    const auto b = unwrap(fillet_geometry(m, s, top, r, FilletOptions<T>{}));
    ASSERT_TRUE(b.exact);
    ASSERT_NE(dynamic_cast<const ToroidalSurface<T> *>(b.surface.get()), nullptr);
    const auto [a, x] = corner(r);
    const auto f = unwrap(fillet(m, s, std::vector{top}, r));
    expect_solid(m, f, 4, pi - 2. * pi * (1. - x) * a, 1e-6);
}

TEST(tests_brep_fillet, fillet_two_closed_bsplines)
{
    // a frustum (r 1.3 to 1) under a cylinder of radius 1, both sides closed rational B-splines:
    // the concave edge between them filleted, the fillet round both faces seen past their seams
    const T r = 0.1, dev = std::atan(0.3); // the turn of the generatrix at the edge
    auto ring = [&](const Model<T> &m, SolidId s) {
        std::vector<EdgeId> es;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate || !ed.curve)
                continue;
            const auto p = ed.curve->value((ed.u1 + ed.u2) / 2.);
            if (std::abs(p[2] - 1.) < 1e-6 && std::abs(std::hypot(p[0], p[1]) - 1.) < 1e-6)
                es.push_back(e);
        }
        return es;
    };
    auto solid = [&](Model<T> &m, T turn, T twist) {
        return unwrap(fuse(m, ShapeId{frustum(m, 1.3, 1., twist)}, ShapeId{nurbs_cylinder(m, 1., 1., 1., turn)})).solids.front();
    };
    // the seams in the same plane: one closed edge, the loop of the fillet from the seams
    T aligned{};
    {
        Model<T> m;
        const auto s = solid(m, 0., 0.);
        const auto e = ring(m, s);
        ASSERT_EQ(e.size(), 1u);
        ASSERT_EQ(m.edge(e[0]).v1, m.edge(e[0]).v2);
        const T v0 = volume(m, s);
        const auto f = unwrap(fillet(m, s, e, r));
        // the corner r^2 (tan(dev / 2) - dev / 2) added round the edge, its centroid about 1 from the axis
        expect_solid(m, f, 5, v0 + 2. * pi * r * r * (std::tan(dev / 2.) - dev / 2.), 1e-6);
        aligned = volume(m, f);
    }
    // the seams apart: the edge split at both, a chain of two edges, the same fillet
    {
        Model<T> m;
        const auto s = solid(m, 1., 0.);
        const auto e = ring(m, s);
        ASSERT_EQ(e.size(), 2u);
        const auto f = unwrap(fillet(m, s, e, r));
        expect_solid(m, f, 6, aligned, 1e-8);
    }
    // the seam of the frustum oblique: the seams meet at the vertex of the closed edge but not in
    // the same section of the fillet, the marching stops at each (not handled)
    {
        Model<T> m;
        const auto s = solid(m, 0.3, 0.3);
        const auto e = ring(m, s);
        ASSERT_EQ(e.size(), 1u);
        const auto f = fillet(m, s, e, r);
        ASSERT_FALSE(f.has_value());
        ASSERT_NE(f.error().message.find("seams do not meet in the same section"), std::string::npos);
    }
}

TEST(tests_brep_fillet, fillet_sharp_corner)
{
    // the apex of a tetrahedron whose faces are 150 degrees apart round it, its three edges
    // filleted: the patch of the corner nearly fills a hemisphere, every pole of the sphere is
    // a few degrees from its arcs; their pcurves on the sphere stay within the edges' tolerance
    // between the points interpolated (33 of them left 3e-5 for an edge at 1e-6), clear of the
    // seam and of the poles
    const T zn = 0.05, fan = 150. * pi / 180., L = 3., r = 0.02;
    const auto a = unit_of({1., 0., zn}), b = unit_of({std::cos(fan), std::sin(fan), zn}), c = unit_of({std::cos(fan), -std::sin(fan), zn});
    ASSERT_LT(a * (a + b + c), 0.); // a foot more than 90 degrees from the sum of the feet
    auto edge_dir = [](const point<T, 3> &p, const point<T, 3> &q, const point<T, 3> &o) {
        const auto e = unit_of(cross(p, q));
        return e * o < 0. ? e : T(-1) * e;
    };
    const point<T, 3> O{};
    const auto Pab = L * edge_dir(a, b, c), Pbc = L * edge_dir(b, c, a), Pca = L * edge_dir(c, a, b);
    Model<T> m;
    auto face = [&](const std::vector<point<T, 3>> &p) { return polygon(m, p, unit_of(cross(p[1] - p[0], p[2] - p[0])), unit_of(p[1] - p[0])); };
    const std::vector<FaceId> faces{face({O, Pab, Pca}), face({O, Pbc, Pab}), face({O, Pca, Pbc}), face({Pab, Pbc, Pca})};
    const auto s = unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()));
    const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, Pab / 2.), edge_near(m, s, Pbc / 2.), edge_near(m, s, Pca / 2.)}, r));
    const T v0 = volume(m, s);
    expect_solid(m, f, 8, v0, 0.05 * v0);
    ASSERT_LT(volume(m, f), v0);
    int arcs = 0;
    for (auto fc : explore<FaceId>(m, ShapeId{f}))
        if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(m.face(fc).surface.get()))
            for (auto w : m.face(fc).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    ++arcs;
                    const auto &ed = m.edge(ce.edge);
                    T gap{};
                    for (int i = 0; i <= 1024; ++i)
                    {
                        const T t = ed.u1 + (ed.u2 - ed.u1) * T(i) / 1024.;
                        const auto q = ce.pcurve->value(t);
                        ASSERT_LT(std::abs(q[0]), pi);
                        ASSERT_LT(std::abs(q[1]), pi / 2.);
                        gap = std::max(gap, norm(sp->value(q[0], q[1]) - ed.curve->value(t)));
                    }
                    ASSERT_LE(gap, ed.tol);
                }
    ASSERT_EQ(arcs, 3);
}
