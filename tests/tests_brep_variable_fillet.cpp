// Fillets whose radius follows a law along their edge (stage 5, PR 1): the law of the radius
// (radius_law.h) and the geometry of the fillet (fillet_geometry with a RadiusLaw).
#include <doctest_gtest.hpp>
#include <gbs-brep/fillet.h>
#include <gbs-brep/sweep.h>
#include <gbs/bscinterp.h>

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <numbers>
#include <vector>

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

    FaceId polygon(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    }

    SolidId box(Model<T> &m, T a, T b, T c)
    {
        return solid_of(m, make_prism(m, polygon(m, {{0., 0., 0.}, {a, 0., 0.}, {a, b, 0.}, {0., b, 0.}}), point<T, 3>{0., 0., c}));
    }

    SolidId cylinder(Model<T> &m, T r, T h)
    {
        const ax2<T, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

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

    point<T, 3> outward(const Model<T> &m, SolidId s, FaceId f, const Surface<T, 3> &support, const point<T, 2> &uv)
    {
        auto n = normal_direction(support, uv[0], uv[1]);
        for (auto sh : explore<ShellId>(m, ShapeId{s}))
            for (const auto &fu : m.shell(sh).faces)
                if (fu.face == f && fu.orient == Orientation::Reversed)
                    n = T(-1) * n;
        return n;
    }

    /// The radius of the law at the parameter t of the edge.
    std::function<T(T)> radius_of(const Model<T> &m, EdgeId e, const RadiusLaw<T> &law)
    {
        const auto R = fillet_detail::edge_radius(m, e, law, 1e-6, "test");
        REQUIRE(R.has_value());
        return R->r;
    }

    /// The surface of the blend at t and the contacts on it: the blend's, or that of its piece holding t.
    std::pair<const Surface<T, 3> *, std::array<const Curve<T, 2> *, 2>> surface_at(const EdgeBlend<T> &b, T t)
    {
        if (b.surface)
            return {b.surface.get(), {b.contacts[0].on_blend.get(), b.contacts[1].on_blend.get()}};
        for (const auto &p : b.pieces)
            if (t >= p.a && t <= p.b)
                return {p.surface.get(), {p.on_blend[0].get(), p.on_blend[1].get()}};
        return {nullptr, {nullptr, nullptr}};
    }

    /**
     * A fillet of varying radius r(t): along the edge, its centre in the plane normal to the edge,
     * its contacts at r from it on the faces, the surface through them tangent to the faces;
     * across, the surface on the ball of the station (at r from the centre) and on its
     * characteristic circle, in the plane (q - c) . c' = -r r' (c' and r' by differences, away
     * from the kinks of the law). At a kink, the lune of the ball between the pieces: centred on
     * the centre line, of the radius there, between the contacts, tangent to the pieces.
     */
    void expect_variable(const Model<T> &m, SolidId s, const EdgeBlend<T> &b, const std::function<T(T)> &r_of, T gap, T angle)
    {
        ASSERT_TRUE(b.surface || !b.pieces.empty());
        ASSERT_LE(b.deviation, gap);
        const auto &ed = m.edge(b.edge);
        const T h = 1e-5 * (b.b - b.a);
        auto near_kink = [&](T t) { return std::ranges::any_of(b.lunes, [&](const auto &l) { return std::abs(l.t - t) < 2. * h; }); };
        for (int i = 0; i <= 64; ++i)
        {
            const T t = b.t1 + (b.t2 - b.t1) * T(i) / 64.;
            const T r = r_of(t);
            const auto c = b.spine->value(t);
            const auto e = ed.curve->value(t), te = ed.curve->value(t, 1) / norm(ed.curve->value(t, 1));
            ASSERT_LE(std::abs((c - e) * te), gap);
            const auto [srf, ob] = surface_at(b, t);
            REQUIRE(srf);
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &ct = b.contacts[k];
                const auto p = ct.curve->value(t);
                const auto q = ob[k]->value(t);
                ASSERT_LE(norm(srf->value(q[0], q[1]) - p), gap);
                ASSERT_NEAR(norm(p - c), r, gap);
                const auto nf = outward(m, s, ct.face, *ct.support, ct.pcurve->value(t));
                ASSERT_LE(norm(cross(nf, normal_direction(*srf, q[0], q[1]))), angle);
                // the centre along the normal of the face, inside a convex edge
                ASSERT_LE(norm(c - (p + (b.convex ? -r : r) * nf)), gap);
            }
            if (near_kink(t))
                continue;
            const T t0 = std::max(b.a, t - h), t1 = std::min(b.b, t + h);
            const auto dc = (b.spine->value(t1) - b.spine->value(t0)) / (t1 - t0);
            const T dr = (r_of(t1) - r_of(t0)) / (t1 - t0);
            const auto q1 = ob[0]->value(t), q2 = ob[1]->value(t);
            for (int j = 1; j < 8; ++j)
            {
                const auto w = q1 + (j / 8.) * (q2 - q1);
                const auto q = srf->value(w[0], w[1]);
                ASSERT_NEAR(norm(q - c), r, gap);
                ASSERT_NEAR((q - c) * dc / norm(dc), -r * dr / norm(dc), 10. * gap);
            }
        }
        for (const auto &l : b.lunes)
        {
            const auto c = b.spine->value(l.t);
            const T r = r_of(l.t);
            ASSERT_LE(norm(l.sphere->frame().o - c), gap);
            ASSERT_NEAR(l.sphere->radius(), r, gap);
            // the pcurves on the sphere: within tol of the arcs between the points interpolated as
            // at them (`settle` and `check` sample at the points)
            ASSERT_LE(l.deviation, 1e-6);
            for (std::size_t j = 0; j < 2; ++j)
                for (int k = 0; k <= 1024; ++k)
                {
                    const auto uv = l.on_sphere[j]->value(k / 1024.);
                    ASSERT_LE(norm(l.sphere->value(uv[0], uv[1]) - l.arcs[j]->value(k / 1024.)), 1e-6);
                }
            T apart{};
            for (std::size_t j = 0; j < 2; ++j)
            {
                ASSERT_LE(norm(l.arcs[j]->value(0.) - b.contacts[0].curve->value(l.t)), gap);
                ASSERT_LE(norm(l.arcs[j]->value(1.) - b.contacts[1].curve->value(l.t)), gap);
                for (int k = 0; k <= 16; ++k)
                {
                    const T w = k / 16.;
                    const auto q = l.arcs[j]->value(w);
                    ASSERT_NEAR(norm(q - c), r, gap);
                    apart = std::max(apart, norm(l.arcs[0]->value(w) - l.arcs[1]->value(w)));
                    // tangent to the piece along its end (G1): the normal of the sphere there
                    const auto [srf, ob] = surface_at(b, j == 0 ? l.t - 1e-12 : l.t + 1e-12);
                    (void)ob;
                    const auto *cs = dynamic_cast<const CurveOnSurface<T, 3> *>(l.arcs[j].get());
                    REQUIRE(cs);
                    const auto pq = cs->p_basisCurve()->value(w);
                    if (k > 0 && k < 16)
                        ASSERT_LE(norm(cross(normal_direction(*cs->p_basisSurface(), pq[0], pq[1]), (q - c) / r)), angle);
                    (void)srf;
                }
            }
            // a lune of some width
            ASSERT_GT(apart, 1e-3);
        }
    }
}

TEST(tests_brep_variable_fillet, radius_law)
{
    using radius_detail::radius_function;
    const std::vector<std::pair<T, T>> pts{{0., 0.1}, {1., 0.3}, {1.5, 0.3}, {3., 0.05}};
    // linear: through the points, linear between them
    const auto lin = radius_function(pts, RadiusVariation::Linear);
    for (const auto &[s, r] : pts)
        ASSERT_NEAR(lin(s), r, 1e-14);
    ASSERT_NEAR(lin(0.5), 0.2, 1e-14);
    ASSERT_NEAR(lin(2.25), 0.175, 1e-14);
    // smooth: through the points, C1, monotone between them (flat between equal radii), its
    // slopes 0 at the ends
    const auto sm = radius_function(pts, RadiusVariation::Smooth);
    for (const auto &[s, r] : pts)
        ASSERT_NEAR(sm(s), r, 1e-14);
    ASSERT_NEAR(sm(0., 1), 0., 1e-12);
    ASSERT_NEAR(sm(3., 1), 0., 1e-12);
    for (std::size_t i = 0; i + 1 < pts.size(); ++i)
    {
        const auto [s0, r0] = pts[i];
        const auto [s1, r1] = pts[i + 1];
        for (int k = 0; k <= 50; ++k)
        {
            const T s = s0 + (s1 - s0) * k / 50.;
            ASSERT_GE(sm(s), std::min(r0, r1) - 1e-14);
            ASSERT_LE(sm(s), std::max(r0, r1) + 1e-14);
            if (k > 0)
                ASSERT_GE((sm(s) - sm(s - (s1 - s0) / 50.)) * (r1 - r0), -1e-15);
        }
        if (i > 0)
            ASSERT_NEAR(sm(s0 - 1e-9, 1), sm(s0 + 1e-9, 1), 1e-6);
    }
    // closed: the slope across the closure
    const auto cl = radius_function<T>({{0., 0.1}, {1., 0.3}, {2., 0.2}}, RadiusVariation::Smooth, true, 3.);
    ASSERT_NEAR(cl(0.), 0.1, 1e-14);
    ASSERT_NEAR(cl(3.), 0.1, 1e-14);
    ASSERT_NEAR(cl(0., 1), cl(3., 1), 1e-12);
    // one point: constant
    const auto one = radius_function<T>({{0.5, 0.2}}, RadiusVariation::Smooth);
    ASSERT_NEAR(one(0.5), 0.2, 1e-14);
    ASSERT_NEAR(one(1.2), 0.2, 1e-14);
    // refused
    ASSERT_THROW(radius_function<T>({}, RadiusVariation::Linear), std::invalid_argument);
    ASSERT_THROW(radius_function<T>({{0., 0.1}, {0., 0.2}}, RadiusVariation::Linear), std::invalid_argument);
    ASSERT_THROW(radius_function<T>({{0., 0.1}, {1., -0.2}}, RadiusVariation::Linear), std::invalid_argument);
}

TEST(tests_brep_variable_fillet, box_edge)
{
    // an edge of a box between two planes at right angles: the ball of radius r(x) centred at
    // (x, r, 1 - r) touching y = 0 and z = 1
    for (const auto variation : {RadiusVariation::Linear, RadiusVariation::Smooth})
    {
        INFO("smooth " << (variation == RadiusVariation::Smooth));
        Model<T> m;
        const auto s = box(m, 2., 1., 1.);
        const auto e = edge_near(m, s, {1., 0., 1.});
        const RadiusLaw<T> law{{{e, 0., 0.1}, {e, 0.4, 0.3}, {e, 1., 0.15}}, variation};
        FilletOptions<T> o;
        const auto b = unwrap(fillet_geometry(m, s, e, law, o));
        ASSERT_TRUE(b.convex);
        ASSERT_FALSE(b.closed);
        ASSERT_FALSE(b.exact);
        // a linear law: two pieces and the lune between them at 40 %; a smooth one: one surface
        ASSERT_EQ(b.lunes.size(), variation == RadiusVariation::Linear ? 1u : 0u);
        ASSERT_EQ(b.pieces.size(), variation == RadiusVariation::Linear ? 2u : 0u);
        const auto r_of = radius_of(m, e, law);
        expect_variable(m, s, b, r_of, o.approx_tol, o.angular_tol);
        const auto &ed = m.edge(e);
        for (int i = 0; i <= 32; ++i)
        {
            const T t = ed.u1 + (ed.u2 - ed.u1) * i / 32.;
            const T r = r_of(t);
            const auto c = b.spine->value(t);
            ASSERT_LE(norm(c - point<T, 3>{ed.curve->value(t)[0], r, 1. - r}), o.approx_tol);
        }
        // the law along the edge from its first vertex: 0.1 there, 0.3 at 40 % of its length
        const auto v1 = m.vertex(ed.v1).pnt;
        T t04 = ed.u1;
        for (int i = 0; i <= 1000; ++i)
        {
            const T t = ed.u1 + (ed.u2 - ed.u1) * i / 1000.;
            if (std::abs(norm(ed.curve->value(t) - v1) - 0.8) < std::abs(norm(ed.curve->value(t04) - v1) - 0.8))
                t04 = t;
        }
        ASSERT_NEAR(r_of(ed.u1), 0.1, 1e-12);
        ASSERT_NEAR(r_of(ed.u2), 0.15, 1e-12);
        ASSERT_NEAR(r_of(t04), 0.3, 1e-3);
        // continued past the vertices: at constant radius (smooth: its slope 0 there), with its
        // slope (linear: no lune there)
        ASSERT_LT(b.a, ed.u1);
        ASSERT_GT(b.b, ed.u2);
        ASSERT_NEAR(norm(b.contacts[0].curve->value(b.a) - b.spine->value(b.a)), r_of(b.a), o.approx_tol);
        const T past = (b.a - ed.u1) * norm(ed.curve->value(ed.u1, 1)); // (negative)
        ASSERT_NEAR(r_of(b.a), variation == RadiusVariation::Linear ? 0.1 + 0.2 / 0.8 * past : 0.1, 1e-9);
    }
}

TEST(tests_brep_variable_fillet, constant_law)
{
    // the same radius everywhere: the fillet of constant radius, a cylinder
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto b = unwrap(fillet_geometry(m, s, e, RadiusLaw<T>{{{e, 0., 0.2}, {e, 1., 0.2}}}));
    ASSERT_TRUE(b.exact);
    ASSERT_NE(dynamic_cast<const CylindricalSurface<T> *>(b.surface.get()), nullptr);
}

TEST(tests_brep_variable_fillet, cylinder_rim)
{
    // the top rim of a cylinder, a closed edge between a plane and a cylinder: a smooth law
    // closed round the rim
    Model<T> m;
    const auto s = cylinder(m, 1., 1.5);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.5});
    ASSERT_EQ(m.edge(e).v1, m.edge(e).v2);
    const RadiusLaw<T> law{{{e, 0., 0.1}, {e, 0.25, 0.25}, {e, 0.6, 0.15}}};
    FilletOptions<T> o;
    const auto b = unwrap(fillet_geometry(m, s, e, law, o));
    ASSERT_TRUE(b.closed);
    expect_variable(m, s, b, radius_of(m, e, law), o.approx_tol, o.angular_tol);
    // a loop: its ends together
    ASSERT_LE(norm(b.spine->value(b.a) - b.spine->value(b.b)), o.approx_tol);
    ASSERT_TRUE(b.lunes.empty());
    // linear: a lune at each of its three points, the one at the vertex between the last piece
    // and the first
    const RadiusLaw<T> lin{law.points, RadiusVariation::Linear};
    const auto l = unwrap(fillet_geometry(m, s, e, lin, o));
    ASSERT_EQ(l.pieces.size(), 3u);
    ASSERT_EQ(l.lunes.size(), 3u);
    expect_variable(m, s, l, radius_of(m, e, lin), o.approx_tol, o.angular_tol);
}

TEST(tests_brep_variable_fillet, bspline_face)
{
    // the top edge of the prism of a B-spline bump, square at its ends: a plane and a B-spline face
    Model<T> m;
    const auto bump = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {0.4, 0., 0.}, {0.8, -0.3, 0.}, {1.2, 0.2, 0.}, {1.6, 0., 0.}, {2., 0., 0.}},
                                                      std::vector<T>{0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.}, std::size_t{3});
    const auto v0 = unwrap(make_vertex(m, point<T, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<T, 3>{2., 0., 0.}));
    const auto v2 = unwrap(make_vertex(m, point<T, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<T, 3>{0., 1., 0.}));
    const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                unwrap(make_edge(m, v3, v0))};
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                    unwrap(make_wire(m, w))));
    const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., 0.8}));
    const auto e = edge_near(m, s, bump->value(0.5) + point<T, 3>{0., 0., 0.8});
    for (const auto variation : {RadiusVariation::Linear, RadiusVariation::Smooth})
    {
        INFO("smooth " << (variation == RadiusVariation::Smooth));
        const RadiusLaw<T> law{{{e, 0., 0.05}, {e, 0.5, 0.12}, {e, 1., 0.08}}, variation};
        FilletOptions<T> o;
        const auto b = unwrap(fillet_geometry(m, s, e, law, o));
        expect_variable(m, s, b, radius_of(m, e, law), o.approx_tol, o.angular_tol);
    }
}

TEST(tests_brep_variable_fillet, refusals)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto other = edge_near(m, s, {1., 1., 1.});
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    // no radius on the edge, two radii at a point, a point out of the edge, a radius not positive
    ASSERT_EQ(code(fillet_geometry(m, s, e, RadiusLaw<T>{{{other, 0., 0.1}}})), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet_geometry(m, s, e, RadiusLaw<T>{{{e, 0.5, 0.1}, {e, 0.5, 0.2}}})), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet_geometry(m, s, e, RadiusLaw<T>{{{e, 1.5, 0.1}}})), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet_geometry(m, s, e, RadiusLaw<T>{{{e, 0., 0.1}, {e, 1., 0.}}})), BuildErrc::InvalidBounds);
    // wider than the faces next to it
    ASSERT_EQ(code(fillet_geometry(m, s, e, RadiusLaw<T>::linear(e, 0.1, 1.5))), BuildErrc::NonGeneric);
    // a closed edge with two radii at its vertex
    const auto c = cylinder(m, 1., 1.5);
    const auto rim = edge_near(m, c, {std::cos(1.), std::sin(1.), 1.5});
    ASSERT_EQ(code(fillet_geometry(m, c, rim, RadiusLaw<T>::linear(rim, 0.1, 0.2))), BuildErrc::InvalidBounds);
}

TEST(tests_brep_variable_fillet, law_past_its_points)
{
    // points inside the edge: past them the law goes on with its slope, into the edge (linear)
    // or at constant radius (smooth); a linear law no longer positive at an end of the blend is
    // refused
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto &ed = m.edge(e);
    FilletOptions<T> o;
    {
        const RadiusLaw<T> law{{{e, 0.3, 0.15}, {e, 0.7, 0.2}}, RadiusVariation::Linear};
        const auto r_of = radius_of(m, e, law);
        // 0.15 at 0.6 from v1, slope 0.05 / 0.8: 0.1125 at v1, 0.2375 at v2
        ASSERT_NEAR(r_of(ed.u1), 0.1125, 1e-9);
        ASSERT_NEAR(r_of(ed.u2), 0.2375, 1e-9);
        const auto b = unwrap(fillet_geometry(m, s, e, law, o));
        ASSERT_TRUE(b.lunes.empty()); // no kink at the first and last points
        expect_variable(m, s, b, r_of, o.approx_tol, o.angular_tol);
    }
    {
        const RadiusLaw<T> law{{{e, 0.3, 0.15}, {e, 0.7, 0.2}}, RadiusVariation::Smooth};
        const auto r_of = radius_of(m, e, law);
        ASSERT_NEAR(r_of(ed.u1), 0.15, 1e-12);
        ASSERT_NEAR(r_of(ed.u2), 0.2, 1e-12);
        expect_variable(m, s, unwrap(fillet_geometry(m, s, e, law, o)), r_of, o.approx_tol, o.angular_tol);
    }
    {
        // 0.05 at 0.6 from v1, slope 0.25 / 0.8: below 0 at v1
        const auto r = fillet_geometry(m, s, e, RadiusLaw<T>{{{e, 0.3, 0.05}, {e, 0.7, 0.3}}, RadiusVariation::Linear}, o);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
        ASSERT_NE(r.error().message.find("not positive at an end of the blend"), std::string::npos);
    }
}

TEST(tests_brep_variable_fillet, extension_and_largest_radius)
{
    // a linear law steep at its first vertex, growing past it, the blend continued by three
    // times the larger radius at the vertices: the balls at the ends of the blend, larger than
    // any radius on the edge (the supports and the offsets checked for them), on the faces
    Model<T> m;
    const auto bump = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {0.4, 0., 0.}, {0.8, -0.3, 0.}, {1.2, 0.2, 0.}, {1.6, 0., 0.}, {2., 0., 0.}},
                                                      std::vector<T>{0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.}, std::size_t{3});
    const auto v0 = unwrap(make_vertex(m, point<T, 3>{0., 0., 0.})), v1 = unwrap(make_vertex(m, point<T, 3>{2., 0., 0.}));
    const auto v2 = unwrap(make_vertex(m, point<T, 3>{2., 1., 0.})), v3 = unwrap(make_vertex(m, point<T, 3>{0., 1., 0.}));
    const std::vector<EdgeId> w{unwrap(make_edge(m, bump, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v2)), unwrap(make_edge(m, v2, v3)),
                                unwrap(make_edge(m, v3, v0))};
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                    unwrap(make_wire(m, w))));
    const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., 0.8}));
    const auto e = edge_near(m, s, bump->value(0.5) + point<T, 3>{0., 0., 0.8});
    const RadiusLaw<T> law{{{e, 0., 0.05}, {e, 0.02, 0.03}, {e, 1., 0.08}}, RadiusVariation::Linear};
    FilletOptions<T> o;
    o.extension = 3.;
    const auto r_of = radius_of(m, e, law);
    const auto b = unwrap(fillet_geometry(m, s, e, law, o));
    ASSERT_GT(r_of(b.a), 0.08); // larger than every radius on the edge
    ASSERT_NEAR(norm(b.contacts[0].curve->value(b.a) - b.spine->value(b.a)), r_of(b.a), o.approx_tol);
    expect_variable(m, s, b, r_of, o.approx_tol, o.angular_tol);
}

TEST(tests_brep_variable_fillet, closed_law_closing_inside)
{
    // a smooth law round the rim of a cylinder whose first point is not at the vertex: it closes
    // at that point, C1 only there, a knot of the law (a station)
    Model<T> m;
    const auto s = cylinder(m, 1., 1.5);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.5});
    const auto &ed = m.edge(e);
    const RadiusLaw<T> law{{{e, 0.25, 0.1}, {e, 0.5, 0.2}, {e, 0.8, 0.12}}};
    const auto R = fillet_detail::edge_radius(m, e, law, 1e-6, "test");
    REQUIRE(R.has_value());
    const T closes = ed.u1 + 0.25 * (ed.u2 - ed.u1); // (a circle: the parameter along the length)
    ASSERT_TRUE(std::ranges::any_of(R->knots, [&](T t) { return std::abs(t - closes) < 1e-9; }));
    FilletOptions<T> o;
    const auto b = unwrap(fillet_geometry(m, s, e, law, o));
    expect_variable(m, s, b, R->r, o.approx_tol, o.angular_tol);
}

TEST(tests_brep_variable_fillet, degenerate_edge)
{
    // the edge at a pole of a ball (made by revolution): no curve, refused as such
    Model<T> m;
    const auto e1 = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., ax2<T, 3>{{{0., 0., 0.}, cross(X, Z), X}}, -pi / 2, pi / 2)));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., 1.}, point<T, 3>{0., 0., -1.}));
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                    unwrap(make_wire(m, std::vector<EdgeId>{e1, e2}))));
    const auto ball = solid_of(m, make_revolution(m, f, ax1<T, 3>{{{0., 0., 0.}, Z}}));
    std::optional<EdgeId> pole;
    for (auto e : explore<EdgeId>(m, ShapeId{ball}))
        if (m.edge(e).degenerate)
            pole = e;
    REQUIRE(pole.has_value());
    const auto r = fillet_geometry(m, ball, *pole, RadiusLaw<T>{{{*pole, 0., 0.1}, {*pole, 1., 0.2}}});
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::InvalidId);
}
