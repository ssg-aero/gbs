// Fillets of continuous curvature (stage 7, PR 3, docs/sources/design/brep_palier7_blend_sections.md,
// § 4.4, docs/sources/design/brep7_pr03_g2_sections.md): the quartic G2 section of rule C, exact
// extrusions along straight edges (between planes, between a plane and a cylinder along the edge),
// skinned round curved ones; its contacts, its curvature at them and its volume against the
// quartic worked out here on its own; corners in vertex blends, chains, a radius law; refusals.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs/bssurf.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>

#include <array>
#include <cmath>
#include <functional>
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

    SolidId prism(Model<T> &m, const std::vector<point<T, 3>> &p, T h)
    {
        return solid_of(m, make_prism(m, polygon(m, p), point<T, 3>{0., 0., h}));
    }

    SolidId box(Model<T> &m, T a, T b, T c)
    {
        return prism(m, {{0., 0., 0.}, {a, 0., 0.}, {a, b, 0.}, {0., b, 0.}}, c);
    }

    SolidId cylinder(Model<T> &m, T r, T h)
    {
        const ax2<T, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    /// The prism of height h on the half disc of radius R (y >= 0): its vertical edges at (+-R, 0) between a plane and the cylinder, at right angles.
    SolidId half_disc(Model<T> &m, T R, T h)
    {
        std::vector<EdgeId> e{unwrap(make_edge(m, point<T, 3>{-R, 0., 0.}, point<T, 3>{R, 0., 0.})),
                              unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(R, R, ax2<T, 3>{{{0., 0., 0.}, Z, X}}, 0., pi)))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    /// The slab on a stadium: straight sides of length L, half circles of radius R, height h.
    SolidId oblong(Model<T> &m, T L, T R, T h)
    {
        auto arc = [](T r, const point<T, 3> &c, T a1, T a2) { return std::make_shared<EllipseCurve<T>>(r, r, ax2<T, 3>{{c, Z, X}}, a1, a2); };
        std::vector<EdgeId> e{unwrap(make_edge(m, point<T, 3>{0., -R, 0.}, point<T, 3>{L, -R, 0.})), unwrap(make_edge(m, arc(R, {L, 0., 0.}, -pi / 2., pi / 2.))),
                              unwrap(make_edge(m, point<T, 3>{L, R, 0.}, point<T, 3>{0., R, 0.})), unwrap(make_edge(m, arc(R, {0., 0., 0.}, pi / 2., 3. * pi / 2.)))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, e))));
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

    T volume(const Model<T> &m, SolidId s)
    {
        return signed_volume(m, m.solid(s).outer, 512);
    }

    std::size_t n_faces(const Model<T> &m, SolidId s)
    {
        return explore<FaceId>(m, ShapeId{s}).size();
    }

    std::size_t n_spheres(const Model<T> &m, SolidId s)
    {
        std::size_t n{};
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (dynamic_cast<const SphericalSurface<T> *>(m.face(f).surface.get()))
                ++n;
        return n;
    }

    void expect_valid(const Model<T> &m, SolidId s)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
    }

    /// Gauss-Legendre on [a, b], 10 points each side of the middle of each of n panels.
    template <typename F>
    T gauss(F &&f, T a, T b, int n = 32)
    {
        static const std::array<T, 10> x{0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
                                         0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949};
        static const std::array<T, 10> w{0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
                                         0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521};
        T sum{};
        for (int k = 0; k < n; ++k)
        {
            const T a1 = a + (b - a) * T(k) / T(n), b1 = a + (b - a) * T(k + 1) / T(n);
            const T c = (a1 + b1) / 2., h = (b1 - a1) / 2.;
            for (std::size_t i = 0; i < 10; ++i)
                sum += h * w[i] * (f(c - h * x[i]) + f(c + h * x[i]));
        }
        return sum;
    }

    using P2 = std::array<T, 2>;

    /// The planar quartic Bezier of poles P (or its derivative of order k) at t.
    P2 bez(std::array<P2, 5> P, T t, int k = 0)
    {
        std::size_t n = 4;
        for (int j = 0; j < k; ++j, --n)
            for (std::size_t i = 0; i < n; ++i)
                P[i] = {T(n) * (P[i + 1][0] - P[i][0]), T(n) * (P[i + 1][1] - P[i][1])};
        for (std::size_t j = n; j > 0; --j)
            for (std::size_t i = 0; i < j; ++i)
                P[i] = {(1. - t) * P[i][0] + t * P[i + 1][0], (1. - t) * P[i][1] + t * P[i + 1][1]};
        return P[0];
    }

    T curvature(const std::array<P2, 5> &P, T t)
    {
        const auto d1 = bez(P, t, 1), d2 = bez(P, t, 2);
        return (d1[0] * d2[1] - d1[1] * d2[0]) / std::pow(std::hypot(d1[0], d1[1]), 3.);
    }

    /**
     * The G2 section of rule C between two planes at the opening angle a, the edge at S = 0, the
     * faces along x and along (cos a, sin a), worked out here (§ 4.4): the poles at m from S,
     * lambda of the smoothest shape (golden section on the integral of (dk/ds)^2), m for the sag
     * of the arc of radius r. Its poles, lambda and m.
     */
    struct Rule
    {
        std::array<P2, 5> P;
        T lambda, m;
    };
    Rule rule_c(T a, T r)
    {
        auto poles = [&](T lam, T m) {
            const P2 d1{1., 0.}, d2{std::cos(a), std::sin(a)};
            return std::array<P2, 5>{P2{m * d1[0], m * d1[1]}, P2{(1. - lam) * m * d1[0], (1. - lam) * m * d1[1]}, P2{0., 0.},
                                     P2{(1. - lam) * m * d2[0], (1. - lam) * m * d2[1]}, P2{m * d2[0], m * d2[1]}};
        };
        // dk/dt = (B' x B''') / |B'|^3 - 3 (B' x B'') (B' . B'') / |B'|^5
        auto fairness = [&](T lam) {
            const auto P = poles(lam, 1.);
            return gauss(
                [&](T t) {
                    const auto d1 = bez(P, t, 1), d2 = bez(P, t, 2), d3 = bez(P, t, 3);
                    const T s = std::hypot(d1[0], d1[1]);
                    const T dk = (d1[0] * d3[1] - d1[1] * d3[0]) / std::pow(s, 3.) - 3. * (d1[0] * d2[1] - d1[1] * d2[0]) * (d1[0] * d2[0] + d1[1] * d2[1]) / std::pow(s, 5.);
                    return dk * dk / s;
                },
                0., 1., 8);
        };
        T lo = 0.02, hi = 0.98;
        const T g = (std::sqrt(5.) - 1.) / 2.;
        while (hi - lo > 1e-10)
        {
            const T x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
            if (fairness(x1) < fairness(x2))
                hi = x2;
            else
                lo = x1;
        }
        // (a flat minimum: then Newton on the derivative by 5 points, beyond sqrt(epsilon))
        T lam = (lo + hi) / 2.;
        for (int it = 0; it < 4; ++it)
        {
            const T d = 1e-3, e2 = fairness(lam - 2. * d), e1 = fairness(lam - d), e0 = fairness(lam), f1 = fairness(lam + d), f2 = fairness(lam + 2. * d);
            lam -= ((e2 - 8. * e1 + 8. * f1 - f2) / (12. * d)) / ((-e2 + 16. * e1 - 30. * e0 + 16. * f1 - f2) / (12. * d * d));
        }
        const auto mid = bez(poles(lam, 1.), 0.5);
        const T sag_arc = r / std::sin(a / 2.) - r; // S to the middle of the arc
        const T m = sag_arc / std::hypot(mid[0], mid[1]);
        return {poles(lam, m), lam, m};
    }

    /// The area between the quartic and S = 0 (Green: the straight parts pass through S).
    T area(const std::array<P2, 5> &P)
    {
        return 0.5 * std::abs(gauss([&](T t) { const auto p = bez(P, t), d = bez(P, t, 1); return p[0] * d[1] - p[1] * d[0]; }, 0., 1.));
    }

    /// The edge along x at y = 1, z = 1 of the box [0, 3] x [0, 1] x [0, 1].
    struct BoxEdge
    {
        Model<T> m;
        SolidId s;
        EdgeId e;
        BoxEdge()
        {
            s = box(m, 3., 1., 1.);
            e = edge_near(m, s, {1.5, 1., 1.});
        }
    };

    FilletOptions<T> g2()
    {
        FilletOptions<T> o;
        o.section = BlendSection<T>::g2();
        return o;
    }

    /// The curvature in the plane of normal n of the section u -> s(u, v) of a surface at (u, v), signed along n.
    T across(const Surface<T, 3> &s, T u, T v, const point<T, 3> &n)
    {
        const auto d1 = s.value(u, v, 1, 0), d2 = s.value(u, v, 2, 0);
        return (cross(d1, d2) * n) / std::pow(norm(d1), 3.);
    }
}

TEST(tests_brep_g2_sections, rule_c_at_right_angles)
{
    // the numbers of the note (§ 4.4): lambda 0.352, contacts at 1.305 r, curvature peak 1.234 / r
    const auto c = rule_c(pi / 2., 1.);
    ASSERT_NEAR(c.lambda, 0.35247, 1e-4);
    ASSERT_NEAR(c.m, 1.30533, 1e-4);
    T peak{};
    for (int i = 0; i <= 2000; ++i)
        peak = std::max(peak, std::abs(curvature(c.P, i / 2000.)));
    ASSERT_NEAR(peak, 1.2341, 1e-3);
    ASSERT_NEAR(curvature(c.P, 0.), 0., 1e-12); // G2 with the planes
}

TEST(tests_brep_g2_sections, box_edge)
{
    // r = 0.3 on the edge of the box: the extrusion of the quartic, exact (degree (4, 1)); its
    // contacts at m r from the edge, its curvature 0 at them, the curvature of the section the
    // one worked out here; the volume 3 less the area of the section along 3
    BoxEdge b;
    const T r = 0.3;
    const auto g = unwrap(fillet_geometry(b.m, b.s, b.e, r, g2()));
    ASSERT_TRUE(g.exact);
    ASSERT_FALSE(g.circular);
    const auto rule = rule_c(pi / 2., r);
    const auto *bs = dynamic_cast<const BSSurfaceRational<T, 3> *>(g.surface.get());
    ASSERT_NE(bs, nullptr);
    ASSERT_EQ(bs->degreeU(), 4u);
    for (const T t : {g.a, (g.a + g.b) / 2., g.b})
    {
        // the edge at (x, 1, 1): the contacts on the top (z = 1) and on the side (y = 1)
        const auto p1 = g.contacts[0].curve->value(t), p2 = g.contacts[1].curve->value(t);
        const T d1 = std::hypot(p1[1] - 1., p1[2] - 1.), d2 = std::hypot(p2[1] - 1., p2[2] - 1.);
        ASSERT_NEAR(d1, rule.m, 1e-9);
        ASSERT_NEAR(d2, rule.m, 1e-9);
        const auto w0 = blend_detail::across(g, t, 0.), w1 = blend_detail::across(g, t, 1.), wm = blend_detail::across(g, t, 0.5);
        ASSERT_NEAR(across(*g.surface, w0[0], w0[1], X), 0., 1e-9);
        ASSERT_NEAR(across(*g.surface, w1[0], w1[1], X), 0., 1e-9);
        ASSERT_NEAR(std::abs(across(*g.surface, wm[0], wm[1], X)), std::abs(curvature(rule.P, 0.5)), 1e-6);
        // the middle of the section as far from the edge as that of the arc of radius r
        const auto mid = g.surface->value(wm[0], wm[1]);
        ASSERT_NEAR(std::hypot(mid[1] - 1., mid[2] - 1.), r * std::sqrt(2.) - r, 1e-9);
    }
    const auto f = unwrap(fillet(b.m, b.s, std::vector{b.e}, r, g2()));
    expect_valid(b.m, f);
    ASSERT_EQ(n_faces(b.m, f), 7u);
    ASSERT_NEAR(volume(b.m, f), 3. - 3. * area(rule.P), 1e-9);
}

TEST(tests_brep_g2_sections, acute_edge)
{
    // the prism on the trapezium (0, 0), (2, 0), (1.5, 1), (0, 1): the edge at (2, 0) between faces
    // at an opening of atan2(1, 0.5) (63.4 degrees), r = 0.2: the section of rule C at that angle
    Model<T> m;
    const auto s = prism(m, {{0., 0., 0.}, {2., 0., 0.}, {1.5, 1., 0.}, {0., 1., 0.}}, 1.);
    const auto e = edge_near(m, s, {2., 0., 0.5});
    const T a = std::atan2(1., 0.5), r = 0.2;
    const auto g = unwrap(fillet_geometry(m, s, e, r, g2()));
    ASSERT_TRUE(g.exact);
    const auto rule = rule_c(a, r);
    const T t = (g.a + g.b) / 2.;
    for (std::size_t k = 0; k < 2; ++k)
    {
        const auto p = g.contacts[k].curve->value(t);
        ASSERT_NEAR(std::hypot(p[0] - 2., p[1]), rule.m, 1e-9);
    }
    const auto f = unwrap(fillet(m, s, std::vector{e}, r, g2()));
    expect_valid(m, f);
    ASSERT_NEAR(volume(m, f), 1.75 - area(rule.P), 1e-9);
    // the edge at (1.5, 1), obtuse (116.6 degrees)
    {
        Model<T> m2;
        const auto s2 = prism(m2, {{0., 0., 0.}, {2., 0., 0.}, {1.5, 1., 0.}, {0., 1., 0.}}, 1.);
        const auto e2 = edge_near(m2, s2, {1.5, 1., 0.5});
        const auto g2b = unwrap(fillet_geometry(m2, s2, e2, r, g2()));
        ASSERT_TRUE(g2b.exact);
        const auto rule2 = rule_c(pi - a, r);
        for (std::size_t k = 0; k < 2; ++k)
        {
            const auto p = g2b.contacts[k].curve->value((g2b.a + g2b.b) / 2.);
            ASSERT_NEAR(std::hypot(p[0] - 1.5, p[1] - 1.), rule2.m, 1e-9);
        }
        const auto f2 = unwrap(fillet(m2, s2, std::vector{e2}, r, g2()));
        expect_valid(m2, f2);
        ASSERT_NEAR(volume(m2, f2), 1.75 - area(rule2.P), 1e-9);
    }
}

TEST(tests_brep_g2_sections, plane_and_cylinder_along_the_edge)
{
    // the vertical edge at (1, 0) of the prism on the half disc of radius 1: a plane and the
    // cylinder at right angles, the cylinder curving across the edge (1 / R in the section). The
    // centre line straight, the section the same all along: the extrusion, exact; at its contact on
    // the cylinder the curvature of the quartic that of the circle, 0 at the plane
    Model<T> m;
    const auto s = half_disc(m, 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 0.5});
    const auto g = unwrap(fillet_geometry(m, s, e, 0.2, g2()));
    ASSERT_TRUE(g.exact);
    for (const T t : {g.a, (g.a + g.b) / 2., g.b})
        for (const T u : {0., 1.})
        {
            const auto w = blend_detail::across(g, t, u);
            const auto p = g.surface->value(w[0], w[1]);
            const bool on_cylinder = std::abs(p[1]) > 1e-9;
            const T k = std::abs(across(*g.surface, w[0], w[1], Z));
            ASSERT_NEAR(k, on_cylinder ? 1. : 0., 1e-7);
            if (on_cylinder)
                ASSERT_NEAR(std::hypot(p[0], p[1]), 1., 1e-12);
        }
    const auto f = unwrap(fillet(m, s, std::vector{e}, 0.2, g2()));
    expect_valid(m, f);
    ASSERT_LT(volume(m, f), pi / 2.);
}

TEST(tests_brep_g2_sections, round_a_cylinder)
{
    // the top rim of a cylinder of radius 1, height 1 (skinned): in each meridian plane the section
    // of the box edge, the corner at x = 1: Pappus, the volume pi less 2 pi times the first moment
    // of the area about the axis; the curvature across at the contacts 0 (both faces straight in
    // the meridian planes), within the approximation
    Model<T> m;
    const auto s = cylinder(m, 1., 1.);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.});
    const T r = 0.2;
    const auto g = unwrap(fillet_geometry(m, s, e, r, g2()));
    ASSERT_FALSE(g.exact);
    const auto rule = rule_c(pi / 2., r);
    // the first moment about the axis of the section at the corner (1, 1): x' = 1 - x, z' = 1 - y
    const T M = std::abs(gauss([&](T t) {
        const auto p = bez(rule.P, t), d = bez(rule.P, t, 1);
        const T x = 1. - p[0];
        return x * x * (-d[1]) / 2.;
    }, 0., 1.) + 0.5 * rule.m);
    for (int i = 0; i <= 8; ++i)
    {
        const T t = g.a + (g.b - g.a) * i / 8.;
        // the contacts at m r from the rim, on the top and down the side
        const auto p1 = g.contacts[0].curve->value(t), p2 = g.contacts[1].curve->value(t);
        ASSERT_NEAR(std::abs(std::hypot(p1[0], p1[1]) - 1.) + std::abs(p1[2] - 1.), rule.m, 1e-6);
        ASSERT_NEAR(std::abs(std::hypot(p2[0], p2[1]) - 1.) + std::abs(p2[2] - 1.), rule.m, 1e-6);
        for (const T u : {0., 1.})
        {
            const auto w = blend_detail::across(g, t, u);
            const auto p = g.surface->value(w[0], w[1]);
            auto n = cross(Z, point<T, 3>{p[0], p[1], 0.});
            n = n / norm(n);
            ASSERT_LT(std::abs(across(*g.surface, w[0], w[1], n)), 1e-3 / r);
        }
    }
    const auto f = unwrap(fillet(m, s, std::vector{e}, r, g2()));
    expect_valid(m, f);
    ASSERT_EQ(n_faces(m, f), 4u);
    ASSERT_NEAR(volume(m, f), pi - 2. * pi * M, 1e-5);
}

TEST(tests_brep_g2_sections, corner_and_chain)
{
    // the three edges at a corner of a box, G2 fillets: no sphere, a vertex blend; the rim of a
    // stadium slab, a chain of G2 fillets
    {
        Model<T> m;
        const auto s = box(m, 2., 2., 2.);
        const std::vector<EdgeId> at{edge_near(m, s, {1., 2., 2.}), edge_near(m, s, {2., 1., 2.}), edge_near(m, s, {2., 2., 1.})};
        const auto r = fillet(m, s, at, 0.3, g2());
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        expect_valid(m, *r);
        ASSERT_EQ(n_faces(m, *r), 10u);
        ASSERT_EQ(n_spheres(m, *r), 0u);
    }
    {
        const T L = 3., R = 1., h = 0.8;
        Model<T> m;
        const auto s = oblong(m, L, R, h);
        const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
        const auto r = fillet(m, s, rim, 0.2, g2());
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        expect_valid(m, *r);
        ASSERT_EQ(n_faces(m, *r), 10u);
    }
}

TEST(tests_brep_g2_sections, radius_law)
{
    // a G2 fillet whose radius follows a smooth law along the edge of the box (0.2 to 0.35): the
    // contacts of rule C for the radius at each station, a valid solid, less material taken away
    // than with the largest radius all along and more than with the smallest
    BoxEdge b;
    RadiusLaw<T> law{{{b.e, 0., 0.2}, {b.e, 1., 0.35}}, RadiusVariation::Smooth};
    const auto f = fillet(b.m, b.s, std::vector{b.e}, law, g2());
    INFO((f ? std::string{} : f.error().message));
    ASSERT_TRUE(f.has_value());
    expect_valid(b.m, *f);
    const T v = volume(b.m, *f);
    ASSERT_LT(v, 3. - 3. * area(rule_c(pi / 2., 0.2).P));
    ASSERT_GT(v, 3. - 3. * area(rule_c(pi / 2., 0.35).P));
}

TEST(tests_brep_g2_sections, radius_law_on_a_cylinder)
{
    // the vertical edge of the half disc, the radius from 0.15 to 0.25 (skinned): at the contact on
    // the cylinder the curvature across that of the circle (1), 0 at the plane, within 1e-3 of the
    // approximated surface
    Model<T> m;
    const auto s = half_disc(m, 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 0.5});
    const RadiusLaw<T> law{{{e, 0., 0.15}, {e, 1., 0.25}}, RadiusVariation::Smooth};
    const auto g = unwrap(fillet_geometry(m, s, e, law, g2()));
    ASSERT_FALSE(g.exact);
    for (int i = 0; i <= 8; ++i)
    {
        const T t = g.a + (g.b - g.a) * i / 8.;
        for (const T u : {0., 1.})
        {
            const auto w = blend_detail::across(g, t, u);
            const auto p = g.surface->value(w[0], w[1]);
            const bool on_cylinder = std::abs(p[1]) > 1e-6;
            // the curvature in the horizontal plane (the plane of the station, the edge vertical)
            const auto d1 = g.surface->value(w[0], w[1], 1, 0), d2 = g.surface->value(w[0], w[1], 2, 0);
            const auto h1 = point<T, 3>{d1[0], d1[1], 0.}, h2 = point<T, 3>{d2[0], d2[1], 0.};
            const T k = std::abs(cross(h1, h2)[2]) / std::pow(norm(h1), 3.);
            ASSERT_NEAR(k, on_cylinder ? 1. : 0., 1e-3);
        }
    }
    const auto f = unwrap(fillet(m, s, std::vector{e}, law, g2()));
    expect_valid(m, f);
}

TEST(tests_brep_g2_sections, refusals)
{
    auto says = [](const auto &r, const std::string &what) { return !r.has_value() && r.error().message.find(what) != std::string::npos; };
    {
        // constant chord and G2: the contacts of the section are not those of its ball
        BoxEdge b;
        const auto before = n_faces(b.m, b.s);
        const auto r = fillet_chordal(b.m, b.s, std::vector{b.e}, 0.3, g2());
        ASSERT_TRUE(says(r, "a fillet of constant chord with a G2 section"));
        ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
        ASSERT_EQ(n_faces(b.m, b.s), before);
    }
    {
        // a linear radius law with a kink: the lune of the ball meets arcs of circle only
        BoxEdge b;
        RadiusLaw<T> law{{{b.e, 0., 0.2}, {b.e, 0.5, 0.3}, {b.e, 1., 0.2}}, RadiusVariation::Linear};
        const auto r = fillet(b.m, b.s, std::vector{b.e}, law, g2());
        ASSERT_TRUE(says(r, "a G2 fillet whose radius law kinks"));
    }
}
