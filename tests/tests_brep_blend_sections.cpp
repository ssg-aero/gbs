// Fillets whose section is not an arc of circle (stage 7, PR 2,
// docs/sources/design/brep_palier7_blend_sections.md, § 4): conic sections of rho constant or
// along a law, exact extrusions along straight edges, skinned round curved ones; corners in vertex
// blends; a law of rho with a kink (pieces), fillets of constant chord, along a chain too;
// refusals with their messages; the volumes against the areas of the conics, worked out on their
// own.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs/bssurf.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>

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

    /// The planar face of the solid through p whose normal is along n (either way).
    FaceId face_at(const Model<T> &m, SolidId s, const point<T, 3> &p, const point<T, 3> &n)
    {
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (const auto *pl = dynamic_cast<const Plane<T> *>(m.face(f).surface.get()))
            {
                const auto uv = pl->parameters(p);
                const auto q = pl->value(uv[0], uv[1]);
                if (norm(q - p) < 1e-9 && std::abs(normal_direction(*pl, uv[0], uv[1]) * n) > 1. - 1e-9)
                    return f;
            }
        return FaceId{};
    }

    T volume(const Model<T> &m, SolidId s)
    {
        return signed_volume(m, m.solid(s).outer, 512);
    }

    std::size_t n_faces(const Model<T> &m, SolidId s)
    {
        return explore<FaceId>(m, ShapeId{s}).size();
    }

    void expect_solid(const Model<T> &m, SolidId s, std::size_t faces, T v, T eps)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
        ASSERT_EQ(n_faces(m, s), faces);
        ASSERT_NEAR(volume(m, s), v, eps);
    }

    /// Gauss-Legendre on [a, b], 10 points each side of the middle of each of n panels.
    template <typename F>
    T gauss(F &&f, T a, T b, int n = 64)
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

    /// The law on an edge of length L as a function of the arc length from its first vertex (`radius_function`).
    std::function<T(T)> law_of(const BlendLaw<T> &law, T L)
    {
        std::vector<std::pair<T, T>> pts;
        for (const auto &p : law.points)
            pts.emplace_back(p.at * L, p.value);
        if (pts.size() == 1)
            return [v = pts.front().second](T) { return v; };
        const auto f = radius_detail::radius_function(pts, law.variation);
        return [f](T s) { return f(s); };
    }

    /// The edge along x at y = 1, z = 1 of the box [0, 3] x [0, 1] x [0, 1], and its top face (z = 1).
    struct BoxEdge
    {
        Model<T> m;
        SolidId s;
        EdgeId e;
        FaceId top, side;
        BoxEdge()
        {
            s = box(m, 3., 1., 1.);
            e = edge_near(m, s, {1.5, 1., 1.});
            top = face_at(m, s, {1.5, 0.5, 1.}, Z);
            side = face_at(m, s, {1.5, 1., 0.5}, Y);
        }
    };

    /**
     * The conic of weight w between p1 = (r, 0) and p2 = (0, r) with S = (0, 0) its middle pole
     * (a fillet of radius r at a right-angled edge, the edge at S): the area between it and S, and
     * its first moment about the line x = 0 (Green, the rational form differentiated exactly).
     */
    std::array<T, 2> conic_area(T r, T w, T axis_x = 0., T sign = 1.)
    {
        // P(t) = N(t) / D(t), N = r ((1 - t)^2, t^2) (S = 0), D = (1 - t)^2 + 2 w t (1 - t) + t^2
        auto P = [&](T t) {
            const T a = (1. - t) * (1. - t), c = t * t, d = a + 2. * w * t * (1. - t) + c;
            return std::array<T, 2>{r * a / d, r * c / d};
        };
        auto dP = [&](T t) {
            const T h = 1e-6;
            const auto p = P(t + h), q = P(t - h);
            return std::array<T, 2>{(p[0] - q[0]) / (2. * h), (p[1] - q[1]) / (2. * h)};
        };
        // the closed path S -> p1 (straight), the conic p1 -> p2, p2 -> S (straight): only the
        // conic contributes to x dy - y dx (the straight parts pass through S); first moment of
        // x' = axis_x + sign x about the axis: x'^2 dy / 2 along the conic and p2 -> S
        const T area = 0.5 * std::abs(gauss([&](T t) { const auto p = P(t), d = dP(t); return p[0] * d[1] - p[1] * d[0]; }, 0., 1.));
        const T m_conic = gauss([&](T t) { const auto p = P(t), d = dP(t); const T x = axis_x + sign * p[0]; return x * x * d[1] / 2.; }, 0., 1.);
        const T m_line = axis_x * axis_x / 2. * (0. - r); // p2 -> S along x = 0: y from r to 0
        return {area, std::abs(m_conic + m_line)};
    }

    T rho_circle()
    {
        const T c = std::cos(pi / 4.);
        return c / (1. + c);
    }

    std::size_t n_spheres(const Model<T> &m, SolidId s)
    {
        std::size_t n{};
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (dynamic_cast<const SphericalSurface<T> *>(m.face(f).surface.get()))
                ++n;
        return n;
    }
}

TEST(tests_brep_blend_sections, conic_at_the_circle_is_the_circle)
{
    // rho = cos 45 / (1 + cos 45): the conic of the arc, the same solid as the circular fillet
    BoxEdge b;
    FilletOptions<T> o;
    o.section = BlendSection<T>::conic(rho_circle());
    const auto g = unwrap(fillet_geometry(b.m, b.s, b.e, 0.3, o));
    ASSERT_TRUE(g.exact);
    ASSERT_FALSE(g.circular);
    const auto c = unwrap(fillet_geometry(b.m, b.s, b.e, 0.3));
    for (int i = 0; i <= 8; ++i)
        for (int j = 0; j <= 8; ++j)
        {
            const T t = g.a + (g.b - g.a) * i / 8., u = j / 8.;
            const auto w = blend_detail::across(g, t, u);
            const auto p = g.surface->value(w[0], w[1]);
            const auto ax = c.spine->value(t);
            ASSERT_NEAR(norm(p - ax), 0.3, 1e-12);
        }
    const auto f = unwrap(fillet(b.m, b.s, std::vector{b.e}, 0.3, o));
    const auto fc = unwrap(fillet(b.m, b.s, std::vector{b.e}, 0.3));
    ASSERT_NEAR(volume(b.m, f), volume(b.m, fc), 1e-9);
}

TEST(tests_brep_blend_sections, conics_on_a_box_edge)
{
    // rho 0.3, 0.5, 0.7 at r = 0.3: the extrusion of the conic, exact; the volume 3 less the area
    // between the conic and the edge along 3, falling as rho grows (fuller sections)
    std::vector<T> vols;
    for (const T rho : {0.3, 0.5, 0.7})
    {
        BoxEdge b;
        FilletOptions<T> o;
        o.section = BlendSection<T>::conic(rho);
        const auto g = unwrap(fillet_geometry(b.m, b.s, b.e, 0.3, o));
        ASSERT_TRUE(g.exact);
        const auto f = unwrap(fillet(b.m, b.s, std::vector{b.e}, 0.3, o));
        const T A = conic_area(0.3, rho / (1. - rho))[0];
        expect_solid(b.m, f, 7, 3. - 3. * A, 1e-9);
        vols.push_back(volume(b.m, f));
    }
    ASSERT_LT(vols[0], vols[1]);
    ASSERT_LT(vols[1], vols[2]);
    ASSERT_NEAR(conic_area(1., 1.)[0], 1. / 6., 1e-12);                       // the parabola
    ASSERT_NEAR(conic_area(1., std::cos(pi / 4.))[0], 1. - pi / 4., 1e-12); // the circle
}

TEST(tests_brep_blend_sections, conic_law_along_an_edge)
{
    // rho smooth from 0.3 to 0.7 along the edge: sections skinned (not exact), the volume within
    // approx_tol of 3 less the area of the conic of rho(s) integrated along the edge
    BoxEdge b;
    FilletOptions<T> o;
    o.section = BlendSection<T>::conic(BlendLaw<T>{{{b.e, 0., 0.3}, {b.e, 1., 0.7}}});
    const auto g = unwrap(fillet_geometry(b.m, b.s, b.e, 0.3, o));
    ASSERT_FALSE(g.exact);
    ASSERT_LE(g.deviation, o.approx_tol);
    const auto f = unwrap(fillet(b.m, b.s, std::vector{b.e}, 0.3, o));
    const auto rho = law_of(*o.section.rho_law, 3.);
    const T off = gauss([&](T s) { const T r = rho(s); return conic_area(0.3, r / (1. - r))[0]; }, 0., 3., 32);
    expect_solid(b.m, f, 7, 3. - off, 1e-5);
}

TEST(tests_brep_blend_sections, conic_law_with_a_kink)
{
    // rho linear 0.3, 0.7, 0.3 along a box edge: the surface creased at the kink (the sections
    // sweep on otherwise), in two pieces sharing the section there; the volume 3 less the area of
    // the conic of rho(s) along the edge; along a fillet whose radius varies, refused
    BoxEdge b;
    FilletOptions<T> o;
    o.section = BlendSection<T>::conic(BlendLaw<T>{{{b.e, 0., 0.3}, {b.e, 0.5, 0.7}, {b.e, 1., 0.3}}, BlendVariation::Linear});
    const auto g = unwrap(fillet_geometry(b.m, b.s, b.e, 0.3, o));
    ASSERT_EQ(g.pieces.size(), 2u);
    ASSERT_LE(g.deviation, o.approx_tol);
    const auto f = unwrap(fillet(b.m, b.s, std::vector{b.e}, 0.3, o));
    const auto rho = law_of(*o.section.rho_law, 3.);
    T off{};
    for (const auto &[s0, s1] : {std::pair{0., 1.5}, std::pair{1.5, 3.}})
        off += gauss([&](T s) { const T r = rho(s); return conic_area(0.3, r / (1. - r))[0]; }, s0, s1, 16);
    expect_solid(b.m, f, 8, 3. - off, 1e-5);
    const auto refused = fillet(b.m, b.s, std::vector{b.e}, BlendLaw<T>::linear(b.e, 0.2, 0.3), o);
    ASSERT_FALSE(refused.has_value());
    ASSERT_NE(refused.error().message.find("a conic fillet whose law of rho kinks along a fillet whose radius varies"), std::string::npos);
}

TEST(tests_brep_blend_sections, conic_round_a_cylinder)
{
    // the top rim of a cylinder of radius 1, height 1: a conic of rho 0.6 and radius 0.2 round a
    // circular centre line (skinned); in each meridian plane the section of the box edge, the
    // corner at x = 1: Pappus, 2 pi times the first moment of the area about the axis
    Model<T> m;
    const auto s = cylinder(m, 1., 1.);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.});
    FilletOptions<T> o;
    o.section = BlendSection<T>::conic(0.6);
    const auto f = unwrap(fillet(m, s, std::vector{e}, 0.2, o));
    const T M = conic_area(0.2, 0.6 / 0.4, 1., -1.)[1];
    expect_solid(m, f, 4, pi - 2. * pi * M, 1e-5);
}

TEST(tests_brep_blend_sections, conic_corner_and_chain)
{
    // the three edges at a corner of a box, conic fillets of rho 0.6: no sphere, a vertex blend;
    // the rim of a stadium slab, a chain of conic fillets
    {
        Model<T> m;
        const auto s = box(m, 2., 2., 2.);
        const std::vector<EdgeId> at{edge_near(m, s, {1., 2., 2.}), edge_near(m, s, {2., 1., 2.}), edge_near(m, s, {2., 2., 1.})};
        FilletOptions<T> o;
        o.section = BlendSection<T>::conic(0.6);
        const auto r = fillet(m, s, at, 0.3, o);
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        ASSERT_TRUE(check(m, ShapeId{*r}).ok());
        ASSERT_EQ(n_faces(m, *r), 10u);
        ASSERT_EQ(n_spheres(m, *r), 0u);
    }
    {
        const T L = 3., R = 1., h = 0.8;
        Model<T> m;
        const auto s = oblong(m, L, R, h);
        const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
        FilletOptions<T> o;
        o.section = BlendSection<T>::conic(0.6);
        const auto r = fillet(m, s, rim, 0.2, o);
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        ASSERT_TRUE(check(m, ShapeId{*r}).ok());
        ASSERT_EQ(n_faces(m, *r), 10u);
    }
}

TEST(tests_brep_blend_sections, chordal)
{
    // between two planes at right angles, the chord c: the cylinder of radius c / (2 sin 45),
    // the same solid as the fillet of that radius; at 120 degrees (the prism on a quadrilateral)
    // c / (2 sin 30)
    {
        BoxEdge b;
        const T c = 0.4, r = c / (2. * std::sin(pi / 4.));
        const auto g = unwrap(fillet_chordal_geometry(b.m, b.s, b.e, c));
        // the contacts c apart across the edge (lines parallel to it, each on its own parameters)
        const auto &ed = b.m.edge(b.e);
        const auto z = ed.curve->value(ed.u1, 1) / norm(ed.curve->value(ed.u1, 1));
        const auto [a1, b1] = g.contacts[1].curve->bounds();
        const auto q = g.contacts[1].curve->value((a1 + b1) / 2.);
        const auto [a0, b0] = g.contacts[0].curve->bounds();
        for (int i = 0; i <= 8; ++i)
        {
            const auto d = q - g.contacts[0].curve->value(a0 + (b0 - a0) * i / 8.);
            ASSERT_NEAR(norm(d - (d * z) * z), c, 1e-7);
        }
        // the exact cylinder of the fillet of that radius (it was a rational B-spline on it)
        ASSERT_TRUE(g.exact);
        ASSERT_FALSE(g.variable);
        const auto *cyl = dynamic_cast<const CylindricalSurface<T> *>(g.surface.get());
        ASSERT_TRUE(cyl != nullptr);
        ASSERT_NEAR(cyl->radius(), r, 1e-14);
        const auto f = unwrap(fillet_chordal(b.m, b.s, std::vector{b.e}, c));
        expect_solid(b.m, f, 7, 3. - 3. * (1. - pi / 4.) * r * r, 1e-7);        // the chord and a conic section together: the conic of radius c / (2 sin 45)
        FilletOptions<T> o;
        o.section = BlendSection<T>::conic(0.6);
        const auto fc = unwrap(fillet_chordal(b.m, b.s, std::vector{b.e}, c, o));
        expect_solid(b.m, fc, 7, 3. - 3. * conic_area(r, 0.6 / 0.4)[0], 1e-6);
    }
    {
        Model<T> m;
        const T h = 2., cs = std::cos(2. * pi / 3.), s3 = std::sin(2. * pi / 3.);
        const auto s = prism(m, {{0., 0., 0.}, {3., 0., 0.}, {3., 2. * s3, 0.}, {2. * cs, 2. * s3, 0.}}, h);
        const auto e = edge_near(m, s, {0., 0., 1.});
        const T c = 0.4, r = c / (2. * std::sin(pi / 6.));
        const auto f = unwrap(fillet_chordal(m, s, std::vector{e}, c));
        const auto fr = unwrap(fillet(m, s, std::vector{e}, r));
        ASSERT_NEAR(volume(m, f), volume(m, fr), 1e-7);
        const auto g = unwrap(fillet_chordal_geometry(m, s, e, c));
        const auto *cyl = dynamic_cast<const CylindricalSurface<T> *>(g.surface.get());
        ASSERT_TRUE(g.exact && cyl != nullptr);
        ASSERT_NEAR(cyl->radius(), r, 1e-14);
    }
}

TEST(tests_brep_blend_sections, chordal_along_a_chain)
{
    // the rim of a stadium slab, a chain of two lines and two half circles, the faces at right
    // angles all round: the radius c / (2 sin 45) everywhere, the fillet of that radius (the
    // fillet continued by 1.5 radius at the vertices of the chain: refused before, the radius of
    // the chord not known that far)
    const T L = 3., R = 1., h = 0.8, c = 0.2;
    Model<T> m;
    const auto s = oblong(m, L, R, h);
    const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
    const auto f = fillet_chordal(m, s, rim, c);
    INFO((f ? std::string{} : f.error().message));
    ASSERT_TRUE(f.has_value());
    ASSERT_TRUE(check(m, ShapeId{*f}).ok());
    ASSERT_NEAR(volume(m, *f), volume(m, unwrap(fillet(m, s, rim, c / (2. * std::sin(pi / 4.))))), 1e-6);
}

TEST(tests_brep_blend_sections, chordal_between_plane_and_cylinder)
{
    // the top rim of a cylinder: the faces at right angles all round, the chord constant; the
    // radius c / (2 sin 45) everywhere, the torus of the fillet of that radius (by volume)
    Model<T> m;
    const auto s = cylinder(m, 1., 1.);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.});
    const T c = 0.2;
    const auto f = unwrap(fillet_chordal(m, s, std::vector{e}, c));
    const auto fr = unwrap(fillet(m, s, std::vector{e}, c / (2. * std::sin(pi / 4.))));
    ASSERT_TRUE(check(m, ShapeId{f}).ok());
    ASSERT_NEAR(volume(m, f), volume(m, fr), 1e-6);
}

TEST(tests_brep_blend_sections, chordal_angle_varying)
{
    // a cylinder of radius 1 cut by a plane tilted by 30 degrees: along the elliptic rim the angle
    // of the faces varies, the radius with it; the contacts c apart at every station
    Model<T> m;
    const auto cyl = cylinder(m, 1., 2.);
    const T a = pi / 6.;
    const point<T, 3> z{std::sin(a), 0., std::cos(a)}, x{std::cos(a), 0., -std::sin(a)};
    const auto big = box(m, 6., 6., 3.);
    const auto cutter = std::get<SolidId>(copy_shape(m, ShapeId{big}, Placement<T>::frame(point<T, 3>{0., 0., 1.5} - 3. * x - 3. * Y, x, Y, z)));
    const auto cu = unwrap(cut(m, ShapeId{cyl}, ShapeId{cutter}));
    REQUIRE(cu.solids.size() == 1u);
    const auto s = cu.solids.front();
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.5 - std::cos(1.) * std::tan(a)});
    const T c = 0.2;
    const auto g = unwrap(fillet_chordal_geometry(m, s, e, c));
    const auto &ed = m.edge(e);
    T rmin = 1e300, rmax = 0.;
    for (int i = 0; i <= 32; ++i)
    {
        const T t = ed.u1 + (ed.u2 - ed.u1) * i / 32.;
        const auto p1 = g.contacts[0].curve->value(t), p2 = g.contacts[1].curve->value(t);
        ASSERT_NEAR(norm(p1 - p2), c, 1e-6);
        const T r = norm(p1 - g.spine->value(t));
        rmin = std::min(rmin, r), rmax = std::max(rmax, r);
    }
    ASSERT_GT(rmax - rmin, 0.01);
    // (the solid: the rim of this boolean is not seen as a chain, its pieces not tangent at their
    // vertices within 1e-6: the fillet of constant radius refused as well)
}

TEST(tests_brep_blend_sections, refusals)
{
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    auto says = [](const auto &r, const std::string &what) { return !r && r.error().message.find(what) != std::string::npos; };
    BoxEdge b;
    const auto n = b.m.ids<FaceId>().size();
    for (const T rho : {0., 1., -0.2, 1.5})
    {
        FilletOptions<T> o;
        o.section = BlendSection<T>::conic(rho);
        const auto r = fillet(b.m, b.s, std::vector{b.e}, 0.3, o);
        ASSERT_EQ(code(r), BuildErrc::InvalidBounds);
        ASSERT_TRUE(says(r, "the rho of a conic section must be within ]0, 1["));
    }
    // a conic whose radius law kinks (lunes of the ball), a law of rho on an edge not filleted
    FilletOptions<T> o;
    o.section = BlendSection<T>::conic(0.6);
    const auto kinked = fillet(b.m, b.s, std::vector{b.e}, BlendLaw<T>{{{b.e, 0., 0.1}, {b.e, 0.5, 0.3}, {b.e, 1., 0.15}}, BlendVariation::Linear}, o);
    ASSERT_EQ(code(kinked), BuildErrc::NonGeneric);
    ASSERT_TRUE(says(kinked, "a conic fillet whose radius law kinks"));
    o.section = BlendSection<T>::conic(BlendLaw<T>{{{b.e, 0., 0.3}, {edge_near(b.m, b.s, {1.5, 0., 1.}), 0., 0.4}}});
    const auto elsewhere = fillet(b.m, b.s, std::vector{b.e}, 0.3, o);
    ASSERT_EQ(code(elsewhere), BuildErrc::InvalidBounds);
    ASSERT_TRUE(says(elsewhere, "a point of the law on an edge not filleted"));
    ASSERT_EQ(b.m.ids<FaceId>().size(), n);
    // conic fillets of a chain whose radius law kinks at a vertex of the chain (the rim of a
    // stadium slab, a hollow of the law there): the lune of the ball, refused
    {
        const T L = 3., R = 1., h = 0.8;
        Model<T> m;
        const auto s = oblong(m, L, R, h);
        const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
        const std::vector<point<T, 3>> corners{{0., -R, h}, {L, -R, h}, {L, R, h}, {0., R, h}};
        const std::vector<T> value{0.2, 0.1, 0.2, 0.2};
        BlendLaw<T> law{{}, BlendVariation::Linear};
        for (auto e : rim)
            for (const T at : {0., 1.})
            {
                const auto p = m.edge(e).curve->value(at == 0. ? m.edge(e).u1 : m.edge(e).u2);
                for (std::size_t i = 0; i < 4; ++i)
                    if (norm(p - corners[i]) < 1e-9)
                        law.points.push_back({e, at, value[i]});
            }
        FilletOptions<T> oc;
        oc.section = BlendSection<T>::conic(0.6);
        const auto r = fillet(m, s, rim, law, oc);
        ASSERT_EQ(code(r), BuildErrc::NonGeneric);
        ASSERT_TRUE(says(r, "conic fillets of a chain whose radius law kinks at their vertex"));
    }
}
