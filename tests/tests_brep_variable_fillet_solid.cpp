// Fillets whose radius follows a law along their edges, on the solid (stage 5, PR 2): fillet
// with a RadiusLaw, its pieces and lunes as faces, chains, loops, corners; the volume against
// the union of the balls, integrated on its own.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/sweep.h>

#include <blade_fixture.h>

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
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

    /// A prism of height h on the oblong of two half circles of radius R, centres (0, 0) and (L, 0).
    SolidId oblong(Model<T> &m, T L, T R, T h)
    {
        auto arc = [](T r, const point<T, 3> &c, T a1, T a2) {
            return std::make_shared<EllipseCurve<T>>(r, r, ax2<T, 3>{{c, Z, X}}, a1, a2);
        };
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
        std::size_t n = 0;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            n += dynamic_cast<const SphericalSurface<T> *>(m.face(f).surface.get()) != nullptr;
        return n;
    }

    void expect_solid(const Model<T> &m, SolidId s, std::size_t faces)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
        ASSERT_EQ(n_faces(m, s), faces);
    }

    /// Gauss-Legendre on [a, b], 20 points.
    template <typename F>
    T gauss(F &&f, T a, T b)
    {
        static const std::array<T, 10> x{0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
                                         0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949};
        static const std::array<T, 10> w{0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
                                         0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521};
        const T c = (a + b) / 2., h = (b - a) / 2.;
        T sum{};
        for (std::size_t i = 0; i < 10; ++i)
            sum += w[i] * (f(c - h * x[i]) + f(c + h * x[i]));
        return sum * h;
    }

    /**
     * The volume a convex fillet of radius r(x) takes off the edge y = 0, z = 1 of a box, from
     * x = 0 to x = L, its balls centred at (s, r(s), 1 - r(s)) for s in [s1, s2] (the blend
     * continued past the ends): in the plane x, the corner between y = 0, z = 1 and the top of the
     * union of the discs the balls leave there, from y = 0 to the contact y = r(x). The top at y:
     * the largest over s (sampled, then golden section) of 1 - r(s) + sqrt(r(s)^2 - (x - s)^2 -
     * (y - r(s))^2); y = r(x) w^2 (the top grows as sqrt(y) from y = 0); x in panels cut at the
     * kinks of r.
     */
    T removed_volume(const std::function<T(T)> &r, T L, T s1, T s2, std::vector<T> kinks)
    {
        auto top = [&](T x, T y) {
            auto g = [&](T s) {
                const T rs = r(s), d = rs * rs - (x - s) * (x - s) - (y - rs) * (y - rs);
                return d >= 0. ? 1. - rs + std::sqrt(d) : -1e300;
            };
            const T R = 1.5;
            T best = x, gb = g(x);
            const T lo = std::max(s1, x - R), hi = std::min(s2, x + R);
            for (int i = 0; i <= 400; ++i)
                if (const T s = lo + (hi - lo) * i / 400., v = g(s); v > gb)
                    gb = v, best = s;
            T a = std::max(lo, best - (hi - lo) / 400.), b = std::min(hi, best + (hi - lo) / 400.);
            const T phi = (std::sqrt(5.) - 1.) / 2.;
            for (int i = 0; i < 80; ++i)
            {
                const T c = b - phi * (b - a), d = a + phi * (b - a);
                (g(c) > g(d) ? b : a) = g(c) > g(d) ? d : c;
            }
            return std::max(gb, g((a + b) / 2.));
        };
        auto area = [&](T x) {
            const T rx = r(x);
            T a{};
            for (int k = 0; k < 8; ++k)
                a += gauss([&](T w) { return (1. - std::min(1., top(x, rx * w * w))) * 2. * rx * w; }, k / 8., (k + 1) / 8.);
            return a;
        };
        kinks.push_back(0.);
        kinks.push_back(L);
        std::ranges::sort(kinks);
        T v{};
        for (std::size_t i = 0; i + 1 < kinks.size(); ++i)
            for (int k = 0; k < 16; ++k)
            {
                const T a = kinks[i] + (kinks[i + 1] - kinks[i]) * k / 16., b = kinks[i] + (kinks[i + 1] - kinks[i]) * (k + 1) / 16.;
                v += gauss(area, a, b);
            }
        return v;
    }
}

TEST(tests_brep_variable_fillet_solid, box_edge_volume)
{
    // the edge y = 0, z = 1 of a box 2 x 1 x 1 along x, its radius r0 at x = 0, r1 at x = 0.8,
    // r2 at x = 2: the volume taken off against the union of the balls. Linear, a peak (0.1,
    // 0.3, 0.15): two pieces meeting along a crease; a hollow (0.3, 0.1, 0.25): two pieces and
    // the lune between them
    struct Case
    {
        RadiusVariation variation;
        std::array<T, 3> r;
        std::size_t faces, spheres;
    };
    for (const auto &c : {Case{RadiusVariation::Smooth, {0.1, 0.3, 0.15}, 7, 0}, Case{RadiusVariation::Linear, {0.1, 0.3, 0.15}, 8, 0},
                          Case{RadiusVariation::Smooth, {0.3, 0.1, 0.25}, 7, 0}, Case{RadiusVariation::Linear, {0.3, 0.1, 0.25}, 9, 1}})
    {
        INFO("smooth " << (c.variation == RadiusVariation::Smooth) << " r1 " << c.r[1]);
        Model<T> m;
        const auto s = box(m, 2., 1., 1.);
        const auto e = edge_near(m, s, {1., 0., 1.});
        const auto &ed = m.edge(e);
        const bool along = ed.curve->value(ed.u1)[0] < ed.curve->value(ed.u2)[0];
        // the law from the first vertex of the edge
        const RadiusLaw<T> law{{{e, along ? 0. : 1., c.r[0]}, {e, along ? 0.4 : 0.6, c.r[1]}, {e, along ? 1. : 0., c.r[2]}}, c.variation};
        FilletOptions<T> o;
        o.approx_tol = 1e-7;
        const auto f = unwrap(fillet(m, s, std::vector{e}, law, o));
        expect_solid(m, f, c.faces);
        ASSERT_EQ(n_spheres(m, f), c.spheres);
        // r(x), as the law
        const auto F = radius_detail::radius_function<T>({{0., c.r[0]}, {0.8, c.r[1]}, {2., c.r[2]}}, c.variation);
        const T d1 = F(0., 1), d2 = F(2., 1);
        auto r = [&](T x) { return x < 0. ? c.r[0] + d1 * x : x > 2. ? c.r[2] + d2 * (x - 2.) : F(x); };
        const T v = removed_volume(r, 2., -0.3, 2.3, {0.8});
        INFO(std::setprecision(12) << "volume " << volume(m, f) << " expected " << 2. - v);
        ASSERT_NEAR(volume(m, f), 2. - v, 2e-7);
    }
}

TEST(tests_brep_variable_fillet_solid, constant_law_is_constant_fillet)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto f = unwrap(fillet(m, s, std::vector{e}, RadiusLaw<T>{{{e, 0., 0.2}, {e, 1., 0.2}}}));
    const auto g = unwrap(fillet(m, s, std::vector{e}, 0.2));
    expect_solid(m, f, 7);
    ASSERT_NEAR(volume(m, f), volume(m, g), 1e-12);
}

TEST(tests_brep_variable_fillet_solid, cylinder_rim)
{
    // the top rim of a cylinder (a closed edge): a smooth law closed round it, one face; a
    // linear one, a piece between each two of its three points
    for (const auto variation : {RadiusVariation::Smooth, RadiusVariation::Linear})
    {
        INFO("smooth " << (variation == RadiusVariation::Smooth));
        Model<T> m;
        const auto s = cylinder(m, 1., 1.5);
        const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.5});
        const RadiusLaw<T> law{{{e, 0., 0.1}, {e, 0.25, 0.25}, {e, 0.6, 0.15}}, variation};
        const auto f = unwrap(fillet(m, s, std::vector{e}, law));
        // linear: three pieces, lunes at the hollows of the law (at the vertex and at 60 %), a
        // crease at its peak (25 %)
        expect_solid(m, f, variation == RadiusVariation::Linear ? 8u : 4u);
        ASSERT_EQ(n_spheres(m, f), variation == RadiusVariation::Linear ? 2u : 0u);
        // between the fillets of constant radius 0.1 and 0.25
        const auto lo = unwrap(fillet(m, s, std::vector{e}, 0.1)), hi = unwrap(fillet(m, s, std::vector{e}, 0.25));
        ASSERT_LT(volume(m, f), volume(m, lo));
        ASSERT_GT(volume(m, f), volume(m, hi));
    }
}

TEST(tests_brep_variable_fillet_solid, oblong_chain)
{
    // the top rim of an oblong prism, a closed chain of four tangent edges: the law along the
    // chain (by the arc length round it); smooth, one face per edge; linear with a point at a
    // vertex of the chain where its slope jumps, a lune between the fillets of two edges
    const T L = 2., R = 1., h = 0.8;
    for (const auto variation : {RadiusVariation::Smooth, RadiusVariation::Linear})
    {
        INFO("smooth " << (variation == RadiusVariation::Smooth));
        Model<T> m;
        const auto s = oblong(m, L, R, h);
        const std::vector<EdgeId> rim{edge_near(m, s, {1., -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1., R, h}), edge_near(m, s, {-R, 0., h})};
        const auto &e0 = m.edge(rim[0]);
        // 0.1 at the start of the first straight edge, 0.25 at its end, 0.15 at the middle of the
        // second straight edge
        const bool along = e0.curve->value(e0.u1)[0] < e0.curve->value(e0.u2)[0];
        const RadiusLaw<T> law{{{rim[0], along ? 0. : 1., 0.1}, {rim[0], along ? 1. : 0., 0.25}, {rim[2], 0.5, 0.15}}, variation};
        const auto f = unwrap(fillet(m, s, rim, law));
        // a lune or a crease wherever the law kinks, and where a straight edge meets a half
        // circle with the radius varying there (the centre line has a corner: its speed along
        // the edge drops from 1 to 1 - r / R); linear: the second straight edge in two pieces
        // smooth: lunes where the half circles meet the second straight edge, the slope of the
        // law not 0 there (0 at the ends of the first one, extrema of the law); linear: also at the
        // start of the first straight edge (a hollow) and in the middle of the second, a crease at
        // the end of the first (a peak)
        const std::size_t spheres = n_spheres(m, f);
        ASSERT_EQ(spheres, variation == RadiusVariation::Smooth ? 2u : 4u);
        expect_solid(m, f, (variation == RadiusVariation::Smooth ? 10u : 11u) + spheres);
        // without extension past the vertices: the blends continued past a vertex of the chain
        // all the same, room for the crease there
        FilletOptions<T> o0;
        o0.extension = 0.;
        const auto f0 = unwrap(fillet(m, s, rim, law, o0));
        expect_solid(m, f0, n_faces(m, f));
        ASSERT_NEAR(volume(m, f0), volume(m, f), 1e-6);
        // between the fillets of constant radius 0.1 and 0.25
        const auto lo = unwrap(fillet(m, s, rim, 0.1)), hi = unwrap(fillet(m, s, rim, 0.25));
        ASSERT_LT(volume(m, f), volume(m, lo));
        ASSERT_GT(volume(m, f), volume(m, hi));
    }
}

TEST(tests_brep_variable_fillet_solid, corner)
{
    // the three edges at a corner of a box, each smooth from 0.15 at its far vertex to 0.25 at
    // the three edges at a corner of a box, each smooth from 0.15 at its far vertex to 0.25 at
    // the corner: the ends of the fillets circles of the ball of the corner shifted away from it,
    // the face of the sphere bounded by them; the radius dropping toward the corner (0.25 to
    // 0.15): the balls next to the corner cover its sphere, a vertex blend (stage 6)

    for (const bool rising : {true, false})
    {
        INFO("rising " << rising);
        Model<T> m;
        const auto s = box(m, 1., 1., 1.);
        const std::vector<EdgeId> edges{edge_near(m, s, {0.5, 1., 1.}), edge_near(m, s, {1., 0.5, 1.}), edge_near(m, s, {1., 1., 0.5})};
        auto at_corner = [&](EdgeId e, T at) {
            const auto &ed = m.edge(e);
            return norm(m.vertex(ed.v2).pnt - point<T, 3>{1., 1., 1.}) < 1e-9 ? at : 1. - at;
        };
        RadiusLaw<T> law;
        for (auto e : edges)
        {
            law.points.push_back({e, at_corner(e, 0.), rising ? 0.15 : 0.25});
            law.points.push_back({e, at_corner(e, 1.), rising ? 0.25 : 0.15});
        }
        const auto f = fillet(m, s, edges, law);
        if (!rising)
        {
            ASSERT_TRUE(f.has_value());
            expect_solid(m, *f, 10);
            ASSERT_EQ(n_spheres(m, *f), 0u);
            continue;
        }
        expect_solid(m, *f, 10);
        ASSERT_EQ(n_spheres(m, *f), 1u);
        // less than the corner of constant radius 0.25, more than that of 0.15
        ASSERT_GT(volume(m, *f), volume(m, unwrap(fillet(m, s, edges, 0.25))));
        ASSERT_LT(volume(m, *f), volume(m, unwrap(fillet(m, s, edges, 0.15))));
    }
}

TEST(tests_brep_variable_fillet_solid, refusals)
{
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    const auto e = edge_near(m, s, {1., 0., 1.});
    const auto other = edge_near(m, s, {1., 1., 1.});
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    const auto n = m.ids<FaceId>().size();
    // a point on an edge not filleted, no point on an edge filleted
    ASSERT_EQ(code(fillet(m, s, std::vector{e}, RadiusLaw<T>{{{e, 0., 0.1}, {other, 0., 0.2}}})), BuildErrc::InvalidBounds);
    ASSERT_EQ(code(fillet(m, s, std::vector{e, other}, RadiusLaw<T>{{{e, 0., 0.1}, {e, 1., 0.2}}})), BuildErrc::InvalidBounds);
    // wider than a face
    ASSERT_EQ(code(fillet(m, s, std::vector{e}, RadiusLaw<T>::linear(e, 0.1, 1.5))), BuildErrc::NonGeneric);
    ASSERT_EQ(m.ids<FaceId>().size(), n);
}

TEST(tests_brep_variable_fillet_solid, blade_root)
{
    // the root of a blade fused to its hub, a closed edge round a closed B-spline face that is
    // not periodic: the loop of the fillet from the seam of the face; a smooth law closed round
    // the root, from 0.02 at its vertex to 0.04 halfway round
    Model<T> m;
    const auto [s, root] = blade_fixture::blade_on_hub(m);
    REQUIRE(root.size() == 1u);
    const auto e = root.front();
    const auto f = unwrap(fillet(m, s, root, RadiusLaw<T>{{{e, 0., 0.02}, {e, 0.5, 0.04}}}));
    // (the root closes at its vertex with a corner of 1.2e-3 rad, the radius varying there: the
    // balls on either side leave a thin lune, the fillet in two pieces)
    ASSERT_EQ(n_spheres(m, f), 1u);
    expect_solid(m, f, n_faces(m, s) + 3);
    // between the fillets of constant radius 0.02 and 0.04 (concave: adding material)
    const auto lo = unwrap(fillet(m, s, root, 0.02)), hi = unwrap(fillet(m, s, root, 0.04));
    ASSERT_GT(volume(m, f), volume(m, lo));
    ASSERT_LT(volume(m, f), volume(m, hi));
}

TEST(tests_brep_variable_fillet_solid, peak_on_an_acute_edge)
{
    // the edge of a wedge of 30 degrees, a linear law with a peak in the middle: the pieces
    // continued past the peak with their slopes reach 0.25 (1 + 1.5 x 0.23) = 0.34, beyond the
    // largest radius of the law; the supports and the reach for that radius (refused before:
    // "the balls of the fillet cannot be followed along the edge")
    Model<T> m;
    const T a = pi / 6.;
    const std::vector<point<T, 3>> p{{0., 0., 0.}, {0., 4., 0.}, {0., 4. * std::cos(a), 4. * std::sin(a)}};
    std::vector<EdgeId> es;
    for (std::size_t i = 0; i < 3; ++i)
        es.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % 3])));
    const auto face = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, es))));
    const auto s = solid_of(m, make_prism(m, face, point<T, 3>{2., 0., 0.}));
    const auto e = edge_near(m, s, {1., 0., 0.});
    const RadiusLaw<T> law{{{e, 0., 0.02}, {e, 0.5, 0.25}, {e, 1., 0.02}}, RadiusVariation::Linear};
    const auto g = unwrap(fillet_geometry(m, s, e, law));
    ASSERT_EQ(g.pieces.size(), 2u);
    ASSERT_EQ(g.creases.size(), 1u);
    const auto f = unwrap(fillet(m, s, std::vector{e}, law));
    expect_solid(m, f, 7);
    ASSERT_LT(volume(m, f), volume(m, s));
}

TEST(tests_brep_variable_fillet_solid, largest_radius_per_edge)
{
    // a stadium slab of half circles of radius 0.3, its rim at 0.15 but 0.31 in the middle of a
    // straight edge: the largest radius of each edge its own (the half circles checked for 0.15,
    // not for 0.31, beyond their radius of curvature)
    const T L = 3., R = 0.3, h = 0.8;
    Model<T> m;
    const auto s = oblong(m, L, R, h);
    const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
    const RadiusLaw<T> law{{{rim[0], 0.1, 0.15}, {rim[0], 0.5, 0.31}, {rim[0], 0.9, 0.15}, {rim[1], 0.5, 0.15}, {rim[2], 0.5, 0.15}, {rim[3], 0.5, 0.15}}};
    const auto f = unwrap(fillet(m, s, rim, law));
    expect_solid(m, f, n_faces(m, s) + 4 + n_spheres(m, f));
}

TEST(tests_brep_variable_fillet_solid, degenerate_edge)
{
    // the edge at a pole of a ball: refused as such (InvalidId), not "map::at"
    Model<T> m;
    const auto e1 = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., ax2<T, 3>{{{0., 0., 0.}, cross(X, Z), X}}, -pi / 2, pi / 2)));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., 1.}, point<T, 3>{0., 0., -1.}));
    const auto face = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                       unwrap(make_wire(m, std::vector<EdgeId>{e1, e2}))));
    const auto ball = solid_of(m, make_revolution(m, face, ax1<T, 3>{{{0., 0., 0.}, Z}}));
    std::optional<EdgeId> pole;
    for (auto e : explore<EdgeId>(m, ShapeId{ball}))
        if (m.edge(e).degenerate)
            pole = e;
    REQUIRE(pole.has_value());
    const auto r = fillet(m, ball, std::vector{*pole}, RadiusLaw<T>{{{*pole, 0., 0.1}, {*pole, 1., 0.2}}});
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::InvalidId);
}
