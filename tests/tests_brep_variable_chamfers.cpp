// Chamfers whose distances vary (stage 7, PR 1, docs/sources/design/brep_palier7_blend_sections.md,
// § 2 and § 3; note brep7_pr01_variable_chamfers.md): two laws, a law on both faces, distance and
// angle (between planes and on a face curved across the edge); ruled pieces and creases where the
// slope of a law jumps; chains, loops, the corner of three; the volumes against references worked
// out on their own (the triangles taken off, integrated along the edge); refusals with their
// messages.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <numbers>
#include <ranges>
#include <set>
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
}

TEST(tests_brep_variable_chamfers, two_laws_on_a_box_edge)
{
    // d1 on the top, d2 on the side: the triangle d1 d2 / 2 taken off at each abscissa. Linear
    // d1 with a constant d2: contacts straight, the ruled surface a hyperbolic paraboloid (one
    // piece, exact); a smooth law on both; a linear law with a point at the middle: two pieces, a
    // crease
    struct Case
    {
        BlendLaw<T> d1, d2;
        std::size_t pieces;
    };
    BoxEdge b;
    const std::vector<Case> cases{
        {BlendLaw<T>::linear(b.e, 0.1, 0.4), BlendLaw<T>::constant(b.e, 0.2), 1},
        {BlendLaw<T>{{{b.e, 0., 0.1}, {b.e, 0.5, 0.3}, {b.e, 1., 0.15}}, BlendVariation::Smooth}, BlendLaw<T>{{{b.e, 0., 0.35}, {b.e, 1., 0.15}}}, 1},
        {BlendLaw<T>{{{b.e, 0., 0.1}, {b.e, 0.5, 0.4}, {b.e, 1., 0.18}}, BlendVariation::Linear}, BlendLaw<T>::constant(b.e, 0.2), 2},
    };
    for (const auto &c : cases)
    {
        const auto g = unwrap(chamfer_geometry(b.m, b.s, b.e, b.top, c.d1, c.d2));
        ASSERT_EQ(std::max<std::size_t>(g.pieces.size(), 1), c.pieces);
        ASSERT_EQ(g.creases.size(), c.pieces - 1);
        // the contacts B-splines on the planes: the ruled surface built exactly, piece by piece
        ASSERT_TRUE(g.exact);
        ASSERT_EQ(g.deviation, 0.);
        // the contacts at the distances asked from the edge (within the accuracy of the contacts of a
        // chamfer, that of the intersections)
        const auto &ed = b.m.edge(b.e);
        const auto f1 = law_of(c.d1, 3.), f2 = law_of(c.d2, 3.);
        for (int i = 0; i <= 16; ++i)
        {
            const T t = ed.u1 + (ed.u2 - ed.u1) * i / 16.;
            const auto p = ed.curve->value(t);
            const T s = norm(p - ed.curve->value(ed.u1));
            ASSERT_NEAR(norm(g.contacts[0].curve->value(t) - p), f1(s), 1e-6);
            ASSERT_NEAR(norm(g.contacts[1].curve->value(t) - p), f2(s), 1e-6);
            ASSERT_NEAR(g.contacts[0].curve->value(t)[2], 1., 1e-9);
        }
        const auto r = unwrap(chamfer(b.m, b.s, std::vector{b.e}, b.top, c.d1, c.d2));
        const T off = gauss([&](T s) { return f1(s) * f2(s) / 2.; }, 0., 3.);
        expect_solid(b.m, r, 6 + c.pieces, 3. - off, 1e-9);
    }
}

TEST(tests_brep_variable_chamfers, one_law_on_both_faces)
{
    BoxEdge b;
    const BlendLaw<T> d{{{b.e, 0., 0.1}, {b.e, 1., 0.3}}, BlendVariation::Linear};
    const auto r = unwrap(chamfer(b.m, b.s, std::vector{b.e}, d));
    const auto f = law_of(d, 3.);
    expect_solid(b.m, r, 7, 3. - gauss([&](T s) { return f(s) * f(s) / 2.; }, 0., 3.), 1e-9);
}

TEST(tests_brep_variable_chamfers, distance_and_angle)
{
    // between two planes at right angles, d2 = d1 tan(angle); at 120 degrees (the prism on a
    // quadrilateral, its vertical edge at the origin), by the law of sines d2 = d1 sin(angle) /
    // sin(120 + angle), the triangle d1 d2 sin(120) / 2 taken off along the edge of height 2
    {
        BoxEdge b;
        const T d1 = 0.3, a = pi / 6.;
        const auto g = unwrap(chamfer_angle_geometry(b.m, b.s, b.e, b.top, d1, a));
        ASSERT_TRUE(g.exact); // a plane
        const auto r = unwrap(chamfer_angle(b.m, b.s, std::vector{b.e}, b.top, d1, a));
        expect_solid(b.m, r, 7, 3. - 3. * d1 * d1 * std::tan(a) / 2., 1e-9);
    }
    {
        Model<T> m;
        const T h = 2., c = std::cos(2. * pi / 3.), s3 = std::sin(2. * pi / 3.);
        const auto s = prism(m, {{0., 0., 0.}, {3., 0., 0.}, {3., 2. * s3, 0.}, {2. * c, 2. * s3, 0.}}, h);
        const auto e = edge_near(m, s, {0., 0., 1.});
        const auto f = face_at(m, s, {1., 0., 1.}, Y);
        const T d1 = 0.5, a = 40. * pi / 180.;
        const T d2 = d1 * std::sin(a) / std::sin(2. * pi / 3. + a);
        const auto g = unwrap(chamfer_angle_geometry(m, s, e, f, d1, a));
        ASSERT_NEAR(norm(g.contacts[1].curve->value(1.)[0] * X + g.contacts[1].curve->value(1.)[1] * Y), d2, 1e-9);
        const auto r = unwrap(chamfer_angle(m, s, std::vector{e}, f, d1, a));
        expect_solid(m, r, 7, 2. * s3 * (3. + 3. - 2. * c) / 2. * h - h * d1 * d2 * s3 / 2., 1e-9);
    }
}

TEST(tests_brep_variable_chamfers, distance_and_angle_on_a_curved_face)
{
    // a half cylinder of radius 1 lying on its flat face, its straight edge at x = 1; d1 on the
    // cylinder (curved across the edge), 30 degrees from its tangent plane at the contact (not
    // from the chord to the edge: 36 degrees for d1 = 0.3 so): the line from the contact meets
    // the flat face; the section taken off the triangle (edge, contacts) and the segment of the
    // circle under the chord, along 3
    Model<T> m;
    const ax2<T, 3> ax{{{0., 0., 0.}, -1. * Y, X}};
    const auto arc = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., ax, 0., pi)));
    const auto dia = unwrap(make_edge(m, point<T, 3>{-1., 0., 0.}, point<T, 3>{1., 0., 0.}));
    const auto base = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{arc, dia}))));
    const auto s = solid_of(m, make_prism(m, base, point<T, 3>{0., 3., 0.}));
    const auto e = edge_near(m, s, {1., 1.5, 0.});
    FaceId cyl{};
    for (auto f : explore<FaceId>(m, ShapeId{s}))
        if (!dynamic_cast<const Plane<T> *>(m.face(f).surface.get()))
            cyl = f;
    const T d1 = 0.3, a = pi / 6.;
    const auto g = unwrap(chamfer_angle_geometry(m, s, e, cyl, d1, a));
    const auto &ed = m.edge(e);
    for (int i = 0; i <= 8; ++i)
    {
        const T t = ed.u1 + (ed.u2 - ed.u1) * i / 8.;
        const auto p1 = g.contacts[0].curve->value(t), q = g.contacts[1].curve->value(t);
        const point<T, 3> n{p1[0], 0., p1[2]}; // the normal of the cylinder at p1
        ASSERT_NEAR(norm(p1 - ed.curve->value(t)), d1, 1e-6);
        ASSERT_NEAR(std::asin(std::abs((q - p1) * n) / norm(q - p1)), a, 1e-6);
        ASSERT_NEAR(q[2], 0., 1e-9);
    }
    // the section: p1 at the angle al on the circle (chord d1), q where the line at 30 degrees from
    // the tangent there meets z = 0
    const T al = 2. * std::asin(d1 / 2.);
    const point<T, 3> p1{std::cos(al), 0., std::sin(al)}, n{std::cos(al), 0., std::sin(al)}, u{std::sin(al), 0., -std::cos(al)};
    const auto dir = std::cos(a) * u - std::sin(a) * n;
    const auto q = p1 + (-p1[2] / dir[2]) * dir;
    const point<T, 3> c{1., 0., 0.};
    const T triangle = norm(cross(p1 - c, q - c)) / 2., segment = (al - std::sin(al)) / 2.;
    const auto r = unwrap(chamfer_angle(m, s, std::vector{e}, cyl, d1, a));
    expect_solid(m, r, 5, 1.5 * pi - 3. * (triangle + segment), 1e-6);
}

TEST(tests_brep_variable_chamfers, loop_of_a_cylinder)
{
    // the top rim of a cylinder of radius 1, the distance smooth round the loop (0.1, then 0.3
    // half way round): in each meridian plane the triangle d^2 / 2, its centroid d / 3 in from the
    // rim (Pappus slice by slice)
    Model<T> m;
    const auto s = cylinder(m, 1., 1.);
    const auto e = edge_near(m, s, {std::cos(1.), std::sin(1.), 1.});
    const BlendLaw<T> d{{{e, 0., 0.1}, {e, 0.5, 0.3}}};
    const auto r = unwrap(chamfer(m, s, std::vector{e}, d));
    std::vector<std::pair<T, T>> pts{{0., 0.1}, {pi, 0.3}};
    const auto f = radius_detail::radius_function(pts, BlendVariation::Smooth, true, 2. * pi);
    const T off = gauss([&](T th) { return f(th) * f(th) / 2. * (1. - f(th) / 3.); }, 0., 2. * pi);
    expect_solid(m, r, 4, pi - off, 1e-6);
}

TEST(tests_brep_variable_chamfers, corner_of_three)
{
    // the three edges at a corner of a box, a linear law on each: three chamfers meeting at a
    // point (mitre); a valid solid, nine faces, the three chamfers sharing one vertex; the
    // material taken off at most the sum of the three chamfers alone, more than the largest
    Model<T> m;
    const auto s = box(m, 2., 2., 2.);
    const std::vector<EdgeId> at{edge_near(m, s, {1., 2., 2.}), edge_near(m, s, {2., 1., 2.}), edge_near(m, s, {2., 2., 1.})};
    BlendLaw<T> d{{}, BlendVariation::Linear};
    for (auto e : at)
        d.points.push_back({e, 0., 0.15}), d.points.push_back({e, 1., 0.25});
    const auto r = chamfer(m, s, at, d);
    INFO((r ? std::string{} : r.error().message));
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(check(m, ShapeId{*r}).ok());
    ASSERT_EQ(n_faces(m, *r), 9u);
    // (the chamfers: the faces whose normal is along no axis; the same law on both faces makes
    // two coplanar contacts, a plane)
    std::map<VertexId, int> uses;
    for (auto f : explore<FaceId>(m, ShapeId{*r}))
    {
        const auto *pl = dynamic_cast<const Plane<T> *>(m.face(f).surface.get());
        const auto nz = pl ? normal_direction(*pl, 0., 0.) : point<T, 3>{};
        if (pl && std::abs(std::abs(nz[0]) + std::abs(nz[1]) + std::abs(nz[2]) - 1.) < 1e-9)
            continue;
        for (auto v : std::set<VertexId>{std::from_range, explore<VertexId>(m, ShapeId{f})})
            ++uses[v];
    }
    ASSERT_EQ(std::ranges::count(uses | std::views::values, 3), 1);
    T alone{}, largest{};
    for (auto e : at)
    {
        BlendLaw<T> one{{}, d.variation}; // the points of the law on e
        for (const auto &p : d.points)
            if (p.edge == e)
                one.points.push_back(p);
        const T off = 8. - volume(m, unwrap(chamfer(m, s, std::vector{e}, one)));
        alone += off, largest = std::max(largest, off);
    }
    ASSERT_GE(volume(m, *r), 8. - alone);
    ASSERT_LT(volume(m, *r), 8. - largest);
}

TEST(tests_brep_variable_chamfers, chain_round_a_slab)
{
    // the top rim of a stadium slab, a loop of two lines and two half circles, the distance smooth
    // along the whole chain (0.1 in the middle of one straight side, 0.3 in the middle of the
    // other): continuous from edge to edge; in each normal plane the triangle d^2 / 2, its centroid
    // d / 3 in from the rim (the half circles: the slice weighted by 1 - d / (3 R))
    const T L = 3., R = 1., h = 0.8;
    Model<T> m;
    const auto s = oblong(m, L, R, h);
    const std::vector<EdgeId> rim{edge_near(m, s, {1.5, -R, h}), edge_near(m, s, {L + R, 0., h}), edge_near(m, s, {1.5, R, h}), edge_near(m, s, {-R, 0., h})};
    const BlendLaw<T> d{{{rim[0], 0.5, 0.1}, {rim[2], 0.5, 0.3}}};
    const auto r = unwrap(chamfer(m, s, rim, d));
    // the chain from the middle of the first side: L / 2, a half circle, L, a half circle, L / 2
    const T P = 2. * L + 2. * pi * R;
    std::vector<std::pair<T, T>> pts{{0., 0.1}, {P / 2., 0.3}};
    const auto f = radius_detail::radius_function(pts, BlendVariation::Smooth, true, P);
    auto kappa = [&](T x) { return (x > L / 2. && x < L / 2. + pi * R) || (x > 1.5 * L + pi * R && x < 1.5 * L + 2. * pi * R) ? 1. / R : 0.; };
    T off{};
    const std::vector<T> cuts{0., L / 2., L / 2. + pi * R, 1.5 * L + pi * R, 1.5 * L + 2. * pi * R, P};
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
        off += gauss([&](T x) { return f(x) * f(x) / 2. * (1. - kappa(x) * f(x) / 3.); }, cuts[i], cuts[i + 1]);
    expect_solid(m, r, 10, (2. * R * L + pi * R * R) * h - off, 1e-6);
}

TEST(tests_brep_variable_chamfers, refusals)
{
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    auto says = [](const auto &r, const std::string &what) { return !r && r.error().message.find(what) != std::string::npos; };
    BoxEdge b;
    const auto n = b.m.ids<FaceId>().size();
    const auto other = edge_near(b.m, b.s, {1.5, 0., 1.});
    // an angle out of ]0, 90[, a distance not positive, a law falling under the tolerance past a
    // vertex (0.1 to 0.3 over a fifth of the edge, continued by 0.3: 0 there), a point on an edge
    // not chamfered, two values at one point
    for (const T a : {pi / 2., 0.})
    {
        const auto r = chamfer_angle(b.m, b.s, std::vector{b.e}, b.top, 0.2, a);
        ASSERT_EQ(code(r), BuildErrc::InvalidBounds);
        ASSERT_TRUE(says(r, "the angle of a chamfer must be within ]0, 90[ degrees"));
    }
    const auto negative = chamfer(b.m, b.s, std::vector{b.e}, BlendLaw<T>::linear(b.e, 0.2, -0.1));
    ASSERT_EQ(code(negative), BuildErrc::InvalidBounds);
    ASSERT_TRUE(says(negative, "a distance of the law not larger than the tolerance"));
    const BlendLaw<T> steep{{{b.e, 0., 0.1}, {b.e, 0.2, 0.3}, {b.e, 1., 0.3}}, BlendVariation::Linear};
    const auto past = chamfer(b.m, b.s, std::vector{b.e}, steep);
    ASSERT_EQ(code(past), BuildErrc::InvalidBounds);
    ASSERT_TRUE(says(past, "FilletOptions::extension shorter"));
    {
        // the same law with a shorter extension
        BoxEdge b2;
        FilletOptions<T> shorter;
        shorter.extension = 0.5;
        ASSERT_TRUE(chamfer(b2.m, b2.s, std::vector{b2.e}, steep, shorter).has_value());
    }
    const auto elsewhere = chamfer(b.m, b.s, std::vector{b.e}, BlendLaw<T>{{{b.e, 0., 0.2}, {other, 0., 0.2}}});
    ASSERT_EQ(code(elsewhere), BuildErrc::InvalidBounds);
    ASSERT_TRUE(says(elsewhere, "a point of the law on an edge not chamfered"));
    const auto twice = chamfer(b.m, b.s, std::vector{b.e}, BlendLaw<T>{{{b.e, 0.5, 0.2}, {b.e, 0.5, 0.3}}});
    ASSERT_EQ(code(twice), BuildErrc::InvalidBounds);
    ASSERT_TRUE(says(twice, "two values of the distance at the same point of the edge"));
    // the face of d1 not a face of the edge
    const auto wrong = chamfer(b.m, b.s, std::vector{b.e}, face_at(b.m, b.s, {1.5, 0., 0.5}, Y), BlendLaw<T>::constant(b.e, 0.2), BlendLaw<T>::constant(b.e, 0.2));
    ASSERT_EQ(code(wrong), BuildErrc::InvalidId);
    ASSERT_TRUE(says(wrong, "the face of d1 is not a face of the edge"));
    // wider than the face at an angle near 90 degrees: the line misses the side within it
    ASSERT_EQ(code(chamfer_angle(b.m, b.s, std::vector{b.e}, b.top, 0.5, 89. * pi / 180.)), BuildErrc::NonGeneric);
    ASSERT_EQ(b.m.ids<FaceId>().size(), n);
}
