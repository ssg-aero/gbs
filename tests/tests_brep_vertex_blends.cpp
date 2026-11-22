// Blends meeting at a vertex (stage 6, PR 1, docs/sources/design/brep_palier6_vertex_blends.md):
// mitres of two blends at a vertex of three edges (fillets and chamfers, convex and concave),
// three chamfers at a corner, a free end at a vertex of more than three edges, the sphere of the
// n fillets at the top of a pyramid; the volumes against references worked out on their own.
// Vertex blends (PR 3): radii that differ at a corner, the limit toward the sphere, the top of a
// pyramid, setbacks (a suitcase corner), convex and concave fillets, a radius dropping toward the
// corner; tangency along every side, the curvature inside, refusals.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/sew.h>
#include <gbs-brep/split.h>
#include <gbs-brep/sweep.h>

#include <array>
#include <cmath>
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

    /// A planar face on the polygon p (its plane through p[0], normal along (p1 - p0) x (p2 - p0)).
    FaceId planar(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        auto x = p[1] - p[0];
        x = x / norm(x);
        auto n = cross(p[1] - p[0], p[2] - p[0]);
        n = n / norm(n);
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n, x}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    }

    SolidId box(Model<T> &m, T a, T b, T c, const point<T, 3> &o = {})
    {
        return solid_of(m, make_prism(m, planar(m, {o, o + point<T, 3>{a, 0., 0.}, o + point<T, 3>{a, b, 0.}, o + point<T, 3>{0., b, 0.}}), point<T, 3>{0., 0., c}));
    }

    /// The pyramid on the square [-1, 1]^2 at z = 0, its top at (0, 0, h).
    SolidId pyramid(Model<T> &m, T h)
    {
        const std::array<point<T, 3>, 4> b{{{-1., -1., 0.}, {1., -1., 0.}, {1., 1., 0.}, {-1., 1., 0.}}};
        const point<T, 3> top{0., 0., h};
        std::vector<FaceId> faces{planar(m, {b[3], b[2], b[1], b[0]})};
        for (std::size_t i = 0; i < 4; ++i)
            faces.push_back(planar(m, {b[i], b[(i + 1) % 4], top}));
        const auto r = unwrap(sew(m, faces, {.tol = 1e-9}));
        REQUIRE(r.shells.size() == 1u);
        auto s = unwrap(make_solid(m, r.shells.front()));
        if (signed_volume(m, m.solid(s).outer, 64) < 0.)
        {
            for (auto &fu : m.shell(r.shells.front()).faces)
                fu.orient = reverse(fu.orient);
        }
        return s;
    }

    /// The solid of the planar faces given (their polygons), sewn, oriented outward.
    SolidId polyhedron(Model<T> &m, const std::vector<std::vector<point<T, 3>>> &polygons)
    {
        std::vector<FaceId> faces;
        for (const auto &p : polygons)
            faces.push_back(planar(m, p));
        const auto r = unwrap(sew(m, faces, {.tol = 1e-9}));
        REQUIRE(r.shells.size() == 1u);
        auto s = unwrap(make_solid(m, r.shells.front()));
        if (signed_volume(m, m.solid(s).outer, 64) < 0.)
            for (auto &fu : m.shell(r.shells.front()).faces)
                fu.orient = reverse(fu.orient);
        return s;
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

    void expect_solid(const Model<T> &m, SolidId s, std::size_t faces, T v, T eps)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
        ASSERT_EQ(n_faces(m, s), faces);
        ASSERT_NEAR(volume(m, s), v, eps);
    }

    /// Along the edges of the face fc: the largest gap between its surface and those of the faces next to it, and the largest sine of the angle of their normals.
    std::array<T, 2> seams(const Model<T> &m, SolidId s, FaceId fc)
    {
        std::array<T, 2> worst{};
        const auto faces = explore<FaceId>(m, ShapeId{s});
        auto uv_of = [&](FaceId f, EdgeId e) -> std::shared_ptr<Curve<T, 2>> {
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (ce.edge == e)
                        return ce.pcurve;
            return nullptr;
        };
        for (auto w : m.face(fc).wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                const auto &ed = m.edge(ce.edge);
                for (auto f2 : faces)
                {
                    if (f2 == fc)
                        continue;
                    const auto p2 = uv_of(f2, ce.edge);
                    if (!p2)
                        continue;
                    const auto &s1 = *m.face(fc).surface, &s2 = *m.face(f2).surface;
                    for (int i = 0; i <= 64; ++i)
                    {
                        const T u = ed.u1 + (ed.u2 - ed.u1) * T(i) / 64.;
                        const auto a = ce.pcurve->value(u), b = p2->value(u);
                        worst[0] = std::max(worst[0], norm(s1.value(a[0], a[1]) - s2.value(b[0], b[1])));
                        const auto n1 = normal_direction(s1, a[0], a[1]), n2 = normal_direction(s2, b[0], b[1]);
                        worst[1] = std::max(worst[1], norm(cross(n1, n2)));
                    }
                }
            }
        return worst;
    }

    /// The face of the solid on a B-spline surface (the vertex blend).
    std::vector<FaceId> spline_faces(const Model<T> &m, SolidId s)
    {
        std::vector<FaceId> r;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (dynamic_cast<const BSSurface<T, 3> *>(m.face(f).surface.get()))
                r.push_back(f);
        return r;
    }

    /// The largest principal curvature of the face fc (its surface, within the polygon of its wire, convex) on a grid of barycentric points, the sides left out (a layer where the curvature jumps, G1 only).
    T largest_curvature(const Model<T> &m, FaceId fc)
    {
        std::vector<point<T, 2>> poly;
        for (auto w : m.face(fc).wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                const auto &ed = m.edge(ce.edge);
                poly.push_back(ce.pcurve->value(ce.orient == Orientation::Forward ? ed.u1 : ed.u2));
            }
        point<T, 2> o{};
        for (const auto &p : poly)
            o = o + p / T(poly.size());
        const auto &srf = *m.face(fc).surface;
        T worst{};
        for (std::size_t i = 0; i < poly.size(); ++i)
            for (int a = 0; a <= 12; ++a)
                for (int b = 0; a + b <= 11; ++b)
                {
                    const T x = T(a) / 12., y = T(b) / 12.;
                    const auto uv = o + x * (poly[i] - o) + y * (poly[(i + 1) % poly.size()] - o);
                    const auto su = srf.value(uv[0], uv[1], 1, 0), sv = srf.value(uv[0], uv[1], 0, 1);
                    const auto suu = srf.value(uv[0], uv[1], 2, 0), suv = srf.value(uv[0], uv[1], 1, 1), svv = srf.value(uv[0], uv[1], 0, 2);
                    auto n = cross(su, sv);
                    n = n / norm(n);
                    const T E = su * su, F = su * sv, G = sv * sv, L = suu * n, M = suv * n, N = svv * n;
                    const T K = (L * N - M * M) / (E * G - F * F), H = (E * N - 2. * F * M + G * L) / (2. * (E * G - F * F));
                    const T kk = std::abs(H) + std::sqrt(std::max(T(0), H * H - K));
                    worst = std::max(worst, kk);
                }
        return worst;
    }

    /// Gauss-Legendre on [a, b], 20 points on each of n panels.
    template <typename F>
    T gauss(F &&f, T a, T b, int n = 16)
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

    /// The depth a fillet of radius r takes off at the distance q from one of its faces (q < r), from the other.
    T depth(T r, T q)
    {
        return r - std::sqrt(std::max(T(0), r * r - (r - q) * (r - q)));
    }

    /**
     * The volume two fillets of radii r1 and r2 both take off at a corner where their edges meet
     * at right angles (each with its own other face, the face between them common): at the
     * distance q from the common face, the rectangle of the depths of both; q = u^2 (the depths
     * grow as sqrt(q) from 0).
     */
    T both(T r1, T r2)
    {
        const T q = std::min(r1, r2);
        return gauss([&](T u) { return depth(r1, u * u) * depth(r2, u * u) * 2. * u; }, 0., std::sqrt(q), 32);
    }

    T quarter(T r)
    {
        return (1. - pi / 4.) * r * r;
    }
}

TEST(tests_brep_vertex_blends, mitre_two_fillets)
{
    // the box [0, 2] x [0, 1] x [0, 1.5]; the edges along x (y = 0, z = 1.5) and along z (x = 2,
    // y = 0) at the corner (2, 0, 1.5), the edge along y there left sharp: the fillets meet in a
    // mitre on the face y = 0, the wider one cut by the face beyond the sharp edge
    for (const auto [r1, r2] : {std::array{0.2, 0.2}, std::array{0.15, 0.3}, std::array{0.3, 0.15}, std::array{0.25, 0.1}})
    {
        INFO("r1 = " << r1 << ", r2 = " << r2);
        Model<T> m;
        const auto s = box(m, 2., 1., 1.5);
        const auto ex = edge_near(m, s, {1., 0., 1.5}), ez = edge_near(m, s, {2., 0., 0.75});
        RadiusLaw<T> law{{{ex, 0., r1}, {ex, 1., r1}, {ez, 0., r2}, {ez, 1., r2}}, RadiusVariation::Linear};
        const auto f = unwrap(fillet(m, s, std::vector{ex, ez}, law));
        expect_solid(m, f, 8, 3. - quarter(r1) * 2. - quarter(r2) * 1.5 + both(r1, r2), 1e-9);
    }
}

TEST(tests_brep_vertex_blends, mitre_two_chamfers)
{
    // two chamfers of the distance d at the corner (1, 0, 1) of the unit cube: the triangles
    // q + w < d and q + p < d (q from the common face) overlap on (d - q)^2 at q
    const T d = 0.2;
    Model<T> m;
    const auto s = box(m, 1., 1., 1.);
    const auto f = unwrap(chamfer(m, s, std::vector{edge_near(m, s, {0.5, 0., 1.}), edge_near(m, s, {1., 0., 0.5})}, d));
    expect_solid(m, f, 8, 1. - d * d + d * d * d / 3., 1e-12);
}

TEST(tests_brep_vertex_blends, mitre_three_chamfers)
{
    // three chamfers of the distance d at a corner of the unit cube: they meet at (1 - d/2)^3;
    // the three prisms d^2 / 2, less their overlaps two by two d^3 / 3, plus the overlap of the
    // three d^3 / 4
    const T d = 0.2;
    Model<T> m;
    const auto s = box(m, 1., 1., 1.);
    const std::vector<EdgeId> edges{edge_near(m, s, {0.5, 1., 1.}), edge_near(m, s, {1., 0.5, 1.}), edge_near(m, s, {1., 1., 0.5})};
    const auto f = unwrap(chamfer(m, s, edges, d));
    expect_solid(m, f, 9, 1. - 1.5 * d * d + 0.75 * d * d * d, 1e-12);
    // the three faces of the chamfers meet at one vertex
    bool found = false;
    for (auto v : explore<VertexId>(m, ShapeId{f}))
        found = found || norm(m.vertex(v).pnt - point<T, 3>{1. - d / 2., 1. - d / 2., 1. - d / 2.}) < 1e-9;
    ASSERT_TRUE(found);
}

TEST(tests_brep_vertex_blends, mitre_concave)
{
    // a pocket [1, 2]^3 in the block [0, 2]^3: two of the three concave edges at its inner
    // corner filleted; the fillets add what they take off at a convex corner
    const T r = 0.2;
    Model<T> m;
    const auto block = box(m, 2., 2., 2.);
    const auto pocket = box(m, 1., 1., 1., {1., 1., 1.});
    const auto s = unwrap(cut(m, ShapeId{block}, ShapeId{pocket})).solids.front();
    const auto n0 = n_faces(m, s);
    const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {1.5, 1., 1.}), edge_near(m, s, {1., 1.5, 1.})}, r));
    expect_solid(m, f, n0 + 2, 7. + 2. * quarter(r) - both(r, r), 1e-9);
}

TEST(tests_brep_vertex_blends, free_end_at_four_edges)
{
    // the top of the box [0, 2] x [0, 1] x [0, 1] cut along its diagonal from (2, 0, 1): the vertical
    // edge at (2, 0) ends at a vertex of four edges, cut by both halves of the top (one plane)
    const T r = 0.15;
    Model<T> m;
    const auto s = box(m, 2., 1., 1.);
    FaceId top{};
    for (auto fc : explore<FaceId>(m, ShapeId{s}))
    {
        const auto *pl = dynamic_cast<const Plane<T> *>(m.face(fc).surface.get());
        if (pl && std::abs(pl->value(0., 0.)[2] - 1.) < 1e-12 && std::abs(std::abs(normal_direction(*pl, 0., 0.)[2]) - 1.) < 1e-12)
            top = fc;
    }
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
    (void)unwrap(split_face(m, top, std::vector<FaceCut<T>>{{d, pc->pcurve}}));
    REQUIRE(check(m, ShapeId{s}).ok());
    const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {2., 0., 0.5})}, r));
    expect_solid(m, f, 8, 2. - quarter(r), 1e-9);
}

TEST(tests_brep_vertex_blends, rounded_pyramid)
{
    // the eight edges of a pyramid on a square: four spheres of three arcs at the base, one of
    // four arcs at the top (the four centre lines meet on the axis). The pyramid has a sphere
    // touching its five faces (radius rho): the faces moved in by r bound the pyramid shrunk by
    // (rho - r) / rho about its centre, and the rounded pyramid is that one grown by r (Steiner:
    // volume, area, edges times their outer angles, the ball)
    const T h = 1.5, r = 0.1;
    Model<T> m;
    const auto s = pyramid(m, h);
    REQUIRE(check(m, ShapeId{s}).ok());
    ASSERT_NEAR(volume(m, s), 4. * h / 3., 1e-12);
    const auto all = explore<EdgeId>(m, ShapeId{s});
    ASSERT_EQ(all.size(), 8u);
    const auto f = unwrap(fillet(m, s, all, r));
    const T l = std::sqrt(h * h + 1.), rho = h / (1. + l), k = (rho - r) / rho;
    const T V = 4. * h / 3., A = 4. + 4. * l;
    const T base = 4. * 2. * std::acos(-1. / l), side = 4. * std::sqrt(2. + h * h) * std::acos(1. / (h * h + 1.));
    const T v = k * k * k * V + r * k * k * A + r * r * k * (base + side) / 2. + 4. / 3. * pi * r * r * r;
    expect_solid(m, f, 18, v, 1e-8);
    // the sphere at the top: of radius r, centred on the axis at rho from the faces of the pyramid shrunk
    int spheres = 0;
    for (auto fc : explore<FaceId>(m, ShapeId{f}))
        if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(m.face(fc).surface.get()))
        {
            ++spheres;
            ASSERT_NEAR(sp->radius(), r, 1e-12);
        }
    ASSERT_EQ(spheres, 5);
}

TEST(tests_brep_vertex_blends, free_end_at_the_top_of_a_pyramid)
{
    // one side edge of the pyramid filleted: at the top, a vertex of four edges, the fillet cut by
    // the two faces across (meeting on the edge between them, cut short there)
    const T h = 1.5, r = 0.05;
    Model<T> m;
    const auto s = pyramid(m, h);
    const auto f = unwrap(fillet(m, s, std::vector{edge_near(m, s, {0.5, 0.5, h / 2.})}, r));
    const auto report = check(m, ShapeId{f});
    ASSERT_TRUE(report.ok());
    ASSERT_EQ(n_faces(m, f), 6u);
    ASSERT_LT(volume(m, f), 4. * h / 3.);
    ASSERT_GT(volume(m, f), 4. * h / 3. - quarter(r) * std::sqrt(2. + h * h));
}

TEST(tests_brep_vertex_blends, refusals)
{
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    // two side edges of a pyramid at its top, a vertex of four edges: no mitre there
    {
        Model<T> mp;
        const auto py = pyramid(mp, 1.5);
        const auto np = mp.ids<FaceId>().size();
        const auto rp = fillet(mp, py, std::vector{edge_near(mp, py, {0.5, 0.5, 0.75}), edge_near(mp, py, {-0.5, 0.5, 0.75})}, 0.05);
        ASSERT_EQ(code(rp), BuildErrc::NonGeneric);
        ASSERT_NE(rp.error().message.find("more than three edges"), std::string::npos);
        ASSERT_EQ(mp.ids<FaceId>().size(), np);
    }
    // a pyramid on a dart (its base point (0.2, -0.1) inside): the side edge to (-1, 1) ends at
    // the top, a vertex of four edges, where its fillet meets the edges round it out of turn
    {
        Model<T> md;
        const std::array<point<T, 3>, 4> b{point<T, 3>{-1., -1., 0.}, {1., -1., 0.}, {0.2, -0.1, 0.}, {-1., 1., 0.}};
        const point<T, 3> top{0., 0.4, 0.6};
        std::vector<std::vector<point<T, 3>>> faces{{b[3], b[2], b[1], b[0]}};
        for (std::size_t i = 0; i < 4; ++i)
            faces.push_back({b[i], b[(i + 1) % 4], top});
        const auto dart = polyhedron(md, faces);
        const auto nd = md.ids<FaceId>().size();
        const auto rd = fillet(md, dart, std::vector{edge_near(md, dart, (b[3] + top) / 2.)}, 0.1);
        ASSERT_EQ(code(rd), BuildErrc::NonGeneric);
        ASSERT_NE(rd.error().message.find("out of turn"), std::string::npos);
        ASSERT_EQ(md.ids<FaceId>().size(), nd);
    }
    // three chamfers at a corner of a box, 0.05 along x and y, 0.6 along z (not by `chamfer`,
    // one distance for all: its geometry per edge): their common point (0.7, 0.7, 1.25) lies off
    // the chamfers along x and y
    {
        Model<T> mc;
        const auto cube = box(mc, 1., 1., 1.);
        const std::vector<EdgeId> at{edge_near(mc, cube, {0.5, 1., 1.}), edge_near(mc, cube, {1., 0.5, 1.}), edge_near(mc, cube, {1., 1., 0.5})};
        const std::map<EdgeId, T> dist{{at[0], 0.05}, {at[1], 0.05}, {at[2], 0.6}};
        const auto nc = mc.ids<FaceId>().size();
        std::map<EdgeId, EdgeId> original;
        const auto rc = blend_detail::guarded(mc, "chamfer", [&] {
            return blend_detail::blend(
                mc, cube, std::span<const EdgeId>{at},
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    if (original.empty())
                        for (auto x : at)
                            original[c.copy(x)] = x;
                    const T d = dist.at(original.at(e));
                    return fillet_detail::chamfer(mc, ns, e, shell_edge_uses(mc, mc.solid(ns).outer).at(e).front().face, d, d, oe);
                },
                FilletOptions<T>{}, "chamfer", false);
        });
        ASSERT_EQ(code(rc), BuildErrc::NonGeneric);
        ASSERT_NE(rc.error().message.find("common point is not within each"), std::string::npos);
        ASSERT_EQ(mc.ids<FaceId>().size(), nc);
    }
    // a convex and a concave edge at the inner top corner of an L: no mitre of the two
    {
        Model<T> ml;
        const auto l = solid_of(ml, make_prism(ml, planar(ml, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}), point<T, 3>{0., 0., 1.}));
        const auto nl = ml.ids<FaceId>().size();
        const auto rl = fillet(ml, l, std::vector{edge_near(ml, l, {1., 1., 0.5}), edge_near(ml, l, {1.5, 1., 1.})}, 0.2);
        ASSERT_EQ(code(rl), BuildErrc::NonGeneric);
        ASSERT_NE(rl.error().message.find("convex and concave"), std::string::npos);
        ASSERT_EQ(ml.ids<FaceId>().size(), nl);
    }
}

TEST(tests_brep_vertex_blends, mitre_obtuse)
{
    // a prism on a quadrilateral whose angle at B is 120 or 150 degrees: the top edge AB and the
    // vertical edge at B filleted, the other top edge at B left sharp. At 150 degrees the other
    // contact of the top fillet meets that edge 1.73 r past B, beyond the blend continued by one
    // radius: continued further as at a free end. The same solid whatever the extension asked.
    for (const T deg : {120., 150.})
    {
        INFO("angle " << deg);
        const T a = (180. - deg) * pi / 180.;
        Model<T> m;
        const auto s = solid_of(m, make_prism(m, planar(m, {{0., 0., 0.}, {2., 0., 0.}, {2. + std::cos(a), std::sin(a), 0.}, {0., 1.5, 0.}}), point<T, 3>{0., 0., 1.}));
        const std::vector<EdgeId> edges{edge_near(m, s, {1., 0., 1.}), edge_near(m, s, {2., 0., 0.5})};
        const auto f = unwrap(fillet(m, s, edges, 0.1));
        FilletOptions<T> o;
        o.extension = 3.;
        const auto g = unwrap(fillet(m, s, edges, 0.1, o));
        expect_solid(m, f, 8, volume(m, g), 1e-9);
        ASSERT_LT(volume(m, f), volume(m, s));
    }
}

TEST(tests_brep_vertex_blends, chain_through_a_vertex)
{
    // two edges going on one into the other at a vertex (a chain through it) among the blended
    // edges there: found whatever their order (`chain_through`, refused by `blend`). With planar
    // faces no blendable edge meets such a vertex (the faces on either side of a straight chain
    // hold the line of the chain: an edge between two of them is flat), hence bare edges here.
    Model<T> m;
    const auto v = unwrap(make_vertex(m, {0., 0., 0.}));
    const auto a = unwrap(make_vertex(m, {-1., 0., 0.})), b = unwrap(make_vertex(m, {1., 0., 0.})), c = unwrap(make_vertex(m, {0., 1., 0.}));
    const auto e1 = unwrap(make_edge(m, a, v)), e2 = unwrap(make_edge(m, v, b)), e2r = unwrap(make_edge(m, b, v)), e3 = unwrap(make_edge(m, v, c));
    for (const auto &edges : {std::vector{e1, e2, e3}, std::vector{e3, e1, e2}, std::vector{e3, e2, e1}, std::vector{e2r, e3, e1}})
        ASSERT_TRUE(blend_detail::chain_through(m, std::span<const EdgeId>{edges}, v));
    for (const auto &edges : {std::vector{e1, e3}, std::vector{e3, e2}, std::vector{e3}})
        ASSERT_FALSE(blend_detail::chain_through(m, std::span<const EdgeId>{edges}, v));
}

TEST(tests_brep_vertex_blends, vertex_blend_three_radii)
{
    // the corner (3, 3, 3) of the box [0, 3]^3, its three edges filleted at 0.5, 0.7 and 0.9: no
    // sphere, a vertex blend: tangent to the fillets and the faces round it within the
    // tolerances, its curvature bounded (no wave: a few times that of the smallest fillet), the
    // same volume whichever edge takes which radius (up to the grid of the B-spline)
    const std::array<T, 3> radii{0.5, 0.7, 0.9};
    std::vector<T> volumes;
    for (int turn = 0; turn < 3; ++turn)
    {
        Model<T> m;
        const auto s = box(m, 3., 3., 3.);
        const auto ex = edge_near(m, s, {1.5, 3., 3.}), ey = edge_near(m, s, {3., 1.5, 3.}), ez = edge_near(m, s, {3., 3., 1.5});
        const RadiusLaw<T> law{{{ex, 0.5, radii[turn % 3]}, {ey, 0.5, radii[(turn + 1) % 3]}, {ez, 0.5, radii[(turn + 2) % 3]}}};
        const auto r = fillet(m, s, std::vector{ex, ey, ez}, law);
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        const auto report = check(m, ShapeId{*r});
        ASSERT_TRUE(report.ok());
        ASSERT_EQ(n_faces(m, *r), 10u);
        const auto vb = spline_faces(m, *r);
        ASSERT_EQ(vb.size(), 1u);
        const auto [gap, sine] = seams(m, *r, vb.front());
        ASSERT_LE(gap, 1e-5);
        ASSERT_LE(sine, 1e-5);
        ASSERT_LE(largest_curvature(m, vb.front()), 4. / 0.5);
        volumes.push_back(volume(m, *r));
    }
    // between the box less its three fillets taken apart and that with the corner less the sphere of the largest
    const T apart = 27. - 3. * (quarter(0.5) + quarter(0.7) + quarter(0.9));
    for (auto v : volumes)
    {
        ASSERT_GT(v, apart);
        ASSERT_NEAR(v, volumes.front(), 1e-5);
    }
}

TEST(tests_brep_vertex_blends, vertex_blend_near_sphere)
{
    // radii (r, r, r (1 + eps)) at the corner of a box: the vertex blend tends to the sphere as
    // eps goes to 0 (its volume, less the part eps takes off along the edges, to that of the
    // sphere), by ends at the crossings of the contacts where they would end close to one another
    const T r = 0.5;
    auto run = [&](T eps) {
        Model<T> m;
        const auto s = box(m, 3., 3., 3.);
        const auto ex = edge_near(m, s, {1.5, 3., 3.}), ey = edge_near(m, s, {3., 1.5, 3.}), ez = edge_near(m, s, {3., 3., 1.5});
        const RadiusLaw<T> law{{{ex, 0.5, r}, {ey, 0.5, r}, {ez, 0.5, r * (1. + eps)}}};
        const auto res = fillet(m, s, std::vector{ex, ey, ez}, law);
        INFO((res ? std::string{} : res.error().message));
        REQUIRE(res.has_value());
        REQUIRE(check(m, ShapeId{*res}).ok());
        const auto vb = spline_faces(m, *res);
        if (eps > 0.)
        {
            REQUIRE(vb.size() == 1u);
            const auto [gap, sine] = seams(m, *res, vb.front());
            CHECK(gap <= 1e-5);
            CHECK(sine <= 1e-5);
        }
        return volume(m, *res);
    };
    const T sphere = run(0.);
    // (the slope of the volume in eps, from the larger eps, where the blend is far from the sphere)
    const T slope = (run(0.03) - run(0.02)) / 0.01;
    for (const T eps : {1e-4, 1e-3, 3e-3, 6e-3, 1e-2, 0.1, 0.3})
    {
        const T v = run(eps);
        if (eps <= 1e-2)
            ASSERT_NEAR(v, sphere + slope * eps, 2e-4);
        ASSERT_LT(v, sphere);
    }
    ASSERT_NEAR(run(1e-4), sphere, 5e-5);
}

TEST(tests_brep_vertex_blends, vertex_blend_top_of_a_pyramid)
{
    // the four side edges of a pyramid at its top, filleted at four radii: a vertex blend of
    // four fillets (a vertex of four edges); at the base, free ends
    const T h = 1.5;
    Model<T> m;
    const auto s = pyramid(m, h);
    std::vector<EdgeId> sides;
    for (const auto &p : {point<T, 3>{0.5, 0.5, h / 2.}, {-0.5, 0.5, h / 2.}, {-0.5, -0.5, h / 2.}, {0.5, -0.5, h / 2.}})
        sides.push_back(edge_near(m, s, p));
    const RadiusLaw<T> law{{{sides[0], 0.5, 0.05}, {sides[1], 0.5, 0.07}, {sides[2], 0.5, 0.09}, {sides[3], 0.5, 0.11}}};
    const auto r = fillet(m, s, sides, law);
    INFO((r ? std::string{} : r.error().message));
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(check(m, ShapeId{*r}).ok());
    ASSERT_EQ(n_faces(m, *r), 10u);
    const auto vb = spline_faces(m, *r);
    ASSERT_EQ(vb.size(), 1u);
    const auto [gap, sine] = seams(m, *r, vb.front());
    ASSERT_LE(gap, 1e-5);
    ASSERT_LE(sine, 1e-5);
    ASSERT_LE(largest_curvature(m, vb.front()), 4. / 0.05);
    ASSERT_LT(volume(m, *r), 4. * h / 3.);
}

TEST(tests_brep_vertex_blends, vertex_blend_setbacks)
{
    // the corner of a box, its three edges filleted at 0.5, setbacks of 1.2 along each: a
    // suitcase corner (the ends on sections at 1.2 from the corner, joining curves on the faces
    // between them, the vertex blend of six sides); the same with the setbacks on two edges only,
    // the third ending where its contacts cross the others (a smaller blend: more material kept)
    std::vector<T> volumes;
    for (const int asked : {3, 2})
    {
        Model<T> m;
        const auto s = box(m, 3., 3., 3.);
        const point<T, 3> corner{3., 3., 3.};
        const std::vector<EdgeId> at{edge_near(m, s, {1.5, 3., 3.}), edge_near(m, s, {3., 1.5, 3.}), edge_near(m, s, {3., 3., 1.5})};
        VertexId v{};
        for (auto x : explore<VertexId>(m, ShapeId{s}))
            if (norm(m.vertex(x).pnt - corner) < 1e-9)
                v = x;
        FilletOptions<T> o;
        for (int i = 0; i < asked; ++i)
            o.setbacks.push_back({v, at[std::size_t(i)], 1.2});
        const auto r = fillet(m, s, at, 0.5, o);
        INFO((r ? std::string{} : r.error().message));
        ASSERT_TRUE(r.has_value());
        ASSERT_TRUE(check(m, ShapeId{*r}).ok());
        ASSERT_EQ(n_faces(m, *r), 10u);
        const auto vb = spline_faces(m, *r);
        ASSERT_EQ(vb.size(), 1u);
        const auto [gap, sine] = seams(m, *r, vb.front());
        ASSERT_LE(gap, 1e-5);
        ASSERT_LE(sine, 1e-5);
        ASSERT_LE(largest_curvature(m, vb.front()), 4. / 0.5);
        volumes.push_back(volume(m, *r));
    }
    ASSERT_LT(volumes[0], volumes[1]);
    ASSERT_LT(volumes[1], 27.);
}

TEST(tests_brep_vertex_blends, vertex_blend_mixed)
{
    // the inner top corner of an L (a prism): the concave vertical edge and the two convex top
    // edges at it filleted: a vertex blend of convex and concave fillets
    Model<T> m;
    const auto l = solid_of(m, make_prism(m, planar(m, {{0., 0., 0.}, {2., 0., 0.}, {2., 1., 0.}, {1., 1., 0.}, {1., 2., 0.}, {0., 2., 0.}}), point<T, 3>{0., 0., 1.}));
    const std::vector<EdgeId> at{edge_near(m, l, {1., 1., 0.5}), edge_near(m, l, {1.5, 1., 1.}), edge_near(m, l, {1., 1.5, 1.})};
    const RadiusLaw<T> law{{{at[0], 0.5, 0.2}, {at[1], 0.5, 0.1}, {at[2], 0.5, 0.1}}};
    const auto r = fillet(m, l, at, law);
    INFO((r ? std::string{} : r.error().message));
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(check(m, ShapeId{*r}).ok());
    const auto vb = spline_faces(m, *r);
    ASSERT_EQ(vb.size(), 1u);
    const auto [gap, sine] = seams(m, *r, vb.front());
    ASSERT_LE(gap, 1e-5);
    ASSERT_LE(sine, 1e-5);
    ASSERT_EQ(n_faces(m, *r), 12u);
}

TEST(tests_brep_vertex_blends, vertex_blend_radius_dropping)
{
    // the three edges at a corner of a box, their radius dropping linearly from 0.4 at their
    // other vertex to 0.15 at the corner (refused at stage 5: the balls next to the corner
    // cover its sphere): a vertex blend
    Model<T> m;
    const auto s = box(m, 3., 3., 3.);
    const point<T, 3> corner{3., 3., 3.};
    const std::vector<EdgeId> at{edge_near(m, s, {1.5, 3., 3.}), edge_near(m, s, {3., 1.5, 3.}), edge_near(m, s, {3., 3., 1.5})};
    RadiusLaw<T> law{{}, RadiusVariation::Linear};
    for (auto e : at)
    {
        const bool v1_at_corner = norm(m.vertex(m.edge(e).v1).pnt - corner) < 1e-9;
        law.points.push_back({e, 0., v1_at_corner ? 0.15 : 0.4});
        law.points.push_back({e, 1., v1_at_corner ? 0.4 : 0.15});
    }
    const auto r = fillet(m, s, at, law);
    INFO((r ? std::string{} : r.error().message));
    ASSERT_TRUE(r.has_value());
    ASSERT_TRUE(check(m, ShapeId{*r}).ok());
    ASSERT_EQ(n_faces(m, *r), 10u);
    const auto vb = spline_faces(m, *r);
    ASSERT_EQ(vb.size(), 1u);
    const auto [gap, sine] = seams(m, *r, vb.front());
    ASSERT_LE(gap, 1e-5);
    ASSERT_LE(sine, 1e-5);
    // less than the corner of constant radius 0.15 (the radius at the corner), more than that of 0.4
    ASSERT_LT(volume(m, *r), volume(m, unwrap(fillet(m, s, at, 0.15))));
    ASSERT_GT(volume(m, *r), volume(m, unwrap(fillet(m, s, at, 0.4))));
}

TEST(tests_brep_vertex_blends, vertex_blend_refusals)
{
    auto code = [](const auto &r) { return r ? BuildErrc{} : r.error().code; };
    Model<T> m;
    const auto s = box(m, 3., 3., 3.);
    const point<T, 3> corner{3., 3., 3.};
    const std::vector<EdgeId> at{edge_near(m, s, {1.5, 3., 3.}), edge_near(m, s, {3., 1.5, 3.}), edge_near(m, s, {3., 3., 1.5})};
    VertexId v{}, other{};
    for (auto x : explore<VertexId>(m, ShapeId{s}))
        if (norm(m.vertex(x).pnt - corner) < 1e-9)
            v = x;
        else if (norm(m.vertex(x).pnt - point<T, 3>{0., 3., 3.}) < 1e-9)
            other = x;
    const auto n = m.ids<FaceId>().size();
    auto with = [&](std::vector<Setback<T>> sb) {
        FilletOptions<T> o;
        o.setbacks = std::move(sb);
        return o;
    };
    // shorter than where the contacts cross, longer than the edge, not positive
    for (const T d : {0.2, 3.5})
    {
        const auto r = fillet(m, s, at, 0.5, with({{v, at[0], d}}));
        ASSERT_EQ(code(r), BuildErrc::NonGeneric);
        ASSERT_EQ(m.ids<FaceId>().size(), n);
    }
    ASSERT_EQ(code(fillet(m, s, at, 0.5, with({{v, at[0], -1.}}))), BuildErrc::InvalidBounds);
    // at the other vertex of the edge (a free end there), on an edge not filleted, on a chamfer
    ASSERT_EQ(code(fillet(m, s, at, 0.5, with({{other, at[0], 1.}}))), BuildErrc::NonGeneric);
    ASSERT_EQ(code(fillet(m, s, std::vector{at[1], at[2]}, 0.5, with({{v, at[0], 1.}}))), BuildErrc::InvalidId);
    ASSERT_EQ(code(chamfer(m, s, at, 0.3, with({{v, at[0], 1.}}))), BuildErrc::InvalidBounds);
    ASSERT_EQ(m.ids<FaceId>().size(), n);
}
