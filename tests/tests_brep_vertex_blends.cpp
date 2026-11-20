// Blends meeting at a vertex (stage 6, PR 1, docs/sources/design/brep_palier6_vertex_blends.md):
// mitres of two blends at a vertex of three edges (fillets and chamfers, convex and concave),
// three chamfers at a corner, a free end at a vertex of more than three edges, the sphere of the
// n fillets at the top of a pyramid; the volumes against references worked out on their own.
#include <doctest_gtest.hpp>
#include <gbs-brep/blend.h>
#include <gbs-brep/boolean.h>
#include <gbs-brep/sew.h>
#include <gbs-brep/split.h>
#include <gbs-brep/sweep.h>

#include <array>
#include <cmath>
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
    Model<T> m;
    const auto s = box(m, 1., 1., 1.);
    const auto n = m.ids<FaceId>().size();
    // three fillets of different radii at a corner: a vertex blend (PR 3)
    const auto ex = edge_near(m, s, {0.5, 1., 1.}), ey = edge_near(m, s, {1., 0.5, 1.}), ez = edge_near(m, s, {1., 1., 0.5});
    RadiusLaw<T> law{{{ex, 0., 0.1}, {ex, 1., 0.1}, {ey, 0., 0.2}, {ey, 1., 0.2}, {ez, 0., 0.3}, {ez, 1., 0.3}}, RadiusVariation::Linear};
    const auto r = fillet(m, s, std::vector{ex, ey, ez}, law);
    ASSERT_EQ(code(r), BuildErrc::NonGeneric);
    ASSERT_EQ(m.ids<FaceId>().size(), n);
}
