#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include "tests_brep_geom.h"
#include <gbs/bscbuild.h>

#include <cmath>
#include <numbers>
#include <set>

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using namespace brep_tests;

    /// The largest distance of the curve to the chords of its parameters, sampled.
    T deviation(const Curve<T, 3> &c, const std::vector<T> &t)
    {
        T d{};
        for (std::size_t k = 0; k + 1 < t.size(); ++k)
            for (int i = 1; i < 20; ++i)
                d = std::max(d, intersect_detail::distance_to_segment(c.value(t[k] + (t[k + 1] - t[k]) * i / 20.), c.value(t[k]), c.value(t[k + 1])));
        return d;
    }
}

TEST(tests_brep_mesh, curve_parameters)
{
    // a segment: its ends only, unless the length is bounded
    const BSCurve<T, 3> line{points_vector<T, 3>{{0., 0., 0.}, {4., 0., 0.}}, std::vector<T>{0., 0., 1., 1.}, 1};
    ASSERT_EQ(curve_parameters(line, 0., 1., MeshOptions<T>{.deflection = 1e-3}).size(), 2u);
    const auto bounded = curve_parameters(line, 0., 1., MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.3});
    ASSERT_EQ(bounded.size(), 17u); // halved down to 0.25
    for (std::size_t k = 0; k + 1 < bounded.size(); ++k)
        ASSERT_LE(norm(line.value(bounded[k + 1]) - line.value(bounded[k])), 0.3);
    // a rational circle of radius 2: within the deflection and the angle, its knots kept
    const auto circle = build_circle<T, 3>(2.);
    const T a = circle.bounds()[0], b = circle.bounds()[1];
    for (T deflection : {1e-2, 1e-3, 1e-5})
    {
        INFO("deflection " << deflection);
        const MeshOptions<T> o{.deflection = deflection, .angle = 0.5};
        const auto t = curve_parameters(circle, a, b, o);
        ASSERT_TRUE(std::ranges::is_sorted(t));
        ASSERT_EQ(t.front(), a);
        ASSERT_EQ(t.back(), b);
        ASSERT_LE(deviation(circle, t), deflection);
        for (std::size_t k = 0; k + 1 < t.size(); ++k)
            ASSERT_LE(mesh_detail::angle_between(circle.value(t[k], 1), circle.value(t[k + 1], 1)), 0.5 + 1e-12);
        for (T knot : circle.knots())
            ASSERT_TRUE(std::ranges::any_of(t, [&](T x) { return x == knot; }));
        // the chordal deviation of a circle is r (1 - cos(phi / 2)): about as many points as needed
        const T phi = 2. * std::acos(1. - deflection / 2.);
        ASSERT_LE(t.size(), 4. * (2. * std::numbers::pi / phi) + 8.);
    }
    // the length bounded on the circle
    {
        const auto t = curve_parameters(circle, a, b, MeshOptions<T>{.deflection = 1e-2, .max_edge = 0.1});
        for (std::size_t k = 0; k + 1 < t.size(); ++k)
            ASSERT_LE(norm(circle.value(t[k + 1]) - circle.value(t[k])), 0.1);
    }
    // a wave of one quintic span through its chord at the quarter points, y = c t (t - 1/4)
    // (t - 1/2) (t - 3/4) (t - 1): three samples see no deviation, the poles bound it
    {
        const std::array<T, 6> power{0., 3. / 32., -25. / 32., 35. / 16., -5. / 2., 1.}; // t (t - 1/4) (t - 1/2) (t - 3/4) (t - 1)
        const T c = 3.;
        points_vector<T, 3> poles;
        for (int i = 0; i <= 5; ++i)
        {
            // Bernstein coefficients of the power basis: b_i = sum_j C(i, j) / C(5, j) a_j
            T y{};
            for (int j = 0; j <= i; ++j)
                y += std::tgamma(i + 1.) / (std::tgamma(j + 1.) * std::tgamma(i - j + 1.)) /
                     (std::tgamma(6.) / (std::tgamma(j + 1.) * std::tgamma(6. - j))) * c * power[j];
            poles.push_back({i / 5., y, 0.});
        }
        const BSCurve<T, 3> wave{poles, std::vector<T>{0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1.}, 5};
        for (T x : {0.25, 0.5, 0.75})
            REQUIRE(std::abs(wave.value(x)[1]) <= 1e-15);
        REQUIRE(std::abs(wave.value(0.1)[1]) > 5e-3);
        // (the angle out of play: the deflection alone decides)
        const auto t = curve_parameters(wave, 0., 1., MeshOptions<T>{.deflection = 1e-3, .angle = 10.});
        ASSERT_LE(deviation(wave, t), 1e-3);
    }
    // a C0 knot (multiplicity 3 in a cubic): a kink, never inside an interval
    {
        const BSCurve<T, 3> kink{points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {3., 0., 0.}, {3., 1., 0.}, {3., 2., 0.}, {3., 3., 0.}},
                                 std::vector<T>{0., 0., 0., 0., 1., 1., 1., 2., 2., 2., 2.}, 3};
        const auto t = curve_parameters(kink, 0., 2., MeshOptions<T>{.deflection = 1e-3});
        ASSERT_TRUE(std::ranges::any_of(t, [](T x) { return x == 1.; }));
        ASSERT_EQ(t.size(), 3u); // two straight pieces
    }
    // a deflection or an angle that is not positive is refused
    ASSERT_THROW(std::ignore = curve_parameters(circle, a, b, MeshOptions<T>{}), std::invalid_argument);
    ASSERT_THROW(std::ignore = curve_parameters(circle, a, b, MeshOptions<T>{.deflection = 1e-3, .angle = 0.}), std::invalid_argument);
}

TEST(tests_brep_mesh, face_boundaries)
{
    const MeshOptions<T> o{.deflection = 1e-3, .angle = 0.5};
    // a planar quad: one loop of its four corners
    {
        Model<T> m;
        const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.})));
        MeshBuilder<T> mb{m, o};
        const auto loops = unwrap(mb.boundary(f));
        ASSERT_EQ(loops.size(), 1u);
        ASSERT_EQ(loops[0].size(), 4u);
        ASSERT_EQ(mb.vertices().size(), 4u);
    }
    // a cylinder: the seam twice (at u1 and at u2) with the same vertices, the circles shared
    // between the loop and the vertices of the mesh
    for (auto srf : {cylinder(), nurbs_cylinder()})
    {
        Model<T> m;
        const auto f = unwrap(make_face(m, srf));
        MeshBuilder<T> mb{m, o};
        const auto loops = unwrap(mb.boundary(f));
        ASSERT_EQ(loops.size(), 1u);
        const auto &loop = loops[0];
        // every vertex of the mesh is a point of the loop, the seam ones twice
        std::multiset<std::uint32_t> uses;
        for (const auto &p : loop)
            uses.insert(p.vertex);
        ASSERT_EQ(std::set<std::uint32_t>(uses.begin(), uses.end()).size(), mb.vertices().size());
        std::size_t twice{};
        for (std::uint32_t v = 0; v < mb.vertices().size(); ++v)
            twice += uses.count(v) == 2;
        ASSERT_EQ(twice, 2u); // the two ends of the seam
        // each point is the surface at its parameters
        for (const auto &p : loop)
            ASSERT_LE(norm(srf->value(p.uv[0], p.uv[1]) - mb.vertices()[p.vertex]), 1e-9);
        // two points of the loop at the same vertex are on both sides of the seam
        for (const auto &p : loop)
            for (const auto &q : loop)
                if (p.vertex == q.vertex && p.uv[0] != q.uv[0])
                    ASSERT_NEAR(std::abs(p.uv[0] - q.uv[0]), srf->bounds()[1] - srf->bounds()[0], 1e-9);
    }
    // a sphere: the points of each pole are one vertex
    {
        Model<T> m;
        const auto f = unwrap(make_face(m, sphere()));
        MeshBuilder<T> mb{m, o};
        const auto loops = unwrap(mb.boundary(f));
        ASSERT_EQ(loops.size(), 1u);
        std::size_t at_poles{};
        for (const auto &p : loops[0])
            if (std::abs(std::abs(p.uv[1]) - std::numbers::pi / 2.) < 1e-12)
            {
                ++at_poles;
                ASSERT_NEAR(std::abs(mb.vertices()[p.vertex][2]), 1., 1e-12);
            }
        // a full turn at each pole: 2 pi / 0.5 rounded up, 13 segments
        ASSERT_GE(at_poles, 2u * 13u);
        std::set<std::uint32_t> pole_vertices;
        for (const auto &p : loops[0])
            if (std::abs(std::abs(p.uv[1]) - std::numbers::pi / 2.) < 1e-12)
                pole_vertices.insert(p.vertex);
        ASSERT_EQ(pole_vertices.size(), 2u);
    }
}

TEST(tests_brep_mesh, holes_and_options)
{
    // a square face with a circular hole: the outer loop first, counter-clockwise in (u, v), the
    // hole clockwise; the hole is one closed edge, its loop closes on its vertex
    Model<T> m;
    auto segment = [&](point<T, 3> p, point<T, 3> q) {
        return unwrap(make_edge(m, std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{p, q}, std::vector<T>{0., 0., 1., 1.}, 1)));
    };
    const auto outer = unwrap(make_wire(m, std::vector{segment({-1., -1., 0.}, {1., -1., 0.}), segment({1., -1., 0.}, {1., 1., 0.}),
                                                         segment({1., 1., 0.}, {-1., 1., 0.}), segment({-1., 1., 0.}, {-1., -1., 0.})}));
    const auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(0.5));
    const auto hole = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle))}));
    const auto f = unwrap(make_face(m, quad({-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}), outer, std::vector{hole}));
    MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3}};
    const auto loops = unwrap(mb.boundary(f));
    ASSERT_EQ(loops.size(), 2u);
    auto area = [](const std::vector<BoundaryPoint<T>> &l) {
        T a{};
        for (std::size_t k = 0; k < l.size(); ++k)
        {
            const auto &p = l[k].uv, &q = l[(k + 1) % l.size()].uv;
            a += p[0] * q[1] - q[0] * p[1];
        }
        return a / 2.;
    };
    ASSERT_GT(area(loops[0]), 0.);
    ASSERT_LT(area(loops[1]), 0.);
    ASSERT_EQ(loops[0].size(), 4u);
    // the hole: its vertex once, at the start of the loop, the other points inside the edge
    const auto v = mb.vertex(m.edge(m.wire(hole).coedges.front().edge).v1);
    ASSERT_EQ(std::ranges::count_if(loops[1], [&](const auto &p) { return p.vertex == v; }), 1);
    ASSERT_GE(loops[1].size(), 8u);
    // the default deflection: 1e-3 of the diagonal of the shape (2 sqrt 2), from the shape
    MeshBuilder<T> by_default{m, ShapeId{f}, MeshOptions<T>{}};
    ASSERT_NEAR(by_default.options().deflection, 2e-3 * std::sqrt(2.), 1e-5);
    ASSERT_EQ(unwrap(by_default.boundary(f)).size(), 2u);
    // options as given must be positive
    ASSERT_THROW((MeshBuilder<T>{m, MeshOptions<T>{}}), std::invalid_argument);
    // a co-edge without pcurve, an edge not parametrized like its pcurves: refused
    {
        Model<T> m2 = m;
        m2.wire(hole).coedges.front().pcurve = nullptr;
        MeshBuilder<T> mb2{m2, MeshOptions<T>{.deflection = 1e-3}};
        const auto r = mb2.boundary(f);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::MissingPCurve);
    }
    {
        Model<T> m2 = m;
        m2.edge(m2.wire(hole).coedges.front().edge).same_parameter = false;
        MeshBuilder<T> mb2{m2, MeshOptions<T>{.deflection = 1e-3}};
        const auto r = mb2.boundary(f);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::NotSameParameter);
    }
}

TEST(tests_brep_mesh, shared_edges)
{
    // the faces of a cube: an edge shared by two faces gives both the same vertices, in opposite
    // orders; 8 corners and no other vertex (straight edges)
    Model<T> m;
    const auto solid = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front()));
    MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.5}};
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> sides;
    for (auto f : explore<FaceId>(m, ShapeId{solid}))
    {
        const auto loops = unwrap(mb.boundary(f));
        ASSERT_EQ(loops.size(), 1u);
        ASSERT_EQ(loops[0].size(), 16u); // four edges of length 2 in segments of 0.5
        for (std::size_t k = 0; k < loops[0].size(); ++k)
        {
            const auto a = loops[0][k].vertex, b = loops[0][(k + 1) % loops[0].size()].vertex;
            ++sides[{std::min(a, b), std::max(a, b)}];
        }
    }
    // every segment of the boundaries is used by two faces
    ASSERT_EQ(sides.size(), 12u * 4u);
    for (const auto &[s, n] : sides)
        ASSERT_EQ(n, 2);
    ASSERT_EQ(mb.vertices().size(), 8u + 12u * 3u);
}
