#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include "tests_brep_geom.h"
#include <gbs/bscbuild.h>
#include <gbs-io/mesh_writer.h>
#include <gbs-render/vtkMeshRender.h>
#include <gbs-io/step/reader.h>

#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
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

    /// What a face mesh is: its triangles counter-clockwise in (u, v), the deviation of the
    /// surface from them (sampled), the uses of its sides by vertex, its area in 3D.
    struct FaceStats
    {
        T deviation{};
        T area{};
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> sides; ///< directed, by vertex of the mesh
        std::size_t vertices{};
    };
    FaceStats stats(const Surface<T, 3> &srf, const MeshBuilder<T> &mb, const FaceMesh<T> &fm)
    {
        FaceStats r;
        std::set<std::uint32_t> used;
        for (const auto &t : fm.triangles)
        {
            const auto &a = fm.uv[t[0]], &b = fm.uv[t[1]], &c = fm.uv[t[2]];
            REQUIRE((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0.);
            const auto &p0 = mb.vertices()[fm.vertices[t[0]]], &p1 = mb.vertices()[fm.vertices[t[1]]], &p2 = mb.vertices()[fm.vertices[t[2]]];
            r.area += norm(cross(p1 - p0, p2 - p0)) / 2.;
            constexpr int n = 6;
            for (int i = 0; i <= n; ++i)
                for (int j = 0; i + j <= n; ++j)
                {
                    const T x = T(i) / n, y = T(j) / n;
                    const auto uv = a * (1. - x - y) + b * x + c * y;
                    r.deviation = std::max(r.deviation, mesh_detail::distance_to_triangle(srf.value(uv[0], uv[1]), p0, p1, p2));
                }
            for (int k = 0; k < 3; ++k)
            {
                ++r.sides[{fm.vertices[t[k]], fm.vertices[t[(k + 1) % 3]]}];
                used.insert(fm.vertices[t[k]]);
            }
        }
        r.vertices = used.size();
        return r;
    }
    /// Every side used once in each sense: a closed surface.
    bool closed(const FaceStats &s)
    {
        return std::ranges::all_of(s.sides, [&](const auto &e) {
            const auto back = s.sides.find({e.first.second, e.first.first});
            return e.second == 1 && back != s.sides.end() && back->second == 1;
        });
    }
    /// V - E + F.
    long euler(const FaceStats &s)
    {
        std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
        std::size_t corners{};
        for (const auto &[e, n] : s.sides)
        {
            edges.insert({std::min(e.first, e.second), std::max(e.first, e.second)});
            corners += n;
        }
        return long(s.vertices) - long(edges.size()) + long(corners / 3);
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
    // a rational circle of radius 2: within the deflection and the angle
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
    // a cubic of 199 spans along a straight line (an approximated intersection curve on a plane):
    // its pieces merged, its ends only
    {
        points_vector<T, 3> poles;
        std::vector<T> knots{0., 0., 0., 0.};
        for (int k = 0; k < 202; ++k)
            poles.push_back({3. * k / 201., 0., 0.});
        // (202 poles, degree 3: 206 flat knots, 199 inside)
        for (int k = 1; k < 199; ++k)
            knots.push_back(k / 199.);
        knots.insert(knots.end(), {1., 1., 1., 1.});
        const BSCurve<T, 3> line{poles, knots, 3};
        ASSERT_EQ(line.knots().size(), 200u);
        ASSERT_EQ(curve_parameters(line, 0., 1., MeshOptions<T>{.deflection = 1e-3}).size(), 2u);
        // the length still bounds the merged pieces
        const auto t = curve_parameters(line, 0., 1., MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.25});
        ASSERT_LE(t.size(), 14u); // (3 / 0.25 = 12 segments, at the knots: 13)
        for (std::size_t k = 0; k + 1 < t.size(); ++k)
            ASSERT_LE(norm(line.value(t[k + 1]) - line.value(t[k])), 0.25);
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

TEST(tests_brep_mesh, planar_faces)
{
    // a unit square: two triangles; its inner sides bounded by max_edge
    Model<T> m;
    const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.})));
    const auto srf = m.face(f).surface;
    {
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3}};
        const auto fm = unwrap(mb.face(f));
        ASSERT_TRUE(fm.complete);
        ASSERT_EQ(fm.triangles.size(), 2u);
        const auto s = stats(*srf, mb, fm);
        ASSERT_NEAR(s.area, 1., 1e-12);
        ASSERT_EQ(euler(s), 1);
    }
    {
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.2}};
        const auto fm = unwrap(mb.face(f));
        ASSERT_TRUE(fm.complete);
        const auto s = stats(*srf, mb, fm);
        ASSERT_NEAR(s.area, 1., 1e-12);
        ASSERT_EQ(euler(s), 1);
        for (const auto &[e, n] : s.sides)
            ASSERT_LE(norm(mb.vertices()[e.first] - mb.vertices()[e.second]), 0.2 + 1e-12);
        // the points are on the surface
        for (std::size_t k = 0; k < fm.uv.size(); ++k)
            ASSERT_LE(norm(srf->value(fm.uv[k][0], fm.uv[k][1]) - mb.vertices()[fm.vertices[k]]), 1e-12);
    }
    // a square with a circular hole: the hole left out
    {
        Model<T> m2;
        auto segment = [&](point<T, 3> p, point<T, 3> q) {
            return unwrap(make_edge(m2, std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{p, q}, std::vector<T>{0., 0., 1., 1.}, 1)));
        };
        const auto outer = unwrap(make_wire(m2, std::vector{segment({-1., -1., 0.}, {1., -1., 0.}), segment({1., -1., 0.}, {1., 1., 0.}),
                                                              segment({1., 1., 0.}, {-1., 1., 0.}), segment({-1., 1., 0.}, {-1., -1., 0.})}));
        const auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(0.5));
        const auto hole = unwrap(make_wire(m2, std::vector{unwrap(make_edge(m2, circle))}));
        const auto f2 = unwrap(make_face(m2, quad({-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}), outer, std::vector{hole}));
        MeshBuilder<T> mb{m2, MeshOptions<T>{.deflection = 1e-4}};
        const auto fm = unwrap(mb.face(f2));
        ASSERT_TRUE(fm.complete);
        const auto s = stats(*m2.face(f2).surface, mb, fm);
        ASSERT_NEAR(s.area, 4. - pi * 0.25, 1e-3); // (the hole a polygon within 1e-4 of the circle)
        ASSERT_EQ(euler(s), 0);                      // an annulus
        for (const auto &t : fm.triangles)
        {
            const auto g = (mb.vertices()[fm.vertices[t[0]]] + mb.vertices()[fm.vertices[t[1]]] + mb.vertices()[fm.vertices[t[2]]]) / 3.;
            ASSERT_GT(norm(g), 0.45);
        }
    }
}

TEST(tests_brep_mesh, curved_faces)
{
    // surfaces of revolution and a NURBS: within the deflection (sampled), along S_u x S_v, the
    // seams and the poles closing the mesh
    struct Case
    {
        std::string name;
        std::shared_ptr<Surface<T, 3>> srf;
        T area;      // exact
        long euler;  // V - E + F of the face
        bool closed; // every side used in both senses
        std::array<std::size_t, 2> triangles; // at 1e-2 and 1e-3, as measured (the tests bound them by 25 % more)
    };
    const std::vector<Case> cases{
        {"cylinder", cylinder(), 4. * pi, 0, false, {64, 256}},
        {"nurbs cylinder", nurbs_cylinder(), 4. * pi, 0, false, {64, 256}},
        {"sphere", sphere(), 4. * pi, 2, true, {1276, 13540}},
        {"torus", torus(), 4. * pi * pi * 3., 0, true, {4126, 42760}},
        {"cone", cone(), pi * std::sqrt(2.), 1, false, {414, 4948}},
        {"half revolution", revolution(pi), pi * 2. * std::sqrt(5.), 1, false, {398, 4760}}, // (half a frustum: radii 1 and 3, slant sqrt 5)
    };
    for (const auto &c : cases)
        for (T deflection : {1e-2, 1e-3})
        {
            INFO(c.name << ", deflection " << deflection);
            Model<T> m;
            const auto f = unwrap(make_face(m, c.srf));
            MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = deflection, .angle = 0.5}};
            const auto fm = unwrap(mb.face(f));
            ASSERT_TRUE(fm.complete);
            const auto s = stats(*c.srf, mb, fm);
            MESSAGE(c.name << " " << deflection << ": " << fm.triangles.size() << " triangles, deviation " << s.deviation << ", area " << s.area);
            ASSERT_EQ(euler(s), c.euler);
            ASSERT_EQ(closed(s), c.closed);
            ASSERT_LE(s.deviation, 1.05 * deflection); // (four samples per triangle in the criterion: 4 % over at most)
            ASSERT_LE(fm.triangles.size(), (deflection > 5e-3 ? c.triangles[0] : c.triangles[1]) * 5 / 4);
            ASSERT_NEAR(s.area, c.area, 2. * c.area * deflection);
            // along S_u x S_v
            for (const auto &t : fm.triangles)
            {
                const auto &p0 = mb.vertices()[fm.vertices[t[0]]], &p1 = mb.vertices()[fm.vertices[t[1]]], &p2 = mb.vertices()[fm.vertices[t[2]]];
                const auto g = (fm.uv[t[0]] + fm.uv[t[1]] + fm.uv[t[2]]) / 3.;
                INFO(fm.uv[t[0]][0] << "," << fm.uv[t[0]][1] << " " << fm.uv[t[1]][0] << "," << fm.uv[t[1]][1] << " " << fm.uv[t[2]][0] << "," << fm.uv[t[2]][1] << " v " << fm.vertices[t[0]] << " " << fm.vertices[t[1]] << " " << fm.vertices[t[2]]);
                ASSERT_GT(cross(p1 - p0, p2 - p0) * cross(c.srf->value(g[0], g[1], 1, 0), c.srf->value(g[0], g[1], 0, 1)), 0.);
            }
        }
}

TEST(tests_brep_mesh, faces_of_a_solid)
{
    // the faces of a cube: they share their edges, the mesh of the solid is closed
    Model<T> m;
    const auto solid = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front()));
    MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.5}};
    FaceStats all;
    for (const auto &fu : m.shell(m.solid(solid).outer).faces)
    {
        const auto fm = unwrap(mb.face(fu.face));
        ASSERT_TRUE(fm.complete);
        const auto s = stats(*m.face(fu.face).surface, mb, fm);
        for (const auto &[e, n] : s.sides)
            all.sides[fu.orient == Orientation::Forward ? e : std::pair{e.second, e.first}] += n; // in the sense of the shell
    }
    std::set<std::uint32_t> used;
    for (const auto &[e, n] : all.sides)
        used.insert(e.first);
    all.vertices = used.size();
    ASSERT_TRUE(closed(all));
    ASSERT_EQ(euler(all), 2);
}

TEST(tests_brep_mesh, faces_of_a_boolean)
{
    // a box 0..4 minus the sphere R 3 centred at its corner (4, 4, 4): trimmed planes, a trimmed
    // sphere, its pole and seam inside the faces; the mesh of the solid is closed, its volume
    // that of the solid
    Model<T> m;
    const auto bx = ShapeId{unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 4.))).shells.front()))};
    const auto sph = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/sphere.stp")).root; // R 3 at the origin
    const point<T, 3> z{1. / std::sqrt(3.), 1. / std::sqrt(3.), 1. / std::sqrt(3.)}, x{1. / std::sqrt(2.), -1. / std::sqrt(2.), 0.};
    const auto s = copy_shape(m, sph, Placement<T>::frame({4., 4., 4.}, x, cross(z, x), z));
    const auto d = unwrap(cut(m, bx, s));
    ASSERT_EQ(d.solids.size(), 1u);
    const T deflection = 1e-3;
    MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = deflection}};
    FaceStats all;
    T volume{}, area{};
    for (const auto &fu : m.shell(m.solid(d.solids.front()).outer).faces)
    {
        const auto fm = unwrap(mb.face(fu.face));
        ASSERT_TRUE(fm.complete);
        const auto st = stats(*m.face(fu.face).surface, mb, fm);
        ASSERT_LE(st.deviation, 1.05 * deflection);
        area += st.area;
        for (const auto &[e, n] : st.sides)
            all.sides[fu.orient == Orientation::Forward ? e : std::pair{e.second, e.first}] += n;
        for (const auto &t : fm.triangles)
        {
            const auto &p0 = mb.vertices()[fm.vertices[t[0]]], &p1 = mb.vertices()[fm.vertices[t[1]]], &p2 = mb.vertices()[fm.vertices[t[2]]];
            volume += (fu.orient == Orientation::Forward ? 1. : -1.) * (p0 * cross(p1, p2)) / 6.;
        }
    }
    std::set<std::uint32_t> used;
    for (const auto &[e, n] : all.sides)
        used.insert(e.first);
    all.vertices = used.size();
    ASSERT_TRUE(closed(all));
    ASSERT_EQ(euler(all), 2);
    const T eighth = 4. / 3. * pi * 27. / 8.;
    ASSERT_NEAR(volume, 64. - eighth, area * deflection);
}

TEST(tests_brep_mesh, face_options_and_failures)
{
    // the angle alone (a loose deflection): the normals at the corners of each triangle within it
    {
        Model<T> m;
        const auto srf = cylinder();
        const auto f = unwrap(make_face(m, srf));
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1., .angle = 0.1}};
        const auto fm = unwrap(mb.face(f));
        ASSERT_TRUE(fm.complete);
        for (const auto &t : fm.triangles)
            for (int k = 0; k < 3; ++k)
            {
                auto n = [&](std::uint32_t i) { return cross(srf->value(fm.uv[i][0], fm.uv[i][1], 1, 0), srf->value(fm.uv[i][0], fm.uv[i][1], 0, 1)); };
                ASSERT_LE(mesh_detail::angle_between(n(t[k]), n(t[(k + 1) % 3])), 0.1 + 1e-9);
            }
        ASSERT_GE(fm.triangles.size(), 2u * 63u); // (2 pi / 0.1 around)
    }
    // the length on a curved face: every side, of the boundary and inside
    {
        Model<T> m;
        const auto srf = torus();
        const auto f = unwrap(make_face(m, srf));
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-2, .max_edge = 0.2}};
        const auto fm = unwrap(mb.face(f));
        ASSERT_TRUE(fm.complete);
        const auto s = stats(*srf, mb, fm);
        ASSERT_TRUE(closed(s));
        for (const auto &[e, n] : s.sides)
            ASSERT_LE(norm(mb.vertices()[e.first] - mb.vertices()[e.second]), 0.2 + 1e-12);
    }
    // the limit of points: the face left as it is, a valid closed mesh, not complete
    {
        Model<T> m;
        const auto srf = sphere();
        const auto f = unwrap(make_face(m, srf));
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3, .max_points = 300}};
        const auto fm = unwrap(mb.face(f));
        ASSERT_FALSE(fm.complete);
        ASSERT_LE(fm.uv.size(), 300u);
        const auto s = stats(*srf, mb, fm);
        ASSERT_TRUE(closed(s));
        ASSERT_EQ(euler(s), 2);
    }
    // a boundary crossing itself in (u, v) (the pcurve of the hole across the outer wire):
    // NonGeneric, the interior points of the face dropped
    {
        Model<T> m;
        auto segment = [&](point<T, 3> p, point<T, 3> q) {
            return unwrap(make_edge(m, std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{p, q}, std::vector<T>{0., 0., 1., 1.}, 1)));
        };
        const auto outer = unwrap(make_wire(m, std::vector{segment({-1., -1., 0.}, {1., -1., 0.}), segment({1., -1., 0.}, {1., 1., 0.}),
                                                             segment({1., 1., 0.}, {-1., 1., 0.}), segment({-1., 1., 0.}, {-1., -1., 0.})}));
        const auto hole = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(0.5))))}));
        const auto f = unwrap(make_face(m, quad({-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}), outer, std::vector{hole}));
        // (the face maps (u, v) in [0, 1] onto [-2, 2]: the outer wire is the square 0.25..0.75)
        m.wire(hole).coedges.front().pcurve = std::make_shared<BSCurveRational<T, 2>>(build_circle<T, 2>(0.3, {0.5, 0.5}));
        MeshBuilder<T> mb{m, MeshOptions<T>{.deflection = 1e-3}};
        std::ignore = unwrap(mb.boundary(f));
        const auto before = mb.vertices().size();
        const auto r = mb.face(f);
        ASSERT_FALSE(r.has_value());
        ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
        ASSERT_EQ(mb.vertices().size(), before);
    }
}

TEST(tests_brep_mesh, mesh_of_shapes)
{
    // solids read from STEP files: closed meshes facing outwards, their volumes those of the
    // solids (sampled by signed_volume), the normals of the corners outwards and unit
    for (const char *name : {"box.stp", "cylinder.stp", "sphere.stp", "cone.stp", "occt_partial_torus.stp"})
    {
        INFO(name);
        Model<T> m;
        const auto root = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/" + name)).root;
        const auto solids = explore<SolidId>(m, root);
        const auto mesh = unwrap(brep::mesh(m, root, MeshOptions<T>{}));
        ASSERT_TRUE(mesh.incomplete.empty());
        ASSERT_TRUE(mesh.folded.empty());
        ASSERT_EQ(mesh.uv.size(), 3 * mesh.triangles.size());
        ASSERT_EQ(mesh.normals.size(), 3 * mesh.triangles.size());
        ASSERT_EQ(mesh.faces.size(), explore<FaceId>(m, root).size());
        std::uint32_t next{};
        for (const auto &fr : mesh.faces)
        {
            ASSERT_EQ(fr.first, next);
            next += fr.count;
        }
        ASSERT_EQ(next, mesh.triangles.size());
        // every vertex used
        std::vector<int> used(mesh.vertices.size());
        for (const auto &t : mesh.triangles)
            for (auto v : t)
                used[v] = 1;
        ASSERT_TRUE(std::ranges::all_of(used, [](int u) { return u == 1; }));
        for (std::size_t k = 0; k < mesh.triangles.size(); ++k)
        {
            const auto &t = mesh.triangles[k];
            const auto nt = cross(mesh.vertices[t[1]] - mesh.vertices[t[0]], mesh.vertices[t[2]] - mesh.vertices[t[0]]);
            for (int c = 0; c < 3; ++c)
            {
                ASSERT_NEAR(norm(mesh.normals[3 * k + c]), 1., 1e-12);
                ASSERT_GT(mesh.normals[3 * k + c] * nt, 0.);
            }
        }
        MESSAGE(std::string{name} << ": " << mesh.triangles.size() << " triangles, volume " << volume(mesh) << ", area " << area(mesh));
        if (solids.size() != 1)
            continue;
        ASSERT_EQ(open_sides(mesh), 0u);
        const auto b = bounding_box(m, root);
        const T deflection = 1e-3 * norm(b.max - b.min);
        const T exact = signed_volume(m, m.solid(solids.front()).outer, 160);
        ASSERT_GT(volume(mesh), 0.);
        ASSERT_NEAR(volume(mesh), exact, area(mesh) * deflection + 1e-4 * std::abs(exact));
    }
    // a face alone: open, its boundary the open sides
    {
        Model<T> m;
        const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.})));
        const auto mesh = unwrap(brep::mesh(m, ShapeId{f}, MeshOptions<T>{.deflection = 1e-3, .max_edge = 0.25}));
        ASSERT_EQ(open_sides(mesh), 16u);
        ASSERT_NEAR(area(mesh), 1., 1e-12);
        for (const auto &n : mesh.normals)
            ASSERT_NEAR(n[2], 1., 1e-12);
    }
    // the uses of the faces in their shell orient the mesh: its volume has the sign of the
    // shell's (a sewn shell may face inwards), and turns over with them
    {
        Model<T> m;
        const auto shell = unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front();
        const T before = volume(unwrap(brep::mesh(m, ShapeId{shell}, MeshOptions<T>{.deflection = 1e-3})));
        ASSERT_NEAR(std::abs(before), 8., 1e-12);
        ASSERT_EQ(before > 0., signed_volume(m, shell) > 0.);
        for (auto &fu : m.shell(shell).faces)
            fu.orient = fu.orient == Orientation::Forward ? Orientation::Reversed : Orientation::Forward;
        const auto mesh = unwrap(brep::mesh(m, ShapeId{shell}, MeshOptions<T>{.deflection = 1e-3}));
        ASSERT_EQ(open_sides(mesh), 0u);
        ASSERT_NEAR(volume(mesh), -before, 1e-12);
    }
}

TEST(tests_brep_mesh, stl_and_obj)
{
    Model<T> m;
    const auto root = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/cylinder.stp")).root;
    const auto mesh = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-2}));
    const auto n = mesh.triangles.size();
    // binary STL: 80 bytes, the count, 50 bytes per triangle, little endian; read back
    {
        const auto stl = io::write_stl(mesh);
        ASSERT_EQ(stl.size(), 84 + 50 * n);
        auto u32 = [&](std::size_t at) {
            std::uint32_t x{};
            for (int k = 3; k >= 0; --k)
                x = (x << 8) | static_cast<unsigned char>(stl[at + k]);
            return x;
        };
        auto f32 = [&](std::size_t at) {
            const auto x = u32(at);
            float f;
            std::memcpy(&f, &x, 4);
            return f;
        };
        ASSERT_EQ(u32(80), n);
        T v{};
        for (std::size_t k = 0; k < n; ++k)
        {
            const std::size_t at = 84 + 50 * k;
            point<T, 3> nrm{f32(at), f32(at + 4), f32(at + 8)}, p[3];
            for (int c = 0; c < 3; ++c)
                p[c] = {f32(at + 12 + 12 * c), f32(at + 16 + 12 * c), f32(at + 20 + 12 * c)};
            ASSERT_NEAR(norm(nrm), 1., 1e-6);
            ASSERT_NEAR(norm(p[0] - mesh.vertices[mesh.triangles[k][0]]), 0., 1e-5);
            v += p[0] * cross(p[1], p[2]) / 6.;
        }
        ASSERT_NEAR(v, volume(mesh), 1e-4 * volume(mesh)); // (single precision)
    }
    // ASCII STL: one facet per triangle, the vertices exact
    {
        const auto stl = io::write_stl(mesh, false, "cylinder");
        ASSERT_TRUE(stl.starts_with("solid cylinder\n"));
        ASSERT_TRUE(stl.ends_with("endsolid cylinder\n"));
        std::istringstream in{stl};
        std::string line;
        std::size_t facets{}, vertices{};
        T v{};
        point<T, 3> p[3];
        while (std::getline(in, line))
        {
            std::istringstream l{line};
            std::string word;
            l >> word;
            if (word == "facet")
                ++facets;
            else if (word == "vertex")
            {
                auto &q = p[vertices % 3];
                l >> q[0] >> q[1] >> q[2];
                if (++vertices % 3 == 0)
                    v += p[0] * cross(p[1], p[2]) / 6.;
            }
        }
        ASSERT_EQ(facets, n);
        ASSERT_EQ(vertices, 3 * n);
        ASSERT_NEAR(v, volume(mesh), 1e-9 * volume(mesh));
    }
    // OBJ: the vertices once, a normal per corner, a group per face
    {
        const auto obj = io::write_obj(mesh);
        std::istringstream in{obj};
        std::string line;
        std::size_t v{}, vn{}, f{}, g{};
        while (std::getline(in, line))
        {
            v += line.starts_with("v ");
            vn += line.starts_with("vn ");
            f += line.starts_with("f ");
            g += line.starts_with("g ");
        }
        ASSERT_EQ(v, mesh.vertices.size());
        ASSERT_EQ(vn, 3 * n);
        ASSERT_EQ(f, n);
        ASSERT_EQ(g, mesh.faces.size());
    }
    // files
    {
        const auto dir = std::filesystem::temp_directory_path();
        const auto stl = dir / "gbs_tests_mesh.stl", obj = dir / "gbs_tests_mesh.obj";
        ASSERT_TRUE(io::write_stl_file(mesh, stl).has_value());
        ASSERT_EQ(std::filesystem::file_size(stl), 84 + 50 * n);
        ASSERT_TRUE(io::write_obj_file(mesh, obj).has_value());
        ASSERT_EQ(std::filesystem::file_size(obj), io::write_obj(mesh).size());
        std::filesystem::remove(stl);
        std::filesystem::remove(obj);
        const auto bad = io::write_stl_file(mesh, dir / "no_such_directory" / "x.stl");
        ASSERT_FALSE(bad.has_value());
        ASSERT_EQ(bad.error().code, BuildErrc::InvalidFile);
    }
}

TEST(tests_brep_mesh, vtk)
{
    // the polydata of a mesh: three points per triangle, the normals of the corners
    Model<T> m;
    const auto root = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/sphere.stp")).root;
    const auto mesh = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-2}));
    const auto poly = make_polydata(mesh);
    ASSERT_EQ(poly->GetNumberOfCells(), static_cast<vtkIdType>(mesh.triangles.size()));
    ASSERT_EQ(poly->GetNumberOfPoints(), static_cast<vtkIdType>(3 * mesh.triangles.size()));
    auto *normals = poly->GetPointData()->GetNormals();
    ASSERT_NE(normals, nullptr);
    for (vtkIdType i = 0; i < poly->GetNumberOfPoints(); ++i)
    {
        double p[3], n[3];
        poly->GetPoint(i, p);
        normals->GetTuple(i, n);
        // (the sphere R 3 at the origin: the normal is the point over 3, at the poles too)
        INFO(i << " p " << p[0] << " " << p[1] << " " << p[2] << " n " << n[0] << " " << n[1] << " " << n[2]);
        for (int c = 0; c < 3; ++c)
            ASSERT_NEAR(n[c], p[c] / norm(point<T, 3>{p[0], p[1], p[2]}), 1e-2);
    }
    ASSERT_NE(make_actor(mesh, {0.2, 0.6, 0.8})->GetMapper(), nullptr);
}

TEST(tests_brep_mesh, options_of_shapes)
{
    Model<T> m;
    const auto root = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/cone.stp")).root;
    // without the smallest angle: fewer triangles, some under 20 degrees, still within the
    // deflection and closed
    {
        const auto with = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-3}));
        const auto without = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-3, .min_angle = 0.}));
        ASSERT_LT(without.triangles.size(), with.triangles.size());
        ASSERT_EQ(open_sides(without), 0u);
        ASSERT_TRUE(without.folded.empty());
        MESSAGE("cone: " << with.triangles.size() << " triangles at 20 degrees, " << without.triangles.size() << " without");
    }
    // the limit of points reached: the face reported, the mesh still closed
    {
        const auto mesh = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-4, .max_points = 100}));
        ASSERT_FALSE(mesh.incomplete.empty());
        ASSERT_EQ(open_sides(mesh), 0u);
    }
}

TEST(tests_brep_mesh, writers_header_and_locale)
{
    Model<T> m;
    const auto root = unwrap(step::read_step_file(m, std::string(GBS_TESTS_IN_DIR) + "/step/cylinder.stp")).root;
    const auto mesh = unwrap(brep::mesh(m, root, MeshOptions<T>{.deflection = 1e-2}));
    // a binary STL never starts with "solid" (the ASCII signature), whatever its name
    {
        const auto stl = io::write_stl(mesh, true, "solid_blade");
        ASSERT_FALSE(stl.starts_with("solid"));
        ASSERT_TRUE(stl.starts_with("binary STL solid_blade"));
        ASSERT_EQ(stl.size(), 84 + 50 * mesh.triangles.size());
    }
    // a point as decimal separator under a locale with a comma
    const char *locales[] = {"fr_FR.UTF-8", "fr_FR.utf8", "de_DE.UTF-8", "de_DE.utf8", "French_France.1252"};
    const char *found = nullptr;
    for (const char *l : locales)
        if (std::setlocale(LC_NUMERIC, l))
        {
            found = l;
            break;
        }
    if (!found)
    {
        MESSAGE("no locale with a decimal comma on this host");
        return;
    }
    char probe[16];
    std::snprintf(probe, sizeof probe, "%.1f", 1.5);
    const auto stl = io::write_stl(mesh, false), obj = io::write_obj(mesh);
    std::setlocale(LC_NUMERIC, "C");
    INFO(std::string{found} << " formats 1.5 as " << probe);
    ASSERT_EQ(std::string{probe}, "1,5");
    ASSERT_EQ(stl.find(','), std::string::npos);
    ASSERT_EQ(obj.find(','), std::string::npos);
    ASSERT_NE(stl.find('.'), std::string::npos);
}
