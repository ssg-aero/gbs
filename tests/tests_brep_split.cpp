#include <doctest_gtest.hpp>
#include <gbs-brep/split.h>
#include <gbs-brep/check.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <sstream>
#include <string>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using brep_tests::pi;
    using T = double;

    std::string file(const std::string &name) { return std::string(GBS_TESTS_IN_DIR) + "/step/" + name; }

    ShapeId read(Model<T> &m, const std::string &name)
    {
        auto r = step::read_step_file(m, file(name));
        REQUIRE(r.has_value());
        return r->root;
    }

    std::shared_ptr<Plane<T>> plane(const point<T, 3> &o, const point<T, 3> &n, T half = 100.)
    {
        const auto z = intersect_detail::unit<T>(n);
        return std::make_shared<Plane<T>>(intersect_detail::frame_of<T>(o, z, intersect_detail::orthogonal<T>(z)), std::array<T, 4>{-half, half, -half, half});
    }

    void expect_valid(const Model<T> &m, const ShapeId &s)
    {
        const auto rep = check(m, s);
        std::ostringstream os;
        for (const auto &e : rep.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        INFO(os.str());
        REQUIRE(rep.ok());
    }

    // volume by the sampling of the faces (signed_volume): about 1e-4 relative on trimmed faces
    T volume(const Model<T> &m, const ShapeId &s)
    {
        T v{};
        for (auto so : explore<SolidId>(m, s))
            v += signed_volume(m, m.solid(so).outer, 256);
        return v;
    }

    // area of the faces by the surface element on their interior points, summed over a grid
    T area(const Model<T> &m, FaceId f)
    {
        const auto &s = *m.face(f).surface;
        const auto b = classify_detail::face_uv_box(m, f);
        constexpr int n = 60;
        T a{};
        const T du = (b[1] - b[0]) / n, dv = (b[3] - b[2]) / n;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
            {
                const point<T, 2> uv{b[0] + (i + 0.5) * du, b[2] + (j + 0.5) * dv};
                if (classify(m, f, uv) == Location::In)
                    a += norm(cross(s.value(uv[0], uv[1], 1, 0), s.value(uv[0], uv[1], 0, 1))) * du * dv;
            }
        return a;
    }
}

TEST(tests_brep_split, plate_through_its_hole)
{
    // the plate 4 x 2 in with a hole of radius 0.5 in at (1, 1), cut by the plane x = 1 in: two
    // faces, the hole cut in two halves
    Model<T> m;
    const auto plate = read(m, "plate_inch.stp");
    const auto f0 = explore<FaceId>(m, plate).front();
    const T a0 = area(m, f0);
    auto r = unwrap(split(m, plate, plane({25.4, 0., 0.}, {1., 0., 0.})));
    ASSERT_EQ(r.faces.size(), 1u);
    ASSERT_EQ(r.faces[0].second.size(), 2u);
    ASSERT_TRUE(r.dangling.empty());
    const auto faces = explore<FaceId>(m, plate);
    ASSERT_EQ(faces.size(), 2u);
    T a{};
    for (auto f : faces)
    {
        ASSERT_EQ(m.face(f).wires.size(), 1u); // no hole left in either piece
        a += area(m, f);
    }
    ASSERT_NEAR(a, a0, 1e-2 * a0);
    expect_valid(m, plate);
}

TEST(tests_brep_split, solids_stay_closed)
{
    // a box cut by an oblique plane: the crossed faces are split, the shell stays closed
    Model<T> m;
    const auto box = read(m, "box.stp"); // 2 x 3 x 4
    const T v0 = volume(m, box);
    auto r = unwrap(split(m, box, plane({1., 1.5, 2.}, {0.3, 0.2, 1.})));
    ASSERT_EQ(r.faces.size(), 4u); // the four sides
    ASSERT_EQ(explore<FaceId>(m, box).size(), 10u);
    expect_valid(m, box);
    ASSERT_NEAR(volume(m, box), v0, 1e-3 * v0);

    // the cylinder cut across: the lateral face in two, across its seam
    const auto cyl = read(m, "cylinder.stp");
    const T vc = volume(m, cyl);
    auto c = unwrap(split(m, cyl, plane({0., 0., 5.}, {0., 0., 1.})));
    ASSERT_EQ(c.faces.size(), 1u);
    ASSERT_EQ(c.faces[0].second.size(), 2u);
    ASSERT_EQ(explore<FaceId>(m, cyl).size(), 4u);
    expect_valid(m, cyl);
    ASSERT_NEAR(volume(m, cyl), vc, 1e-3 * vc);
    // obliquely: an ellipse cut at the seam into two edges
    const auto cyl2 = read(m, "cylinder.stp");
    unwrap(split(m, cyl2, plane({0., 0., 5.}, {0., 0.5, 1.})));
    ASSERT_EQ(explore<FaceId>(m, cyl2).size(), 4u);
    expect_valid(m, cyl2);
    ASSERT_NEAR(volume(m, cyl2), vc, 1e-3 * vc);
}

TEST(tests_brep_split, holes_and_dangling_cuts)
{
    // a box pierced by a thin cylinder: the top and the bottom get a disc and a holed face
    Model<T> m;
    const auto box = read(m, "box.stp");
    const T v0 = volume(m, box);
    const auto tool = std::make_shared<CylindricalSurface<T>>(0.5, ax2<T, 3>{{{1., 1.5, 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., 2 * pi, -1., 5.});
    auto r = unwrap(split(m, box, tool));
    ASSERT_EQ(r.faces.size(), 2u);
    for (const auto &[f, pieces] : r.faces)
    {
        ASSERT_EQ(pieces.size(), 2u);
        std::size_t holed{};
        for (auto p : pieces)
            holed += m.face(p).wires.size() == 2;
        ASSERT_EQ(holed, 1u);
    }
    ASSERT_EQ(explore<FaceId>(m, box).size(), 8u);
    expect_valid(m, box);
    ASSERT_NEAR(volume(m, box), v0, 1e-3 * v0);

    // a cut from the boundary into the middle of a face: dangling, the face is kept
    using namespace brep_tests;
    Model<T> mp;
    const auto f = unwrap(make_face(mp, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
    const auto e = unwrap(make_edge(mp, point<T, 3>{1., 0., 0.}, point<T, 3>{1., 1., 0.}));
    const auto &ed = mp.edge(e);
    auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0.5, 0.}, {0.5, 0.5}}, std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1);
    auto d = unwrap(split_face(mp, f, std::vector<FaceCut<T>>{{e, pc}}));
    ASSERT_EQ(d.faces.size(), 1u);
    ASSERT_EQ(d.faces[0], f);
    ASSERT_EQ(d.dangling.size(), 1u);
}

TEST(tests_brep_split, by_another_shape)
{
    // boxes 0..2 and 1..3: the faces of the first one crossed by the second one are split
    using namespace brep_tests;
    Model<T> m;
    auto a = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front()));
    auto b = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {1., 1., 1.}, 2.))).shells.front()));
    const T va = volume(m, a);
    auto r = unwrap(split(m, a, b));
    ASSERT_EQ(r.faces.size(), 3u); // the faces x = 2, y = 2, z = 2
    ASSERT_EQ(explore<FaceId>(m, a).size(), 9u);
    expect_valid(m, a);
    expect_valid(m, b);
    ASSERT_NEAR(volume(m, a), va, 1e-3 * va);
    // errors
    ASSERT_FALSE(split_face(m, FaceId{999}, std::vector<FaceCut<T>>{}).has_value());
}

TEST(tests_brep_split, hole_against_pinched_regions)
{
    // two diamonds touching at their common vertex (1, 1) inside a square face, turned around it:
    // the hole of the outer piece runs along both of them, through the vertex twice, and must not
    // be given to a diamond whose boundary its first point lies on (an even-odd test counts the
    // leftmost vertex of a diamond inside it)
    using namespace brep_tests;
    for (int turn = 0; turn < 16; ++turn)
    {
        INFO("turned by " << turn << " pi / 8");
        const T a = turn * pi / 8., ca = std::cos(a), sa = std::sin(a);
        auto at = [&](T x, T y) { return point<T, 3>{1. + ca * (x - 1.) - sa * (y - 1.), 1. + sa * (x - 1.) + ca * (y - 1.), 0.}; };
        Model<T> m;
        const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
        const auto square_edges = explore<EdgeId>(m, ShapeId{f});
        const std::vector<point<T, 3>> left{at(1., 1.), at(0.8, 1.2), at(0.6, 1.), at(0.8, 0.8)};
        const std::vector<point<T, 3>> right{at(1., 1.), at(1.2, 0.8), at(1.4, 1.), at(1.2, 1.2)};
        const auto centre = unwrap(make_vertex(m, left[0]));
        std::vector<FaceCut<T>> cuts;
        for (const auto *d : {&left, &right})
        {
            std::vector<VertexId> vs{centre};
            for (std::size_t k = 1; k < d->size(); ++k)
                vs.push_back(unwrap(make_vertex(m, (*d)[k])));
            for (std::size_t k = 0; k < d->size(); ++k)
            {
                const auto va = vs[k], vb = vs[(k + 1) % vs.size()];
                const auto e = unwrap(make_edge(m, va, vb));
                const auto &ed = m.edge(e);
                const auto pa = m.vertex(va).pnt, pb = m.vertex(vb).pnt;
                auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{pa[0] / 2., pa[1] / 2.}, {pb[0] / 2., pb[1] / 2.}},
                                                          std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1);
                cuts.push_back({e, pc});
            }
        }
        auto r = unwrap(split_face(m, f, cuts));
        REQUIRE(r.faces.size() == 3u);
        std::size_t holed{};
        for (auto p : r.faces)
            if (m.face(p).wires.size() == 2)
            {
                ++holed;
                for (const auto &ce : m.wire(m.face(p).wires[0]).coedges) // the square of the face
                    REQUIRE(std::ranges::find(square_edges, ce.edge) != square_edges.end());
                REQUIRE(m.wire(m.face(p).wires[1]).coedges.size() == 8u); // around both diamonds
            }
            else
                REQUIRE(m.wire(m.face(p).wires[0]).coedges.size() == 4u);
        REQUIRE(holed == 1u);
    }
}

namespace
{
    /// Closed polygons of segment edges as cuts of the square face (x, y in [0, 2], uv = xy / 2);
    /// each polygon starts at its first point, the first edge of the first polygon is the first cut
    auto polygon_cuts(Model<T> &m, const std::vector<std::vector<point<T, 3>>> &polys) -> std::vector<FaceCut<T>>
    {
        std::vector<FaceCut<T>> cuts;
        for (const auto &d : polys)
        {
            std::vector<VertexId> vs;
            for (const auto &q : d)
                vs.push_back(unwrap(make_vertex(m, q)));
            for (std::size_t k = 0; k < d.size(); ++k)
            {
                const auto va = vs[k], vb = vs[(k + 1) % vs.size()];
                const auto e = unwrap(make_edge(m, va, vb));
                const auto &ed = m.edge(e);
                const auto pa = m.vertex(va).pnt, pb = m.vertex(vb).pnt;
                auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{pa[0] / 2., pa[1] / 2.}, {pb[0] / 2., pb[1] / 2.}},
                                                          std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1);
                cuts.push_back({e, pc});
            }
        }
        return cuts;
    }
}

TEST(tests_brep_split, hole_from_a_leftmost_vertex)
{
    // the hole around two diamonds touching at (1, 1) starts at the leftmost vertex (0.6, 1) of the
    // left one, which an even-odd test counts inside that diamond: the first cut is the edge
    // (0.8, 1.2) -> (0.6, 1), whose reverse half-edge opens the hole
    using namespace brep_tests;
    Model<T> m;
    const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
    const auto square_edges = explore<EdgeId>(m, ShapeId{f});
    const auto centre = unwrap(make_vertex(m, point<T, 3>{1., 1., 0.}));
    // the two diamonds share the vertex at the centre
    std::vector<FaceCut<T>> cuts;
    auto add = [&](const std::vector<point<T, 3>> &d) {
        std::vector<VertexId> vs;
        for (const auto &q : d)
            vs.push_back(norm(q - point<T, 3>{1., 1., 0.}) < 1e-12 ? centre : unwrap(make_vertex(m, q)));
        for (std::size_t k = 0; k < d.size(); ++k)
        {
            const auto va = vs[k], vb = vs[(k + 1) % vs.size()];
            const auto e = unwrap(make_edge(m, va, vb));
            const auto &ed = m.edge(e);
            const auto pa = m.vertex(va).pnt, pb = m.vertex(vb).pnt;
            auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{pa[0] / 2., pa[1] / 2.}, {pb[0] / 2., pb[1] / 2.}},
                                                      std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1);
            cuts.push_back({e, pc});
        }
    };
    add({{0.8, 1.2, 0.}, {0.6, 1., 0.}, {0.8, 0.8, 0.}, {1., 1., 0.}});
    add({{1., 1., 0.}, {1.2, 0.8, 0.}, {1.4, 1., 0.}, {1.2, 1.2, 0.}});
    auto r = unwrap(split_face(m, f, cuts));
    REQUIRE(r.faces.size() == 3u);
    std::size_t holed{};
    for (auto p : r.faces)
        if (m.face(p).wires.size() == 2)
        {
            ++holed;
            // the outer piece: its outer wire is the square of the face
            for (const auto &ce : m.wire(m.face(p).wires[0]).coedges)
                REQUIRE(std::ranges::find(square_edges, ce.edge) != square_edges.end());
            REQUIRE(m.wire(m.face(p).wires[1]).coedges.size() == 8u);
        }
        else
            REQUIRE(m.wire(m.face(p).wires[0]).coedges.size() == 4u);
    REQUIRE(holed == 1u);
}

TEST(tests_brep_split, hole_of_a_thin_ring)
{
    // an island [0.5, 1.5]² in a ring 1e-7 thick (thinner than a fixed offset of 1e-6 of the
    // diagonal of the face): the hole of the ring must stay with the ring
    using namespace brep_tests;
    Model<T> m;
    const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
    const T d = 1e-7;
    auto square = [](T a, T b) { return std::vector<point<T, 3>>{{a, a, 0.}, {b, a, 0.}, {b, b, 0.}, {a, b, 0.}}; };
    auto cuts = polygon_cuts(m, {square(0.5, 1.5), square(0.5 - d, 1.5 + d)});
    auto r = unwrap(split_face(m, f, cuts));
    REQUIRE(r.faces.size() == 3u);
    std::size_t one{}, two{};
    for (auto p : r.faces)
        (m.face(p).wires.size() == 1 ? one : two) += 1;
    REQUIRE(one == 1u); // the island
    REQUIRE(two == 2u); // the ring and the outer piece, each with its hole
}

TEST(tests_brep_split, cuts_crossing_each_other)
{
    // a square cut by a cross: the cuts were not split where they cross, the face stayed one face
    // of area 4 with slits; a third cut through the centre (0, 0.5) - (2, 1.5) shares its vertex
    using namespace brep_tests;
    for (bool diagonal : {false, true})
    {
        INFO("diagonal " << diagonal);
        Model<T> m;
        const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
        std::vector<FaceCut<T>> cuts;
        auto cut = [&](point<T, 3> a, point<T, 3> b) {
            const auto e = unwrap(make_edge(m, a, b));
            const auto &ed = m.edge(e);
            cuts.push_back({e, std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{a[0] / 2., a[1] / 2.}, {b[0] / 2., b[1] / 2.}},
                                                              std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1)});
        };
        cut({1., 0., 0.}, {1., 2., 0.});
        cut({0., 1., 0.}, {2., 1., 0.});
        if (diagonal)
            cut({0., 0.5, 0.}, {2., 1.5, 0.});
        auto r = unwrap(split_face(m, f, cuts));
        ASSERT_EQ(r.faces.size(), diagonal ? 6u : 4u);
        ASSERT_TRUE(r.dangling.empty());
        T sum{};
        for (auto p : r.faces)
        {
            const T a = area(m, p);
            // quarters, two of them cut in 0.75 and 0.25
            ASSERT_TRUE(std::abs(a - 1.) < 2e-2 || (diagonal && (std::abs(a - 0.75) < 2e-2 || std::abs(a - 0.25) < 2e-2)));
            sum += a;
            expect_valid(m, p);
        }
        ASSERT_NEAR(sum, 4., 1e-2);
    }

    // a box split by a compound of two crossing planar faces: the top and the bottom in four
    // (they were kept whole, with slits), the sides in two
    Model<T> m;
    const auto box = read(m, "box.stp"); // 2 x 3 x 4
    const auto p1 = unwrap(make_face(m, quad({1., -1., -1.}, {1., 4., -1.}, {1., -1., 5.}, {1., 4., 5.})));
    const auto p2 = unwrap(make_face(m, quad({-1., 1.5, -1.}, {3., 1.5, -1.}, {-1., 1.5, 5.}, {3., 1.5, 5.})));
    unwrap(split(m, box, ShapeId{unwrap(make_compound(m, std::vector<ShapeId>{p1, p2}))}));
    ASSERT_EQ(explore<FaceId>(m, box).size(), 16u);
    expect_valid(m, box);
}

TEST(tests_brep_split, face_of_a_compound)
{
    // a face held by a compound (not by a shell) was erased and lost from the compound; its pieces
    // also keep its name and its external id
    using namespace brep_tests;
    Model<T> m;
    const auto f = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
    m.setName(f, "plate");
    m.setExternalId(f, 42);
    const auto e = unwrap(make_edge(m, point<T, 3>{5., 0., 0.}, point<T, 3>{6., 0., 0.}));
    const auto c = unwrap(make_compound(m, std::vector<ShapeId>{e, f}));
    auto r = unwrap(split(m, ShapeId{c}, plane({1., 0., 0.}, {1., 0., 0.})));
    ASSERT_EQ(r.faces.size(), 1u);
    const auto &pieces = r.faces[0].second;
    ASSERT_EQ(pieces.size(), 2u);
    const auto &shapes = m.compound(c).shapes;
    ASSERT_EQ(shapes.size(), 3u);
    ASSERT_EQ(shapes[0], ShapeId{e});
    ASSERT_EQ(shapes[1], ShapeId{pieces[0]});
    ASSERT_EQ(shapes[2], ShapeId{pieces[1]});
    ASSERT_EQ(explore<FaceId>(m, ShapeId{c}).size(), 2u);
    for (auto p : pieces)
    {
        ASSERT_EQ(m.name(p), "plate");
        ASSERT_EQ(m.externalId(p), std::optional<std::int64_t>{42});
    }
}
