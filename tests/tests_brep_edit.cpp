#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/bselementary.h>
#include "tests_brep_geom.h"

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using namespace brep_tests;
    constexpr auto Fw = Orientation::Forward;
    constexpr auto Rv = Orientation::Reversed;
}

TEST(tests_brep_edit, split_edge_in_free_wires)
{
    Model<T> m;
    const auto a = unwrap(make_vertex(m, {0., 0., 0.})), b = unwrap(make_vertex(m, {2., 0., 0.})), c = unwrap(make_vertex(m, {2., 1., 0.}));
    const auto ab = unwrap(make_edge(m, a, b)), bc = unwrap(make_edge(m, b, c)), ca = unwrap(make_edge(m, c, a));
    m.setName(ab, "AB");
    m.setExternalId(ab, 42);
    const auto w1 = unwrap(make_wire_ordered(m, {{ab, Fw}, {bc, Fw}, {ca, Fw}}));
    const auto w2 = unwrap(make_wire_ordered(m, {{ca, Rv}, {bc, Rv}, {ab, Rv}}));

    auto s = unwrap(split_edge(m, ab, 0.5));
    ASSERT_EQ(s.first, ab);
    ASSERT_LT(distance(m.vertex(s.vertex).pnt, point<T, 3>{0.5, 0., 0.}), 1e-15);
    ASSERT_EQ(m.edge(ab).v2, s.vertex);
    ASSERT_EQ(m.edge(s.second).v1, s.vertex);
    ASSERT_EQ(m.edge(s.second).v2, b);
    ASSERT_EQ(m.edge(s.second).curve, m.edge(ab).curve); // shared curve
    ASSERT_EQ(m.name(s.second), "AB");
    ASSERT_EQ(*m.externalId(s.second), 42);
    // both wires updated in their own sense, still closed
    ASSERT_EQ(m.wire(w1).coedges.size(), 4u);
    ASSERT_EQ(m.wire(w1).coedges[0].edge, ab);
    ASSERT_EQ(m.wire(w1).coedges[1].edge, s.second);
    ASSERT_EQ(m.wire(w2).coedges[2].edge, s.second);
    ASSERT_EQ(m.wire(w2).coedges[3].edge, ab);
    ASSERT_TRUE(is_closed(m, w1));
    ASSERT_TRUE(is_closed(m, w2));

    // errors leave the model unchanged
    const auto n = m.ids<EdgeId>().size();
    ASSERT_EQ(split_edge(m, ab, 0.).error().code, BuildErrc::InvalidBounds);
    ASSERT_EQ(split_edge(m, ab, 0.7).error().code, BuildErrc::InvalidBounds); // beyond the new u2
    ASSERT_EQ(split_edge(m, EdgeId{999}, 0.5).error().code, BuildErrc::InvalidId);
    ASSERT_EQ(split_edge(m, unwrap(make_degenerate_edge(m, a)), 0.5).error().code, BuildErrc::DegenerateEdge);
    ASSERT_EQ(m.ids<EdgeId>().size(), n + 1);
}

TEST(tests_brep_edit, split_edge_of_a_face_keeps_it_valid)
{
    // a disk bounded by a closed circle edge, split in two arcs after the face is built
    Model<T> m;
    const ax2<T, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(1., 0., 2 * pi, ax));
    const auto v = unwrap(make_vertex(m, {1., 0., 0.}));
    const auto e = unwrap(make_edge(m, circle, 0., 2 * pi, v, v));
    auto plane = quad({-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.});
    const auto f = unwrap(make_face(m, plane, unwrap(make_wire_ordered(m, {{e, Fw}}))));
    const T area = uv_signed_area(m, m.face(f).wires[0], 512);
    auto s = unwrap(split_edge(m, e, pi));
    ASSERT_LT(distance(m.vertex(s.vertex).pnt, point<T, 3>{-1., 0., 0.}), 1e-12);
    const auto w = m.face(f).wires[0];
    ASSERT_EQ(m.wire(w).coedges.size(), 2u);
    ASSERT_EQ(m.wire(w).coedges[0].pcurve, m.wire(w).coedges[1].pcurve); // the pcurve is shared
    ASSERT_NEAR(uv_signed_area(m, w, 256), area, 1e-4);
    auto c = check(m, f);
    INFO(c.entries.size());
    ASSERT_TRUE(c.ok());
}

TEST(tests_brep_edit, projection_near_the_seam)
{
    // a point of a closed surface just before its seam: the bracketed search alone stops on the
    // other border of the parametric rectangle
    const double R = 5., us = pi / 8;
    const ax2<double, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    const auto s = build_cylinder<double>(R, ax, -1., 11., us, us + 2 * pi);
    const point<double, 3> p{R * std::cos(0.3), R * std::sin(0.3), 4.5};
    const auto r = gbs::brep::detail::project_point_global(s, p, 1e-8);
    ASSERT_LT(r[2], 1e-10);
    ASSERT_NEAR(r[1], 4.5, 1e-10);
}
