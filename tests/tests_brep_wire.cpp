#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/curves>
#include <gbs/bscbuild.h>
#include <gbs/bscanalysis.h>

#include <cmath>
#include <limits>
#include <numbers>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;
    const T tol = 1e-6;

    std::shared_ptr<Curve<T, 3>> segment(const point<T, 3> &a, const point<T, 3> &b)
    {
        return std::make_shared<BSCurve<T, 3>>(build_segment(a, b));
    }

    // the chain of co-edges is consistent: end of i == start of i+1 (cyclically if closed)
    bool chained(const Model<T> &m, WireId wid)
    {
        return m.wire(wid).closed ? is_closed(m, wid) : is_chained(m, wid);
    }
}

// ---- vertices ---------------------------------------------------------------

TEST(tests_brep_wire, vertex)
{
    Model<T> m;
    auto v = make_vertex(m, {1., 2., 3.}, 1e-5);
    ASSERT_TRUE(v.has_value());
    ASSERT_NEAR(distance(m.vertex(*v).pnt, point<T, 3>{1., 2., 3.}), 0., 1e-15);
    ASSERT_DOUBLE_EQ(m.vertex(*v).tol, 1e-5);
    ASSERT_DOUBLE_EQ(m.vertex(unwrap(make_vertex(m, point<T, 3>{}))).tol, brep_default_tolerance<T>);

    auto bad_tol = make_vertex(m, {0., 0., 0.}, 0.);
    ASSERT_FALSE(bad_tol.has_value());
    ASSERT_TRUE(bad_tol.error().code == BuildErrc::InvalidTolerance);
    auto bad_pnt = make_vertex(m, {std::numeric_limits<T>::quiet_NaN(), 0., 0.});
    ASSERT_TRUE(bad_pnt.error().code == BuildErrc::InvalidPoint);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 2); // failures left the model unchanged

    ASSERT_THROW((void)unwrap(make_vertex(m, {0., 0., 0.}, -1.)), BRepError);
}

// ---- edges ------------------------------------------------------------------

TEST(tests_brep_wire, edge_from_curve)
{
    Model<T> m;

    // segment between points: arc-length parametrization, two vertices
    auto e = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{3., 4., 0.}));
    const auto &ed = m.edge(e);
    ASSERT_DOUBLE_EQ(ed.u1, 0.);
    ASSERT_DOUBLE_EQ(ed.u2, 5.);
    ASSERT_TRUE(ed.v1 != ed.v2);
    ASSERT_FALSE(ed.degenerate);
    ASSERT_NEAR(distance(m.vertex(ed.v2).pnt, point<T, 3>{3., 4., 0.}), 0., 1e-14);
    ASSERT_GE(m.vertex(ed.v1).tol, ed.tol); // vertex >= edge

    // closed curve: a single vertex
    auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(2.));
    auto ec = unwrap(make_edge(m, circle));
    ASSERT_TRUE(m.edge(ec).v1 == m.edge(ec).v2);

    // part of the circle, explicit bounds
    auto [c1, c2] = circle->bounds();
    auto ep = unwrap(make_edge(m, circle, c1, 0.5 * (c1 + c2)));
    ASSERT_TRUE(m.edge(ep).v1 != m.edge(ep).v2);

    // unbounded line: needs explicit bounds
    auto line = std::make_shared<Line<T, 3>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.});
    ASSERT_TRUE(make_edge(m, line).error().code == BuildErrc::UnboundedCurve);
    auto el = unwrap(make_edge(m, line, T(0.), T(2.)));
    ASSERT_NEAR(distance(m.vertex(m.edge(el).v2).pnt, point<T, 3>{2., 0., 0.}), 0., 1e-14);

    const auto n_vtx = m.count<ShapeType::Vertex>();
    const auto n_edg = m.count<ShapeType::Edge>();
    ASSERT_TRUE(make_edge(m, std::shared_ptr<Curve<T, 3>>{}).error().code == BuildErrc::NullCurve);
    ASSERT_TRUE(make_edge(m, circle, c2, c1).error().code == BuildErrc::InvalidBounds);
    ASSERT_TRUE(make_edge(m, circle, c1, c2 + 1.).error().code == BuildErrc::InvalidBounds);
    ASSERT_TRUE(make_edge(m, circle, c1, c2, 0.).error().code == BuildErrc::InvalidTolerance);
    ASSERT_TRUE(make_edge(m, point<T, 3>{1., 1., 1.}, point<T, 3>{1., 1., 1. + 1e-8}).error().code == BuildErrc::DegenerateEdge);
    // a curve collapsed on a point
    auto dot = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{1., 1., 1.}, {1., 1., 1.}}, std::vector<T>{0., 0., 1., 1.}, 1);
    ASSERT_TRUE(make_edge(m, dot).error().code == BuildErrc::DegenerateEdge);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), n_vtx);
    ASSERT_EQ(m.count<ShapeType::Edge>(), n_edg);
}

TEST(tests_brep_wire, edge_on_vertices)
{
    Model<T> m;
    auto va = unwrap(make_vertex(m, {0., 0., 0.}));
    auto vb = unwrap(make_vertex(m, {1., 0., 0.}));

    // straight edge between existing vertices
    auto e = unwrap(make_edge(m, va, vb));
    ASSERT_TRUE(m.edge(e).v1 == va);
    ASSERT_TRUE(m.edge(e).v2 == vb);

    // curve ending 0.5 tol away from vb: accepted, vb's tolerance contains the end
    auto crv = segment({0., 0., 0.}, {1., 0.5 * tol, 0.});
    auto [u1, u2] = crv->bounds();
    auto e2 = unwrap(make_edge(m, crv, u1, u2, va, vb, tol));
    ASSERT_TRUE(m.edge(e2).v2 == vb);
    ASSERT_GE(m.vertex(vb).tol, 0.5 * tol);
    ASSERT_GE(m.vertex(vb).tol, m.edge(e2).tol);

    // vertex tolerance larger than the edge tolerance also admits the end
    m.vertex(vb).tol = 1e-3;
    auto crv3 = segment({0., 0., 0.}, {1., 5e-4, 0.});
    auto [w1, w2] = crv3->bounds();
    ASSERT_TRUE(make_edge(m, crv3, w1, w2, va, vb, tol).has_value());

    // too far: rejected, model unchanged
    const auto n = m.count<ShapeType::Edge>();
    auto far = segment({0., 0., 0.}, {1., 0.1, 0.});
    auto [f1, f2] = far->bounds();
    auto r = make_edge(m, far, f1, f2, va, vb, tol);
    ASSERT_FALSE(r.has_value());
    ASSERT_TRUE(r.error().code == BuildErrc::VertexOffCurve);
    ASSERT_EQ(r.error().shapes.size(), 1);
    ASSERT_EQ(m.count<ShapeType::Edge>(), n);

    // same vertex at both ends, dead vertex
    ASSERT_TRUE(make_edge(m, va, va).error().code == BuildErrc::DegenerateEdge);
    auto vc = unwrap(make_vertex(m, {5., 5., 5.}));
    m.erase(vc);
    ASSERT_TRUE(make_edge(m, va, vc).error().code == BuildErrc::InvalidId);

    // closed edge on one given vertex
    auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1.));
    auto vs = unwrap(make_vertex(m, circle->begin()));
    auto [k1, k2] = circle->bounds();
    auto ec = unwrap(make_edge(m, circle, k1, k2, vs, vs));
    ASSERT_TRUE(m.edge(ec).v1 == vs && m.edge(ec).v2 == vs);
}

TEST(tests_brep_wire, degenerate_edge)
{
    Model<T> m;
    auto v = unwrap(make_vertex(m, {0., 0., 1.}, 1e-5));
    auto e = unwrap(make_degenerate_edge(m, v, T(0.), 2. * std::numbers::pi));
    const auto &ed = m.edge(e);
    ASSERT_TRUE(ed.degenerate);
    ASSERT_FALSE(ed.curve);
    ASSERT_TRUE(ed.v1 == v && ed.v2 == v);
    ASSERT_DOUBLE_EQ(ed.tol, 1e-5);
    ASSERT_NEAR(distance(edge_point(m, e, T(1.)), point<T, 3>{0., 0., 1.}), 0., 1e-15);
    ASSERT_TRUE(make_degenerate_edge(m, v, T(1.), T(1.)).error().code == BuildErrc::InvalidBounds);
    ASSERT_TRUE(make_degenerate_edge(m, VertexId{}).error().code == BuildErrc::InvalidId);
}

// ---- wires ------------------------------------------------------------------

TEST(tests_brep_wire, wire_square_any_order)
{
    // port of the legacy tests_topo.edge_wire: four independent segments, given
    // shuffled and with mixed senses; vertices are merged, the wire closes.
    Model<T> m;
    point<T, 3> p1{0., 0., 0.}, p2{1., 0., 0.}, p3{1., 1., 0.}, p4{0., 1., 0.};
    auto e12 = unwrap(make_edge(m, p1, p2));
    auto e32 = unwrap(make_edge(m, p3, p2)); // reversed sense
    auto e34 = unwrap(make_edge(m, p3, p4));
    auto e14 = unwrap(make_edge(m, p1, p4)); // reversed sense
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 8);

    auto w = unwrap(make_wire(m, std::vector{e12, e34, e14, e32}));
    const auto &wire = m.wire(w);
    ASSERT_TRUE(wire.closed);
    ASSERT_EQ(wire.coedges.size(), 4);
    ASSERT_TRUE(wire.coedges[0].edge == e12); // the first input edge leads…
    ASSERT_TRUE(wire.coedges[0].orient == Orientation::Forward); // …Forward
    ASSERT_TRUE(chained(m, w));
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 4); // 4 merged away
    ASSERT_EQ(explore<VertexId>(m, w).size(), 4);
    for (const auto &ce : wire.coedges)
        ASSERT_FALSE(ce.pcurve); // free wire

    // the reversed segments are used Reversed
    for (const auto &ce : wire.coedges)
        if (ce.edge == e32 || ce.edge == e14)
            ASSERT_TRUE(ce.orient == Orientation::Reversed);
        else
            ASSERT_TRUE(ce.orient == Orientation::Forward);
}

TEST(tests_brep_wire, wire_open)
{
    Model<T> m;
    auto a = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}));
    auto b = unwrap(make_edge(m, point<T, 3>{2., 0., 0.}, point<T, 3>{1., 0., 0.}));
    auto c = unwrap(make_edge(m, point<T, 3>{2., 0., 0.}, point<T, 3>{2., 1., 0.}));

    // middle edge given first: the chain is reoriented so that it stays Forward
    auto w = unwrap(make_wire(m, std::vector{b, c, a}));
    const auto &wire = m.wire(w);
    ASSERT_FALSE(wire.closed);
    ASSERT_TRUE(chained(m, w));
    ASSERT_FALSE(is_closed(m, w));
    ASSERT_EQ(wire.coedges.size(), 3);
    for (const auto &ce : wire.coedges)
        if (ce.edge == b)
            ASSERT_TRUE(ce.orient == Orientation::Forward);
    // b runs from (2,0,0) to (1,0,0), so the chain is c(rev) b a(rev)
    ASSERT_TRUE(wire.coedges[0].edge == c);
    ASSERT_TRUE(wire.coedges[2].edge == a);
    ASSERT_NEAR(distance(m.vertex(coedge_start(m, wire.coedges[0])).pnt, point<T, 3>{2., 1., 0.}), 0., 1e-14);

    // a single edge is an open wire
    auto w1 = unwrap(make_wire(m, std::vector{a}));
    ASSERT_FALSE(m.wire(w1).closed);
}

TEST(tests_brep_wire, wire_closed_single_edge)
{
    Model<T> m;
    auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1.));
    auto e = unwrap(make_edge(m, circle));
    auto w = unwrap(make_wire(m, std::vector{e}));
    ASSERT_TRUE(m.wire(w).closed);
    ASSERT_TRUE(is_closed(m, w));
    ASSERT_EQ(m.wire(w).coedges.size(), 1);
}

TEST(tests_brep_wire, wire_vertex_merge)
{
    Model<T> m;
    // two segments meeting with a gap of 0.4 tol: merged, survivor tolerance contains both ends
    auto a = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}, tol));
    auto b = unwrap(make_edge(m, point<T, 3>{1., 0.4 * tol, 0.}, point<T, 3>{1., 1., 0.}, tol));
    const auto absorbed = m.edge(b).v1;
    // an edge outside the wire also uses the vertex that will be absorbed
    auto other = unwrap(make_edge(m, m.edge(b).v1, unwrap(make_vertex(m, {3., 3., 3.}))));
    ASSERT_TRUE(m.edge(other).v1 == absorbed);

    auto w = unwrap(make_wire(m, std::vector{a, b}, tol));
    ASSERT_TRUE(chained(m, w));
    ASSERT_FALSE(m.alive(absorbed));
    const auto survivor = m.edge(a).v2;
    ASSERT_TRUE(m.edge(b).v1 == survivor);
    ASSERT_TRUE(m.edge(other).v1 == survivor); // redirected too
    ASSERT_GE(m.vertex(survivor).tol, 0.4 * tol);
    ASSERT_GE(m.vertex(survivor).tol, m.edge(b).tol);
    // survivor point not moved
    ASSERT_NEAR(distance(m.vertex(survivor).pnt, point<T, 3>{1., 0., 0.}), 0., 1e-15);
    // every curve end is inside its vertex ball
    for (auto eid : {a, b, other})
    {
        const auto &e = m.edge(eid);
        ASSERT_LE(distance(e.curve->value(e.u1), m.vertex(e.v1).pnt), m.vertex(e.v1).tol + 1e-15);
        ASSERT_LE(distance(e.curve->value(e.u2), m.vertex(e.v2).pnt), m.vertex(e.v2).tol + 1e-15);
    }

    // a larger fusion tolerance merges a gap the default would not
    Model<T> m2;
    auto c = unwrap(make_edge(m2, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}));
    auto d = unwrap(make_edge(m2, point<T, 3>{1., 1e-4, 0.}, point<T, 3>{1., 1., 0.}));
    ASSERT_TRUE(make_wire(m2, std::vector{c, d}).error().code == BuildErrc::Disconnected);
    ASSERT_EQ(m2.count<ShapeType::Vertex>(), 4); // unchanged after the failure
    auto w2 = unwrap(make_wire(m2, std::vector{c, d}, T(1e-3)));
    ASSERT_TRUE(chained(m2, w2));
    ASSERT_EQ(m2.count<ShapeType::Vertex>(), 3);
}

TEST(tests_brep_wire, wire_errors)
{
    Model<T> m;
    auto o = point<T, 3>{0., 0., 0.};
    auto a = unwrap(make_edge(m, o, point<T, 3>{1., 0., 0.}));
    auto b = unwrap(make_edge(m, o, point<T, 3>{0., 1., 0.}));
    auto c = unwrap(make_edge(m, o, point<T, 3>{0., 0., 1.}));
    const auto n_vtx = m.count<ShapeType::Vertex>();

    auto r = make_wire(m, std::vector{a, b, c}); // three edges meet at the origin
    ASSERT_TRUE(r.error().code == BuildErrc::Branching);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), n_vtx); // no merge applied
    ASSERT_EQ(m.count<ShapeType::Wire>(), 0);

    ASSERT_TRUE(make_wire(m, std::vector<EdgeId>{}).error().code == BuildErrc::EmptyInput);
    ASSERT_TRUE(make_wire(m, std::vector{a, a}).error().code == BuildErrc::DuplicateEdge);
    ASSERT_TRUE(make_wire(m, std::vector{a, EdgeId{999}}).error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(make_wire(m, std::vector{a}, T(0.)).error().code == BuildErrc::InvalidTolerance);
    auto vd = unwrap(make_vertex(m, {5., 5., 5.}));
    auto ed = unwrap(make_degenerate_edge(m, vd));
    ASSERT_TRUE(make_wire(m, std::vector{ed}).error().code == BuildErrc::DegenerateEdgeInWire);

    // two disjoint closed loops
    auto c1 = unwrap(make_edge(m, std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1.))));
    auto c2 = unwrap(make_edge(m, std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1., {10., 0., 0.}))));
    ASSERT_TRUE(make_wire(m, std::vector{c1, c2}).error().code == BuildErrc::Disconnected);

    // a path plus a separate segment
    auto far = unwrap(make_edge(m, point<T, 3>{7., 7., 7.}, point<T, 3>{8., 8., 8.}));
    ASSERT_TRUE(make_wire(m, std::vector{a, far}).error().code == BuildErrc::Disconnected);
}

TEST(tests_brep_wire, discretize_wire)
{
    // port of the legacy tests_topo.mesh_wire: boundary points of a square wire,
    // one point per vertex plus the interior points of each edge.
    Model<T> m;
    point<T, 3> p1{0., 0., 0.}, p2{1., 0., 0.}, p3{1., 1., 0.}, p4{0., 1., 0.};
    auto w = unwrap(make_wire(m, std::vector{
                                     unwrap(make_edge(m, p1, p2)), unwrap(make_edge(m, p2, p3)),
                                     unwrap(make_edge(m, p3, p4)), unwrap(make_edge(m, p4, p1))}));
    const T dm = 0.1;
    std::vector<point<T, 3>> boundary;
    for (const auto &ce : m.wire(w).coedges)
    {
        const auto &e = m.edge(ce.edge);
        const auto l = length(*e.curve, e.u1, e.u2);
        const auto n = static_cast<size_t>(std::round(l / dm)) + 1;
        auto pts = make_points(*e.curve, uniform_distrib_params(*e.curve, e.u1, e.u2, n));
        if (ce.orient == Orientation::Reversed)
            std::ranges::reverse(pts);
        boundary.push_back(m.vertex(coedge_start(m, ce)).pnt);
        boundary.insert(boundary.end(), std::next(pts.begin()), std::prev(pts.end()));
    }
    ASSERT_EQ(boundary.size(), 40);
    for (size_t i = 0; i < boundary.size(); ++i)
        ASSERT_NEAR(distance(boundary[i], boundary[(i + 1) % boundary.size()]), dm, 1e-9);
}
