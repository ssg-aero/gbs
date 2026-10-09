#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include "tests_brep_helpers.h"
#include "tests_brep_geom.h"

#include <limits>
#include <sstream>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using brep_tests::T;
    using brep_tests::pi;

    std::string dump(const CheckReport &r)
    {
        std::ostringstream os;
        for (const auto &e : r.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        return os.str();
    }

    std::shared_ptr<Curve<T, 2>> seg2(point<T, 2> a, point<T, 2> b)
    {
        return std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>(a, b));
    }

    std::vector<std::shared_ptr<Curve<T, 2>>> polygon2(std::vector<point<T, 2>> p)
    {
        std::vector<std::shared_ptr<Curve<T, 2>>> res;
        for (size_t i = 0; i < p.size(); ++i)
            res.push_back(seg2(p[i], p[(i + 1) % p.size()]));
        return res;
    }

    std::shared_ptr<Surface<T, 3>> xy_plane()
    {
        points_vector<T, 3> poles{{-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}};
        std::vector<T> k{-2., -2., 2., 2.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);
    }

    // first face wire of a face
    Wire<T> &outer(Model<T> &m, FaceId f) { return m.wire(m.face(f).wires.front()); }
}

// ---- valid models from every builder ------------------------------------------

TEST(tests_brep_check, builders_produce_valid_shapes)
{
    // hand-made box (fixture of the model tests)
    auto b = brep_tests::make_box();
    auto r = check(b.m, b.solid);
    INFO(dump(r));
    ASSERT_TRUE(r.ok());

    Model<T> m;
    std::vector<ShapeId> all;
    // natural faces, including seams, poles and apex
    for (auto srf : {brep_tests::cylinder(), brep_tests::sphere(), brep_tests::torus(), brep_tests::cone(),
                     brep_tests::nurbs_cylinder(), brep_tests::plane(), brep_tests::revolution(2. * pi)})
        all.push_back(unwrap(make_face(m, srf)));
    // closed single-face shells: sphere, torus
    all.push_back(m.add(Shell{{FaceUse{std::get<FaceId>(all[1])}}, true}));
    all.push_back(m.add(Shell{{FaceUse{std::get<FaceId>(all[2])}}, true}));
    // trimmed faces: disk with a projected pcurve, 2D loops with a hole
    auto disk = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1.))))}));
    all.push_back(unwrap(make_face(m, xy_plane(), disk)));
    all.push_back(unwrap(make_face(m, xy_plane(), polygon2({{-1.5, -1.}, {1.5, -1.}, {1.5, 1.}, {-1.5, 1.}}),
                                   {polygon2({{0., 0.5}, {0.5, -0.5}, {-0.5, -0.5}})})));
    // free wires, open and closed
    auto e1 = unwrap(make_edge(m, point<T, 3>{5., 0., 0.}, point<T, 3>{6., 0., 0.}));
    all.push_back(unwrap(make_wire(m, std::vector{e1})));
    auto c = m.add(Compound{all});

    r = check(m, c);
    INFO(dump(r));
    ASSERT_TRUE(r.ok());

    // geometry off: only topology
    ASSERT_TRUE(check(m, c, CheckOptions{.geometry = false}).ok());
}

// ---- one invalid case per issue -----------------------------------------------

TEST(tests_brep_check, vertex_and_edge_issues)
{
    {
        auto b = brep_tests::make_box();
        b.m.vertex(b.v[0]).tol = 0.;
        b.m.vertex(b.v[1]).pnt[0] = std::numeric_limits<T>::quiet_NaN();
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::InvalidTolerance, b.v[0]));
        ASSERT_TRUE(r.has(Issue::InvalidPoint, b.v[1]));
        ASSERT_GE(r.count(Issue::ToleranceOrderViolated), 3); // the three edges of v[0]
    }
    {
        auto b = brep_tests::make_box();
        auto &e0 = b.m.edge(b.e[0]);
        e0.curve.reset();
        std::swap(b.m.edge(b.e[1]).u1, b.m.edge(b.e[1]).u2);
        b.m.vertex(b.m.edge(b.e[2]).v2).pnt[2] += 1e-3; // off its curves
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::EdgeWithoutCurve, b.e[0]));
        ASSERT_TRUE(r.has(Issue::EdgeBoundsInverted, b.e[1]));
        ASSERT_TRUE(r.has(Issue::VertexOffCurve, b.e[2]));
    }
    {
        Model<T> m;
        auto va = unwrap(make_vertex(m, {0., 0., 0.}));
        auto vb = unwrap(make_vertex(m, {1., 0., 0.}));
        auto bad = m.add(Edge<T>{nullptr, 0., 1., va, vb, 1e-7, true}); // degenerate between two vertices
        ASSERT_TRUE(check(m, bad).has(Issue::InvalidDegenerateEdge, bad));
        auto ok = unwrap(make_degenerate_edge(m, va));
        ASSERT_TRUE(check(m, ok).ok());
    }
    {
        auto b = brep_tests::make_box();
        b.m.erase(b.v[3]);
        auto r = check(b.m, b.solid);
        ASSERT_EQ(r.count(Issue::DeadReference), 3); // the three edges of the corner
    }
}

TEST(tests_brep_check, tolerance_order)
{
    auto b = brep_tests::make_box();
    b.m.edge(b.e[4]).tol = 1e-5; // above its vertices (1e-7)
    b.m.face(b.f[0]).tol = 1e-3; // above its edges
    auto r = check(b.m, b.solid);
    ASSERT_TRUE(r.has(Issue::ToleranceOrderViolated, b.e[4]));
    ASSERT_EQ(r.count(Issue::ToleranceOrderViolated), 2 + 4); // two vertices of e[4], four edges of f[0]
}

TEST(tests_brep_check, wire_issues)
{
    {
        auto b = brep_tests::make_box();
        auto &w = outer(b.m, b.f[0]);
        std::swap(w.coedges[1], w.coedges[2]);
        ASSERT_TRUE(check(b.m, b.solid).has(Issue::WireNotChained, b.m.face(b.f[0]).wires.front()));
    }
    {
        auto b = brep_tests::make_box();
        outer(b.m, b.f[1]).coedges.pop_back(); // chained but open, as a face boundary
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::WireNotClosed, b.m.face(b.f[1]).wires.front()));
        ASSERT_TRUE(r.has(Issue::ShellNotClosed, b.shell));
    }
    {
        Model<T> m;
        auto e = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}));
        auto w = unwrap(make_wire(m, std::vector{e}));
        ASSERT_TRUE(check(m, w).ok()); // a free open wire is fine…
        m.wire(w).closed = true;       // …unless it claims to be closed
        ASSERT_TRUE(check(m, w).has(Issue::WireNotClosed, w));
    }
}

TEST(tests_brep_check, face_issues)
{
    {
        auto b = brep_tests::make_box();
        b.m.face(b.f[0]).surface.reset();
        b.m.face(b.f[1]).wires.clear();
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::FaceWithoutSurface, b.f[0]));
        ASSERT_TRUE(r.has(Issue::FaceWithoutBoundary, b.f[1]));
    }
    {
        auto b = brep_tests::make_box();
        auto &w0 = outer(b.m, b.f[0]);
        w0.coedges[0].pcurve.reset();
        auto &w1 = outer(b.m, b.f[1]);
        w1.coedges[0].pcurve = seg2({0., 0.}, {0.5, 0.}); // parameter range [0, 0.5] < edge range [0, 1]
        auto &w2 = outer(b.m, b.f[2]);
        auto pc = w2.coedges[0].pcurve;
        w2.coedges[0].pcurve = std::make_shared<BSCurve<T, 2>>(
            points_vector<T, 2>{pc->value(0.) + point<T, 2>{0., 0.01}, pc->value(1.) + point<T, 2>{0., 0.01}}, std::vector<T>{0., 0., 1., 1.}, 1);
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::CoEdgeWithoutPCurve, b.m.face(b.f[0]).wires.front()));
        ASSERT_TRUE(r.has(Issue::PCurveRangeMismatch, b.m.face(b.f[1]).wires.front()));
        ASSERT_TRUE(r.has(Issue::SameParameterViolated, b.m.face(b.f[2]).wires.front()));
    }
    {
        // outer boundary turned clockwise
        auto b = brep_tests::make_box();
        auto &w = outer(b.m, b.f[3]);
        std::ranges::reverse(w.coedges);
        for (auto &ce : w.coedges)
            ce.orient = reverse(ce.orient);
        auto r = check(b.m, b.solid);
        ASSERT_TRUE(r.has(Issue::OuterWireNotDirect, b.m.face(b.f[3]).wires.front()));
        ASSERT_TRUE(r.has(Issue::ShellNotOrientable, b.shell)); // the face boundary now runs against its neighbours
    }
}

TEST(tests_brep_check, hole_issues)
{
    Model<T> m;
    auto srf = xy_plane();
    const std::vector<point<T, 2>> big{{-1.9, -1.9}, {1.9, -1.9}, {1.9, 1.9}, {-1.9, 1.9}};
    auto tri = [](T x, T y) { return polygon2({{x, y + 0.3}, {x + 0.3, y - 0.3}, {x - 0.3, y - 0.3}}); };
    auto fa = unwrap(make_face(m, srf, polygon2(big), {tri(0., 0.)}));
    auto fb = unwrap(make_face(m, srf, polygon2(big), {tri(0.2, 0.)}));     // overlaps fa's hole
    auto fc = unwrap(make_face(m, srf, polygon2({{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}}), {tri(0., 0.)}));
    auto fd = unwrap(make_face(m, srf, polygon2(big), {tri(1.5, 1.5)}));    // a hole outside fc's outer
    ASSERT_TRUE(check(m, fa).ok());

    // hole turned counter-clockwise
    {
        auto &h = m.wire(m.face(fa).wires[1]);
        std::ranges::reverse(h.coedges);
        for (auto &ce : h.coedges)
            ce.orient = reverse(ce.orient);
        ASSERT_TRUE(check(m, fa).has(Issue::InnerWireNotIndirect, m.face(fa).wires[1]));
        std::ranges::reverse(h.coedges);
        for (auto &ce : h.coedges)
            ce.orient = reverse(ce.orient);
        ASSERT_TRUE(check(m, fa).ok());
    }
    // a second hole overlapping the first one
    m.face(fa).wires.push_back(m.face(fb).wires[1]);
    ASSERT_TRUE(check(m, fa).has(Issue::InnerWiresOverlap, m.face(fb).wires[1]));
    m.face(fa).wires.pop_back();
    // a hole outside the outer boundary
    m.face(fc).wires.push_back(m.face(fd).wires[1]);
    ASSERT_TRUE(check(m, fc).has(Issue::InnerWireOutsideOuter, m.face(fd).wires[1]));
}

TEST(tests_brep_check, shell_issues)
{
    {
        // non-manifold: a seventh face reusing the wire of face 0
        auto b = brep_tests::make_box();
        auto f7 = b.m.add(Face<T>{b.m.face(b.f[0]).surface, {b.m.face(b.f[0]).wires.front()}, 1e-7, true});
        b.m.shell(b.shell).faces.push_back(FaceUse{f7, Orientation::Forward});
        auto r = check(b.m, b.solid);
        ASSERT_EQ(r.count(Issue::NonManifoldEdge), 4);
    }
    {
        // open shell bounding a solid, open shell flagged closed, open free shell
        auto b = brep_tests::make_box();
        b.m.shell(b.shell).faces.pop_back();
        b.m.shell(b.shell).closed = false;
        ASSERT_TRUE(check(b.m, b.solid).has(Issue::ShellNotClosed, b.shell)); // bounds a solid
        ASSERT_TRUE(check(b.m, b.shell).ok());                                // a free open shell is fine
        b.m.shell(b.shell).closed = true;
        ASSERT_TRUE(check(b.m, b.shell).has(Issue::ShellNotClosed, b.shell));
    }
    {
        // one flipped face use
        auto b = brep_tests::make_box();
        b.m.shell(b.shell).faces[2].orient = Orientation::Reversed;
        auto r = check(b.m, b.solid);
        ASSERT_EQ(r.count(Issue::ShellNotOrientable), 4); // its four edges
    }
    {
        // dead face in a shell, dead shell in a solid, dead member of a compound
        auto b = brep_tests::make_box();
        auto c = b.m.add(Compound{{ShapeId{b.solid}, ShapeId{FaceId{99}}}});
        b.m.erase(b.f[5]);
        auto r = check(b.m, c);
        ASSERT_TRUE(r.has(Issue::DeadReference, b.shell));
        ASSERT_TRUE(r.has(Issue::DeadReference, c));
        b.m.erase(b.shell);
        ASSERT_TRUE(check(b.m, b.solid).has(Issue::DeadReference, b.solid));
        ASSERT_TRUE(check(b.m, b.shell).has(Issue::DeadReference, b.shell)); // dead root
    }
}
