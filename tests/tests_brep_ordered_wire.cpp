#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/bselementary.h>
#include "tests_brep_geom.h"

#include <map>
#include <sstream>

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

    std::string dump(const CheckReport &r)
    {
        std::ostringstream os;
        for (const auto &e : r.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        return os.str();
    }

    const ax2<T, 3> ax_z{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};

    // Closed edge on a full circle of radius r at height z, starting and ending at v.
    EdgeId circle_edge(Model<T> &m, T r, T z, VertexId v)
    {
        auto c = std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(r, 0., 2 * pi, ax2<T, 3>{{{0., 0., z}, {0., 0., 1.}, {1., 0., 0.}}}));
        return unwrap(make_edge(m, c, 0., 2 * pi, v, v));
    }

    // Plane z = h with normal +z, large enough for the disks.
    std::shared_ptr<Surface<T, 3>> plane_z(T z, T s = 4.)
    {
        return quad({-s, -s, z}, {s, -s, z}, {-s, s, z}, {s, s, z});
    }

    std::size_t degenerate_coedges(const Model<T> &m, WireId w)
    {
        return static_cast<std::size_t>(std::ranges::count_if(m.wire(w).coedges, [&](const CoEdge<T> &ce) { return m.edge(ce.edge).degenerate; }));
    }

    // Shell of the given uses, checked closed and orientable; returns its solid.
    SolidId solid_of(Model<T> &m, const std::vector<FaceUse> &uses)
    {
        const auto sh = m.add(Shell{uses});
        REQUIRE(is_closed(m, sh));
        REQUIRE(is_orientable(m, sh));
        return unwrap(make_solid(m, sh));
    }
}

TEST(tests_brep_ordered_wire, make_wire_ordered)
{
    Model<T> m;
    const auto a = unwrap(make_vertex(m, {0., 0., 0.})), b = unwrap(make_vertex(m, {1., 0., 0.})), c = unwrap(make_vertex(m, {1., 1., 0.}));
    const auto ab = unwrap(make_edge(m, a, b)), bc = unwrap(make_edge(m, b, c)), ca = unwrap(make_edge(m, c, a));

    // order and senses kept, closed
    auto w = unwrap(make_wire_ordered(m, {{ab, Fw}, {bc, Fw}, {ca, Fw}}));
    ASSERT_TRUE(m.wire(w).closed);
    ASSERT_TRUE(is_closed(m, w));
    // reversed loop, open chain
    auto r = unwrap(make_wire_ordered(m, {{ca, Rv}, {bc, Rv}}));
    ASSERT_FALSE(m.wire(r).closed);
    ASSERT_EQ(m.wire(r).coedges[0].orient, Rv);
    // a seam: the same edge twice in opposite senses
    auto s = unwrap(make_wire_ordered(m, {{ab, Fw}, {ab, Rv}}));
    ASSERT_TRUE(m.wire(s).closed);

    const auto n_edges = m.ids<EdgeId>().size();
    const auto n_wires = m.ids<WireId>().size();
    auto code = [&](std::vector<OrientedEdge> ce) { return make_wire_ordered(m, ce).error().code; };
    ASSERT_EQ(code({{ab, Fw}, {ca, Fw}}), BuildErrc::Disconnected);           // b is not c
    ASSERT_EQ(code({{ab, Fw}, {ab, Fw}}), BuildErrc::DuplicateEdge);          // same sense twice
    ASSERT_EQ(code({{ab, Fw}, {ab, Rv}, {ab, Fw}}), BuildErrc::DuplicateEdge); // three uses
    ASSERT_EQ(code({}), BuildErrc::EmptyInput);
    ASSERT_EQ(code({{EdgeId{999}, Fw}}), BuildErrc::InvalidId);
    const auto d = unwrap(make_degenerate_edge(m, a));
    ASSERT_EQ(code({{d, Fw}}), BuildErrc::DegenerateEdgeInWire);
    ASSERT_EQ(m.ids<EdgeId>().size(), n_edges + 1);
    ASSERT_EQ(m.ids<WireId>().size(), n_wires);

    // vertices within tolerance are merged at the junctions
    const auto p = unwrap(make_vertex(m, {5., 0., 0.})), q = unwrap(make_vertex(m, {6., 0., 0.}));
    const auto q2 = unwrap(make_vertex(m, {6., 1e-8, 0.})), p2 = unwrap(make_vertex(m, {5., 1e-8, 0.}));
    const auto pq = unwrap(make_edge(m, p, q)), qp = unwrap(make_edge(m, std::make_shared<BSCurve<T, 3>>(build_segment<T, 3>({6., 1e-8, 0.}, {5.5, 1., 0.}, true)), 0., 1., q2, unwrap(make_vertex(m, {5.5, 1., 0.}))));
    auto back = unwrap(make_edge(m, m.edge(qp).v2, p2));
    auto merged = unwrap(make_wire_ordered(m, {{pq, Fw}, {qp, Fw}, {back, Fw}}));
    ASSERT_TRUE(m.wire(merged).closed);
    ASSERT_TRUE(is_closed(m, merged));
    ASSERT_FALSE(m.alive(q2)); // q (smaller index) survives
    ASSERT_FALSE(m.alive(p2));
}

TEST(tests_brep_ordered_wire, cylinder_with_seam)
{
    const T R = 1.5, h = 2.;
    for (int variant = 0; variant < 3; ++variant)
    {
        Model<T> m;
        auto lateral = std::make_shared<BSSurfaceRational<T, 3>>(build_cylinder<T>(R, ax_z, 0., h));
        const auto vb = unwrap(make_vertex(m, {R, 0., 0.})), vt = unwrap(make_vertex(m, {R, 0., h}));
        const auto bottom = circle_edge(m, R, 0., vb), top = circle_edge(m, R, h, vt);
        const auto seam = unwrap(make_edge(m, vb, vt));
        // loop counter-clockwise seen from outside (the surface normal), as STEP writes it; the
        // same loop starting on the seam; and the loop in the other sense (face normal inwards)
        const std::vector<std::vector<OrientedEdge>> loops{
            {{bottom, Fw}, {seam, Fw}, {top, Rv}, {seam, Rv}},
            {{seam, Fw}, {top, Rv}, {seam, Rv}, {bottom, Fw}},
            {{seam, Fw}, {top, Fw}, {seam, Rv}, {bottom, Rv}},
        };
        const auto w = unwrap(make_wire_ordered(m, loops[variant]));
        const auto use = unwrap(make_face_use(m, lateral, w));
        ASSERT_EQ(use.orient, variant == 2 ? Rv : Fw);
        auto c = check(m, use.face);
        INFO(dump(c));
    ASSERT_TRUE(c.ok());
        // the two seam co-edges are on both sides of the parametric rectangle
        std::vector<T> seam_u;
        for (const auto &ce : m.wire(w).coedges)
            if (ce.edge == seam)
                seam_u.push_back(ce.pcurve->value(0.5 * (m.edge(seam).u1 + m.edge(seam).u2))[0]);
        ASSERT_EQ(seam_u.size(), 2u);
        ASSERT_NEAR(std::abs(seam_u[0] - seam_u[1]), 2 * pi, 1e-9);
        ASSERT_NEAR(uv_signed_area(m, w), 2 * pi * h, 1e-6);

        if (variant == 2)
            continue;
        // closed by two disks: bottom seen from below (circle reversed), top from above
        const auto wb = unwrap(make_wire_ordered(m, {{bottom, Rv}}));
        const auto wt = unwrap(make_wire_ordered(m, {{top, Fw}}));
        const auto ub = unwrap(make_face_use(m, plane_z(0.), wb));
        const auto ut = unwrap(make_face_use(m, plane_z(h), wt));
        ASSERT_EQ(ub.orient, Rv); // the plane normal is +z, the bottom face looks down
        ASSERT_EQ(ut.orient, Fw);
        const auto sh = m.add(Shell{{use, ub, ut}});
        ASSERT_TRUE(is_closed(m, sh));
        ASSERT_TRUE(is_orientable(m, sh));
        ASSERT_NEAR(signed_volume(m, sh, 256), pi * R * R * h, 1e-3); // outwards without any turn
        const auto so = unwrap(make_solid(m, sh));
        auto cs = check(m, so);
        INFO(dump(cs));
    ASSERT_TRUE(cs.ok());
    }
}

TEST(tests_brep_ordered_wire, sphere_poles)
{
    Model<T> m;
    const T R = 2.;
    auto sphere = std::make_shared<BSSurfaceRational<T, 3>>(build_sphere<T>(R, ax_z));
    const auto S = unwrap(make_vertex(m, {0., 0., -R})), N = unwrap(make_vertex(m, {0., 0., R}));
    // the meridian u = 0, from the south pole to the north pole
    auto meridian = std::make_shared<BSCurveRational<T, 3>>(build_ellipse_arc<T, 3>(R, R, -pi / 2, pi / 2, {0., 0., 0.}, {1., 0., 0.}, {0., 0., 1.}));
    const auto seam = unwrap(make_edge(m, meridian, -pi / 2, pi / 2, S, N));
    // STEP: a sphere bounded by its seam only, both senses; the poles are missing in (u, v)
    const auto w = unwrap(make_wire_ordered(m, {{seam, Fw}, {seam, Rv}}));
    const auto use = unwrap(make_face_use(m, sphere, w));
    ASSERT_EQ(m.wire(w).coedges.size(), 4u);
    ASSERT_EQ(degenerate_coedges(m, w), 2u); // one per pole
    ASSERT_NEAR(std::abs(uv_signed_area(m, w)), 2 * pi * pi, 1e-6);
    auto c = check(m, use.face);
    INFO(dump(c));
    ASSERT_TRUE(c.ok());
    const auto so = solid_of(m, {use});
    ASSERT_NEAR(signed_volume(m, m.solid(so).outer, 256), 4. / 3. * pi * R * R * R, 1e-2);
    auto cs = check(m, so);
    INFO(dump(cs));
    ASSERT_TRUE(cs.ok());
}

TEST(tests_brep_ordered_wire, cone_apex)
{
    Model<T> m;
    const T R = 1., alpha = pi / 6, H = R / std::tan(alpha);
    // apex at z = -H, base circle of radius R at z = 0
    auto cone = std::make_shared<BSSurfaceRational<T, 3>>(build_cone<T>(R, alpha, ax_z, -H, 0.));
    const auto A = unwrap(make_vertex(m, {0., 0., -H})), B = unwrap(make_vertex(m, {R, 0., 0.}));
    const auto base = circle_edge(m, R, 0., B);
    const auto seam = unwrap(make_edge(m, A, B));
    // counter-clockwise from outside: up the seam, base backwards, down the seam; the apex is missing
    const auto w = unwrap(make_wire_ordered(m, {{seam, Fw}, {base, Rv}, {seam, Rv}}));
    const auto use = unwrap(make_face_use(m, cone, w));
    ASSERT_EQ(use.orient, Fw);
    ASSERT_EQ(m.wire(w).coedges.size(), 4u);
    ASSERT_EQ(degenerate_coedges(m, w), 1u);
    auto c = check(m, use.face);
    INFO(dump(c));
    ASSERT_TRUE(c.ok());
    const auto top = unwrap(make_face_use(m, plane_z(0.), unwrap(make_wire_ordered(m, {{base, Fw}}))));
    ASSERT_EQ(top.orient, Fw);
    const auto so = solid_of(m, {use, top});
    ASSERT_NEAR(signed_volume(m, m.solid(so).outer, 256), pi * R * R * H / 3., 1e-3);
    auto cs = check(m, so);
    INFO(dump(cs));
    ASSERT_TRUE(cs.ok());
}

TEST(tests_brep_ordered_wire, torus_two_seams)
{
    Model<T> m;
    const T R = 3., r = 1.;
    auto torus = std::make_shared<BSSurfaceRational<T, 3>>(build_torus<T>(R, r, ax_z));
    const auto V = unwrap(make_vertex(m, {R + r, 0., 0.}));
    const auto major = circle_edge(m, R + r, 0., V); // v = 0
    auto minor_c = std::make_shared<BSCurveRational<T, 3>>(build_ellipse_arc<T, 3>(r, r, 0., 2 * pi, {R, 0., 0.}, {1., 0., 0.}, {0., 0., 1.}));
    const auto minor = unwrap(make_edge(m, minor_c, 0., 2 * pi, V, V)); // u = 0
    // every co-edge lies on a seam: the sides follow from the neighbours
    const auto w = unwrap(make_wire_ordered(m, {{major, Fw}, {minor, Fw}, {major, Rv}, {minor, Rv}}));
    const auto use = unwrap(make_face_use(m, torus, w));
    ASSERT_EQ(use.orient, Fw);
    ASSERT_EQ(degenerate_coedges(m, w), 0u);
    ASSERT_NEAR(uv_signed_area(m, w), 4 * pi * pi, 1e-6);
    auto c = check(m, use.face);
    INFO(dump(c));
    ASSERT_TRUE(c.ok());
    const auto so = solid_of(m, {use});
    ASSERT_NEAR(signed_volume(m, m.solid(so).outer, 256), 2 * pi * pi * R * r * r, 1e-2);
}

TEST(tests_brep_ordered_wire, box_face_senses)
{
    // a box read as STEP would give it: shared vertices and edges, every loop counter-clockwise
    // seen from outside, planes whose normals point inwards on half of the faces
    Model<T> m;
    std::array<VertexId, 8> v;
    for (unsigned i = 0; i < 8; ++i)
        v[i] = unwrap(make_vertex(m, box_corner(i)));
    std::map<std::pair<unsigned, unsigned>, EdgeId> edges;
    for (unsigned i = 0; i < 8; ++i)
        for (unsigned b = 0; b < 3; ++b)
            if (!(i & (1u << b)))
                edges[{i, i | (1u << b)}] = unwrap(make_edge(m, v[i], v[i | (1u << b)]));
    auto coedge = [&](unsigned p, unsigned q) {
        return p < q ? OrientedEdge{edges.at({p, q}), Fw} : OrientedEdge{edges.at({q, p}), Rv};
    };
    std::vector<FaceUse> uses;
    for (unsigned axis = 0; axis < 3; ++axis)
        for (unsigned side = 0; side < 2; ++side)
        {
            const unsigned a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
            auto idx = [&](unsigned u, unsigned w) { return (side << axis) | (u << a1) | (w << a2); };
            // the plane normal is +axis; the face looks towards +axis on side 1, -axis on side 0
            auto srf = quad(box_corner(idx(0, 0)), box_corner(idx(1, 0)), box_corner(idx(0, 1)), box_corner(idx(1, 1)));
            std::vector<unsigned> c{idx(0, 0), idx(1, 0), idx(1, 1), idx(0, 1)};
            if (side == 0)
                std::ranges::reverse(c);
            std::vector<OrientedEdge> loop;
            for (std::size_t k = 0; k < 4; ++k)
                loop.push_back(coedge(c[k], c[(k + 1) % 4]));
            auto use = unwrap(make_face_use(m, srf, unwrap(make_wire_ordered(m, loop))));
            ASSERT_EQ(use.orient, side == 1 ? Fw : Rv);
            uses.push_back(use);
        }
    const auto sh = m.add(Shell{uses});
    ASSERT_TRUE(is_closed(m, sh));
    ASSERT_TRUE(is_orientable(m, sh));
    ASSERT_NEAR(signed_volume(m, sh), 1., 1e-12); // outwards as read, no turn needed
    auto cs = check(m, unwrap(make_solid(m, sh)));
    INFO(dump(cs));
    ASSERT_TRUE(cs.ok());
}

TEST(tests_brep_ordered_wire, failures_leave_the_model_unchanged)
{
    Model<T> m;
    const T R = 1.;
    auto lateral = std::make_shared<BSSurfaceRational<T, 3>>(build_cylinder<T>(R, ax_z, 0., 1.));
    const auto vb = unwrap(make_vertex(m, {R, 0., 0.})), vt = unwrap(make_vertex(m, {R, 0., 1.}));
    const auto bottom = circle_edge(m, R, 0., vb), top = circle_edge(m, R, 1., vt);
    // a band without seam edge: its circles go around the seam
    const auto wb = unwrap(make_wire_ordered(m, {{bottom, Fw}})), wt = unwrap(make_wire_ordered(m, {{top, Rv}}));
    const auto n_edges = m.ids<EdgeId>().size();
    auto r = make_face_use(m, lateral, wb, std::vector<WireId>{wt});
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::CrossesSeam);
    ASSERT_EQ(m.ids<EdgeId>().size(), n_edges);
    ASSERT_FALSE(m.wire(wb).coedges[0].pcurve);

    // a sphere seam whose poles are inserted, then a hole outside: the inserted edges are erased
    auto sphere = std::make_shared<BSSurfaceRational<T, 3>>(build_sphere<T>(R, ax_z));
    const auto S = unwrap(make_vertex(m, {0., 0., -R})), N = unwrap(make_vertex(m, {0., 0., R}));
    auto meridian = std::make_shared<BSCurveRational<T, 3>>(build_ellipse_arc<T, 3>(R, R, -pi / 2, pi / 2, {0., 0., 0.}, {1., 0., 0.}, {0., 0., 1.}));
    const auto seam = unwrap(make_edge(m, meridian, -pi / 2, pi / 2, S, N));
    const auto ws = unwrap(make_wire_ordered(m, {{seam, Fw}, {seam, Rv}}));
    const auto off = unwrap(make_wire_ordered(m, {{circle_edge(m, 0.5, 5., unwrap(make_vertex(m, {0.5, 0., 5.}))), Fw}}));
    const auto n2 = m.ids<EdgeId>().size();
    auto r2 = make_face_use(m, sphere, ws, std::vector<WireId>{off});
    ASSERT_FALSE(r2.has_value());
    ASSERT_EQ(m.ids<EdgeId>().size(), n2);
    ASSERT_EQ(m.wire(ws).coedges.size(), 2u);
}
