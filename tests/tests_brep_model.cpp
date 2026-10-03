#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/bscbuild.h>

#include <unordered_map>
#include <unordered_set>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;

    // A unit box built by hand, the way a STEP reader or a builder would:
    // 8 vertices, 12 edges (degree-1 B-spline segments), 6 planar degree-1
    // B-spline faces with their pcurves, 1 shell, 1 solid.
    struct Box
    {
        Model<T> m;
        std::array<VertexId, 8> v;
        std::vector<EdgeId> e;
        std::array<FaceId, 6> f;
        ShellId shell;
        SolidId solid;
    };

    // Corner indices: bit 0 -> x, bit 1 -> y, bit 2 -> z
    point<T, 3> corner(unsigned i)
    {
        return {T(i & 1), T((i >> 1) & 1), T((i >> 2) & 1)};
    }

    EdgeId find_edge(const Model<T> &m, const std::vector<EdgeId> &edges, VertexId a, VertexId b, Orientation &orient)
    {
        for (auto eid : edges)
        {
            const auto &ed = m.edge(eid);
            if (ed.v1 == a && ed.v2 == b) { orient = Orientation::Forward; return eid; }
            if (ed.v1 == b && ed.v2 == a) { orient = Orientation::Reversed; return eid; }
        }
        throw std::runtime_error("edge not found");
    }

    Box make_box()
    {
        Box b;
        auto &m = b.m;

        for (unsigned i = 0; i < 8; ++i)
            b.v[i] = m.add(Vertex<T>{corner(i), 1e-7});

        // 12 edges along the axes, from the lower to the upper corner
        for (unsigned i = 0; i < 8; ++i)
            for (unsigned axis = 0; axis < 3; ++axis)
                if (!(i & (1u << axis)))
                {
                    unsigned j = i | (1u << axis);
                    auto crv = std::make_shared<BSCurve<T, 3>>(build_segment(corner(i), corner(j), true));
                    b.e.push_back(m.add(Edge<T>{crv, 0., 1., b.v[i], b.v[j], 1e-7}));
                }

        // 6 faces: for each axis, the two planes at coordinate 0 and 1.
        // Face corners in (u,v) order: (0,0) (1,0) (1,1) (0,1), surface S(u,v) = c00 + u (c10 - c00) + v (c01 - c00)
        size_t fi = 0;
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned side = 0; side < 2; ++side)
            {
                // (u, v) axes chosen so that the normal du ^ dv points outwards on both sides
                unsigned a1 = side ? (axis + 1) % 3 : (axis + 2) % 3;
                unsigned a2 = side ? (axis + 2) % 3 : (axis + 1) % 3;
                auto idx = [&](unsigned u, unsigned v) { return (side << axis) | (u << a1) | (v << a2); };
                std::array<unsigned, 4> c{idx(0, 0), idx(1, 0), idx(1, 1), idx(0, 1)};

                // degree-1 bilinear surface, poles u-fastest
                points_vector<T, 3> poles{corner(c[0]), corner(c[1]), corner(c[3]), corner(c[2])};
                std::vector<T> k{0., 0., 1., 1.};
                auto srf = std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);

                std::array<point<T, 2>, 4> uv{{{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}};
                Wire<T> w;
                for (unsigned s = 0; s < 4; ++s)
                {
                    unsigned s2 = (s + 1) % 4;
                    Orientation o;
                    auto eid = find_edge(m, b.e, b.v[c[s]], b.v[c[s2]], o);
                    // pcurve parametrized like the edge: from the edge's v1 to its v2
                    auto p_from = o == Orientation::Forward ? uv[s] : uv[s2];
                    auto p_to = o == Orientation::Forward ? uv[s2] : uv[s];
                    auto pc = std::make_shared<BSCurve<T, 2>>(build_segment(p_from, p_to, true));
                    w.coedges.push_back(CoEdge<T>{eid, o, pc});
                }
                w.closed = true;
                auto wid = m.add(std::move(w));
                b.f[fi++] = m.add(Face<T>{srf, {wid}, 1e-7, true});
            }

        Shell sh;
        for (auto fid : b.f)
            sh.faces.push_back(FaceUse{fid, Orientation::Forward});
        sh.closed = true;
        b.shell = m.add(std::move(sh));
        b.solid = m.add(Solid{b.shell, {}});
        return b;
    }
}

TEST(tests_brep_model, handles)
{
    VertexId v;
    ASSERT_FALSE(v.valid());
    ASSERT_EQ(v.index, npos);

    EdgeId e1{3}, e2{3}, e3{4};
    ASSERT_TRUE(e1 == e2);
    ASSERT_TRUE(e1 != e3);
    ASSERT_TRUE(e1 < e3);

    std::unordered_map<EdgeId, int> map;
    map[e1] = 1;
    map[e3] = 2;
    ASSERT_EQ(map.size(), 2);
    ASSERT_EQ(map[e2], 1);

    ShapeId s = FaceId{7};
    ASSERT_TRUE(shape_type(s) == ShapeType::Face);
    ASSERT_EQ(shape_index(s), 7);
    ASSERT_TRUE(valid(s));
    ASSERT_FALSE(valid(ShapeId{SolidId{}}));

    ASSERT_TRUE(reverse(Orientation::Forward) == Orientation::Reversed);
    ASSERT_TRUE(compose(Orientation::Reversed, Orientation::Reversed) == Orientation::Forward);
    ASSERT_TRUE(compose(Orientation::Forward, Orientation::Reversed) == Orientation::Reversed);
}

TEST(tests_brep_model, box_counts_and_access)
{
    auto b = make_box();
    const auto &m = b.m;

    ASSERT_EQ(m.count<ShapeType::Vertex>(), 8);
    ASSERT_EQ(m.count<ShapeType::Edge>(), 12);
    ASSERT_EQ(m.count<ShapeType::Wire>(), 6);
    ASSERT_EQ(m.count<ShapeType::Face>(), 6);
    ASSERT_EQ(m.count(ShapeType::Shell), 1);
    ASSERT_EQ(m.count(ShapeType::Solid), 1);
    ASSERT_EQ(m.count(ShapeType::Compound), 0);
    ASSERT_FALSE(m.empty());

    ASSERT_EQ(m.ids<EdgeId>().size(), 12);
    ASSERT_EQ(m.shell(b.shell).faces.size(), 6);
    ASSERT_TRUE(m.solid(b.solid).outer == b.shell);

    // every edge is used by exactly two co-edges, in opposite senses
    std::unordered_map<EdgeId, int> uses, forward;
    for (auto fid : b.f)
        for (auto wid : m.face(fid).wires)
            for (const auto &ce : m.wire(wid).coedges)
            {
                uses[ce.edge]++;
                if (ce.orient == Orientation::Forward) forward[ce.edge]++;
            }
    ASSERT_EQ(uses.size(), 12);
    for (auto &[eid, n] : uses)
    {
        ASSERT_EQ(n, 2);
        ASSERT_EQ(forward[eid], 1);
    }

    // wires are chained, and pcurves agree with the 3D curves through the surface
    for (auto fid : b.f)
    {
        const auto &f = m.face(fid);
        const auto &w = m.wire(f.wires.front());
        ASSERT_TRUE(w.closed);
        for (size_t i = 0; i < w.coedges.size(); ++i)
        {
            const auto &ce = w.coedges[i];
            const auto &next = w.coedges[(i + 1) % w.coedges.size()];
            ASSERT_TRUE(coedge_end(m, ce) == coedge_start(m, next));

            const auto &ed = m.edge(ce.edge);
            for (T t : {0., 0.25, 0.5, 1.})
            {
                // coedge_uv follows the co-edge sense, so compare with the edge point at the matching parameter
                T u_edge = ce.orient == Orientation::Forward ? t : ed.u1 + ed.u2 - t;
                auto uv = coedge_uv(m, ce, t);
                auto p_srf = f.surface->value(uv[0], uv[1]);
                auto p_crv = edge_point(m, ce.edge, u_edge);
                ASSERT_NEAR(norm(p_srf - p_crv), 0., ed.tol);
            }
            // extremities on the vertices
            ASSERT_NEAR(norm(edge_point(m, ce.edge, ed.u1) - m.vertex(ed.v1).pnt), 0., m.vertex(ed.v1).tol);
            ASSERT_NEAR(norm(edge_point(m, ce.edge, ed.u2) - m.vertex(ed.v2).pnt), 0., m.vertex(ed.v2).tol);
        }
    }

    // invalid handles throw
    ASSERT_THROW((void)m.vertex(VertexId{}), BRepError);
    ASSERT_THROW((void)m.edge(EdgeId{99}), BRepError);
    ASSERT_THROW((void)m.face(FaceId{6}), BRepError);
}

TEST(tests_brep_model, erase_and_compact)
{
    auto b = make_box();
    auto &m = b.m;

    // drop the solid, the shell and the first face (with its wire): an open box
    m.erase(b.solid);
    m.erase(ShapeId{b.shell});
    m.erase(m.face(b.f[0]).wires.front());
    m.erase(b.f[0]);
    ASSERT_FALSE(m.alive(b.f[0]));
    ASSERT_TRUE(m.alive(b.f[1]));
    ASSERT_EQ(m.count<ShapeType::Face>(), 5);
    ASSERT_EQ(m.count<ShapeType::Wire>(), 5);
    ASSERT_EQ(m.count<ShapeType::Solid>(), 0);
    ASSERT_THROW((void)m.face(b.f[0]), BRepError);      // dead handle
    ASSERT_EQ(m.capacity<ShapeType::Face>(), 6);  // tombstone still occupies its slot

    // erase a vertex nobody references after removing its three edges... here we only
    // check that erase does not cascade: edges still reference a dead vertex until compact
    auto edge_pts_before = std::vector<std::pair<point<T, 3>, point<T, 3>>>{};
    for (auto eid : m.ids<EdgeId>())
        edge_pts_before.emplace_back(m.vertex(m.edge(eid).v1).pnt, m.vertex(m.edge(eid).v2).pnt);

    auto remap = m.compact();
    ASSERT_EQ(m.capacity<ShapeType::Face>(), 5);
    ASSERT_EQ(m.count<ShapeType::Face>(), 5);
    ASSERT_FALSE(remap.map(b.f[0]).valid());
    ASSERT_TRUE(remap.map(b.f[1]).valid());
    ASSERT_EQ(remap.map(b.f[1]).index, 0);
    ASSERT_EQ(remap.map(b.f[5]).index, 4);
    ASSERT_FALSE(remap.map(b.solid).valid());
    ASSERT_FALSE(valid(remap.map(ShapeId{b.shell})));
    // vertices and edges were all alive: identity
    for (unsigned i = 0; i < 8; ++i)
        ASSERT_TRUE(remap.map(b.v[i]) == b.v[i]);

    // the handles stored in the remaining faces still point to live, consistent entities
    for (auto fid : m.ids<FaceId>())
    {
        const auto &f = m.face(fid);
        ASSERT_EQ(f.wires.size(), 1);
        const auto &w = m.wire(f.wires.front());
        ASSERT_EQ(w.coedges.size(), 4);
        for (size_t i = 0; i < 4; ++i)
            ASSERT_TRUE(coedge_end(m, w.coedges[i]) == coedge_start(m, w.coedges[(i + 1) % 4]));
    }
    size_t k = 0;
    for (auto eid : m.ids<EdgeId>())
    {
        ASSERT_NEAR(norm(m.vertex(m.edge(eid).v1).pnt - edge_pts_before[k].first), 0., 1e-12);
        ASSERT_NEAR(norm(m.vertex(m.edge(eid).v2).pnt - edge_pts_before[k].second), 0., 1e-12);
        ++k;
    }

    // compacting again is the identity
    auto remap2 = m.compact();
    for (auto fid : m.ids<FaceId>())
        ASSERT_TRUE(remap2.map(fid) == fid);
}

TEST(tests_brep_model, append)
{
    auto b1 = make_box();
    auto b2 = make_box();
    // make the second box recognizable: translate its vertices
    for (auto vid : b2.m.ids<VertexId>())
        b2.m.vertex(vid).pnt = b2.m.vertex(vid).pnt + point<T, 3>{2., 0., 0.};
    // and leave a tombstone in it, which must not be copied
    auto c = b2.m.add(Compound{{ShapeId{b2.solid}}});
    b2.m.erase(c);

    auto remap = b1.m.append(b2.m);
    const auto &m = b1.m;

    ASSERT_EQ(m.count<ShapeType::Vertex>(), 16);
    ASSERT_EQ(m.count<ShapeType::Edge>(), 24);
    ASSERT_EQ(m.count<ShapeType::Face>(), 12);
    ASSERT_EQ(m.count<ShapeType::Shell>(), 2);
    ASSERT_EQ(m.count<ShapeType::Solid>(), 2);
    ASSERT_EQ(m.count<ShapeType::Compound>(), 0);
    ASSERT_FALSE(remap.map(c).valid());

    // original handles untouched
    ASSERT_TRUE(m.solid(b1.solid).outer == b1.shell);

    // copied solid: its handles were rewritten and point to the copied geometry
    auto s2 = remap.map(b2.solid);
    ASSERT_TRUE(s2.valid());
    ASSERT_TRUE(m.solid(s2).outer == remap.map(b2.shell));
    const auto &sh2 = m.shell(m.solid(s2).outer);
    ASSERT_EQ(sh2.faces.size(), 6);
    std::unordered_set<VertexId> vtx2;
    for (const auto &fu : sh2.faces)
    {
        ASSERT_TRUE(fu.face.index >= 6);
        for (const auto &ce : m.wire(m.face(fu.face).wires.front()).coedges)
        {
            const auto &ed = m.edge(ce.edge);
            ASSERT_TRUE(ed.v1.index >= 8);
            vtx2.insert(ed.v1);
            vtx2.insert(ed.v2);
            ASSERT_GE(m.vertex(ed.v1).pnt[0], 2.);
        }
    }
    ASSERT_EQ(vtx2.size(), 8);

    // geometry is shared, not duplicated
    ASSERT_TRUE(m.face(remap.map(b2.f[0])).surface == b2.m.face(b2.f[0]).surface);
}
