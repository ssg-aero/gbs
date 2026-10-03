#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include "tests_brep_helpers.h"

#include <unordered_map>
#include <unordered_set>

using namespace gbs;
using namespace gbs::brep;
using namespace brep_tests;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

TEST(tests_brep_model, ids)
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

    // invalid ids throw
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
    ASSERT_THROW((void)m.face(b.f[0]), BRepError);      // dead id
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

    // the ids stored in the remaining faces still point to live, consistent entities
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

    // original ids untouched
    ASSERT_TRUE(m.solid(b1.solid).outer == b1.shell);

    // copied solid: its ids were rewritten and point to the copied geometry
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
