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

TEST(tests_brep_model, append_to_itself)
{
    // a model appended to itself is copied once, names and parts with it (the copy loops grew the tables they
    // read: a crash, also from Python)
    auto b = make_box();
    auto &m = b.m;
    m.setName(b.f[2], "pressure side");
    m.setInstance(b.solid, InstanceOf<T>{m.newPart(), {}});
    const auto remap = m.append(m);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 16);
    ASSERT_EQ(m.count<ShapeType::Edge>(), 24);
    ASSERT_EQ(m.count<ShapeType::Face>(), 12);
    ASSERT_EQ(m.count<ShapeType::Solid>(), 2);
    const auto f2 = remap.map(b.f[2]);
    ASSERT_TRUE(f2.valid() && f2 != b.f[2]);
    ASSERT_EQ(m.name(f2), "pressure side");
    ASSERT_EQ(m.name(b.f[2]), "pressure side");
    ASSERT_TRUE(m.solid(remap.map(b.solid)).outer == remap.map(b.shell));
    ASSERT_NE(m.instance(remap.map(b.solid))->part, m.instance(b.solid)->part); // as for another model
}

TEST(tests_brep_model, names_and_const_access)
{
    // a name is a copy: erase() and compact() move the names (a view into them read freed memory)
    auto b = make_box();
    auto &m = b.m;
    m.setName(b.f[3], "trailing edge");
    const auto name = m.name(b.f[3]);
    m.erase(b.f[0]);
    const auto remap = m.compact();
    ASSERT_EQ(name, "trailing edge");
    ASSERT_EQ(m.name(remap.map(b.f[3])), "trailing edge");

    // the accessors are const or not as the model (deducing this)
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::as_const(m).face(FaceId{}))>>);
    static_assert(!std::is_const_v<std::remove_reference_t<decltype(m.face(FaceId{}))>>);
    static_assert(std::is_same_v<decltype(std::as_const(m).vertex(VertexId{})), const Vertex<T> &>);
    static_assert(std::is_same_v<decltype(m.compound(CompoundId{})), Compound &>);
    ASSERT_THROW((void)std::as_const(m).face(b.f[5]), BRepError); // five faces left: index 5 no longer exists
}

TEST(tests_brep_model, attributes)
{
    auto b = make_box();
    auto &m = b.m;

    // names and external ids on any kind of entity; empty / nullopt removes
    m.setName(b.solid, "bracket");
    m.setName(b.f[2], "top");
    m.setName(b.e[0], "x-edge");
    m.setExternalId(b.f[2], 1234);
    m.setExternalId(b.v[0], 7);
    ASSERT_EQ(m.name(b.solid), "bracket");
    ASSERT_EQ(m.name(b.f[2]), "top");
    ASSERT_EQ(m.name(b.f[3]), "");
    ASSERT_EQ(*m.externalId(b.f[2]), 1234);
    ASSERT_FALSE(m.externalId(b.f[3]).has_value());
    m.setName(b.e[0], "");
    ASSERT_EQ(m.name(b.e[0]), "");
    m.setExternalId(b.v[0], std::nullopt);
    ASSERT_FALSE(m.externalId(b.v[0]).has_value());
    // dead or invalid ids are refused
    ASSERT_THROW(m.setName(FaceId{99}, "x"), BRepError);
    ASSERT_THROW((void)m.name(VertexId{}), BRepError);

    // unit
    ASSERT_DOUBLE_EQ(m.unitScale(), 1.);
    m.setUnitScale(25.4);
    ASSERT_DOUBLE_EQ(m.unitScale(), 25.4);
    ASSERT_THROW(m.setUnitScale(0.), BRepError);
    m.setUnitScale(1.);

    // erase drops the attributes of the erased entity, compact renumbers the others
    m.setName(b.f[0], "bottom");
    m.setExternalId(b.f[5], 99);
    m.erase(b.f[0]);
    m.erase(b.f[1]);
    auto remap = m.compact();
    const auto top = remap.map(b.f[2]);
    ASSERT_EQ(top.index, 0);
    ASSERT_EQ(m.name(top), "top");
    ASSERT_EQ(*m.externalId(top), 1234);
    ASSERT_EQ(*m.externalId(remap.map(b.f[5])), 99);
    for (auto f : m.ids<FaceId>())
        ASSERT_NE(m.name(f), "bottom"); // the erased face's name did not move to a survivor

    // append copies the attributes with the ids remapped
    auto b2 = make_box();
    b2.m.setName(b2.solid, "copy");
    b2.m.setExternalId(b2.e[3], 42);
    auto r2 = m.append(b2.m);
    ASSERT_EQ(m.name(r2.map(b2.solid)), "copy");
    ASSERT_EQ(*m.externalId(r2.map(b2.e[3])), 42);
    ASSERT_EQ(m.name(remap.map(b.solid)), "bracket"); // the original ones are untouched

    // models in different units are not merged silently
    auto b3 = make_box();
    b3.m.setUnitScale(25.4);
    ASSERT_THROW((void)m.append(b3.m), BRepError);
}

TEST(tests_brep_model, instances)
{
    auto b = make_box();
    auto &m = b.m;
    const auto k1 = m.newPart(), k2 = m.newPart();
    ASSERT_EQ(k1, 1u);
    ASSERT_EQ(k2, 2u);
    const auto p = Placement<T>::frame({1., 2., 3.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.});
    m.setInstance(b.solid, InstanceOf<T>{k1, p});
    m.setInstance(b.f[0], InstanceOf<T>{k2, {}});
    ASSERT_EQ(m.instance(b.solid)->part, k1);
    ASSERT_LT(distance(m.instance(b.solid)->place.apply({0., 0., 0.}), point<T, 3>{1., 2., 3.}), 1e-15);
    ASSERT_FALSE(m.instance(b.f[1]).has_value());
    // keys come from newPart()
    ASSERT_THROW(m.setInstance(b.f[1], InstanceOf<T>{0, {}}), BRepError);
    ASSERT_THROW(m.setInstance(b.f[1], InstanceOf<T>{99, {}}), BRepError);
    m.setInstance(b.f[0], std::nullopt);
    ASSERT_FALSE(m.instance(b.f[0]).has_value());

    // erase drops, compact re-keys
    m.setInstance(b.f[2], InstanceOf<T>{k2, {}});
    m.erase(b.f[1]);
    auto r = m.compact();
    const auto f2 = r.map(b.f[2]);
    ASSERT_EQ(m.instance(f2)->part, k2);
    ASSERT_EQ(m.instance(r.map(b.solid))->part, k1);

    // append shifts the part keys of the other model after ours
    auto o = make_box();
    const auto ko = o.m.newPart();
    o.m.setInstance(o.solid, InstanceOf<T>{ko, p});
    auto ra = m.append(o.m);
    const auto shifted = m.instance(ra.map(o.solid))->part;
    ASSERT_NE(shifted, k1);
    ASSERT_NE(shifted, k2);
    ASSERT_EQ(m.newPart(), shifted + 1);
    ASSERT_THROW((void)m.instance(SolidId{999}), BRepError);
}
