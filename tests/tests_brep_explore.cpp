#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include "tests_brep_helpers.h"

#include <unordered_set>

using namespace gbs;
using namespace gbs::brep;
using namespace brep_tests;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

TEST(tests_brep_explore, explore_by_type)
{
    auto b = make_box();
    const auto &m = b.m;

    ASSERT_EQ(explore<VertexId>(m, b.solid).size(), 8);
    ASSERT_EQ(explore<EdgeId>(m, b.solid).size(), 12);
    ASSERT_EQ(explore<WireId>(m, b.solid).size(), 6);
    ASSERT_EQ(explore<FaceId>(m, b.solid).size(), 6);
    ASSERT_EQ(explore<ShellId>(m, b.solid).size(), 1);
    ASSERT_EQ(explore<SolidId>(m, b.solid).size(), 1); // the shape itself
    ASSERT_EQ(explore<CompoundId>(m, b.solid).size(), 0);

    // from a face
    ASSERT_EQ(explore<VertexId>(m, b.f[0]).size(), 4);
    ASSERT_EQ(explore<EdgeId>(m, b.f[0]).size(), 4);
    ASSERT_EQ(explore<WireId>(m, b.f[0]).size(), 1);
    ASSERT_EQ(explore<ShellId>(m, b.f[0]).size(), 0); // no upward exploration

    // from an edge, a vertex
    ASSERT_EQ(explore<VertexId>(m, b.e[0]).size(), 2);
    ASSERT_EQ(explore<VertexId>(m, b.v[0]).size(), 1);
    ASSERT_EQ(explore<EdgeId>(m, b.v[0]).size(), 0);

    // discovery order of the vertices of a face follows its wire
    auto vtx = explore<VertexId>(m, b.f[0]);
    const auto &w = m.wire(m.face(b.f[0]).wires.front());
    ASSERT_TRUE(vtx[0] == m.edge(w.coedges[0].edge).v1);

    // invalid or dead handles yield nothing
    ASSERT_EQ(explore<VertexId>(m, FaceId{}).size(), 0);
    auto b2 = make_box();
    b2.m.erase(b2.f[0]);
    ASSERT_EQ(explore<FaceId>(b2.m, b2.shell).size(), 5);
    ASSERT_EQ(explore<VertexId>(b2.m, b2.shell).size(), 8); // still reached through the other faces
}

TEST(tests_brep_explore, compound_and_duplicates)
{
    auto b = make_box();
    auto &m = b.m;

    // a compound holding the solid, one of its faces and one of its edges: no duplicates
    auto c = m.add(Compound{{ShapeId{b.solid}, ShapeId{b.f[2]}, ShapeId{b.e[5]}, ShapeId{b.v[7]}}});
    ASSERT_EQ(explore<FaceId>(m, c).size(), 6);
    ASSERT_EQ(explore<EdgeId>(m, c).size(), 12);
    ASSERT_EQ(explore<VertexId>(m, c).size(), 8);
    ASSERT_EQ(explore<CompoundId>(m, c).size(), 1);

    // nested compound
    auto c2 = m.add(Compound{{ShapeId{c}, ShapeId{c}}});
    ASSERT_EQ(explore<CompoundId>(m, c2).size(), 2);
    ASSERT_EQ(explore<SolidId>(m, c2).size(), 1);
}

TEST(tests_brep_explore, topology_index)
{
    auto b = make_box();
    const auto &m = b.m;
    TopologyIndex<T> idx{m, b.solid};

    ASSERT_EQ(idx.edges().size(), 12);
    ASSERT_EQ(idx.faces().size(), 6);

    for (auto eid : b.e)
    {
        ASSERT_EQ(idx.faces_of(eid).size(), 2);
        ASSERT_EQ(idx.coedges_of(eid).size(), 2);
        ASSERT_TRUE(idx.faces_of(eid)[0] != idx.faces_of(eid)[1]);
        // the co-edge refs point to co-edges of that edge
        for (const auto &ref : idx.coedges_of(eid))
        {
            ASSERT_TRUE(m.wire(ref.wire).coedges[ref.index].edge == eid);
            ASSERT_TRUE(idx.face_of(ref.wire).valid());
        }
    }
    for (auto vid : b.v)
        ASSERT_EQ(idx.edges_of(vid).size(), 3);
    for (auto fid : b.f)
    {
        ASSERT_EQ(idx.shells_of(fid).size(), 1);
        ASSERT_TRUE(idx.shells_of(fid)[0] == b.shell);
        ASSERT_TRUE(idx.face_of(m.face(fid).wires.front()) == fid);
    }

    // out-of-range or unrelated handles give empty spans
    ASSERT_EQ(idx.faces_of(EdgeId{999}).size(), 0);
    ASSERT_EQ(idx.edges_of(VertexId{}).size(), 0);
    ASSERT_FALSE(idx.face_of(WireId{999}).valid());

    // index restricted to one face
    TopologyIndex<T> idx_face{m, b.f[0]};
    ASSERT_EQ(idx_face.edges().size(), 4);
    for (auto eid : idx_face.edges())
        ASSERT_EQ(idx_face.faces_of(eid).size(), 1);
    ASSERT_EQ(idx_face.shells_of(b.f[0]).size(), 0);
}

TEST(tests_brep_explore, wire_closure)
{
    auto b = make_box();
    auto &m = b.m;

    for (auto fid : b.f)
    {
        auto wid = m.face(fid).wires.front();
        ASSERT_TRUE(is_closed(m, wid));
        ASSERT_TRUE(is_chained(m, wid));
    }

    // open chain: drop the last co-edge
    auto w_open = m.wire(m.face(b.f[0]).wires.front());
    w_open.coedges.pop_back();
    auto wid_open = m.add(w_open);
    ASSERT_FALSE(is_closed(m, wid_open));
    ASSERT_TRUE(is_chained(m, wid_open));

    // broken chain: swap two co-edges
    auto w_broken = m.wire(m.face(b.f[0]).wires.front());
    std::swap(w_broken.coedges[1], w_broken.coedges[2]);
    auto wid_broken = m.add(w_broken);
    ASSERT_FALSE(is_closed(m, wid_broken));
    ASSERT_FALSE(is_chained(m, wid_broken));

    // flipping one co-edge breaks the chain
    auto w_flip = m.wire(m.face(b.f[0]).wires.front());
    w_flip.coedges[0].orient = reverse(w_flip.coedges[0].orient);
    ASSERT_FALSE(is_closed(m, m.add(w_flip)));

    // empty wire
    ASSERT_FALSE(is_closed(m, m.add(Wire<T>{})));

    // coedge_point follows the co-edge sense: end of i == start of i+1
    auto wid = m.face(b.f[3]).wires.front();
    const auto &w = m.wire(wid);
    for (size_t i = 0; i < 4; ++i)
    {
        const auto &e = m.edge(w.coedges[i].edge);
        const auto &e_next = m.edge(w.coedges[(i + 1) % 4].edge);
        auto p_end = coedge_point(m, wid, i, e.u2);
        auto p_start = coedge_point(m, wid, (i + 1) % 4, e_next.u1);
        ASSERT_NEAR(norm(p_end - p_start), 0., 1e-12);
        ASSERT_NEAR(norm(coedge_point(m, wid, i, e.u1) - m.vertex(coedge_start(m, w.coedges[i])).pnt), 0., 1e-12);
    }
}

TEST(tests_brep_explore, shell_queries)
{
    auto b = make_box();
    auto &m = b.m;

    ASSERT_TRUE(is_closed(m, b.shell));
    ASSERT_TRUE(is_manifold(m, b.shell));
    ASSERT_TRUE(is_orientable(m, b.shell));
    ASSERT_EQ(free_edges(m, b.shell).size(), 0);
    ASSERT_EQ(non_manifold_edges(m, b.shell).size(), 0);
    ASSERT_EQ(shell_edge_uses(m, b.shell).size(), 12);

    // open box: 5 faces -> 4 free edges, still manifold and consistently oriented
    Shell open_sh;
    for (size_t i = 1; i < 6; ++i)
        open_sh.faces.push_back(FaceUse{b.f[i], Orientation::Forward});
    auto open_id = m.add(open_sh);
    ASSERT_FALSE(is_closed(m, open_id));
    ASSERT_TRUE(is_manifold(m, open_id));
    ASSERT_TRUE(is_orientable(m, open_id));
    auto fe = free_edges(m, open_id);
    ASSERT_EQ(fe.size(), 4);
    ASSERT_TRUE(std::is_sorted(fe.begin(), fe.end()));
    // the free edges are exactly those of the removed face
    auto removed = explore<EdgeId>(m, b.f[0]);
    std::sort(removed.begin(), removed.end());
    ASSERT_TRUE(fe == removed);

    // one flipped face use: closed and manifold, but not consistently oriented
    Shell flipped = m.shell(b.shell);
    flipped.faces[2].orient = Orientation::Reversed;
    auto flipped_id = m.add(flipped);
    ASSERT_TRUE(is_closed(m, flipped_id));
    ASSERT_TRUE(is_manifold(m, flipped_id));
    ASSERT_FALSE(is_orientable(m, flipped_id));

    // the whole shell reversed stays consistently oriented
    Shell all_rev = m.shell(b.shell);
    for (auto &fu : all_rev.faces)
        fu.orient = Orientation::Reversed;
    ASSERT_TRUE(is_orientable(m, m.add(all_rev)));

    // a seventh face reusing the wire of face 0: its 4 edges are used three times
    auto f7 = m.add(Face<T>{m.face(b.f[0]).surface, {m.face(b.f[0]).wires.front()}, 1e-7, true});
    Shell nm = m.shell(b.shell);
    nm.faces.push_back(FaceUse{f7, Orientation::Forward});
    auto nm_id = m.add(nm);
    ASSERT_FALSE(is_manifold(m, nm_id));
    ASSERT_FALSE(is_closed(m, nm_id));
    ASSERT_FALSE(is_orientable(m, nm_id));
    ASSERT_EQ(non_manifold_edges(m, nm_id).size(), 4);
    ASSERT_EQ(free_edges(m, nm_id).size(), 0);

    // degenerate edges are ignored by the usage count
    auto vd = m.add(Vertex<T>{{5., 5., 5.}, 1e-7});
    auto ed = m.add(Edge<T>{nullptr, 0., 1., vd, vd, 1e-7, true});
    Wire<T> wd;
    wd.coedges.push_back(CoEdge<T>{ed, Orientation::Forward, nullptr});
    wd.closed = true;
    auto fd = m.add(Face<T>{m.face(b.f[0]).surface, {m.add(wd)}, 1e-7, false});
    Shell with_deg = m.shell(b.shell);
    with_deg.faces.push_back(FaceUse{fd, Orientation::Forward});
    auto with_deg_id = m.add(with_deg);
    ASSERT_TRUE(is_closed(m, with_deg_id));
    ASSERT_EQ(shell_edge_uses(m, with_deg_id).size(), 12);

    // empty shell is not closed
    ASSERT_FALSE(is_closed(m, m.add(Shell{})));
}

TEST(tests_brep_explore, bounding_box)
{
    auto b = make_box();
    const auto &m = b.m;
    const T tol = 1e-7;

    auto box = bounding_box(m, b.solid);
    ASSERT_FALSE(box.empty());
    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NEAR(box.min[i], -tol, 1e-12);
        ASSERT_NEAR(box.max[i], 1. + tol, 1e-12);
    }
    ASSERT_NEAR(box.diagonal(), std::sqrt(3.) * (1. + 2. * tol), 1e-9);
    ASSERT_TRUE(box.contains({0.5, 0.5, 0.5}));
    ASSERT_FALSE(box.contains({1.5, 0.5, 0.5}));

    // one edge along x from the origin
    Orientation o;
    auto ex = find_edge(m, b.e, b.v[0], b.v[1], o);
    auto bx = bounding_box(m, ex);
    ASSERT_NEAR(bx.min[0], -tol, 1e-12);
    ASSERT_NEAR(bx.max[0], 1. + tol, 1e-12);
    ASSERT_NEAR(bx.max[1], tol, 1e-12);
    ASSERT_NEAR(bx.max[2], tol, 1e-12);

    // a vertex
    auto bv = bounding_box(m, b.v[7]);
    ASSERT_NEAR(bv.min[0], 1. - tol, 1e-12);
    ASSERT_NEAR(bv.max[2], 1. + tol, 1e-12);

    // intersection / inflate / add
    auto b1 = bounding_box(m, b.f[0]);
    auto b2 = bounding_box(m, b.f[1]); // opposite face along the same axis: disjoint
    ASSERT_FALSE(b1.intersects(b2));
    b1.inflate(1.);
    ASSERT_TRUE(b1.intersects(b2));
    BoundingBox<T> empty;
    ASSERT_TRUE(empty.empty());
    ASSERT_FALSE(empty.intersects(b1));
    empty.add(b2);
    ASSERT_FALSE(empty.empty());
    ASSERT_TRUE(empty.intersects(b2));

    // a bulging natural-bounds face exceeds its edges: quadratic patch with a raised middle pole
    Model<T> m2;
    points_vector<T, 3> poles{
        {0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.},
        {0., 1., 0.}, {1., 1., 4.}, {2., 1., 0.},
        {0., 2., 0.}, {1., 2., 0.}, {2., 2., 0.}};
    std::vector<T> k{0., 0., 0., 1., 1., 1.};
    auto srf = std::make_shared<BSSurface<T, 3>>(poles, k, k, 2, 2);
    auto fid = m2.add(Face<T>{srf, {}, 1e-7, true});
    auto bb = bounding_box(m2, fid, 21);
    ASSERT_GT(bb.max[2], 0.9); // the middle of the patch reaches z = 1
    auto fid_trimmed = m2.add(Face<T>{srf, {}, 1e-7, false});
    ASSERT_TRUE(bounding_box(m2, fid_trimmed).empty()); // no edges, no natural bounds: nothing to sample
}
