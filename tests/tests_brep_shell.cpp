// Blades as shells: contours cut in pieces skinned one by one (loft_pieces), extended
// together (extrapolated_v_by_lengths), sewn into one shell (sewn_shell), meshed with their
// edges (Mesh::edges).
#include <doctest_gtest.hpp>

#include <gbs-brep/brep>
#include <gbs/bssbuild.h>
#include <gbs/bsstools.h>

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

using gbs::BSSurface;
using gbs::brep::BuildErrc;
using gbs::brep::Model;
using gbs::brep::ShapeId;
using gbs::Surface;
using gbs::operator-;
using gbs::brep::sewn_shell;

// A skinned contour (loft_pieces), extended past its ends (extrapolated_v_by_lengths), as one
// shell (sewn_shell): every join one edge of two faces, the open ends free; its mesh draws
// every edge once, by the points its faces share.
TEST(tests_brep_shell, sewn_shell_of_pieces)
{
    // Contours at heights z, cut in two pieces sharing their ends: a cubic side from (0, 0) to
    // (1 + 0.3 z, 0) and a straight edge back, closing it (a corner at each join).
    auto contour_pieces = [](const std::vector<double> &zs) {
        using C = gbs::BSCurve3d_d;
        const std::vector<double> k{0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
        const std::vector<double> ke{0., 0., 1., 1.};
        std::vector<std::vector<C>> pieces(2);
        for (double z : zs)
        {
            pieces[0].push_back(C({{0., 0., z}, {0.3, 0.2 + 0.1 * z, z}, {0.6, 0.25, z}, {0.9 + 0.2 * z, 0.1, z}, {1. + 0.3 * z, 0., z}}, k, 3));
            pieces[1].push_back(C({{1. + 0.3 * z, 0., z}, {0., 0., z}}, ke, 1));
        }
        return pieces;
    };
    const std::vector<double> zs{0., 0.6, 1.5, 2.}, v{0., 0.3, 0.75, 1.};
    auto surfaces = gbs::extrapolated_v_by_lengths(gbs::loft_pieces(contour_pieces(zs), v), {0.1, 0.1});
    Model<double> m;
    std::vector<std::shared_ptr<Surface<double, 3>>> srfs;
    for (auto &s : surfaces)
        srfs.push_back(std::make_shared<BSSurface<double, 3>>(std::move(s)));
    auto sewn = sewn_shell(m, srfs);
    ASSERT_TRUE(sewn.has_value());
    ASSERT_EQ(sewn->faces.size(), 2u);
    ASSERT_EQ(sewn->report.merged.size(), 2u);     // the two joins
    ASSERT_EQ(sewn->report.free_edges.size(), 4u); // the ends, two per face
    ASSERT_EQ(sewn->report.rejected.size(), 2u);   // the ends of the side and of the edge, apart
    std::map<std::uint32_t, int> uses;
    for (auto f : sewn->faces)
        for (auto w : m.face(f).wires)
            for (const auto &ce : m.wire(w).coedges)
                ++uses[ce.edge.index];
    ASSERT_EQ(std::ranges::count_if(uses, [](const auto &u) { return u.second == 2; }), 2);

    const auto mesh = gbs::brep::mesh(m, ShapeId{sewn->shell});
    ASSERT_TRUE(mesh.has_value());
    ASSERT_EQ(mesh->edges.size(), uses.size());
    for (const auto &line : mesh->edges)
    {
        ASSERT_GE(line.vertices.size(), 2u);
        const auto &e = m.edge(line.edge);
        ASSERT_LT(gbs::norm(mesh->vertices[line.vertices.front()] - m.vertex(e.v1).pnt), 1e-9);
        ASSERT_LT(gbs::norm(mesh->vertices[line.vertices.back()] - m.vertex(e.v2).pnt), 1e-9);
        // each point a vertex of the triangles: the faces share them
        for (auto vtx : line.vertices)
            ASSERT_TRUE(std::ranges::any_of(mesh->triangles, [&](const auto &t) { return std::ranges::contains(t, vtx); }));
    }

    // Surfaces apart: no single shell.
    auto apart = gbs::loft_pieces(contour_pieces(zs), v);
    auto moved = apart[1].poles();
    for (auto &p : moved)
        p[2] += 5.;
    apart[1] = BSSurface<double, 3>(moved, apart[1].knotsFlatsU(), apart[1].knotsFlatsV(), apart[1].degreeU(), apart[1].degreeV());
    Model<double> m2;
    auto fail = sewn_shell(m2, {std::make_shared<BSSurface<double, 3>>(apart[0]), std::make_shared<BSSurface<double, 3>>(apart[1])});
    ASSERT_FALSE(fail.has_value());
    ASSERT_EQ(fail.error().code, BuildErrc::Disconnected);
}
