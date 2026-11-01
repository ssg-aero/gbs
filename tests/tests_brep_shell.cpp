// Blades as shells: contours cut in pieces skinned one by one (loft_pieces), extended
// together (extrapolated_v_by_lengths), sewn into one shell (sewn_shell), meshed with their
// edges (Mesh::edges).
#include <doctest_gtest.hpp>

#include <gbs-brep/brep>
#include <gbs/bssbuild.h>
#include <gbs/bsstools.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <type_traits>
#include <vector>

using gbs::BSSurface;
using gbs::brep::BuildErrc;
using gbs::brep::Model;
using gbs::brep::ShapeId;
using gbs::Surface;
using gbs::operator-;
using gbs::brep::sewn_shell;

namespace
{
    // Contours at heights z, cut in two pieces sharing their ends: a cubic side from (0, 0) to
    // (1 + 0.3 z, 0) and a straight edge back, closing it (a corner at each join).
    auto contour_pieces(const std::vector<double> &zs) -> std::vector<std::vector<gbs::BSCurve3d_d>>
    {
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
    }

    const std::vector<double> zs{0., 0.6, 1.5, 2.}, vs{0., 0.3, 0.75, 1.};

    auto shared(const std::vector<BSSurface<double, 3>> &surfaces)
    {
        std::vector<std::shared_ptr<Surface<double, 3>>> srfs;
        for (const auto &s : surfaces)
            srfs.push_back(std::make_shared<BSSurface<double, 3>>(s));
        return srfs;
    }
}

// loft_pieces: each piece skinned on its own at the same parameters, through its curves,
// the pieces sharing their boundaries exactly, rational or not; bad parameters and pieces
// refused.
TEST(tests_brep_shell, loft_pieces)
{
    const auto pieces = contour_pieces(zs);
    const auto s = gbs::loft_pieces(pieces, vs);
    ASSERT_EQ(s.size(), 2u);
    for (size_t k = 0; k < 2; ++k)
        for (size_t i = 0; i < vs.size(); ++i)
        {
            const auto [c1, c2] = pieces[k][i].bounds();
            for (double t : {0., 0.3, 1.})
                ASSERT_LT(gbs::norm(s[k].value(c1 + (c2 - c1) * t, vs[i]) - pieces[k][i].value(c1 + (c2 - c1) * t)), 1e-12);
        }
    const auto [a1, a2, av1, av2] = s[0].bounds();
    const auto [b1, b2, bv1, bv2] = s[1].bounds();
    for (int j = 0; j <= 40; ++j)
    {
        const double t = av1 + (av2 - av1) * j / 40.;
        ASSERT_LT(gbs::norm(s[0].value(a2, t) - s[1].value(b1, t)), 1e-14);
        ASSERT_LT(gbs::norm(s[0].value(a1, t) - s[1].value(b2, t)), 1e-14);
    }

    // rational pieces (a weight off 1 inside the side): rational surfaces, still sharing
    std::vector<std::vector<gbs::BSCurveRational3d_d>> rational(2);
    for (size_t i = 0; i < zs.size(); ++i)
    {
        const auto &side = pieces[0][i];
        rational[0].emplace_back(side.poles(), side.knotsFlats(), std::vector<double>{1., 0.8, 1.3, 0.9, 1.}, 3);
        rational[1].emplace_back(pieces[1][i]);
    }
    const auto r = gbs::loft_pieces(rational, vs);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(r.front())>, gbs::BSSurfaceRational<double, 3>>);
    for (size_t i = 0; i < vs.size(); ++i)
        for (double t : {0., 0.4, 1.})
            ASSERT_LT(gbs::norm(r[0].value(t, vs[i]) - rational[0][i].value(t)), 1e-12);
    for (int j = 0; j <= 40; ++j)
    {
        const double t = j / 40.;
        ASSERT_LT(gbs::norm(r[0].value(1., t) - r[1].value(0., t)), 1e-14);
        ASSERT_LT(gbs::norm(r[0].value(0., t) - r[1].value(1., t)), 1e-14);
    }

    auto short_piece = pieces;
    short_piece[1].pop_back();
    EXPECT_THROW(gbs::loft_pieces(short_piece, vs), std::invalid_argument);
    EXPECT_THROW(gbs::loft_pieces(pieces, {0., 0.75, 0.3, 1.}), std::invalid_argument);
    EXPECT_THROW(gbs::loft_pieces(pieces, {0., std::numeric_limits<double>::quiet_NaN(), 0.75, 1.}), std::invalid_argument);
    EXPECT_THROW(gbs::loft_pieces(std::vector<std::vector<gbs::BSCurve3d_d>>{{pieces[0][0]}}, {0.}), std::invalid_argument);
    // pieces whose ends part: the edge moved off the end of the side
    auto apart = pieces;
    auto moved = apart[1][2].poles();
    for (auto &p : moved)
        p[1] += 1e-3;
    apart[1][2] = gbs::BSCurve3d_d(moved, apart[1][2].knotsFlats(), 1);
    EXPECT_THROW(gbs::loft_pieces(apart, vs), std::invalid_argument);
}

// A skinned contour (loft_pieces), extended past its ends (extrapolated_v_by_lengths), as one
// shell (sewn_shell): every join one edge of two faces, the open ends free; its mesh draws
// every edge once, by the points its faces share. Faces that do not sew into one shell are
// refused, the model left as it was.
TEST(tests_brep_shell, sewn_shell_of_pieces)
{
    const auto surfaces = gbs::extrapolated_v_by_lengths(gbs::loft_pieces(contour_pieces(zs), vs), {0.1, 0.1});
    Model<double> m;
    auto sewn = sewn_shell(m, shared(surfaces));
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

    // Surfaces apart: no single shell, nothing left in the model.
    auto apart = gbs::loft_pieces(contour_pieces(zs), vs);
    auto moved = apart[1].poles();
    for (auto &p : moved)
        p[2] += 5.;
    apart[1] = BSSurface<double, 3>(moved, apart[1].knotsFlatsU(), apart[1].knotsFlatsV(), apart[1].degreeU(), apart[1].degreeV());
    Model<double> m2;
    auto fail = sewn_shell(m2, shared(apart));
    ASSERT_FALSE(fail.has_value());
    ASSERT_EQ(fail.error().code, BuildErrc::Disconnected);
    ASSERT_TRUE(m2.empty());
}

// Two halves of a Moebius strip sew into one shell that cannot be oriented: refused, nothing
// left in the model.
TEST(tests_brep_shell, sewn_shell_not_orientable)
{
    const double w = 0.3, pi = std::numbers::pi;
    auto half = [&](double s1) {
        std::vector<gbs::BSCurve3d_d> lines;
        std::vector<double> v;
        for (int j = 0; j <= 4; ++j)
        {
            const double s = s1 + pi * j / 4.;
            auto at = [&](double t) { return std::array<double, 3>{(1. + t * std::cos(s / 2.)) * std::cos(s), (1. + t * std::cos(s / 2.)) * std::sin(s), t * std::sin(s / 2.)}; };
            lines.emplace_back(gbs::points_vector_3d_d{at(-w), at(w)}, std::vector<double>{0., 0., 1., 1.}, 1);
            v.push_back(s);
        }
        return gbs::loft(lines, v, 3);
    };
    Model<double> m;
    auto r = sewn_shell(m, shared({half(0.), half(pi)}));
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::ShellNotOrientable);
    ASSERT_TRUE(m.empty());
}

// A sewing tolerance above the width of a thin piece (a blunt trailing edge 2e-5 thick, sewn
// at 1e-4): its face keeps its four sides, the shell its three joins.
TEST(tests_brep_shell, sewn_shell_thin_piece)
{
    using C = gbs::BSCurve3d_d;
    const double thick = 2e-5;
    const std::vector<double> k{0., 0., 0., 1., 1., 1.}, ke{0., 0., 1., 1.};
    std::vector<std::vector<C>> pieces(3);
    for (double z : zs)
    {
        const double c = 1. + 0.2 * z;
        pieces[0].push_back(C({{0., 0., z}, {0.5 * c, 0.3, z}, {c, 0., z}}, k, 2));               // side 1, over
        pieces[1].push_back(C({{c, 0., z}, {c, -thick, z}}, ke, 1));                              // trailing edge
        pieces[2].push_back(C({{c, -thick, z}, {0.5 * c, -0.3, z}, {0., 0., z}}, k, 2));          // side 2, under
    }
    Model<double> m;
    gbs::brep::SewOptions<double> opts;
    opts.tol = 1e-4;
    auto sewn = sewn_shell(m, shared(gbs::loft_pieces(pieces, vs)), opts);
    ASSERT_TRUE(sewn.has_value());
    ASSERT_EQ(sewn->report.merged.size(), 3u);
    for (auto f : sewn->faces)
    {
        ASSERT_LE(m.face(f).tol, gbs::brep_default_tolerance<double>);
        for (auto wid : m.face(f).wires)
            for (const auto &ce : m.wire(wid).coedges)
                ASSERT_FALSE(m.edge(ce.edge).degenerate);
    }
}
