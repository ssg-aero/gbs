// Fillets between two faces (stage 8, docs/sources/design/brep_palier8_face_blends.md; PR 1,
// docs/sources/design/brep8_pr01_face_blends.md; PR 2, brep8_pr02_partial_strips_laws.md): the
// strip between them swallowed, the fillet that of the virtual edge where the faces meet once
// continued; a partial strip (faces that already meet along part of their length) fused with
// their common edge; laws of radius with points by position; the sections of stage 7. Each case
// against the fillet of the edge that existed before the strip (a chamfer, a tapered face): the
// same solid, the same volume; refusals.
#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/face_blend.h>
#include <gbs-brep/sew.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>

#include <array>
#include <chrono>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;
    const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};

    SolidId solid_of(Model<T> &m, const BuildResult<SweepResult<T>> &r)
    {
        return explore<SolidId>(m, unwrap(r).shape).front();
    }

    FaceId polygon(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    }

    SolidId box(Model<T> &m, T a, T b, T c, const point<T, 3> &o = {0., 0., 0.})
    {
        return solid_of(m, make_prism(m, polygon(m, {o, o + point<T, 3>{a, 0., 0.}, o + point<T, 3>{a, b, 0.}, o + point<T, 3>{0., b, 0.}}), point<T, 3>{0., 0., c}));
    }

    /// The prism along x of length L on the profile (y, z) given in the plane x = x0 (counter-clockwise about x).
    SolidId profile(Model<T> &m, const std::vector<std::array<T, 2>> &yz, T L, T x0 = 0.)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < yz.size(); ++i)
        {
            const auto &a = yz[i], &b = yz[(i + 1) % yz.size()];
            e.push_back(unwrap(make_edge(m, point<T, 3>{x0, a[0], a[1]}, point<T, 3>{x0, b[0], b[1]})));
        }
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{x0, 0., 0.}, X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{L, 0., 0.}));
    }

    /// The live entities of the model (vertices, edges, faces): unchanged by a refusal.
    std::array<std::size_t, 3> census(const Model<T> &m)
    {
        return {m.template ids<VertexId>().size(), m.template ids<EdgeId>().size(), m.template ids<FaceId>().size()};
    }

    /// The refusal of r: its code and a part of its message.
    void expect_refused(const BuildResult<SolidId> &r, BuildErrc code, const std::string &why)
    {
        ASSERT_FALSE(r.has_value());
        MESSAGE(r.error().message);
        ASSERT_EQ(r.error().code, code);
        ASSERT_NE(r.error().message.find(why), std::string::npos);
    }

    SolidId cylinder(Model<T> &m, T r, T h)
    {
        const ax2<T, 3> ax{{{0., 0., 0.}, Z, X}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    EdgeId edge_near(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        EdgeId best{};
        T gap = 1e300;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            for (int i = 0; i <= 16; ++i)
            {
                const T d = norm(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / 16.) - p);
                if (d < gap)
                    gap = d, best = e;
            }
        }
        return best;
    }

    /// The face of s nearest p: its surface through p within 1e-7, p inside it.
    FaceId face_near(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        for (auto f : explore<FaceId>(m, ShapeId{s}))
        {
            const auto &srf = *m.face(f).surface;
            const auto y = extrema_surf_pnt(srf, p, 1e-12, std::array<std::size_t, 2>{20, 20});
            if (y[2] < 1e-7 && classify(m, f, point<T, 2>{y[0], y[1]}, 1e-7) == Location::In)
                return f;
        }
        return FaceId{};
    }

    /// The faces of s on a cylinder.
    std::vector<FaceId> cylinders(const Model<T> &m, SolidId s)
    {
        std::vector<FaceId> r;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (dynamic_cast<const CylindricalSurface<T> *>(m.face(f).surface.get()))
                r.push_back(f);
        return r;
    }

    /// The slab on a stadium: straight sides of length L, half circles of radius R, height h.
    SolidId oblong(Model<T> &m, T L, T R, T h)
    {
        auto arc = [](T r, const point<T, 3> &c, T a1, T a2) { return std::make_shared<EllipseCurve<T>>(r, r, ax2<T, 3>{{c, Z, X}}, a1, a2); };
        std::vector<EdgeId> e{unwrap(make_edge(m, point<T, 3>{0., -R, 0.}, point<T, 3>{L, -R, 0.})), unwrap(make_edge(m, arc(R, {L, 0., 0.}, -pi / 2., pi / 2.))),
                              unwrap(make_edge(m, point<T, 3>{L, R, 0.}, point<T, 3>{0., R, 0.})), unwrap(make_edge(m, arc(R, {0., 0., 0.}, pi / 2., 3. * pi / 2.)))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                        unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    SolidId fused(Model<T> &m, SolidId a, SolidId b)
    {
        const auto r = unwrap(fuse(m, ShapeId{a}, ShapeId{b}));
        return r.solids.front();
    }

    /// The edges of s on the plane z = 0 within the disc of radius rho round the point c.
    std::vector<EdgeId> edges_on_floor(const Model<T> &m, SolidId s, const point<T, 3> &c, T rho)
    {
        std::vector<EdgeId> r;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (ed.degenerate)
                continue;
            bool in = true;
            for (int i = 0; i <= 8 && in; ++i)
            {
                const auto p = ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / 8.);
                in = std::abs(p[2]) < 1e-9 && norm(p - c) < rho;
            }
            if (in)
                r.push_back(e);
        }
        return r;
    }

    T volume(const Model<T> &m, SolidId s)
    {
        return signed_volume(m, m.solid(s).outer, 512);
    }

    std::size_t n_faces(const Model<T> &m, SolidId s)
    {
        return explore<FaceId>(m, ShapeId{s}).size();
    }

    void expect_valid(const Model<T> &m, SolidId s)
    {
        const auto report = check(m, ShapeId{s});
        for (const auto &e : report.entries)
            INFO(to_string(e.issue) << " " << e.detail);
        ASSERT_TRUE(report.ok());
    }

    /// The planar face bounded by the polygon p (counter-clockwise about its normal).
    FaceId planar(Model<T> &m, const std::vector<point<T, 3>> &p)
    {
        const auto n = cross(p[1] - p[0], p[2] - p[0]);
        const auto x = p[1] - p[0];
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            e.push_back(unwrap(make_edge(m, p[i], p[(i + 1) % p.size()])));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n / norm(n), x / norm(x)}}, std::array<T, 4>{-10., 10., -10., 10.}),
                                unwrap(make_wire(m, e))));
    }

    /// The tetrahedron a, b, c, d (b - a, c - a, d - a right-handed).
    SolidId tetrahedron(Model<T> &m, const point<T, 3> &a, const point<T, 3> &b, const point<T, 3> &c, const point<T, 3> &d)
    {
        const std::vector<FaceId> f{planar(m, {a, c, b}), planar(m, {a, b, d}), planar(m, {a, d, c}), planar(m, {b, c, d})};
        SewOptions<T> so;
        so.tol = 1e-9;
        const auto sewn = unwrap(sew(m, f, so));
        return unwrap(make_solid(m, sewn.shells.front()));
    }

    // The U channel of PR 2 (§ 6): its section in the plane x = 0, along x over 2; the left wall
    // y = -0.6 and the floor z = 0 meet along the corner edge (y, z) = (-0.6, 0).
    constexpr T Yw = 0.6;
    SolidId u_channel(Model<T> &m)
    {
        return profile(m, {{-1., -0.6}, {1., -0.6}, {1., 1.}, {Yw, 1.}, {Yw, 0.}, {-Yw, 0.}, {-Yw, 1.}, {-1., 1.}}, 2.);
    }

    /// The U with a face tapering to a point on the corner: a triangle of width c at x = 0, its tip at x = 1 (a tetrahedron fused over the corner, overlapping the walls by dl).
    SolidId tapered_u(Model<T> &m, T c = 0.12, T dl = 0.05)
    {
        const auto t = tetrahedron(m, {0., -Yw - dl, -dl}, {0., -Yw - dl, c + dl}, {0., -Yw + c + dl, -dl}, {(c + 2. * dl) / c, -Yw - dl, -dl});
        return fused(m, u_channel(m), t);
    }

    /// The U with the corner chamfered over x in [x0, x0 + L] (a triangular prism fused over it, its hypotenuse 0.17 / sqrt 2 from the corner).
    SolidId chamfered_u(Model<T> &m, T x0, T L)
    {
        return fused(m, u_channel(m), profile(m, {{-0.65, -0.05}, {-0.38, -0.05}, {-0.65, 0.22}}, L, x0));
    }

    FaceId u_wall(const Model<T> &m, SolidId s) { return face_near(m, s, {1.5, -Yw, 0.6}); }
    FaceId u_floor(const Model<T> &m, SolidId s) { return face_near(m, s, {1.5, 0., 0.}); }

    /// The fillet of the corner of the sharp U at r (the reference of the partial strips).
    SolidId u_reference(Model<T> &m, T r, const FilletOptions<T> &o = {})
    {
        const auto u = u_channel(m);
        return unwrap(fillet(m, u, std::vector<EdgeId>{edge_near(m, u, {1., -Yw, 0.})}, r, o));
    }

    T seconds_since(std::chrono::steady_clock::time_point t0)
    {
        return std::chrono::duration<T>(std::chrono::steady_clock::now() - t0).count();
    }
}

TEST(tests_brep_face_blends, chamfered_box_edge)
{
    // the box 2 x 1 x 1, its edge y = 1, z = 1 chamfered at 0.1; the fillet between the top and
    // the side over the chamfer is the fillet of the edge before the chamfer
    for (T r : {0.2, 0.3, 0.5})
    {
        INFO("r " << r);
        Model<T> m;
        const auto b = box(m, 2., 1., 1.);
        const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{edge_near(m, b, {1., 1., 1.})}, 0.1));
        const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
        const auto size = n_faces(m, ch);
        const auto f = unwrap(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, r));
        expect_valid(m, f);
        ASSERT_EQ(n_faces(m, f), size);
        const auto ref = unwrap(fillet(m, b, std::vector<EdgeId>{edge_near(m, b, {1., 1., 1.})}, r));
        ASSERT_EQ(n_faces(m, ref), n_faces(m, f));
        ASSERT_NEAR(volume(m, f), 2. * (1. - (1. - pi / 4.) * r * r), 1e-9);
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
}

TEST(tests_brep_face_blends, chamfered_box_edge_refusals)
{
    // under r = d / (2 - sqrt 2) the ball cuts the chamfer: refused, the model unchanged; parallel
    // faces (the two sides of a slab): no virtual edge; faces with an edge in common: the fillet
    // of their edges
    Model<T> m;
    const auto b = box(m, 2., 1., 1.);
    const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{edge_near(m, b, {1., 1., 1.})}, 0.1));
    const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
    const auto before = n_faces(m, ch);
    const auto spill = fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, 0.13);
    ASSERT_FALSE(spill.has_value());
    MESSAGE(spill.error().message);
    ASSERT_NE(spill.error().message.find("spills over"), std::string::npos);
    ASSERT_EQ(n_faces(m, ch), before);
    const auto bottom = face_near(m, ch, {1., 0.5, 0.}), back = face_near(m, ch, {1., 0., 0.5});
    const auto n0 = census(m);
    expect_refused(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{bottom}, 0.3), BuildErrc::NonGeneric, "do not meet once continued");
    ASSERT_EQ(census(m), n0);
    // faces with an edge in common: the fillet of that edge
    const auto adj = unwrap(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{back}, 0.2));
    expect_valid(m, adj);
    ASSERT_NEAR(volume(m, adj), volume(m, ch) - 2. * (1. - pi / 4.) * 0.04, 1e-9);
}

TEST(tests_brep_face_blends, sides_refused)
{
    // the sides checked before anything: an empty side, a face of another solid, a face on both
    // sides (it is next to itself along its edges: not the fillet of all of them), faces of a side
    // not joined tangentially (InvalidId); faces that meet along an edge in places only, a
    // chamfer between them elsewhere: not the fillet of that edge alone (NonGeneric); the model
    // unchanged each time (the box with its edge chamfered from x = 1 on, which PR 1 refused here,
    // is a partial strip since PR 2: stopped_box_chamfer)
    Model<T> m;
    const auto b = box(m, 2., 1., 1.);
    const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{edge_near(m, b, {1., 1., 1.})}, 0.1));
    const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5}), end = face_near(m, ch, {0., 0.5, 0.5});
    const auto other = face_near(m, b, {1., 0.5, 1.});
    const auto n0 = census(m);
    using F = std::vector<FaceId>;
    expect_refused(fillet_faces(m, ch, F{}, F{side}, 0.3), BuildErrc::InvalidId, "no face on a side");
    expect_refused(fillet_faces(m, ch, F{top}, F{side, other}, 0.3), BuildErrc::InvalidId, "not a face of the solid");
    expect_refused(fillet_faces(m, ch, F{top}, F{top}, 0.03), BuildErrc::InvalidId, "on both sides");
    expect_refused(fillet_faces(m, ch, F{top}, F{side, end}, 0.05), BuildErrc::InvalidId, "not joined tangentially");
    ASSERT_EQ(census(m), n0);
}

TEST(tests_brep_face_blends, chamfer_ending_on_a_chamfer)
{
    // the edges top / side and top / x = 2 of a box chamfered together: past the end of the strip,
    // the second chamfer touches the top only and the face x = 2 the side only, no end face touches
    // both (the faces next to the strip run round the rest of the solid): refused, the model unchanged
    Model<T> m;
    const auto b = box(m, 2., 1., 1.);
    const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{edge_near(m, b, {1., 1., 1.}), edge_near(m, b, {2., 0.5, 1.})}, 0.1));
    const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
    const auto n0 = census(m);
    expect_refused(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, 0.3), BuildErrc::NonGeneric, "run round the rest of the solid");
    ASSERT_EQ(census(m), n0);
}

TEST(tests_brep_face_blends, chamfered_boss_foot)
{
    // a boss on a plate chamfered at its foot (a concave chamfer): the fillet between the plate
    // and the boss over the chamfer adds matter, the strip runs round the boss (no end face, the
    // seam of the cylinder a spoke); it is the fillet of the foot before the chamfer (a torus)
    for (T r : {0.12, 0.2})
    {
        INFO("r " << r);
        Model<T> m;
        auto boss = [&] { return fused(m, box(m, 3., 3., 0.5, {-1.5, -1.5, -0.5}), cylinder(m, 0.5, 0.8)); };
        const auto s0 = boss();
        const auto foot = edges_on_floor(m, s0, {0., 0., 0.}, 0.6);
        ASSERT_EQ(foot.size(), 1u);
        const auto ch = unwrap(chamfer(m, s0, foot, 0.05));
        const auto plate = face_near(m, ch, {1.2, 0., 0.});
        const auto wall = cylinders(m, ch);
        ASSERT_EQ(wall.size(), 1u);
        const auto f = unwrap(fillet_faces(m, ch, std::vector<FaceId>{plate}, wall, r));
        expect_valid(m, f);
        const auto ref = unwrap(fillet(m, s0, foot, r));
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-8);
    }
}

TEST(tests_brep_face_blends, chamfered_stadium_foot)
{
    // a slab on a stadium on a plate, chamfered at its foot: its side is a set of four faces
    // joined tangentially (two planes, two half cylinders), the virtual edge a closed chain
    Model<T> m;
    const auto s0 = fused(m, box(m, 6.6, 3.6, 0.4, {-1.8, -1.8, -0.4}), oblong(m, 3., 1., 0.6));
    const auto foot = edges_on_floor(m, s0, {1.5, 0., 0.}, 2.6);
    ASSERT_EQ(foot.size(), 4u);
    const auto ch = unwrap(chamfer(m, s0, foot, 0.06));
    const auto plate = face_near(m, ch, {4., 1.5, 0.});
    std::vector<FaceId> sides{face_near(m, ch, {1.5, -1., 0.3}), face_near(m, ch, {1.5, 1., 0.3}), face_near(m, ch, {4., 0., 0.3}), face_near(m, ch, {-1., 0., 0.3})};
    const T r = 0.18;
    const auto f = unwrap(fillet_faces(m, ch, std::vector<FaceId>{plate}, sides, r));
    expect_valid(m, f);
    const auto ref = unwrap(fillet(m, s0, foot, r));
    ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
    ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-7);
}

TEST(tests_brep_face_blends, stepped_box_edge)
{
    // the edge y = 1, z = 1 of the box 2 x 1 x 1 notched by a step (a strip of two faces at a
    // concave edge): the fillet over the step is that of the edge before it
    const T r = 0.3;
    Model<T> m;
    const auto s = profile(m, {{0., 0.}, {1., 0.}, {1., 0.9}, {0.95, 0.9}, {0.95, 1.}, {0., 1.}}, 2.);
    const auto top = face_near(m, s, {1., 0.5, 1.}), side = face_near(m, s, {1., 1., 0.5});
    const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{top}, std::vector<FaceId>{side}, r));
    expect_valid(m, f);
    ASSERT_EQ(n_faces(m, f), 7u);
    ASSERT_NEAR(volume(m, f), 2. * (1. - (1. - pi / 4.) * r * r), 1e-9);
}

TEST(tests_brep_face_blends, v_groove)
{
    // a groove of flanks at 25 degrees from the vertical with a flat bottom of width 0.3: the
    // concave fillet between the flanks fills the bottom, as that of the sharp groove (its
    // flanks continued down to their edge); too small a radius: the ball under the bottom
    const T b = 25. * pi / 180., f = 0.15, d = 0.75, apex = -f / std::tan(b);
    auto groove = [&](Model<T> &m, bool sharp) {
        const T w = f + d * std::tan(b);
        std::vector<std::array<T, 2>> yz{{-1., -0.6}, {1., -0.6}, {1., d}, {w, d}};
        if (sharp)
            yz.push_back({0., apex});
        else
            yz.push_back({f, 0.}), yz.push_back({-f, 0.});
        yz.push_back({-w, d}), yz.push_back({-1., d});
        return profile(m, yz, 2.);
    };
    Model<T> m;
    const auto s = groove(m, false), sharp = groove(m, true);
    const point<T, 3> on_right{1., f + 0.5 * d * std::tan(b), 0.5 * d}, on_left{1., -(f + 0.5 * d * std::tan(b)), 0.5 * d};
    const auto l = face_near(m, s, on_left), rr = face_near(m, s, on_right);
    const T r = 0.3;
    const auto x = unwrap(fillet_faces(m, s, std::vector<FaceId>{l}, std::vector<FaceId>{rr}, r));
    expect_valid(m, x);
    const auto ref = unwrap(fillet(m, sharp, std::vector<EdgeId>{edge_near(m, sharp, {1., 0., apex})}, r));
    ASSERT_EQ(n_faces(m, x), n_faces(m, ref));
    ASSERT_NEAR(volume(m, x), volume(m, ref), 1e-9);
    ASSERT_GT(volume(m, x), volume(m, s));
    const auto under = fillet_faces(m, s, std::vector<FaceId>{l}, std::vector<FaceId>{rr}, 0.12);
    ASSERT_FALSE(under.has_value());
    MESSAGE(under.error().message);
}

TEST(tests_brep_face_blends, nothing_joins)
{
    // the facing flanks of two ribs on a plate: the plate joins them, but it bounds the whole
    // solid (no strip between the ribs bounded by end faces): refused, the model unchanged
    Model<T> m;
    const auto s = profile(m, {{-1., -0.3}, {1., -0.3}, {1., 0.}, {0.6, 0.}, {0.5, 0.8}, {0.3, 0.8}, {0.4, 0.}, {-0.4, 0.}, {-0.3, 0.8}, {-0.5, 0.8}, {-0.6, 0.}, {-1., 0.}}, 2.);
    const auto w1 = face_near(m, s, {1., -0.35, 0.4}), w2 = face_near(m, s, {1., 0.35, 0.4});
    const auto before = explore<FaceId>(m, ShapeId{s}).size();
    const auto x = fillet_faces(m, s, std::vector<FaceId>{w1}, std::vector<FaceId>{w2}, 0.2);
    ASSERT_FALSE(x.has_value());
    MESSAGE(x.error().message);
    ASSERT_EQ(explore<FaceId>(m, ShapeId{s}).size(), before);
    expect_valid(m, s);
}

// ---- PR 2: partial strips (§ 4.1 to 4.4)

TEST(tests_brep_face_blends, tapered_face)
{
    // the U with a face tapering from 0.12 wide at x = 0 to a point on the corner at x = 1: the
    // wall and the floor meet along the corner from x = 1 to 2, the tapered face is under its
    // continuation (it touches both, its edges with them end on the tip, a vertex of the corner),
    // the end face x = 0 is an end face; the fillet is that of the corner of the sharp U, a
    // cylinder in one piece
    for (T r : {0.3, 0.5})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = tapered_u(m);
        expect_valid(m, s);
        ASSERT_EQ(n_faces(m, s), 11u);
        const auto t0 = std::chrono::steady_clock::now();
        const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, r));
        const T dt = seconds_since(t0);
        expect_valid(m, f);
        const auto ref = u_reference(m, r);
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        ASSERT_EQ(cylinders(m, f).size(), 1u);
        MESSAGE("tapered face, r " << r << ": volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
}

TEST(tests_brep_face_blends, tapered_face_refusals)
{
    // r = 0.2: the arc of the fillet is r (sqrt 2 - 1) = 0.083 from the corner, the widest point
    // of the face 0.085: it spills over the face, refused, the model unchanged; the fillet of the
    // corner edge alone does not cross the face, its refusal names fillet_faces
    Model<T> m;
    const auto s = tapered_u(m);
    const auto n0 = census(m);
    expect_refused(fillet_faces(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, 0.2), BuildErrc::NonGeneric, "spills over");
    ASSERT_EQ(census(m), n0);
    for (T r : {0.2, 0.3, 0.5})
    {
        INFO("r " << r);
        expect_refused(fillet(m, s, std::vector<EdgeId>{edge_near(m, s, {1.5, -Yw, 0.})}, r), BuildErrc::NonGeneric, "fillet_faces swallows them");
        ASSERT_EQ(census(m), n0);
    }
    // a face tapering to the corner from each end, the corner between them from x = 1 to 1.4:
    // two strips on both sides of a common edge, refused, without naming fillet_faces to its caller
    const T c = 0.12, dl = 0.05;
    const auto t2 = tetrahedron(m, {2., -Yw - dl, -dl}, {2., -Yw + c + dl, -dl}, {2., -Yw - dl, c + dl}, {2. - 0.6 * (c + 2. * dl) / c, -Yw - dl, -dl});
    const auto s2 = fused(m, tapered_u(m), t2);
    expect_valid(m, s2);
    const auto n2 = census(m);
    const auto r2 = fillet_faces(m, s2, std::vector<FaceId>{u_wall(m, s2)}, std::vector<FaceId>{u_floor(m, s2)}, 0.3);
    expect_refused(r2, BuildErrc::NonGeneric, "strips on both sides of an edge where they meet");
    ASSERT_EQ(r2.error().message.find("fillet_faces swallows"), std::string::npos);
    ASSERT_EQ(census(m), n2);
}

TEST(tests_brep_face_blends, stopped_chamfer)
{
    // the corner of the U chamfered from x = 0 to 1, closed there by a triangle whose edges with
    // the wall and the floor end on the vertex of the corner: the chamfer and the triangle are
    // the strip, the end face x = 0 bounds it
    for (T r : {0.3, 0.5})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = chamfered_u(m, 0., 1.);
        expect_valid(m, s);
        const auto t0 = std::chrono::steady_clock::now();
        const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, r));
        const T dt = seconds_since(t0);
        expect_valid(m, f);
        const auto ref = u_reference(m, r);
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        ASSERT_EQ(cylinders(m, f).size(), 1u);
        MESSAGE("stopped chamfer, r " << r << ": volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
    // the fillet of the corner edge alone: refused at the triangle, naming fillet_faces
    Model<T> m;
    const auto s = chamfered_u(m, 0., 1.);
    expect_refused(fillet(m, s, std::vector<EdgeId>{edge_near(m, s, {1.5, -Yw, 0.})}, 0.3), BuildErrc::NonGeneric, "fillet_faces swallows them");
    // a refusal of the fillet itself (a constant chord with a G2 section) is the one given, not
    // that of the triangle alone, a smaller candidate bounded by the chamfer, tried first and
    // failing before its copy is closed
    FilletOptions<T> o;
    o.section = BlendSection<T>::g2();
    expect_refused(fillet_faces_chordal(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, 0.4, o), BuildErrc::NonGeneric,
                   "with a G2 section");
}

TEST(tests_brep_face_blends, middle_chamfer)
{
    // the corner chamfered from x = 0.5 to 1.5: two common edges, the chamfer and its two
    // triangles between them; the corner, the virtual edge and the other corner fused into one
    // edge
    for (T r : {0.3, 0.5})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = chamfered_u(m, 0.5, 1.);
        expect_valid(m, s);
        const auto t0 = std::chrono::steady_clock::now();
        const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, r));
        const T dt = seconds_since(t0);
        expect_valid(m, f);
        const auto ref = u_reference(m, r);
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        ASSERT_EQ(cylinders(m, f).size(), 1u);
        MESSAGE("middle chamfer, r " << r << ": volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
}

TEST(tests_brep_face_blends, stopped_box_chamfer)
{
    // the edge y = 1, z = 1 of the box chamfered from x = 1 on (a wedge cut away): the top and
    // the side meet along it up to x = 1, across the chamfer past it, the triangle at x = 1 of
    // the strip (refused by PR 1): the fillet of the edge of the box
    const T r = 0.3;
    Model<T> m;
    const auto notched = unwrap(cut(m, ShapeId{box(m, 2., 1., 1.)}, ShapeId{profile(m, {{0.9, 1.}, {1., 0.9}, {1.5, 1.5}}, 2., 1.)})).solids.front();
    expect_valid(m, notched);
    const auto t = face_near(m, notched, {0.5, 0.5, 1.}), sd = face_near(m, notched, {0.5, 1., 0.5});
    ASSERT_EQ(face_near(m, notched, {1.5, 0.5, 1.}), t);
    const auto f = unwrap(fillet_faces(m, notched, std::vector<FaceId>{t}, std::vector<FaceId>{sd}, r));
    expect_valid(m, f);
    ASSERT_EQ(n_faces(m, f), 7u);
    ASSERT_NEAR(volume(m, f), 2. * (1. - (1. - pi / 4.) * r * r), 1e-9);
}

TEST(tests_brep_face_blends, boss_foot_chamfered_half_round)
{
    // a boss on a plate, its foot chamfered over half a turn (a half ring fused there, closed by
    // two triangles on the planes x = 0): the plate and the boss meet along the rest of the foot,
    // the chamfer runs on from it; the fillet is that of the foot of the boss without the chamfer,
    // a torus in one piece. The half ring on the side x < 0 (the seam of the boss at x > 0 cuts
    // the rest of the foot in two), or on the side x > 0, over the seam: the rest of the foot
    // one edge from the strip to the strip, under the virtual edge from end to end
    for (const bool over_seam : {false, true})
    for (T r : {0.12, 0.2})
    {
        INFO("over the seam " << over_seam << ", r " << r);
        Model<T> m;
        auto boss = [&] { return fused(m, box(m, 3., 3., 0.5, {-1.5, -1.5, -0.5}), cylinder(m, 0.5, 0.8)); };
        const auto s0 = boss();
        const auto foot = edges_on_floor(m, s0, {0., 0., 0.}, 0.6);
        ASSERT_EQ(foot.size(), 1u);
        const auto tri = planar(m, {{0., 0.48, -0.02}, {0., 0.57, -0.02}, {0., 0.48, 0.07}});
        const auto ring = solid_of(m, make_revolution(m, ShapeId{tri}, ax1<T, 3>{{{0., 0., 0.}, over_seam ? -1. * Z : Z}}, pi));
        const auto s = fused(m, boss(), ring);
        expect_valid(m, s);
        ASSERT_EQ(edges_on_floor(m, s, {0., 0., 0.}, 0.6).size(), over_seam ? 4u : 5u);
        const auto plate = face_near(m, s, {1.2, 0., 0.});
        const auto wall = cylinders(m, s);
        ASSERT_EQ(wall.size(), 1u);
        const auto t0 = std::chrono::steady_clock::now();
        const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{plate}, wall, r));
        const T dt = seconds_since(t0);
        expect_valid(m, f);
        const auto ref = unwrap(fillet(m, s0, foot, r));
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        MESSAGE("boss foot chamfered on half a turn" << std::string(over_seam ? " over its seam" : "") << ", r " << r << ": volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-8);
    }
}

// ---- PR 2: laws by position (§ 4.5), sections (§ 4.6)

TEST(tests_brep_face_blends, law_by_position)
{
    // the chamfered box, a radius from 0.2 at x = 0 to 0.4 at x = 2 along the virtual edge: a
    // point by position, a point on an edge of the strip (its place on the solid given): the
    // solid of the same law on the edge of the box without the chamfer
    for (const auto variation : {RadiusVariation::Linear, RadiusVariation::Smooth})
    {
        Model<T> m;
        const auto b = box(m, 2., 1., 1.);
        const auto e = edge_near(m, b, {1., 1., 1.});
        const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{e}, 0.1));
        const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
        // the edge between the top and the chamfer, its end at x = 2
        const auto ec = edge_near(m, ch, {1., 0.9, 1.});
        const auto &ed = m.edge(ec);
        const T at2 = ed.curve->value(ed.u1)[0] > 1. ? 0. : 1.;
        RadiusLaw<T> law{{RadiusPoint<T>::at_point({0., 1.03, 1.05}, 0.2), RadiusPoint<T>{ec, at2, 0.4}}, variation};
        const auto t0 = std::chrono::steady_clock::now();
        const auto f = unwrap(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, law));
        const T dt = seconds_since(t0);
        expect_valid(m, f);
        const T at0 = m.edge(e).curve->value(m.edge(e).u1)[0] < 1. ? 0. : 1.;
        const auto ref = unwrap(fillet(m, b, std::vector<EdgeId>{e}, RadiusLaw<T>{{{e, at0, 0.2}, {e, 1. - at0, 0.4}}, variation}));
        ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
        MESSAGE("law by position on the chamfered box: volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
}

TEST(tests_brep_face_blends, edge_law_by_position)
{
    // the fillet of an edge of the box with points by position: the solid of the law (edge,
    // fraction) at their feet
    Model<T> m;
    const auto b = box(m, 2., 1., 1.);
    const auto e = edge_near(m, b, {1., 1., 1.});
    const T at0 = m.edge(e).curve->value(m.edge(e).u1)[0] < 1. ? 0. : 1.;
    const auto a = unwrap(fillet(m, b, std::vector<EdgeId>{e},
                                 RadiusLaw<T>{{RadiusPoint<T>::at_point({-0.3, 1.1, 1.2}, 0.15), RadiusPoint<T>::at_point({0.5, 0.95, 1.1}, 0.3),
                                               RadiusPoint<T>::at_point({2.2, 1.05, 0.9}, 0.2)}}));
    const auto ref = unwrap(fillet(m, b, std::vector<EdgeId>{e},
                                   RadiusLaw<T>{{{e, at0, 0.15}, {e, at0 == 0. ? 0.25 : 0.75, 0.3}, {e, 1. - at0, 0.2}}}));
    expect_valid(m, a);
    ASSERT_EQ(n_faces(m, a), n_faces(m, ref));
    MESSAGE("edge law by position: volume error " << std::abs(volume(m, a) - volume(m, ref)));
    // the feet of the points are projections: 2.5e-11 on x86-64 (Linux, Windows), below 1e-12
    // on arm64 (macOS)
    ASSERT_NEAR(volume(m, a), volume(m, ref), 1e-9);
    // the geometry of the edge alone with a law of one value given by position (a radius, the
    // distances of a chamfer, the rho of a conic section): that of the value (its points by
    // position have no edge: the value was read past the points of the law)
    auto same = [&](const EdgeBlend<T> &x, const EdgeBlend<T> &y) {
        for (std::size_t k = 0; k < 2; ++k)
            for (const T w : {0.1, 0.5, 0.9})
            {
                const T t = x.a + (x.b - x.a) * w;
                ASSERT_LT(norm(x.contacts[k].curve->value(t) - y.contacts[k].curve->value(t)), 1e-12);
            }
    };
    const RadiusLaw<T> one{{RadiusPoint<T>::at_point({1., 1.1, 1.1}, 0.2)}};
    same(unwrap(fillet_geometry(m, b, e, one)), unwrap(fillet_geometry(m, b, e, 0.2)));
    const auto top = face_near(m, b, {1., 0.5, 1.});
    const BlendLaw<T> d{{BlendPoint<T>::at_point({1., 1.1, 1.1}, 0.1)}}, d0{{{e, 0., 0.1}}};
    same(unwrap(chamfer_geometry(m, b, e, top, d, d)), unwrap(chamfer_geometry(m, b, e, top, d0, d0)));
    FilletOptions<T> o, o0;
    o.section = BlendSection<T>::conic(BlendLaw<T>{{BlendPoint<T>::at_point({1., 1.1, 1.1}, 0.4)}});
    o0.section = BlendSection<T>::conic(0.4);
    same(unwrap(fillet_geometry(m, b, e, 0.2, o)), unwrap(fillet_geometry(m, b, e, 0.2, o0)));
}

TEST(tests_brep_face_blends, law_by_position_refusals)
{
    // a point away from every edge, two points at the same foot with different values: refused
    // (InvalidBounds), for a fillet of faces and of an edge, the model unchanged
    Model<T> m;
    const auto b = box(m, 2., 1., 1.);
    const auto e = edge_near(m, b, {1., 1., 1.});
    const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{e}, 0.1));
    const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
    const auto n0 = census(m);
    const RadiusLaw<T> far{{RadiusPoint<T>::at_point({0., 1., 1.}, 0.3), RadiusPoint<T>::at_point({1., 10., 10.}, 0.3)}};
    const RadiusLaw<T> twice{{RadiusPoint<T>::at_point({1., 1.1, 1.1}, 0.3), RadiusPoint<T>::at_point({1., 1.2, 1.2}, 0.35)}};
    using F = std::vector<FaceId>;
    expect_refused(fillet_faces(m, ch, F{top}, F{side}, far), BuildErrc::InvalidBounds, "away from the edges");
    expect_refused(fillet_faces(m, ch, F{top}, F{side}, twice), BuildErrc::InvalidBounds, "at the same foot");
    expect_refused(fillet(m, b, std::vector<EdgeId>{e}, far), BuildErrc::InvalidBounds, "away from the edges");
    expect_refused(fillet(m, b, std::vector<EdgeId>{e}, twice), BuildErrc::InvalidBounds, "at the same foot");
    ASSERT_EQ(census(m), n0);
}

TEST(tests_brep_face_blends, sections)
{
    // the sections of stage 7 on the virtual edge: a conic of rho 0.5, the G2 quartic, the
    // constant chord, on the chamfered box and on the tapered U; each the solid of the fillet of
    // the edge without the strip with the same section
    struct Case
    {
        std::string name;
        BlendSection<T> section;
        bool chordal;
    };
    for (const auto &[name, section, chordal] : {Case{"conic rho 0.5", BlendSection<T>::conic(0.5), false}, Case{"G2", BlendSection<T>::g2(), false},
                                                  Case{"constant chord", BlendSection<T>::circular(), true}})
    {
        INFO(name);
        FilletOptions<T> o;
        o.section = section;
        {
            // the chamfered box: r = 0.4, or the chord 0.6 (r 0.42 between the planes at a right angle)
            Model<T> m;
            const auto b = box(m, 2., 1., 1.);
            const auto e = edge_near(m, b, {1., 1., 1.});
            const auto ch = unwrap(chamfer(m, b, std::vector<EdgeId>{e}, 0.1));
            const auto top = face_near(m, ch, {1., 0.5, 1.}), side = face_near(m, ch, {1., 1., 0.5});
            const auto t0 = std::chrono::steady_clock::now();
            const auto f = chordal ? unwrap(fillet_faces_chordal(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, 0.6, o))
                                   : unwrap(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{side}, 0.4, o));
            const T dt = seconds_since(t0);
            expect_valid(m, f);
            const auto ref = chordal ? unwrap(fillet_chordal(m, b, std::vector<EdgeId>{e}, 0.6, o)) : unwrap(fillet(m, b, std::vector<EdgeId>{e}, 0.4, o));
            ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
            MESSAGE(name << " on the chamfered box: volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
            ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-11);
        }
        {
            // the tapered U: r = 0.5, or the chord 0.7
            Model<T> m;
            const auto s = tapered_u(m);
            const auto t0 = std::chrono::steady_clock::now();
            const auto f = chordal ? unwrap(fillet_faces_chordal(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, 0.7, o))
                                   : unwrap(fillet_faces(m, s, std::vector<FaceId>{u_wall(m, s)}, std::vector<FaceId>{u_floor(m, s)}, 0.5, o));
            const T dt = seconds_since(t0);
            expect_valid(m, f);
            const auto u = u_channel(m);
            const auto corner = std::vector<EdgeId>{edge_near(m, u, {1., -Yw, 0.})};
            const auto ref = chordal ? unwrap(fillet_chordal(m, u, corner, 0.7, o)) : unwrap(fillet(m, u, corner, 0.5, o));
            ASSERT_EQ(n_faces(m, f), n_faces(m, ref));
            MESSAGE(name << " on the tapered U: volume error " << std::abs(volume(m, f) - volume(m, ref)) << ", " << dt << " s");
            // (the same geometry as the fillet of the edge: fillet_faces_chordal that of fillet_chordal)
            ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-11);
        }
    }
}

namespace
{
    /// A rib of width w at its top (z = 1 + h) and of draft a (its flanks nearer at the top; a < 0: an undercut, wider at the top) on a plate 4 x 1, along x over L.
    SolidId drafted_rib(Model<T> &m, T w, T h, T a, T L)
    {
        const T d = h * std::tan(a);
        return profile(m, {{-2., 0.}, {2., 0.}, {2., 1.}, {w / 2. + d, 1.}, {w / 2., 1. + h}, {-w / 2., 1. + h}, {-w / 2. - d, 1.}, {-2., 1.}}, L);
    }
}

TEST(tests_brep_face_blends, flanks_meeting_on_the_far_side)
{
    // the flanks of a rib in undercut (-10 degrees: wider at its top) meet, once continued, below
    // it, inside the plate, not above its top: no fillet between them over the top. Refused, the
    // model as it was, at radii about that of the ball touching the top and both flanks (which
    // a full round fillet builds); the plate was taken as the strip, and swallowed (PR 1 and 2)
    const T w = 0.6, h = 1., L = 2., a = -10. * pi / 180., d = h * std::tan(a) / 2.;
    const T rr = (w / 2.) * std::cos(a) / (1. - std::sin(a));
    for (const T k : {0.5, 0.999, 1., 1.001, 1.1, 2.})
    {
        INFO("r = " << k << " r_fr");
        Model<T> m;
        const auto s = drafted_rib(m, w, h, a, L);
        const auto f1 = face_near(m, s, {1., -w / 2. - d, 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2. + d, 1. + h / 2.});
        const auto before = census(m);
        const auto r = fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, k * rr);
        if (r)
            MESSAGE("accepted: " << n_faces(m, *r) << " faces, volume " << volume(m, *r) << " of " << volume(m, s));
        expect_refused(r, BuildErrc::NonGeneric, "on the far side");
        ASSERT_EQ(census(m), before);
    }
}

TEST(tests_brep_face_blends, flanks_meeting_beyond_the_top)
{
    // the same rib with a draft of +10 degrees: its flanks meet above its top, which a fillet of
    // the flanks larger than the full round swallows; the solid of the fillet of the apex of the
    // rib continued to it (the check of the far side does not refuse the good side)
    const T w = 0.6, h = 1., L = 2., a = 10. * pi / 180., d = h * std::tan(a) / 2.;
    const T apex = 1. + h + (w / 2.) / std::tan(a);
    for (const T r : {0.36, 0.4})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = drafted_rib(m, w, h, a, L);
        const auto f1 = face_near(m, s, {1., -w / 2. - d, 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2. + d, 1. + h / 2.});
        const auto f = unwrap(fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, r));
        expect_valid(m, f);
        ASSERT_EQ(n_faces(m, f), n_faces(m, s));
        const T b = w / 2. + 2. * d;
        const auto sharp = profile(m, {{-2., 0.}, {2., 0.}, {2., 1.}, {b, 1.}, {0., apex}, {-b, 1.}, {-2., 1.}}, L);
        const auto ref = unwrap(fillet(m, sharp, std::vector<EdgeId>{edge_near(m, sharp, {1., 0., apex})}, r));
        MESSAGE("fillet_faces over the top: volume " << volume(m, f) << ", the apex filleted " << volume(m, ref));
        ASSERT_NEAR(volume(m, f), volume(m, ref), 1e-9);
    }
}

TEST(tests_brep_face_blends, strip_in_part_beyond_the_sides)
{
    // (review of PR 253) the rib of a draft of +10 degrees, a bump on its top next to a flank,
    // leaning out of it: past the flank continued, out of the sharp rib, the rest of the top in
    // it. The flanks meet above the top, on the near side: refused, not as a strip on the far
    // side, the model as it was (the fillet swallowed the bump, out of the result)
    const T w = 0.6, h = 1., L = 2., a = 10. * pi / 180., d = h * std::tan(a);
    // the flank continued at the height of the top of the bump, 2.8: inside the bump
    ASSERT_LT(w / 2. - d * (2.8 - (1. + h)), 0.18);
    for (const T r : {0.36, 0.4})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = profile(m,
                               {{-2., 0.}, {2., 0.}, {2., 1.}, {w / 2. + d, 1.}, {w / 2., 1. + h}, {0.25, 1. + h}, {0.22, 2.8}, {0.18, 2.8}, {0.15, 1. + h},
                                {-w / 2., 1. + h}, {-w / 2. - d, 1.}, {-2., 1.}},
                               L);
        const auto f1 = face_near(m, s, {1., -w / 2. - d / 2., 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2. + d / 2., 1. + h / 2.});
        const auto before = census(m);
        const auto x = fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, r);
        if (x)
            MESSAGE("accepted: " << n_faces(m, *x) << " faces, volume " << volume(m, *x) << " of " << volume(m, s));
        expect_refused(x, BuildErrc::NonGeneric, "in part beyond those continued");
        ASSERT_EQ(census(m), before);
    }
}

TEST(tests_brep_face_blends, walls_meeting_on_the_far_side)
{
    // the walls of a dovetail slot (wider at its bottom) meet, once continued, above its opening,
    // out of the solid, not under its bottom: refused, the model as it was. (Concave: the strip
    // taken, the top of the block, is out of the copy as a strip under a concave virtual edge is;
    // the copy closed over it is invalid, and refused as such. Were it valid, the strip would not
    // be in the result, which only adds the fillet to the copy: refused as spilling over. A
    // concave fillet cannot swallow a strip on the far side; the convex one could, before.)
    const T w = 0.6, h = 1., L = 2., a = 10. * pi / 180., d = h * std::tan(a);
    for (const T r : {0.2, 0.4})
    {
        INFO("r " << r);
        Model<T> m;
        const auto s = profile(m, {{-1.5, 0.}, {1.5, 0.}, {1.5, 2.}, {w / 2., 2.}, {w / 2. + d, 1.}, {-w / 2. - d, 1.}, {-w / 2., 2.}, {-1.5, 2.}}, L);
        const auto f1 = face_near(m, s, {1., -w / 2. - d / 2., 1.5}), f2 = face_near(m, s, {1., w / 2. + d / 2., 1.5});
        ASSERT_TRUE(f1.valid() && f2.valid());
        const auto before = census(m);
        const auto x = fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, r);
        if (x)
            MESSAGE("accepted: " << n_faces(m, *x) << " faces, volume " << volume(m, *x) << " of " << volume(m, s));
        expect_refused(x, BuildErrc::NonGeneric, "invalid solid");
        ASSERT_EQ(census(m), before);
    }
}
