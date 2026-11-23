// Fillets between two faces that have no edge in common (stage 8, PR 1,
// docs/sources/design/brep_palier8_face_blends.md, docs/sources/design/brep8_pr01_face_blends.md):
// the strip between them swallowed, the fillet that of the virtual edge where the faces meet once
// continued. Each case against the fillet of the edge that existed before the strip (a chamfer):
// the same solid, the same volume; refusals.
#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/face_blend.h>
#include <gbs-brep/sweep.h>
#include <gbs/bselementary.h>

#include <array>
#include <cmath>
#include <numbers>
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

    /// The prism along x of length L on the profile (y, z) given in the plane x = 0 (counter-clockwise about x).
    SolidId profile(Model<T> &m, const std::vector<std::array<T, 2>> &yz, T L)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < yz.size(); ++i)
        {
            const auto &a = yz[i], &b = yz[(i + 1) % yz.size()];
            e.push_back(unwrap(make_edge(m, point<T, 3>{0., a[0], a[1]}, point<T, 3>{0., b[0], b[1]})));
        }
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{L, 0., 0.}));
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
    const auto parallel = fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{bottom}, 0.3);
    ASSERT_FALSE(parallel.has_value());
    MESSAGE(parallel.error().message);
    // faces with an edge in common: the fillet of that edge
    const auto adj = unwrap(fillet_faces(m, ch, std::vector<FaceId>{top}, std::vector<FaceId>{back}, 0.2));
    expect_valid(m, adj);
    ASSERT_NEAR(volume(m, adj), volume(m, ch) - 2. * (1. - pi / 4.) * 0.04, 1e-9);
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
