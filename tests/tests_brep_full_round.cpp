// Full round fillets (stage 8, docs/sources/design/brep_palier8_face_blends.md, § 3; PR 3,
// docs/sources/design/brep8_pr03_full_round.md): the middle face between two faces replaced by
// the fillet tangent to the three. Ribs against closed forms (a half cylinder between parallel
// flanks, a cylinder of radius (w / 2) cos a / (1 - sin a) with a draft a), the bottom of a slot
// (concave), a tapered rib and the blunt trailing edge of a blade (approximated surfaces: contacts
// tangent, radius smooth), refusals.
#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/full_round.h>
#include <gbs-brep/sew.h>
#include <gbs-brep/sweep.h>
#include <gbs/bscinterp.h>
#include <gbs/bssbuild.h>
#include <blade_fixture.h>

#include <array>
#include <iomanip>
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

    /// The face of s through p (its surface within 1e-7, p inside it).
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

    /// The faces of s on a surface of type S.
    template <typename S>
    std::vector<FaceId> faces_on(const Model<T> &m, SolidId s)
    {
        std::vector<FaceId> r;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            if (dynamic_cast<const S *>(m.face(f).surface.get()))
                r.push_back(f);
        return r;
    }

    /// A rib of width w at its top (z = 1 + h) and of draft a (its flanks nearer at the top) on a plate 4 x 1, along x over L.
    SolidId rib(Model<T> &m, T w, T h, T a, T L)
    {
        const T d = h * std::tan(a);
        return profile(m, {{-2., 0.}, {2., 0.}, {2., 1.}, {w / 2. + d, 1.}, {w / 2., 1. + h}, {-w / 2., 1. + h}, {-w / 2. - d, 1.}, {-2., 1.}}, L);
    }

    /// The options of the tests of volume: tol = 1e-7, the volumes then within 1e-9 of the closed forms (the curves at
    /// the free ends are within tol / 10 of their half circles; at the default tol, 1e-6, the volumes within 3e-9)
    FilletOptions<T> tight()
    {
        FilletOptions<T> o;
        o.tol = 1e-7;
        return o;
    }

    auto ms(std::chrono::steady_clock::time_point t0)
    {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
}

TEST(tests_brep_full_round, straight_rib)
{
    // the top of a rib with parallel flanks: a half cylinder of radius w / 2, exact; the area
    // taken off its section w^2 / 2 - pi w^2 / 8
    for (const T w : {0.4, 0.8})
    {
        Model<T> m;
        const T h = 1., L = 3.;
        const auto s = rib(m, w, h, 0., L);
        const auto top = face_near(m, s, {1., 0., 1. + h}), f1 = face_near(m, s, {1., -w / 2., 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2., 1. + h / 2.});
        const auto before = census(m);
        const auto t0 = std::chrono::steady_clock::now();
        const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, tight()));
        const auto dt = ms(t0);
        expect_valid(m, r);
        ASSERT_EQ(n_faces(m, r), n_faces(m, s));
        const auto cyl = faces_on<CylindricalSurface<T>>(m, r);
        ASSERT_EQ(cyl.size(), 1u);
        const auto &c = static_cast<const CylindricalSurface<T> &>(*m.face(cyl.front()).surface);
        ASSERT_NEAR(c.radius(), w / 2., 1e-12);
        const auto bd = c.bounds();
        ASSERT_NEAR(bd[1] - bd[0], pi, 1e-9);
        const T expected = volume(m, s) - L * (w * w / 2. - pi * w * w / 8.);
        MESSAGE("w = " << w << ": volume gap " << volume(m, r) - expected << ", " << dt << " ms");
        ASSERT_NEAR(volume(m, r), expected, 1e-9);
        // the solid given untouched, nothing left over
        expect_valid(m, s);
        const auto after = census(m), mine = std::array<std::size_t, 3>{explore<VertexId>(m, ShapeId{r}).size(), explore<EdgeId>(m, ShapeId{r}).size(), n_faces(m, r)};
        for (std::size_t i = 0; i < 3; ++i)
            ASSERT_EQ(after[i], before[i] + mine[i]);
    }
}

TEST(tests_brep_full_round, drafted_rib)
{
    // flanks at the draft a: a cylinder of radius (w / 2) cos a / (1 - sin a), its arc pi - 2 a;
    // the area taken off: the two corners between the top, a flank and the arc,
    // r^2 (tan(b) - b) with b = pi / 4 - a / 2 the half angle of each corner seen from the centre
    for (const T deg : {5., 15., -10.})
    {
        Model<T> m;
        const T w = 0.6, h = 1., L = 2., a = deg * pi / 180.;
        const auto s = rib(m, w, h, a, L);
        const T d = h * std::tan(a) / 2.;
        const auto top = face_near(m, s, {1., 0., 1. + h}), f1 = face_near(m, s, {1., -w / 2. - d, 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2. + d, 1. + h / 2.});
        const auto t0 = std::chrono::steady_clock::now();
        const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, tight()));
        const auto dt = ms(t0);
        expect_valid(m, r);
        const auto cyl = faces_on<CylindricalSurface<T>>(m, r);
        ASSERT_EQ(cyl.size(), 1u);
        const auto &c = static_cast<const CylindricalSurface<T> &>(*m.face(cyl.front()).surface);
        const T rr = (w / 2.) * std::cos(a) / (1. - std::sin(a));
        ASSERT_NEAR(c.radius(), rr, 1e-12);
        const auto bd = c.bounds();
        ASSERT_NEAR(bd[1] - bd[0], pi - 2. * a, 1e-9);
        const T b = pi / 4. - a / 2.;
        const T expected = volume(m, s) - L * 2. * rr * rr * (std::tan(b) - b);
        MESSAGE("draft " << deg << " deg: r = " << rr << ", volume gap " << volume(m, r) - expected << ", " << dt << " ms");
        ASSERT_NEAR(volume(m, r), expected, 1e-9);
    }
}

namespace
{
    /// The prism along z of height h on the polygon (x, y) (counter-clockwise).
    SolidId tower(Model<T> &m, const std::vector<std::array<T, 2>> &xy, T h)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < xy.size(); ++i)
        {
            const auto &a = xy[i], &b = xy[(i + 1) % xy.size()];
            e.push_back(unwrap(make_edge(m, point<T, 3>{a[0], a[1], 0.}, point<T, 3>{b[0], b[1], 0.})));
        }
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
    }

    /// The area between a right corner and the arc of radius r tangent to both its sides, and the distance of its centroid from the centre along a side.
    constexpr T corner_area(T r) { return r * r * (1. - pi / 4.); }
    constexpr T corner_centroid(T r) { return r / (6. * (1. - pi / 4.)); }

    /// The face of r that is none of those of s (by its surface): the fillet.
    FaceId new_face(const Model<T> &m, SolidId r, SolidId s)
    {
        std::vector<const Surface<T, 3> *> old;
        for (auto f : explore<FaceId>(m, ShapeId{s}))
            old.push_back(m.face(f).surface.get());
        FaceId found{};
        int n = 0;
        for (auto f : explore<FaceId>(m, ShapeId{r}))
        {
            const auto *srf = m.face(f).surface.get();
            // (a face cut keeps its surface, or that surface continued: a plane stays a plane)
            if (std::ranges::contains(old, srf) || dynamic_cast<const Plane<T> *>(srf))
                continue;
            found = f, ++n;
        }
        return n == 1 ? found : FaceId{};
    }
}

TEST(tests_brep_full_round, tolerance)
{
    // the gap to the closed form is that of the curves where the end faces cut the fillet
    // (cubic B-splines within tol / 10 of the half circles, as at the end of any fillet): at
    // tol = 1e-9 the volume is the closed form within 1e-10
    Model<T> m;
    const T w = 0.8, h = 1., L = 3.;
    const auto s = rib(m, w, h, 0., L);
    const auto top = face_near(m, s, {1., 0., 1. + h}), f1 = face_near(m, s, {1., -w / 2., 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2., 1. + h / 2.});
    FilletOptions<T> o;
    o.tol = 1e-9;
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, o));
    expect_valid(m, r);
    const T expected = volume(m, s) - L * 2. * corner_area(w / 2.);
    MESSAGE("tol 1e-9: volume gap " << volume(m, r) - expected);
    ASSERT_NEAR(volume(m, r), expected, 1e-10);
}

TEST(tests_brep_full_round, slot_bottom)
{
    // the bottom of a slot between parallel walls: concave, the ball outside the solid, the
    // fillet adds the two corners between the bottom, a wall and the half circle
    Model<T> m;
    const T w = 0.6, L = 2.;
    const auto s = profile(m, {{-1., 0.}, {1., 0.}, {1., 2.}, {w / 2., 2.}, {w / 2., 1.}, {-w / 2., 1.}, {-w / 2., 2.}, {-1., 2.}}, L);
    const auto bottom = face_near(m, s, {1., 0., 1.}), f1 = face_near(m, s, {1., -w / 2., 1.5}), f2 = face_near(m, s, {1., w / 2., 1.5});
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{bottom}, std::vector<FaceId>{f2}, tight()));
    const auto dt = ms(t0);
    expect_valid(m, r);
    const auto cyl = faces_on<CylindricalSurface<T>>(m, r);
    ASSERT_EQ(cyl.size(), 1u);
    ASSERT_NEAR(static_cast<const CylindricalSurface<T> &>(*m.face(cyl.front()).surface).radius(), w / 2., 1e-12);
    const T expected = volume(m, s) + L * 2. * corner_area(w / 2.);
    MESSAGE("slot: volume gap " << volume(m, r) - expected << ", " << dt << " ms");
    ASSERT_NEAR(volume(m, r), expected, 1e-9);
}

TEST(tests_brep_full_round, oblique_end)
{
    // a rib standing alone whose end is oblique (at 30 and 60 degrees from the square end): the
    // end face cuts the half cylinder along an arc of ellipse, further along one flank than
    // along the other; the corners taken off are prisms cut by that plane, their volume the
    // area of a corner by its length at its centroid, the sum of both by the length at the middle
    for (const T deg : {30., 60.})
    {
        Model<T> m;
        const T w = 0.5, h = 1., L = 2., k = std::tan(deg * pi / 180.);
        const auto s = tower(m, {{0., -w / 2.}, {L, -w / 2.}, {L + w * k, w / 2.}, {0., w / 2.}}, h);
        const auto top = face_near(m, s, {1., 0., h}), f1 = face_near(m, s, {1., -w / 2., h / 2.}), f2 = face_near(m, s, {1., w / 2., h / 2.});
        const auto t0 = std::chrono::steady_clock::now();
        const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, tight()));
        const auto dt = ms(t0);
        expect_valid(m, r);
        ASSERT_EQ(faces_on<CylindricalSurface<T>>(m, r).size(), 1u);
        const T expected = volume(m, s) - 2. * corner_area(w / 2.) * (L + w * k / 2.);
        MESSAGE("end at " << deg << " deg: volume gap " << volume(m, r) - expected << ", " << dt << " ms");
        ASSERT_NEAR(volume(m, r), expected, 1e-9);
    }
}

TEST(tests_brep_full_round, tapered_rib)
{
    // a rib whose top narrows from 0.6 to 0.2 between vertical flanks: three planes meeting at a
    // point A, the balls of radius proportional to their distance from A: the cone of apex A,
    // exact; its contacts tangent to the flanks
    Model<T> m;
    const T h = 1., L = 2., w0 = 0.6, w1 = 0.2;
    const auto s = tower(m, {{0., -w0 / 2.}, {L, -w1 / 2.}, {L, w1 / 2.}, {0., w0 / 2.}}, h);
    const auto top = face_near(m, s, {1., 0., h}), f1 = face_near(m, s, {1., -(w0 + w1) / 4., h / 2.}), f2 = face_near(m, s, {1., (w0 + w1) / 4., h / 2.});
    const auto g = unwrap(full_round_geometry(m, s, f1, top, f2));
    // (three planes, the radius linear: the exact cone)
    ASSERT_TRUE(g.exact);
    ASSERT_NE(dynamic_cast<const ConicalSurface<T> *>(g.surface.get()), nullptr);
    // the apex: where the flanks and the top meet, x = L w0 / (w0 - w1); the centres on the line
    // from it, the radius proportional to the distance from it
    const point<T, 3> A{L * w0 / (w0 - w1), 0., h};
    const auto c0 = g.spine->value(g.t1);
    const T r0 = norm(c0 - g.contacts[0].curve->value(g.t1)), d0 = norm(c0 - A);
    for (int i = 0; i <= 16; ++i)
    {
        const T t = g.t1 + (g.t2 - g.t1) * T(i) / 16.;
        const auto c = g.spine->value(t);
        ASSERT_LE(norm(cross(c - A, c0 - A)) / d0, 1e-5);
        for (const auto &ct : g.contacts)
            ASSERT_NEAR(norm(c - ct.curve->value(t)), r0 * norm(c - A) / d0, 1e-5);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}));
    const auto dt = ms(t0);
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), n_faces(m, s));
    // every point of the fillet on the cone of apex A round the line of the centres
    const auto f = new_face(m, r, s);
    ASSERT_TRUE(f.valid());
    const auto axis = (c0 - A) / d0;
    const T tan_half = r0 / std::sqrt(d0 * d0 - r0 * r0);
    T worst{};
    for (const auto &uv : interior_points(m, f, 1e-7))
    {
        const auto p = m.face(f).surface->value(uv[0], uv[1]) - A;
        const T along = p * axis;
        worst = std::max(worst, std::abs(norm(p - along * axis) - along * tan_half));
    }
    MESSAGE("tapered rib: " << dt << " ms, gap to the cone " << worst << ", volume taken off " << volume(m, s) - volume(m, r));
    ASSERT_LE(worst, 1e-12);
    ASSERT_GT(volume(m, s), volume(m, r));
}

namespace
{
    /// The faces of a blade with a blunt trailing edge: its suction side, its pressure side, the face of its trailing edge.
    struct Blade
    {
        SolidId solid;
        FaceId suction, pressure, trailing;
    };

    /**
     * A blade of span 1.5 along z, its chord about 1 along x: two cubic B-spline sides from a
     * sharp leading edge to a trailing edge of thickness 0.08 at the root and 0.044 at the tip,
     * swept back and twisted (the trailing edge face ruled between the ends of the sides), closed
     * by the planes of the root and of the tip. The three faces have the same poles along their
     * common edges.
     */
    Blade blunt_blade(Model<T> &m)
    {
        constexpr std::size_t nu = 5, nv = 4;
        const std::array<T, nv> z{0., 0.5, 1., 1.5}, e{0.04, 0.035, 0.03, 0.022}, xte{1., 1.04, 1.1, 1.18}, twist{0., 0.06, 0.12, 0.2};
        const std::array<std::array<T, 2>, nu - 1> up{{{0., 0.}, {0.05, 0.08}, {0.4, 0.16}, {0.8, 0.1}}}, down{{{0., 0.}, {0.05, -0.05}, {0.4, -0.04}, {0.8, -0.03}}};
        auto pole = [&](const std::array<T, 2> &p, std::size_t j) {
            const T c = std::cos(twist[j]), s = std::sin(twist[j]);
            return point<T, 3>{c * p[0] - s * p[1], s * p[0] + c * p[1], z[j]};
        };
        points_vector<T, 3> pu, pd, pt;
        for (std::size_t j = 0; j < nv; ++j)
        {
            for (std::size_t i = 0; i + 1 < nu; ++i)
                pu.push_back(pole(up[i], j)), pd.push_back(pole(down[i], j));
            pu.push_back(pole({xte[j], e[j]}, j)), pd.push_back(pole({xte[j], -e[j]}, j));
            pt.push_back(pole({xte[j], -e[j]}, j)), pt.push_back(pole({xte[j], e[j]}, j));
        }
        const std::vector<T> ku{0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, kv{0., 0., 0., 0., 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const auto su = std::make_shared<BSSurface<T, 3>>(pu, ku, kv, std::size_t{3}, std::size_t{3});
        const auto sd = std::make_shared<BSSurface<T, 3>>(pd, ku, kv, std::size_t{3}, std::size_t{3});
        const auto st = std::make_shared<BSSurface<T, 3>>(pt, k1, kv, std::size_t{1}, std::size_t{3});
        std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{su})), unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{sd})),
                                  unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{st}))};
        for (const std::size_t j : {std::size_t{0}, nv - 1})
        {
            // the section there: the first or the last row of poles of each side, and the trailing edge
            auto row = [&](const points_vector<T, 3> &p) {
                return std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>(p.begin() + std::ptrdiff_t(j * nu), p.begin() + std::ptrdiff_t((j + 1) * nu)), ku, 3);
            };
            const auto cu = row(pu), cd = row(pd);
            const auto v0 = unwrap(make_vertex(m, cu->value(0.))), vu = unwrap(make_vertex(m, cu->value(1.))), vd = unwrap(make_vertex(m, cd->value(1.)));
            const std::vector<EdgeId> es{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{cd}, 0., 1., v0, vd)), unwrap(make_edge(m, vd, vu)),
                                         unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{cu}, 0., 1., v0, vu))};
            const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., z[j]}, Z, X}}, std::array<T, 4>{-5., 5., -5., 5.});
            faces.push_back(unwrap(make_face(m, plane, unwrap(make_wire(m, es)))));
        }
        const auto sewn = unwrap(sew(m, faces));
        const auto s = unwrap(make_solid(m, sewn.shells.front()));
        return {s, face_near(m, s, su->value(0.6, 0.5)), face_near(m, s, sd->value(0.6, 0.5)), face_near(m, s, st->value(0.5, 0.5))};
    }
}

TEST(tests_brep_full_round, blade_trailing_edge)
{
    // the blunt trailing edge of a blade between its B-spline sides: the face of the trailing
    // edge replaced by the fillet tangent to it and to both sides, its radius following the
    // thickness; the surface approximated within approx_tol, its contacts tangent to the sides
    Model<T> m;
    const auto bl = blunt_blade(m);
    expect_valid(m, bl.solid);
    const auto g = unwrap(full_round_geometry(m, bl.solid, bl.suction, bl.trailing, bl.pressure));
    ASSERT_FALSE(g.exact);
    ASSERT_LE(g.deviation, 1e-5);
    // the radius: smooth (its second differences small against it), about half the thickness
    std::vector<T> rs;
    for (int i = 0; i <= 32; ++i)
    {
        const T t = g.t1 + (g.t2 - g.t1) * T(i) / 32.;
        const auto c = g.spine->value(t);
        rs.push_back(norm(c - g.contacts[0].curve->value(t)));
        ASSERT_NEAR(norm(c - g.contacts[1].curve->value(t)), rs.back(), 1e-6);
    }
    T kink{};
    for (std::size_t i = 1; i + 1 < rs.size(); ++i)
        kink = std::max(kink, std::abs(rs[i - 1] - 2. * rs[i] + rs[i + 1]));
    MESSAGE("radius from " << rs.front() << " to " << rs.back() << ", largest second difference " << kink);
    ASSERT_GT(std::ranges::min(rs), 0.02);
    ASSERT_LT(std::ranges::max(rs), 0.06);
    ASSERT_LT(kink, 1e-3 * rs.front());
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = unwrap(fillet_full_round(m, bl.solid, std::vector<FaceId>{bl.suction}, std::vector<FaceId>{bl.trailing}, std::vector<FaceId>{bl.pressure}));
    const auto dt = ms(t0);
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), n_faces(m, bl.solid));
    // tangent to the sides along the contacts: the normals of the fillet and of each side there
    T worst{};
    for (std::size_t k = 0; k < 2; ++k)
        for (int i = 1; i < 16; ++i)
        {
            const T t = g.t1 + (g.t2 - g.t1) * T(i) / 16.;
            const auto uv = g.contacts[k].pcurve->value(t), w = g.contacts[k].on_blend->value(t);
            worst = std::max(worst, norm(cross(normal_direction(*g.contacts[k].support, uv[0], uv[1]), normal_direction(*g.surface, w[0], w[1]))));
        }
    const T dv = volume(m, bl.solid) - volume(m, r);
    MESSAGE("blade: " << dt << " ms, deviation " << g.deviation << ", sine of the angle of the normals at the contacts " << worst << ", volume taken off " << dv);
    ASSERT_LE(worst, 1e-5);
    // what is taken off: about two corners of a section of radius the mean one, over the span
    ASSERT_GT(dv, 0.);
    ASSERT_LT(dv, 1.5 * 2. * 0.06 * 0.06);
}

TEST(tests_brep_full_round, refusals)
{
    Model<T> m;
    const T w = 0.8, h = 1., L = 3.;
    const auto s = rib(m, w, h, 0., L);
    const auto top = face_near(m, s, {1., 0., 1. + h}), f1 = face_near(m, s, {1., -w / 2., 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2., 1. + h / 2.});
    const auto plate = face_near(m, s, {1., 1.5, 1.}), under = face_near(m, s, {1., 0., 0.}), end = face_near(m, s, {0., 0., 0.5});
    const auto elsewhere = rib(m, w, h, 0., L); // (another solid of the model)
    const auto before = census(m);
    auto round = [&](std::vector<FaceId> a, std::vector<FaceId> mid, std::vector<FaceId> b, const FilletOptions<T> &o = {}) {
        const auto r = fillet_full_round(m, s, a, mid, b, o);
        // (a refusal leaves the model as it was)
        if (!r)
            for (std::size_t i = 0; i < 3; ++i)
                ASSERT_EQ(census(m)[i], before[i]);
        return r;
    };
    // the faces given
    expect_refused(round({}, {top}, {f2}), BuildErrc::InvalidId, "no face on a side or in the middle");
    expect_refused(round({f1}, {top}, {f1}), BuildErrc::InvalidId, "the same face on both sides, which does not meet the middle face along two edges");
    expect_refused(round({top}, {top}, {f2}), BuildErrc::InvalidId, "a face given both as a side and as the middle face");
    FilletOptions<T> setback;
    setback.setbacks = {Setback<T>{VertexId{}, EdgeId{}, 0.1}};
    expect_refused(round({f1}, {top}, {f2}, setback), BuildErrc::InvalidBounds, "setbacks");
    expect_refused(round({f1}, {top}, {under}), BuildErrc::InvalidId, "the middle face has no edge with a side");
    expect_refused(round({f1, plate}, {top}, {f2}), BuildErrc::NonGeneric, "several faces (not handled yet)");
    expect_refused(round({f1}, {face_near(m, elsewhere, {1., 0., 1. + h})}, {f2}), BuildErrc::InvalidId, "not a face of the solid");
    // the options: the arc of circle only
    FilletOptions<T> conic;
    conic.section = BlendSection<T>::conic(0.5);
    expect_refused(round({f1}, {top}, {f2}, conic), BuildErrc::InvalidBounds, "a conic or G2 section");
    // a middle face convex along a side and concave along the other: the flank between the top and the plate
    expect_refused(round({top}, {f2}, {plate}), BuildErrc::NonGeneric, "convex along a side and concave along the other");
    // the end face as the middle: the plate under it and the top are not both next to it along one edge each
    expect_refused(round({f1}, {end}, {f2}), BuildErrc::NonGeneric, "several faces at an end");
    ASSERT_TRUE(round({f1}, {top}, {f2}).has_value());
}

TEST(tests_brep_full_round, refusals_of_shape)
{
    // a rib lower than half its width: the ball touches the flanks below their foot
    {
        Model<T> m;
        const T w = 0.8, h = 0.3;
        const auto s = rib(m, w, h, 0., 3.);
        const auto before = census(m);
        const auto r = fillet_full_round(m, s, std::vector<FaceId>{face_near(m, s, {1., -w / 2., 1. + h / 2.})}, std::vector<FaceId>{face_near(m, s, {1., 0., 1. + h})},
                                         std::vector<FaceId>{face_near(m, s, {1., w / 2., 1. + h / 2.})});
        expect_refused(r, BuildErrc::NonGeneric, "a contact runs off its face");
        ASSERT_EQ(census(m), before);
    }
    // a top narrowing to a point: its sides meet there
    {
        Model<T> m;
        const T w = 0.5, h = 1.;
        const auto s = tower(m, {{0., -w / 2.}, {2., 0.}, {0., w / 2.}}, h);
        const auto before = census(m);
        const auto r = fillet_full_round(m, s, std::vector<FaceId>{face_near(m, s, {1., -w / 4., h / 2.})}, std::vector<FaceId>{face_near(m, s, {0.5, 0., h})},
                                         std::vector<FaceId>{face_near(m, s, {1., w / 4., h / 2.})});
        expect_refused(r, BuildErrc::NonGeneric, "narrowing to a point");
        ASSERT_EQ(census(m), before);
    }
    // a rib ending in a prow: two faces at that end of its top
    {
        Model<T> m;
        const T w = 0.5, h = 1., L = 2.;
        const auto s = tower(m, {{0., -w / 2.}, {L, -w / 2.}, {L + 0.4, 0.}, {L, w / 2.}, {0., w / 2.}}, h);
        const auto before = census(m);
        const auto r = fillet_full_round(m, s, std::vector<FaceId>{face_near(m, s, {1., -w / 2., h / 2.})}, std::vector<FaceId>{face_near(m, s, {1., 0., h})},
                                         std::vector<FaceId>{face_near(m, s, {1., w / 2., h / 2.})});
        expect_refused(r, BuildErrc::NonGeneric, "several faces at an end");
        ASSERT_EQ(census(m), before);
    }
    // an end at 88 degrees (past the end of the guide by 57 radii): accepted, the fillet
    // continued from the start as far as the other edge runs (refused before the review: the
    // retries of the blends stop at 8 radii)
    {
        Model<T> m;
        const T w = 0.5, h = 1., L = 2., k = std::tan(88. * pi / 180.);
        const auto s = tower(m, {{0., -w / 2.}, {L, -w / 2.}, {L + w * k, w / 2.}, {0., w / 2.}}, h);
        const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{face_near(m, s, {1., -w / 2., h / 2.})}, std::vector<FaceId>{face_near(m, s, {1., 0., h})},
                                                std::vector<FaceId>{face_near(m, s, {1., w / 2., h / 2.})}, tight()));
        expect_valid(m, r);
        ASSERT_NEAR(volume(m, r), volume(m, s) - 2. * corner_area(w / 2.) * (L + w * k / 2.), 1e-9);
    }
}

TEST(tests_brep_full_round, sides_swapped)
{
    // a short rib (0.5) with its end at 60 degrees: the top is a trapezium, its edge with one flank
    // 0.5 long, with the other 1.37; the guide (the edge with faces1) on either: the same solid.
    // Longer ribs whose end is at 80 and 84 degrees (past the end by 11 and 19 radii): the
    // fillet continued from the start by as far as the other edge runs past the guide
    for (const auto &[L, deg] : {std::pair{0.5, 60.}, std::pair{2., 80.}, std::pair{2., 84.}})
    {
        T v[2]{};
        for (const int swap : {0, 1})
        {
            Model<T> m;
            const T w = 0.5, h = 1., k = std::tan(deg * pi / 180.);
            const auto s = tower(m, {{0., -w / 2.}, {L, -w / 2.}, {L + w * k, w / 2.}, {0., w / 2.}}, h);
            auto top = face_near(m, s, {L / 2., 0., h}), f1 = face_near(m, s, {L / 2., -w / 2., h / 2.}), f2 = face_near(m, s, {L / 2., w / 2., h / 2.});
            if (swap)
                std::swap(f1, f2);
            const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, tight()));
            expect_valid(m, r);
            ASSERT_EQ(faces_on<CylindricalSurface<T>>(m, r).size(), 1u);
            v[swap] = volume(m, r);
            const T expected = volume(m, s) - 2. * corner_area(w / 2.) * (L + w * k / 2.);
            MESSAGE("L = " << L << ", end at " << deg << " deg, guide on the " << std::string(swap ? "long" : "short") << " side: volume gap " << v[swap] - expected);
            ASSERT_NEAR(v[swap], expected, 1e-9);
        }
        ASSERT_NEAR(v[0], v[1], 1e-10);
    }
}

TEST(tests_brep_full_round, extension_zero)
{
    // no extension asked: the fillet goes on past the end faces as far as they need (the retries of the blends)
    Model<T> m;
    const T w = 0.6, h = 1., L = 2.;
    const auto s = rib(m, w, h, 0., L);
    const auto top = face_near(m, s, {1., 0., 1. + h}), f1 = face_near(m, s, {1., -w / 2., 1. + h / 2.}), f2 = face_near(m, s, {1., w / 2., 1. + h / 2.});
    auto o = tight();
    o.extension = 0.;
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, o));
    expect_valid(m, r);
    ASSERT_NEAR(volume(m, r), volume(m, s) - L * 2. * corner_area(w / 2.), 1e-9);
}

TEST(tests_brep_full_round, drafted_slot)
{
    // a slot whose walls open by 10 degrees from a bottom 0.5 wide: concave, the cylinder of
    // radius (w / 2) cos a / (1 - sin a), the two corners added, r^2 (tan b - b) each
    Model<T> m;
    const T w = 0.5, a = 10. * pi / 180., d = std::tan(a), L = 2.;
    const auto s = profile(m, {{-1.5, 0.}, {1.5, 0.}, {1.5, 2.}, {w / 2. + d, 2.}, {w / 2., 1.}, {-w / 2., 1.}, {-w / 2. - d, 2.}, {-1.5, 2.}}, L);
    const auto bottom = face_near(m, s, {1., 0., 1.}), f1 = face_near(m, s, {1., -w / 2. - d / 2., 1.5}), f2 = face_near(m, s, {1., w / 2. + d / 2., 1.5});
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{bottom}, std::vector<FaceId>{f2}, tight()));
    expect_valid(m, r);
    const auto cyl = faces_on<CylindricalSurface<T>>(m, r);
    ASSERT_EQ(cyl.size(), 1u);
    const T rr = (w / 2.) * std::cos(a) / (1. - std::sin(a)), b = pi / 4. - a / 2.;
    ASSERT_NEAR(static_cast<const CylindricalSurface<T> &>(*m.face(cyl.front()).surface).radius(), rr, 1e-12);
    const T expected = volume(m, s) + L * 2. * rr * rr * (std::tan(b) - b);
    MESSAGE("drafted slot: volume gap " << volume(m, r) - expected);
    ASSERT_NEAR(volume(m, r), expected, 1e-9);
}

TEST(tests_brep_full_round, curved_rib)
{
    // a rib along a quarter of a circle of radius 2 (its section revolved about z): flanks on
    // coaxial cylinders, the top a plane, one radius w / 2: the quarter of a torus; the
    // volume by Pappus (each corner of area r^2 (1 - pi / 4), its centroid r / (6 (1 - pi / 4))
    // from the centre toward its flank)
    Model<T> m;
    const T w = 0.4, h = 1., R = 2.;
    std::vector<EdgeId> e;
    const std::vector<std::array<T, 2>> rz{{R - w / 2., 0.}, {R + w / 2., 0.}, {R + w / 2., h}, {R - w / 2., h}};
    for (std::size_t i = 0; i < rz.size(); ++i)
    {
        const auto &a = rz[i], &b = rz[(i + 1) % rz.size()];
        e.push_back(unwrap(make_edge(m, point<T, 3>{a[0], 0., a[1]}, point<T, 3>{b[0], 0., b[1]})));
    }
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-20., 20., -20., 20.}), unwrap(make_wire(m, e))));
    const auto s = explore<SolidId>(m, unwrap(make_revolution(m, ShapeId{f}, ax1<T, 3>{point<T, 3>{0., 0., 0.}, Z}, pi / 2.)).shape).front();
    const T c = std::cos(pi / 4.), sn = std::sin(pi / 4.);
    const auto top = face_near(m, s, {R * c, R * sn, h}), f1 = face_near(m, s, {(R - w / 2.) * c, (R - w / 2.) * sn, h / 2.}),
               f2 = face_near(m, s, {(R + w / 2.) * c, (R + w / 2.) * sn, h / 2.});
    // (a torus, exact since PR 4: the centres on a circle, the radius constant, the contacts
    // through the stations; the guide, a curve on a surface from the revolution, is not
    // parametrized by the angle)
    const auto g = unwrap(full_round_geometry(m, s, f1, top, f2, tight()));
    ASSERT_TRUE(g.exact);
    ASSERT_NE(dynamic_cast<const ToroidalSurface<T> *>(g.surface.get()), nullptr);
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}, tight()));
    expect_valid(m, r);
    ASSERT_EQ(faces_on<ToroidalSurface<T>>(m, r).size(), 1u);
    const T rr = w / 2., dc = rr - corner_centroid(rr);
    const T expected = volume(m, s) - (pi / 2.) * corner_area(rr) * ((R - dc) + (R + dc));
    MESSAGE("curved rib: volume gap " << volume(m, r) - expected);
    ASSERT_NEAR(volume(m, r), expected, 1e-9);
}

namespace
{
    /// The U of one face (a B-spline polyline of degree 1: two walls and a floor, w apart, h high) swept along x over L, closed by a flat top: the U, the top.
    /// `knot_on_cut`: the floor in two spans, a knot at u = 1/2 where the face is split.
    std::pair<SolidId, std::array<FaceId, 2>> one_face_u(Model<T> &m, T w, T h, T L, bool knot_on_cut = false)
    {
        const auto u = knot_on_cut ? std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., -w / 2., h}, {0., -w / 2., 0.}, {0., 0., 0.}, {0., w / 2., 0.}, {0., w / 2., h}},
                                                                     std::vector<T>{0., 0., 0.25, 0.5, 0.75, 1., 1.}, 1)
                                   : std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., -w / 2., h}, {0., -w / 2., 0.}, {0., w / 2., 0.}, {0., w / 2., h}},
                                                                     std::vector<T>{0., 0., 1. / 3., 2. / 3., 1., 1.}, 1);
        const auto a = unwrap(make_vertex(m, point<T, 3>{0., -w / 2., h})), b = unwrap(make_vertex(m, point<T, 3>{0., w / 2., h}));
        const std::vector<EdgeId> e{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{u}, 0., 1., a, b)), unwrap(make_edge(m, b, a))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{L, 0., 0.}));
        return {s, {face_near(m, s, {L / 2., -w / 2., h / 2.}), face_near(m, s, {L / 2., 0., h})}};
    }
}

TEST(tests_brep_full_round, one_face_on_both_sides)
{
    // (with a knot of the U on the line where its face is split: the same solid)
    for (const bool knot : {false, true})
    {
    INFO("knot on the cut " << knot);
    // a U of one face (walls and floor of one B-spline of degree 1) closed by a flat top: the top
    // replaced by the full round between the two walls, the same face on both sides: split on a
    // copy between its edges with the top, a piece on each side, merged again. The result has
    // the faces of the solid given (the U, the fillet, the two ends) and the volume of the U
    // with a half cylinder of radius w / 2 (exact: the centres and the contacts on lines along
    // the guide, though the walls are flat pieces of a B-spline)
    Model<T> m;
    const T w = 0.5, h = 1., L = 2.;
    const auto [s, f] = one_face_u(m, w, h, L, knot);
    expect_valid(m, s);
    ASSERT_EQ(n_faces(m, s), 4u);
    const auto before = census(m);
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f[0]}, std::vector<FaceId>{f[1]}, std::vector<FaceId>{f[0]}, tight()));
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), 4u);
    // nothing left over in the model but the result
    const auto after = census(m), mine = std::array<std::size_t, 3>{explore<VertexId>(m, ShapeId{r}).size(), explore<EdgeId>(m, ShapeId{r}).size(), n_faces(m, r)};
    for (std::size_t i = 0; i < 3; ++i)
        ASSERT_EQ(after[i], before[i] + mine[i]);
    const T expected = volume(m, s) - L * 2. * corner_area(w / 2.);
    ASSERT_EQ(faces_on<CylindricalSurface<T>>(m, r).size(), 1u);
    MESSAGE("one face on both sides: volume gap " << volume(m, r) - expected);
    ASSERT_NEAR(volume(m, r), expected, 1e-9);
    }
}

namespace
{
    /// The edges of s at a point p (an end of each).
    std::size_t edges_at(const Model<T> &m, SolidId s, const point<T, 3> &p)
    {
        std::size_t n = 0;
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
            if (const auto &ed = m.edge(e); norm(m.vertex(ed.v1).pnt - p) < 1e-9 || norm(m.vertex(ed.v2).pnt - p) < 1e-9)
                ++n;
        return n;
    }
}

TEST(tests_brep_full_round, one_face_joins_only_its_cut)
{
    // the curve of the U at an end split at a vertex of two edges (a vertex of the solid given,
    // away from the cut): kept in the result; the edges split by the cut joined back
    Model<T> m;
    const T w = 0.5, h = 1., L = 2.;
    const auto [s, f] = one_face_u(m, w, h, L);
    EdgeId end_curve{};
    for (auto e : explore<EdgeId>(m, ShapeId{s}))
        if (const auto &ed = m.edge(e); !ed.degenerate)
            if (const auto q = ed.curve->value((ed.u1 + ed.u2) / 2.); std::abs(q[0]) < 1e-12 && q[2] < 0.5 * h)
                end_curve = e;
    ASSERT_TRUE(end_curve.valid());
    const auto &ec = m.edge(end_curve);
    const T at = ec.u1 + 0.4 * (ec.u2 - ec.u1);
    const auto p = ec.curve->value(at);
    (void)unwrap(split_edge(m, end_curve, at));
    expect_valid(m, s);
    ASSERT_EQ(edges_at(m, s, p), 2u);
    const auto ne = explore<EdgeId>(m, ShapeId{s}).size();
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f[0]}, std::vector<FaceId>{f[1]}, std::vector<FaceId>{f[0]}, tight()));
    expect_valid(m, r);
    // the vertex given kept; the same edges as the solid given (the top and its four edges gone,
    // the fillet and its two contacts and two ends come), none left of the cut
    ASSERT_EQ(edges_at(m, r, p), 2u);
    ASSERT_EQ(explore<EdgeId>(m, ShapeId{r}).size(), ne);
}

TEST(tests_brep_full_round, one_face_cut_obliquely)
{
    // a smooth U of one B-spline (degree 3) cut obliquely at an end by a boolean: its edges with
    // the top stay isoparametric, but the face starts at v = 0.21 in the middle of the U where
    // its box starts at 0.175: the cut halfway runs from contour to contour, not over the box
    Model<T> m;
    const T w = 0.5, h = 1., L = 2.;
    const auto u = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., -w / 2., h}, {0., -w / 2., h - 0.5}, {0., -w / 2., 0.}, {0., w / 2., 0.}, {0., w / 2., h - 0.5}, {0., w / 2., h}},
                                                   std::vector<T>{0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.}, 3);
    const auto a = unwrap(make_vertex(m, point<T, 3>{0., -w / 2., h})), b = unwrap(make_vertex(m, point<T, 3>{0., w / 2., h}));
    const std::vector<EdgeId> e{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{u}, 0., 1., a, b)), unwrap(make_edge(m, b, a))};
    const auto pf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    const auto s0 = solid_of(m, make_prism(m, pf, point<T, 3>{L, 0., 0.}));
    std::vector<EdgeId> k;
    const std::vector<std::array<T, 2>> xz{{-1., -1.}, {0.5, -1.}, {0.2, 3.}, {-1., 3.}};
    for (std::size_t i = 0; i < xz.size(); ++i)
        k.push_back(unwrap(make_edge(m, point<T, 3>{xz[i][0], -2., xz[i][1]}, point<T, 3>{xz[(i + 1) % 4][0], -2., xz[(i + 1) % 4][1]})));
    const auto kf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., -2., 0.}, Y, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, k))));
    const auto knife = solid_of(m, make_prism(m, kf, point<T, 3>{0., 4., 0.}));
    const auto s = unwrap(cut(m, ShapeId{s0}, ShapeId{knife})).solids.front();
    expect_valid(m, s);
    auto p = u->value(0.05);
    p[0] = 1.;
    const auto side = face_near(m, s, p), top = face_near(m, s, {1., 0., h});
    ASSERT_TRUE(side.valid() && top.valid());
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{side}, std::vector<FaceId>{top}, std::vector<FaceId>{side}));
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), n_faces(m, s));
    ASSERT_GT(volume(m, s), volume(m, r));
}


TEST(tests_brep_full_round, refusals_of_ends)
{
    auto sides_of = [](Model<T> &m, SolidId s, T w, T h, T x) {
        return std::array<FaceId, 3>{face_near(m, s, {x, -w / 2., h / 2.}), face_near(m, s, {x, 0., h}), face_near(m, s, {x, w / 2., h / 2.})};
    };
    auto expect = [](Model<T> &m, SolidId s, const std::array<FaceId, 3> &f, const std::string &why) {
        for (const auto &g : f)
            ASSERT_TRUE(g.valid());
        const auto before = census(m);
        const auto r = fillet_full_round(m, s, std::vector<FaceId>{f[0]}, std::vector<FaceId>{f[1]}, std::vector<FaceId>{f[2]});
        expect_refused(r, BuildErrc::NonGeneric, why);
        ASSERT_EQ(census(m), before);
    };
    const T w = 0.5, h = 1., L = 2.;
    // a rib with a rounded end (a half cylinder tangent to both flanks): a quarter of a sphere, not built
    {
        Model<T> m;
        const auto v0 = unwrap(make_vertex(m, point<T, 3>{0., -w / 2., 0.})), v1 = unwrap(make_vertex(m, point<T, 3>{L, -w / 2., 0.})),
                   v2 = unwrap(make_vertex(m, point<T, 3>{L, w / 2., 0.})), v3 = unwrap(make_vertex(m, point<T, 3>{0., w / 2., 0.}));
        const auto arc = std::make_shared<EllipseCurve<T>>(w / 2., w / 2., ax2<T, 3>{{{L, 0., 0.}, Z, X}}, -pi / 2., pi / 2.);
        const std::vector<EdgeId> e{unwrap(make_edge(m, v0, v1)), unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{arc}, -pi / 2., pi / 2., v1, v2)),
                                    unwrap(make_edge(m, v2, v3)), unwrap(make_edge(m, v3, v0))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
        expect(m, s, sides_of(m, s, w, h, 1.), "an end face tangent to a side (the rounded end of a rib)");
    }
    // a top narrowing to 0.005: refused by PR 3 (the fillet did not meet the end face), built
    // since PR 4 (its march divided where the ball moves fast: a cone exact to the end)
    {
        Model<T> m;
        const T w0 = 0.6, w1 = 0.005;
        const auto s = tower(m, {{0., -w0 / 2.}, {L, -w1 / 2.}, {L, w1 / 2.}, {0., w0 / 2.}}, h);
        const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{face_near(m, s, {1., -(w0 + w1) / 4., h / 2.})}, std::vector<FaceId>{face_near(m, s, {1., 0., h})},
                                                std::vector<FaceId>{face_near(m, s, {1., (w0 + w1) / 4., h / 2.})}));
        expect_valid(m, r);
        ASSERT_EQ(n_faces(m, r), n_faces(m, s));
        ASSERT_EQ(faces_on<ConicalSurface<T>>(m, r).size(), 1u);
    }
    // a top with a hole through the rib
    {
        Model<T> m;
        std::vector<EdgeId> outer, inner;
        const std::vector<std::array<T, 2>> o{{0., -w / 2.}, {L, -w / 2.}, {L, w / 2.}, {0., w / 2.}}, i{{0.9, -0.1}, {0.9, 0.1}, {1.1, 0.1}, {1.1, -0.1}};
        for (std::size_t k = 0; k < 4; ++k)
        {
            outer.push_back(unwrap(make_edge(m, point<T, 3>{o[k][0], o[k][1], 0.}, point<T, 3>{o[(k + 1) % 4][0], o[(k + 1) % 4][1], 0.})));
            inner.push_back(unwrap(make_edge(m, point<T, 3>{i[k][0], i[k][1], 0.}, point<T, 3>{i[(k + 1) % 4][0], i[(k + 1) % 4][1], 0.})));
        }
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, outer)),
                                        std::vector<WireId>{unwrap(make_wire(m, inner))}));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
        expect(m, s, sides_of(m, s, w, h, 0.5), "a middle face with holes");
    }
    // a notch into the top from a flank at mid-length: the top meets that flank along two edges
    {
        Model<T> m;
        const auto rib0 = tower(m, {{0., -w / 2.}, {L, -w / 2.}, {L, w / 2.}, {0., w / 2.}}, h);
        // (from half the height up: the flank stays one face)
        std::vector<EdgeId> e;
        const std::vector<std::array<T, 2>> xy{{0.9, -w}, {1.1, -w}, {1.1, 0.}, {0.9, 0.}};
        for (std::size_t k = 0; k < 4; ++k)
            e.push_back(unwrap(make_edge(m, point<T, 3>{xy[k][0], xy[k][1], h / 2.}, point<T, 3>{xy[(k + 1) % 4][0], xy[(k + 1) % 4][1], h / 2.})));
        const auto nf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., h / 2.}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        const auto notch = solid_of(m, make_prism(m, nf, point<T, 3>{0., 0., h}));
        const auto s = unwrap(cut(m, ShapeId{rib0}, ShapeId{notch})).solids.front();
        expect(m, s, sides_of(m, s, w, h, 0.5), "meets a side along several edges");
    }
    // a top on a patch whose end at the tip of a wedge is a degenerate edge (refused as a top
    // narrowing to a point: an end of degenerate edges only)
    {
        Model<T> m;
        const point<T, 3> A{0., -w / 2., h}, B{0., w / 2., h}, C{L, 0., h}, A0{0., -w / 2., 0.}, B0{0., w / 2., 0.}, C0{L, 0., 0.};
        const auto top = std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{A, B, C, C}, std::vector<T>{0., 0., 1., 1.}, std::vector<T>{0., 0., 1., 1.}, 1, 1);
        auto quad = [&](const std::vector<point<T, 3>> &p) {
            std::vector<EdgeId> e;
            for (std::size_t k = 0; k < p.size(); ++k)
                e.push_back(unwrap(make_edge(m, p[k], p[(k + 1) % p.size()])));
            const auto n = cross(p[1] - p[0], p[2] - p[0]), x = p[1] - p[0];
            return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n / norm(n), x / norm(x)}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        };
        const std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{top})), quad({A0, C0, B0}), quad({A0, A, C, C0}), quad({B0, C0, C, B}), quad({A0, B0, B, A})};
        SewOptions<T> so;
        so.tol = 1e-9;
        const auto sewn = unwrap(sew(m, faces, so));
        const auto s = unwrap(make_solid(m, sewn.shells.front()));
        expect_valid(m, s);
        expect(m, s, {face_near(m, s, {0.5, -w / 2. * 0.75, h / 2.}), face_near(m, s, {0.5, 0., h}), face_near(m, s, {0.5, w / 2. * 0.75, h / 2.})}, "narrowing to a point");
    }
}

namespace
{
    /// The blade of `blunt_blade` with its two sides in one patch (u from the trailing edge of the pressure side round the sharp leading edge, at u = 1/2, to that of the suction side): the solid, the sides, the face of the trailing edge.
    /// `smooth`: the leading edge rounded (the knots uniform, the patch C2 there), for a blade fused to a hub (a boolean marches poorly along a crease).
    std::tuple<SolidId, FaceId, FaceId> blunt_blade_one_patch(Model<T> &m, bool smooth = false)
    {
        constexpr std::size_t nv = 4;
        const std::array<T, nv> z{0., 0.5, 1., 1.5}, e{0.04, 0.035, 0.03, 0.022}, xte{1., 1.04, 1.1, 1.18}, twist{0., 0.06, 0.12, 0.2};
        const std::array<std::array<T, 2>, 4> up{{{0., 0.}, {0.05, 0.08}, {0.4, 0.16}, {0.8, 0.1}}}, down{{{0., 0.}, {0.05, -0.05}, {0.4, -0.04}, {0.8, -0.03}}};
        auto pole = [&](const std::array<T, 2> &p, std::size_t j) {
            const T c = std::cos(twist[j]), s = std::sin(twist[j]);
            return point<T, 3>{c * p[0] - s * p[1], s * p[0] + c * p[1], z[j]};
        };
        // a row: the pressure side reversed (from its trailing edge), then the suction side past the leading edge
        auto row = [&](std::size_t j) {
            points_vector<T, 3> r{pole({xte[j], -e[j]}, j)};
            for (std::size_t i = 4; i-- > 1;)
                r.push_back(pole(down[i], j));
            for (std::size_t i = 0; i < 4; ++i)
                r.push_back(pole(up[i], j));
            r.push_back(pole({xte[j], e[j]}, j));
            return r;
        };
        points_vector<T, 3> ps, pt;
        for (std::size_t j = 0; j < nv; ++j)
        {
            const auto r = row(j);
            ps.insert(ps.end(), r.begin(), r.end());
            pt.push_back(r.front()), pt.push_back(r.back());
        }
        const std::vector<T> ku = smooth ? std::vector<T>{0., 0., 0., 0., 1. / 6., 2. / 6., 3. / 6., 4. / 6., 5. / 6., 1., 1., 1., 1.}
                                         : std::vector<T>{0., 0., 0., 0., 0.25, 0.5, 0.5, 0.5, 0.75, 1., 1., 1., 1.};
        const std::vector<T> kv{0., 0., 0., 0., 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const auto ss = std::make_shared<BSSurface<T, 3>>(ps, ku, kv, std::size_t{3}, std::size_t{3});
        const auto st = std::make_shared<BSSurface<T, 3>>(pt, k1, kv, std::size_t{1}, std::size_t{3});
        std::vector<FaceId> faces{unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{ss})), unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{st}))};
        for (const std::size_t j : {std::size_t{0}, nv - 1})
        {
            const auto c = std::make_shared<BSCurve<T, 3>>(row(j), ku, 3);
            const auto v0 = unwrap(make_vertex(m, c->value(0.))), v1 = unwrap(make_vertex(m, c->value(1.)));
            const std::vector<EdgeId> es{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{c}, 0., 1., v0, v1)), unwrap(make_edge(m, v1, v0))};
            const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., z[j]}, Z, X}}, std::array<T, 4>{-5., 5., -5., 5.});
            faces.push_back(unwrap(make_face(m, plane, unwrap(make_wire(m, es)))));
        }
        const auto sewn = unwrap(sew(m, faces));
        const auto s = unwrap(make_solid(m, sewn.shells.front()));
        return {s, face_near(m, s, ss->value(0.7, 0.5)), face_near(m, s, st->value(0.5, 0.5))};
    }
}

TEST(tests_brep_full_round, blade_of_one_patch)
{
    // the blade of blade_trailing_edge with its sides in one patch (as a loft of open sections):
    // the same face on both sides of the trailing edge; the result has the faces of the blade
    // given and the volume of the blade of three patches rounded the same way
    Model<T> m;
    const auto [s, sides, te] = blunt_blade_one_patch(m);
    expect_valid(m, s);
    ASSERT_EQ(n_faces(m, s), 4u);
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{sides}, std::vector<FaceId>{te}, std::vector<FaceId>{sides}));
    const auto dt = ms(t0);
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), 4u);
    const auto bl = blunt_blade(m);
    const auto r3 = unwrap(fillet_full_round(m, bl.solid, std::vector<FaceId>{bl.suction}, std::vector<FaceId>{bl.trailing}, std::vector<FaceId>{bl.pressure}));
    MESSAGE("blade of one patch: " << dt << " ms; volume " << volume(m, r) << ", of three patches " << volume(m, r3) << ", gap " << volume(m, r) - volume(m, r3));
    ASSERT_NEAR(volume(m, s), volume(m, bl.solid), 1e-12);
    ASSERT_NEAR(volume(m, r), volume(m, r3), 1e-8);
}

TEST(tests_brep_full_round, blade_of_one_patch_on_a_hub)
{
    // the blade of one patch fused to a hub whose top is a plane rising along the chord
    // (z = 0.15 + 0.1 x): its root on the hub is not an isoparametric line of the patch; the
    // trailing edge rounded down to the hub
    Model<T> m;
    const auto [b, sides0, te0] = blunt_blade_one_patch(m, true);
    std::vector<EdgeId> k;
    const std::vector<std::array<T, 2>> xz{{-0.5, -0.5}, {1.6, -0.5}, {1.6, 0.15 + 0.16}, {-0.5, 0.15 - 0.05}};
    for (std::size_t i = 0; i < xz.size(); ++i)
        k.push_back(unwrap(make_edge(m, point<T, 3>{xz[i][0], -0.6, xz[i][1]}, point<T, 3>{xz[(i + 1) % 4][0], -0.6, xz[(i + 1) % 4][1]})));
    const auto hf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., -0.6, 0.}, Y, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, k))));
    const auto hub = solid_of(m, make_prism(m, hf, point<T, 3>{0., 1.2, 0.}));
    const auto s = unwrap(fuse(m, ShapeId{hub}, ShapeId{b})).solids.front();
    expect_valid(m, s);
    const auto sides = face_near(m, s, m.face(sides0).surface->value(0.7, 0.5)), te = face_near(m, s, m.face(te0).surface->value(0.5, 0.5));
    ASSERT_TRUE(sides.valid() && te.valid());
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = unwrap(fillet_full_round(m, s, std::vector<FaceId>{sides}, std::vector<FaceId>{te}, std::vector<FaceId>{sides}));
    MESSAGE("blade of one patch on a hub: " << ms(t0) << " ms, volume taken off " << volume(m, s) - volume(m, r));
    expect_valid(m, r);
    ASSERT_EQ(n_faces(m, r), n_faces(m, s));
    ASSERT_GT(volume(m, s), volume(m, r));
}

TEST(tests_brep_full_round, one_face_refusals)
{
    // the U of one face with a pocket cut into its floor where the face is to be split: the line
    // halfway crosses the contour of the face four times (its outer wire, the pocket): refused,
    // saying to split the face; the model as it was
    Model<T> m;
    const T w = 0.5, h = 1., L = 2.;
    const auto [s0, f0] = one_face_u(m, w, h, L);
    std::vector<EdgeId> k;
    const std::vector<std::array<T, 2>> xy{{0.9, -0.05}, {1.1, -0.05}, {1.1, 0.05}, {0.9, 0.05}};
    for (std::size_t i = 0; i < xy.size(); ++i)
        k.push_back(unwrap(make_edge(m, point<T, 3>{xy[i][0], xy[i][1], -0.5}, point<T, 3>{xy[(i + 1) % 4][0], xy[(i + 1) % 4][1], -0.5})));
    const auto kf = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., -0.5}, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, k))));
    const auto pocket = solid_of(m, make_prism(m, kf, point<T, 3>{0., 0., 0.6}));
    const auto s = unwrap(cut(m, ShapeId{s0}, ShapeId{pocket})).solids.front();
    expect_valid(m, s);
    const auto side = face_near(m, s, {L / 2., -w / 2., h / 2.}), top = face_near(m, s, {L / 2., 0., h});
    ASSERT_TRUE(side.valid() && top.valid());
    const auto before = census(m);
    expect_refused(fillet_full_round(m, s, std::vector<FaceId>{side}, std::vector<FaceId>{top}, std::vector<FaceId>{side}), BuildErrc::NonGeneric,
                   "crosses its contour 4 times");
    ASSERT_EQ(census(m), before);
}

TEST(tests_brep_full_round, one_face_nearly_isoparametric)
{
    // the pcurve of an edge of the U with the top moved by a bump of d in u at its middle (the
    // solid still valid): within tol in space (d = 1e-8, some 3e-8 in space) the edge is taken
    // as isoparametric, the same solid as without the bump; far beyond tol (d = 1e-4) refused
    T v0{};
    for (const T d : {0., 1e-8, 1e-4})
    {
        INFO("d = " << d);
        Model<T> m;
        const T w = 0.5, h = 1., L = 2.;
        const auto [s, f] = one_face_u(m, w, h, L);
        if (d > 0.)
        {
            const blend_detail::SolidView<T> sv{m, s};
            bool done = false;
            for (auto wid : m.face(f[0]).wires)
                for (auto &ce : m.wire(wid).coedges)
                    if (!done && std::ranges::contains(sv.users(ce.edge), f[1]))
                    {
                        const auto &ed = m.edge(ce.edge);
                        const auto a = ce.pcurve->value(ed.u1), b = ce.pcurve->value(ed.u2), mid = (a + b) / 2.;
                        // (into the face: from u = 0 toward +u)
                        const T dir = a[0] < 0.5 ? 1. : -1.;
                        ce.pcurve = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{a, {mid[0] + dir * 2. * d, mid[1]}, b}, std::vector<T>{ed.u1, ed.u1, ed.u1, ed.u2, ed.u2, ed.u2}, 2);
                        done = true;
                    }
            ASSERT_TRUE(done);
        }
        const auto r = fillet_full_round(m, s, std::vector<FaceId>{f[0]}, std::vector<FaceId>{f[1]}, std::vector<FaceId>{f[0]});
        if (d > 1e-6)
        {
            expect_refused(r, BuildErrc::NonGeneric, "not isoparametric in one direction");
            continue;
        }
        ASSERT_TRUE(r.has_value());
        expect_valid(m, *r);
        if (d == 0.)
            v0 = volume(m, *r);
        else
            ASSERT_NEAR(volume(m, *r), v0, 1e-9);
    }
}

TEST(tests_brep_full_round, against_fillet_faces)
{
    // the rib with a draft of +10 degrees: fillet_faces of its two flanks over the top at r is the
    // fillet of the virtual edge where the flanks meet; it swallows the top only past the radius
    // of the full round (refused as spilling over at that radius and below), and its volume then
    // tends to that of the full round, linearly in r - r_fr (-2.19 r_fr per unit of r / r_fr - 1)
    const T w = 0.6, h = 1., L = 2., a = 10. * pi / 180., d = h * std::tan(a) / 2.;
    const T rr = (w / 2.) * std::cos(a) / (1. - std::sin(a));
    auto flanks = [&](Model<T> &m, SolidId s) {
        return std::pair{face_near(m, s, {1., -w / 2. - d, 1. + h / 2.}), face_near(m, s, {1., w / 2. + d, 1. + h / 2.})};
    };
    T vfr{};
    {
        Model<T> m;
        const auto s = rib(m, w, h, a, L);
        const auto [f1, f2] = flanks(m, s);
        vfr = volume(m, unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{face_near(m, s, {1., 0., 1. + h})}, std::vector<FaceId>{f2})));
    }
    for (const T k : {1. - 1e-3, 1.})
    {
        Model<T> m;
        const auto s = rib(m, w, h, a, L);
        const auto [f1, f2] = flanks(m, s);
        const auto before = census(m);
        expect_refused(fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, k * rr), BuildErrc::NonGeneric, "spills over");
        ASSERT_EQ(census(m), before);
    }
    for (const T e : {1e-6, 1e-3})
    {
        Model<T> m;
        const auto s = rib(m, w, h, a, L);
        const auto [f1, f2] = flanks(m, s);
        const auto r = unwrap(fillet_faces(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{f2}, (1. + e) * rr));
        expect_valid(m, r);
        ASSERT_EQ(n_faces(m, r), n_faces(m, s));
        const T slope = (volume(m, r) - vfr) / e;
        MESSAGE("fillet_faces at r (1 + " << e << "): volume - full round = " << volume(m, r) - vfr << ", per unit of e " << slope);
        ASSERT_LT(volume(m, r), vfr);
        ASSERT_NEAR(slope, -2.187, 0.01);
    }
}


namespace
{
    /// The sector of a disc of radius R, opening `angle` from the x axis, thickness t (its rim a part of a cylinder, its ends planes through the axis).
    SolidId disc_sector(Model<T> &m, T R, T angle, T t)
    {
        const point<T, 3> o{0., 0., 0.}, a{R, 0., 0.}, b{R * std::cos(angle), R * std::sin(angle), 0.};
        const auto vo = unwrap(make_vertex(m, o)), va = unwrap(make_vertex(m, a)), vb = unwrap(make_vertex(m, b));
        const auto arc = std::make_shared<EllipseCurve<T>>(R, R, ax2<T, 3>{{o, Z, X}}, 0., angle);
        const std::vector<EdgeId> e{unwrap(make_edge(m, vo, va)), unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{arc}, 0., angle, va, vb)), unwrap(make_edge(m, vb, vo))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{o, Z, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        return solid_of(m, make_prism(m, f, point<T, 3>{0., 0., t}));
    }
}

TEST(tests_brep_full_round, curved_open_guide)
{
    // the rim of a sector of a disc (a third of a turn, R = 1, t = 0.3): an open guide whose edge
    // with the top is an arc of circle on a plane, its pcurve not straight in (u, v): the guide
    // put into the top in its place (PR 3 put the chord). The fillet is a part of a torus of
    // radii R - t / 2 and t / 2; the volume by Pappus over the angle
    for (const T angle : {5. * pi / 3., 2. * pi / 3.})
    {
    INFO("angle " << angle);
    Model<T> m;
    const T R = 1., t = 0.3, r = t / 2.;
    const auto s = disc_sector(m, R, angle, t);
    const T c = std::cos(angle / 2.), sn = std::sin(angle / 2.);
    FaceId rim{};
    for (auto g : explore<FaceId>(m, ShapeId{s}))
        if (dynamic_cast<const CylindricalSurface<T> *>(m.face(g).surface.get()))
            rim = g;
    const auto bottom = face_near(m, s, {0.5 * c, 0.5 * sn, 0.}), top = face_near(m, s, {0.5 * c, 0.5 * sn, t});
    ASSERT_TRUE(rim.valid() && bottom.valid() && top.valid());
    const auto f = fillet_full_round(m, s, std::vector<FaceId>{bottom}, std::vector<FaceId>{rim}, std::vector<FaceId>{top}, tight());
    if (!f)
        MESSAGE(f.error().message);
    ASSERT_TRUE(f.has_value());
    expect_valid(m, *f);
    ASSERT_EQ(n_faces(m, *f), n_faces(m, s));
    // (the volume of the sector given exactly: its rim on a cylinder integrates within 5e-8)
    const T expected = angle * R * R / 2. * t - angle * 2. * corner_area(r) * (R - r + corner_centroid(r));
    MESSAGE("sector of a disc: volume gap " << volume(m, *f) - expected);
    ASSERT_NEAR(volume(m, *f), expected, 1e-9);
    }
}


// ---- PR 4: a middle face round a loop, exact surfaces from the stations, the march, folding

namespace
{
    /// The volume taken off (or added) by a full round round a loop of radius rc of the centres: two corners of area r^2 (1 - pi / 4) revolved, at rc + e1 dc and rc + e2 dc (e = +-1, the side of each corner).
    T pappus(T rc, T r, T e1, T e2)
    {
        return 2. * pi * corner_area(r) * ((rc + e1 * corner_centroid(r)) + (rc + e2 * corner_centroid(r)));
    }

    /// The solid of revolution about z of the polygon (radius, z) given (counter-clockwise in the plane y = 0 seen from -y), over a whole turn.
    SolidId revolved(Model<T> &m, const std::vector<std::array<T, 2>> &rz)
    {
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < rz.size(); ++i)
        {
            const auto &a = rz[i], &b = rz[(i + 1) % rz.size()];
            e.push_back(unwrap(make_edge(m, point<T, 3>{a[0], 0., a[1]}, point<T, 3>{b[0], 0., b[1]})));
        }
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-20., 20., -20., 20.}), unwrap(make_wire(m, e))));
        return explore<SolidId>(m, unwrap(make_revolution(m, ShapeId{f}, ax1<T, 3>{point<T, 3>{0., 0., 0.}, Z}, 2. * pi)).shape).front();
    }
}

TEST(tests_brep_full_round, disc_rim)
{
    // the rim of a disc (R = 1, thickness 0.3) between its two faces: the middle face round a
    // loop, no end face; the half of a torus of radii R - t / 2 and t / 2, exact; Pappus
    Model<T> m;
    const T R = 1., t = 0.3, r = t / 2.;
    const ax2<T, 3> ax{{{0., 0., 0.}, Z, X}};
    const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(R, R, ax)));
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2., 2., -2., 2.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
    const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., t}));
    const auto rim = faces_on<CylindricalSurface<T>>(m, s).front(), bottom = face_near(m, s, {0., 0., 0.}), top = face_near(m, s, {0., 0., t});
    const auto t0 = std::chrono::steady_clock::now();
    const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{bottom}, std::vector<FaceId>{rim}, std::vector<FaceId>{top}, tight()));
    const auto dt = ms(t0);
    expect_valid(m, x);
    ASSERT_EQ(n_faces(m, x), n_faces(m, s));
    const auto tor = faces_on<ToroidalSurface<T>>(m, x);
    ASSERT_EQ(tor.size(), 1u);
    const auto &ts = static_cast<const ToroidalSurface<T> &>(*m.face(tor.front()).surface);
    ASSERT_NEAR(ts.major_radius(), R - r, 1e-12);
    ASSERT_NEAR(ts.minor_radius(), r, 1e-12);
    // (the volume of the disc given exactly: its rim on a cylinder integrates within 5e-8)
    const T expected = pi * R * R * t - pappus(R - r, r, 1., 1.);
    MESSAGE("disc rim: volume gap " << volume(m, x) - expected << ", " << dt << " ms");
    ASSERT_NEAR(volume(m, x), expected, 1e-9);
}

TEST(tests_brep_full_round, tube_end)
{
    // the end of a tube (radii 0.7 and 1, a ring between its two cylinders): the middle face a
    // ring of two wires; the half of a torus of radii 0.85 and 0.15
    Model<T> m;
    const auto s = revolved(m, {{0.7, 0.}, {1., 0.}, {1., 1.}, {0.7, 1.}});
    FaceId top = face_near(m, s, {0., 0.85, 1.}), inner{}, outer{};
    for (auto g : faces_on<CylindricalSurface<T>>(m, s))
        (static_cast<const CylindricalSurface<T> &>(*m.face(g).surface).radius() < 0.8 ? inner : outer) = g;
    ASSERT_TRUE(top.valid() && inner.valid() && outer.valid());
    const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{inner}, std::vector<FaceId>{top}, std::vector<FaceId>{outer}, tight()));
    expect_valid(m, x);
    ASSERT_EQ(n_faces(m, x), n_faces(m, s));
    ASSERT_EQ(faces_on<ToroidalSurface<T>>(m, x).size(), 1u);
    const T expected = pi * (1. - 0.49) * 1. - pappus(0.85, 0.15, 1., -1.);
    MESSAGE("tube end: volume gap " << volume(m, x) - expected);
    ASSERT_NEAR(volume(m, x), expected, 1e-9);
}

TEST(tests_brep_full_round, shaft_groove)
{
    // the bottom of a groove round a shaft (radius 1, the groove 0.4 wide, 0.3 deep): concave,
    // round a loop; the half of a torus of radii 0.9 and 0.2 adds the two corners at the bottom
    Model<T> m;
    const auto s = revolved(m, {{0., 0.}, {1., 0.}, {1., 0.8}, {0.7, 0.8}, {0.7, 1.2}, {1., 1.2}, {1., 2.}, {0., 2.}});
    const auto bottom = face_near(m, s, {0., 0.7, 1.}), w1 = face_near(m, s, {0., 0.85, 0.8}), w2 = face_near(m, s, {0., 0.85, 1.2});
    ASSERT_TRUE(bottom.valid() && w1.valid() && w2.valid());
    const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{w1}, std::vector<FaceId>{bottom}, std::vector<FaceId>{w2}, tight()));
    expect_valid(m, x);
    ASSERT_EQ(n_faces(m, x), n_faces(m, s));
    ASSERT_EQ(faces_on<ToroidalSurface<T>>(m, x).size(), 1u);
    const T expected = pi * 2. - pi * (1. - 0.49) * 0.4 + pappus(0.9, 0.2, -1., -1.);
    MESSAGE("shaft groove: volume gap " << volume(m, x) - expected);
    ASSERT_NEAR(volume(m, x), expected, 1e-9);
}

namespace
{
    /// The path of a hairpin at the distance rho from the centre (L, 0) of its turn: y = -rho from x = 0 to L, a half circle, y = +rho back to 0 (a cubic interpolation of 401 points, at the height H).
    std::shared_ptr<BSCurve<T, 3>> hairpin(T rho, T L, T H)
    {
        points_vector<T, 3> p;
        const int n = 400;
        const T total = 2. * L + pi * rho;
        for (int i = 0; i <= n; ++i)
        {
            const T s = total * i / n;
            if (s <= L)
                p.push_back({s, -rho, H});
            else if (s <= L + pi * rho)
            {
                const T a = -pi / 2. + (s - L) / rho;
                p.push_back({L + rho * std::cos(a), rho * std::sin(a), H});
            }
            else
                p.push_back({L - (s - L - pi * rho), rho, H});
        }
        return std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(p, 3, KnotsCalcMode::CHORD_LENGTH));
    }
}

TEST(tests_brep_full_round, hairpin_wall)
{
    // the top of a thin wall (0.004) bent into a hairpin of inner radius a (the probe of the
    // review of the design of PR 4): the plane of a station across one branch near the turn
    // meets the other branch too, where another ball touches the three faces. PR 3 stepped over a
    // turn shorter than its step (a = 0.004: refused, "approx_tol"); the divided march of PR 4
    // does not. The volume taken off: two corners along the middle line (Pappus), within the
    // skin (approx_tol)
    for (const T a : {0.2, 0.02, 0.004})
    {
        INFO("a " << a);
        Model<T> m;
        const T w = 0.004, h = 0.1, L = 1.;
        const auto ci = hairpin(a, L, 0.), co = hairpin(a + w, L, 0.);
        const auto [i0, i1] = ci->bounds();
        const auto [o0, o1] = co->bounds();
        const auto vi0 = unwrap(make_vertex(m, ci->value(i0))), vi1 = unwrap(make_vertex(m, ci->value(i1))), vo0 = unwrap(make_vertex(m, co->value(o0))),
                   vo1 = unwrap(make_vertex(m, co->value(o1)));
        const std::vector<EdgeId> e{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{ci}, i0, i1, vi0, vi1)), unwrap(make_edge(m, vi1, vo1)),
                                    unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{co}, o0, o1, vo0, vo1)), unwrap(make_edge(m, vo0, vi0))};
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-5., 5., -5., 5.}), unwrap(make_wire(m, e))));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., h}));
        const auto top = face_near(m, s, {0.5, -a - w / 2., h}), f1 = face_near(m, s, {0.5, -a, h / 2.}), f2 = face_near(m, s, {0.5, -a - w, h / 2.});
        const auto t0 = std::chrono::steady_clock::now();
        const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}));
        const auto dt = ms(t0);
        expect_valid(m, x);
        const T rr = w / 2., mid = 2. * L + pi * (a + w / 2.);
        const T removed = volume(m, s) - volume(m, x), expected = 2. * corner_area(rr) * mid;
        MESSAGE("hairpin a = " << a << ": taken off " << removed << " for " << expected << " (relative gap " << removed / expected - 1. << "), " << dt << " ms");
        ASSERT_NEAR(removed / expected, 1., 1e-3);
    }
}

namespace
{
    /// The fin of thickness t: the closed profile of blade_fixture (a closed B-spline, not periodic, its seam at the trailing edge) extruded along z.
    std::tuple<SolidId, FaceId, FaceId, FaceId> fin(Model<T> &m, T t)
    {
        const auto c = blade_fixture::section(blade_fixture::foil(), 0., 0.);
        const auto [u1, u2] = c->bounds();
        const auto v = unwrap(make_vertex(m, c->value(u1)));
        const auto e = unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{c}, u1, u2, v, v));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, Z, X}}, std::array<T, 4>{-5., 5., -5., 5.}), unwrap(make_wire(m, std::vector<EdgeId>{e}))));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{0., 0., t}));
        FaceId rim{};
        for (auto g : explore<FaceId>(m, ShapeId{s}))
            if (!dynamic_cast<const Plane<T> *>(m.face(g).surface.get()))
                rim = g;
        return {s, face_near(m, s, {0.5, 0., 0.}), rim, face_near(m, s, {0.5, 0., t})};
    }
}

TEST(tests_brep_full_round, fin_rim)
{
    // the rim of a fin (the profile of blade_fixture extruded, thickness 0.004: r = 0.002, under
    // the radius of curvature of its trailing edge, 0.0035): a loop round a closed B-spline seen
    // past its seam, the march across the thin noses of its edges; a skin within approx_tol
    Model<T> m;
    const auto [s, bottom, rim, top] = fin(m, 0.004);
    ASSERT_TRUE(bottom.valid() && rim.valid() && top.valid());
    const auto t0 = std::chrono::steady_clock::now();
    const auto g = unwrap(full_round_geometry(m, s, bottom, rim, top));
    ASSERT_LE(g.deviation, 1e-5);
    const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{bottom}, std::vector<FaceId>{rim}, std::vector<FaceId>{top}));
    const auto dt = ms(t0);
    expect_valid(m, x);
    ASSERT_EQ(n_faces(m, x), n_faces(m, s));
    // the radius all round: half the thickness
    for (int i = 0; i <= 16; ++i)
    {
        const T t = g.a + (g.b - g.a) * T(i) / 16.;
        ASSERT_NEAR(norm(g.spine->value(t) - g.contacts[0].curve->value(t)), 0.002, 1e-6);
    }
    MESSAGE("fin rim: " << dt << " ms, deviation " << g.deviation << ", volume taken off " << volume(m, s) - volume(m, x));
    ASSERT_GT(volume(m, s), volume(m, x));
}

TEST(tests_brep_full_round, folding_refused)
{
    // (b) along, refused when the march stops on it: the fin 0.01 thick (r = 0.005 beyond the radius of curvature of its trailing
    // edge, 0.0035): the centre line would fold back; refused, naming the middle face
    {
        Model<T> m;
        const auto [s, bottom, rim, top] = fin(m, 0.01);
        const auto before = census(m);
        expect_refused(fillet_full_round(m, s, std::vector<FaceId>{bottom}, std::vector<FaceId>{rim}, std::vector<FaceId>{top}), BuildErrc::NonGeneric,
                       "exceeds a radius of curvature of the middle face along the guide");
        ASSERT_EQ(census(m), before);
    }
    // not (a): a rib with parallel flanks 0.6 apart (r = 0.3) whose top is a B-spline with a narrow
    // ridge at its middle (radius of curvature 0.033 at its crest): the ball does not touch the
    // crest but the shoulders of the ridge, where the top is less curved than the ball; the crest
    // is out of the ball, taken off with the top. Built (a refusal across needs a contact where
    // the face is more curved than the ball: `bending`, tested on its own below)
    {
        Model<T> m;
        const T w = 0.6;
        auto pt = [](T y, T z) { return point<T, 3>{0., y, z}; };
        const std::vector<point<T, 3>> p{pt(-2., 0.), pt(2., 0.), pt(2., 1.), pt(w / 2., 1.), pt(w / 2., 2.), pt(-w / 2., 2.), pt(-w / 2., 1.), pt(-2., 1.)};
        std::vector<VertexId> v;
        for (const auto &q : p)
            v.push_back(unwrap(make_vertex(m, q)));
        std::vector<EdgeId> e;
        for (std::size_t i = 0; i < p.size(); ++i)
            if (i == 4)
            {
                const auto ridge = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{pt(w / 2., 2.), pt(0.1, 2.), pt(0., 2.15), pt(-0.1, 2.), pt(-w / 2., 2.)},
                                                                   std::vector<T>{0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1.}, 2);
                e.push_back(unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{ridge}, 0., 1., v[4], v[5])));
            }
            else
                e.push_back(unwrap(make_edge(m, v[i], v[(i + 1) % p.size()])));
        const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, X, Y}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
        const auto s = solid_of(m, make_prism(m, f, point<T, 3>{2., 0., 0.}));
        expect_valid(m, s);
        const auto f1 = face_near(m, s, {1., -w / 2., 1.5}), f2 = face_near(m, s, {1., w / 2., 1.5});
        FaceId top{};
        for (auto g : explore<FaceId>(m, ShapeId{s}))
            if (!dynamic_cast<const Plane<T> *>(m.face(g).surface.get()))
                top = g;
        ASSERT_TRUE(f1.valid() && f2.valid() && top.valid());
        const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{top}, std::vector<FaceId>{f2}));
        expect_valid(m, x);
        ASSERT_EQ(n_faces(m, x), n_faces(m, s));
        ASSERT_LT(volume(m, x), volume(m, s));
    }
}

TEST(tests_brep_full_round, bending)
{
    // the normal curvature behind the refusals of folding (question 18): a cylinder of radius
    // 0.2 along x, seen from inside (toward its axis) across it: 1 / 0.2; along it: 0; seen
    // from outside: -1 / 0.2
    const CylindricalSurface<T> c{0.2, Frame3<T>{ax2<T, 3>{{{0., 0., 0.}, X, Y}}}, std::array<T, 4>{-pi, pi, -1., 1.}};
    const point<T, 2> uv{0.3, 0.5};
    const auto p = c.value(uv[0], uv[1]);
    const auto inward = (point<T, 3>{p[0], 0., 0.} - p) / norm(point<T, 3>{p[0], 0., 0.} - p);
    const auto across = cross(X, inward);
    ASSERT_NEAR(full_round_detail::bending(c, uv, across, inward), 5., 1e-9);
    ASSERT_NEAR(full_round_detail::bending(c, uv, X, inward), 0., 1e-12);
    ASSERT_NEAR(full_round_detail::bending(c, uv, across, T(-1) * inward), -5., 1e-9);
}


TEST(tests_brep_full_round, not_a_primitive)
{
    // a straight rib whose top is a B-spline with a bump of 2e-5 between two of the 17 stations
    // where the surface is recognized: the stations alone would read a cylinder; the checks
    // between them keep the skin
    Model<T> m;
    const T w = 0.6, L = 2., H = 2.;
    constexpr std::size_t nu = 80;
    points_vector<T, 3> poles;
    for (std::size_t j = 0; j < 2; ++j)
        for (std::size_t i = 0; i < nu; ++i)
        {
            const T x = L * T(i) / T(nu - 1);
            poles.push_back({x, j == 0 ? -w / 2. : w / 2., H + (std::abs(x - 1.08) < 0.02 ? 2e-5 : 0.)});
        }
    std::vector<T> ku{0., 0., 0., 0.};
    for (std::size_t i = 1; i + 3 < nu; ++i)
        ku.push_back(T(i) / T(nu - 3));
    for (int i = 0; i < 4; ++i)
        ku.push_back(1.);
    const auto top_s = std::make_shared<BSSurface<T, 3>>(poles, ku, std::vector<T>{0., 0., 1., 1.}, std::size_t{3}, std::size_t{1});
    const auto top = unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{top_s}));
    auto side = [&](T y) {
        // the flank y: the boundary of the top at y, the vertical lines and the bottom
        points_vector<T, 3> q;
        for (std::size_t i = 0; i < nu; ++i)
            q.push_back(poles[i + (y < 0. ? 0 : nu)]);
        const auto c = std::make_shared<BSCurve<T, 3>>(q, ku, 3);
        const auto a = unwrap(make_vertex(m, c->value(0.))), b = unwrap(make_vertex(m, c->value(1.)));
        const auto a0 = unwrap(make_vertex(m, point<T, 3>{0., y, 0.})), b0 = unwrap(make_vertex(m, point<T, 3>{L, y, 0.}));
        const std::vector<EdgeId> e{unwrap(make_edge(m, a0, b0)), unwrap(make_edge(m, b0, b)), unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{c}, 0., 1., a, b)),
                                    unwrap(make_edge(m, a, a0))};
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., y, 0.}, Y, X}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    };
    auto quad = [&](const std::vector<point<T, 3>> &p) {
        std::vector<EdgeId> e;
        for (std::size_t k = 0; k < p.size(); ++k)
            e.push_back(unwrap(make_edge(m, p[k], p[(k + 1) % p.size()])));
        const auto n = cross(p[1] - p[0], p[2] - p[0]), x = p[1] - p[0];
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{p[0], n / norm(n), x / norm(x)}}, std::array<T, 4>{-10., 10., -10., 10.}), unwrap(make_wire(m, e))));
    };
    const std::vector<FaceId> faces{top, side(-w / 2.), side(w / 2.), quad({{0., -w / 2., 0.}, {L, -w / 2., 0.}, {L, w / 2., 0.}, {0., w / 2., 0.}}),
                                    quad({{0., -w / 2., 0.}, {0., w / 2., 0.}, {0., w / 2., H}, {0., -w / 2., H}}), quad({{L, -w / 2., 0.}, {L, -w / 2., H}, {L, w / 2., H}, {L, w / 2., 0.}})};
    SewOptions<T> so;
    so.tol = 1e-9;
    const auto s = unwrap(make_solid(m, unwrap(sew(m, faces, so)).shells.front()));
    expect_valid(m, s);
    const auto f1 = face_near(m, s, {1., -w / 2., 1.}), f2 = face_near(m, s, {1., w / 2., 1.}), mid = face_near(m, s, {0.5, 0., H});
    ASSERT_TRUE(f1.valid() && f2.valid() && mid.valid());
    const auto g = unwrap(full_round_geometry(m, s, f1, mid, f2));
    MESSAGE("bumped top: exact " << g.exact << ", deviation " << g.deviation);
    ASSERT_FALSE(g.exact);
    const auto x = unwrap(fillet_full_round(m, s, std::vector<FaceId>{f1}, std::vector<FaceId>{mid}, std::vector<FaceId>{f2}));
    expect_valid(m, x);
}

TEST(tests_brep_full_round, one_face_round_a_loop)
{
    // a trough round z (a U of one B-spline revolved, closed by a flat ring): the same face on
    // both sides of the ring, round a loop: refused, saying to split it
    Model<T> m;
    const T R = 2., w = 0.5, h = 1.;
    const auto u = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{R - w / 2., 0., h}, {R - w / 2., 0., 0.}, {R + w / 2., 0., 0.}, {R + w / 2., 0., h}},
                                                   std::vector<T>{0., 0., 1. / 3., 2. / 3., 1., 1.}, 1);
    const auto a = unwrap(make_vertex(m, point<T, 3>{R - w / 2., 0., h})), b = unwrap(make_vertex(m, point<T, 3>{R + w / 2., 0., h}));
    const std::vector<EdgeId> e{unwrap(make_edge(m, std::shared_ptr<Curve<T, 3>>{u}, 0., 1., a, b)), unwrap(make_edge(m, b, a))};
    const auto f = unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, X}}, std::array<T, 4>{-20., 20., -20., 20.}), unwrap(make_wire(m, e))));
    const auto s = explore<SolidId>(m, unwrap(make_revolution(m, ShapeId{f}, ax1<T, 3>{point<T, 3>{0., 0., 0.}, Z}, 2. * pi)).shape).front();
    expect_valid(m, s);
    const auto trough = face_near(m, s, {0., R, 0.}), ring = face_near(m, s, {0., R, h});
    ASSERT_TRUE(trough.valid() && ring.valid());
    const auto before = census(m);
    expect_refused(fillet_full_round(m, s, std::vector<FaceId>{trough}, std::vector<FaceId>{ring}, std::vector<FaceId>{trough}), BuildErrc::NonGeneric, "split it");
    ASSERT_EQ(census(m), before);
}


