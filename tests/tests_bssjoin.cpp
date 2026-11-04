// Two B-spline surfaces joined along a side they share (gbs/bssjoin.h): the pieces of one surface
// joined back into it, whatever the way their parameters run; knots and degrees made the same
// along the side; a kink kept as a knot; refusals; parameters of a few hundreds; a side collapsed
// to a point; the smoothest of two common sides; the deviation within tol; curves end to end.
// Design: docs/sources/design/brep4_pr01ter_join_bsplines.md.
#include <doctest_gtest.hpp>
#include <gbs/bssjoin.h>
#include <gbs/bsstools.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <random>
#include <utility>
#include <vector>

using namespace gbs;

namespace
{
    using T = double;
    using S = BSSurface<T, 3>;

    /// A bicubic surface of 6 x 5 poles, wavy, on [0, 3] x [0, 2].
    auto wavy() -> S
    {
        const std::vector<T> ku{0, 0, 0, 0, 1, 2, 3, 3, 3, 3}, kv{0, 0, 0, 0, 1, 2, 2, 2, 2};
        points_vector<T, 3> poles;
        for (std::size_t j{}; j < 5; ++j)
            for (std::size_t i{}; i < 6; ++i)
                poles.push_back({T(i), T(j), 0.3 * std::sin(1.1 * T(i)) * std::cos(0.7 * T(j)) + 0.05 * T(i * j)});
        return {poles, ku, kv, 3, 3};
    }

    /// The largest distance between f(u, v) and g(u, v) on an n x n grid of the box of f.
    template <typename F, typename G>
    auto gap(const std::array<T, 4> &box, F &&f, G &&g, int n = 12) -> T
    {
        T d{};
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
            {
                const T u = std::min(box[0] + (box[1] - box[0]) * T(i) / T(n), box[1]);
                const T v = std::min(box[2] + (box[3] - box[2]) * T(j) / T(n), box[3]);
                d = std::max(d, norm(f(u, v) - g(u, v)));
            }
        return d;
    }

    auto same_knots(const std::vector<T> &a, const std::vector<T> &b) -> bool
    {
        return a.size() == b.size() && std::ranges::equal(a, b, [](T x, T y) { return std::abs(x - y) <= 1e-12; });
    }

    /// The parameters of b mapped onto the joined surface, put onto its bounds where the rounding
    /// of the map takes them out (a few ulps).
    auto mapped(const SurfaceJoin<T, 3> &j, T u, T v) -> std::array<T, 2>
    {
        const auto &m = j.map;
        const auto b = j.surface.bounds();
        return {std::clamp(m[0] * u + m[1] * v + m[2], b[0], b[1]), std::clamp(m[3] * u + m[4] * v + m[5], b[2], b[3])};
    }

    /// The point of the joined surface at the parameters of b, by the map.
    auto at_b(const SurfaceJoin<T, 3> &j, T u, T v) -> point<T, 3>
    {
        const auto [x, y] = mapped(j, u, v);
        return j.surface.value(x, y);
    }

    /// A dome: u the latitude (u = u_min its apex, a row of equal poles), v the longitude.
    auto dome() -> S
    {
        points_vector<T, 3> poles;
        for (std::size_t j{}; j < 7; ++j)
            for (std::size_t i{}; i < 5; ++i)
            {
                const T r = T(i) / 4., th = 1.5 * T(j) / 6.;
                poles.push_back({r * std::cos(th), r * std::sin(th), 1 - r * r});
            }
        return {poles, {0, 0, 0, 0, 1, 2, 2, 2, 2}, {0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4}, 3, 3};
    }

    /// Both pieces on the joined surface: a at its own parameters, b through the map.
    void expect_on(const SurfaceJoin<T, 3> &j, const S &a, const S &b, T tol)
    {
        EXPECT_LE(gap(a.bounds(), [&](T u, T v) { return j.surface.value(u, v); }, [&](T u, T v) { return a.value(u, v); }), tol);
        EXPECT_LE(gap(b.bounds(), [&](T u, T v) { return at_b(j, u, v); }, [&](T u, T v) { return b.value(u, v); }), tol);
    }
}

TEST(tests_bssjoin, pieces_joined_back)
{
    const auto s = wavy();
    // split between two knots: the knot of the junction removed three times, the surface itself
    {
        const auto [a, b] = split_u(s, 1.5);
        const auto j = join_surfaces(a, b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 3u);
        EXPECT_LE(j->deviation, 1e-12);
        EXPECT_TRUE(same_knots(j->surface.knotsFlatsU(), s.knotsFlatsU()));
        EXPECT_LE(gap(s.bounds(), [&](T u, T v) { return j->surface.value(u, v); }, [&](T u, T v) { return s.value(u, v); }), 1e-12);
        expect_on(*j, a, b, 1e-12);
        const std::array<T, 6> identity{1, 0, 0, 0, 1, 0};
        for (std::size_t k{}; k < 6; ++k)
            EXPECT_NEAR(j->map[k], identity[k], 1e-14);
    }
    // split at a knot (multiplicity 1): removed twice, back to it
    {
        const auto [a, b] = split_u(s, 1.0);
        const auto j = join_surfaces(a, b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 2u);
        EXPECT_TRUE(same_knots(j->surface.knotsFlatsU(), s.knotsFlatsU()));
    }
    // along v
    {
        const auto [a, b] = split_v(s, 0.6);
        const auto j = join_surfaces(a, b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 3u);
        expect_on(*j, a, b, 1e-12);
        EXPECT_TRUE(same_knots(j->surface.knotsFlatsV(), s.knotsFlatsV()));
    }
    // a's side at u = u_min: a keeps its parameters, b goes before it
    {
        const auto [b, a] = split_u(s, 1.5);
        const auto j = join_surfaces(a, b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 3u);
        expect_on(*j, a, b, 1e-12);
    }
}

TEST(tests_bssjoin, the_parameters_of_b_turned)
{
    const auto s = wavy();
    const auto [a, b0] = split_u(s, 1.5);
    // b with its parameters swapped and one reversed, or both reversed: the same surface, facing
    // the same way; joined, through the map
    auto swapped = b0;
    swapped.invertUV();
    swapped.reverseV();
    auto both = b0;
    both.reverseU();
    both.reverseV();
    for (const auto *b : {&swapped, &both})
    {
        const auto j = join_surfaces(a, *b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 3u);
        expect_on(*j, a, *b, 1e-12);
        const auto &m = j->map;
        EXPECT_GT(m[0] * m[4] - m[1] * m[3], 0.);
    }
    // one reversed: b faces the other way, refused
    auto one = b0;
    one.reverseU();
    EXPECT_FALSE(join_surfaces(a, one, 1e-9));
}

TEST(tests_bssjoin, other_knots_degrees_and_speed)
{
    const auto s = wavy();
    const auto [a, b0] = split_u(s, 1.5);
    // b of another degree along v and across, a knot of its own along v, and run 4 times slower
    auto b = b0;
    b.increaseDegreeV();
    b.increaseDegreeU();
    b.insertKnotV(0.37);
    b.changeUBounds(0., 6.);
    const auto j = join_surfaces(a, b, 1e-9);
    ASSERT_TRUE(j);
    EXPECT_EQ(j->surface.degreeU(), 4u);
    EXPECT_EQ(j->surface.degreeV(), 4u);
    EXPECT_GE(j->knots_removed, 3u); // smooth: at least C2 across
    expect_on(*j, a, b, 1e-9);
    EXPECT_NEAR(j->map[0], 0.25, 1e-12); // b's u scaled to a's speed
}

TEST(tests_bssjoin, a_kink_stays_a_knot)
{
    // two flat pieces at right angles along x = 1: joined C0 only, every point kept
    const std::vector<T> k{0, 0, 1, 1};
    const S a{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, k, k, 1, 1};
    const S b{{{1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}, k, k, 1, 1};
    const auto j = join_surfaces(a, b, 1e-9);
    ASSERT_TRUE(j);
    EXPECT_EQ(j->knots_removed, 0u);
    EXPECT_EQ(j->deviation, 0.);
    expect_on(*j, a, b, 1e-14);
}

TEST(tests_bssjoin, sides_apart_or_close)
{
    const auto s = wavy();
    const auto [a, b] = split_u(s, 1.5);
    // b moved away: no side in common
    auto far = b;
    auto poles = far.poles();
    for (auto &p : poles)
        p[2] += 1e-3;
    far.copyPoles(poles);
    EXPECT_FALSE(join_surfaces(a, far, 1e-6));
    // within tol: the common row their mean, the deviation half the gap
    auto near = b;
    poles = near.poles();
    for (auto &p : poles)
        p[2] += 1e-9;
    near.copyPoles(poles);
    const auto j = join_surfaces(a, near, 1e-7);
    ASSERT_TRUE(j);
    EXPECT_GE(j->deviation, 0.5e-9 * (1 - 1e-6));
    EXPECT_LE(j->deviation, 1e-7);
    expect_on(*j, a, near, 1e-7);
}

// Parameters of a few hundreds (millimetres), whatever the side of a and the way b runs: nothing
// evaluated out of the bounds, the bounds of a kept bit for bit, those of b mapped onto the others
// to a few ulps.
TEST(tests_bssjoin, large_parameters)
{
    std::mt19937 g{3};
    std::uniform_real_distribution<T> d{0, 1};
    const T ulp = std::nextafter(T(2048), T(4096)) - T(2048); // (the parameters stay below 2048)
    for (int k = 0; k < 120; ++k)
    {
        const T u0 = 100 + 900 * d(g), L = 50 + 500 * d(g), v0 = 100 + 900 * d(g), H = 50 + 500 * d(g);
        auto s = wavy();
        s.changeUBounds(u0, u0 + L);
        s.changeVBounds(v0, v0 + H);
        auto [a, b] = split_u(s, u0 + L * (0.3 + 0.4 * d(g)));
        const bool a_first = k % 2 == 0; // a's side at u_max, else at u_min
        if (!a_first)
            std::swap(a, b);
        if (k % 3 == 1)
        {
            b.reverseU();
            b.reverseV();
        }
        if (k % 3 == 2)
        {
            b.reverseV();
            b.invertUV();
        }
        const auto j = join_surfaces(a, b, 1e-9);
        ASSERT_TRUE(j);
        EXPECT_EQ(j->knots_removed, 3u);
        const auto jb = j->surface.bounds(), ab = a.bounds();
        EXPECT_EQ(jb[2], ab[2]);
        EXPECT_EQ(jb[3], ab[3]);
        EXPECT_EQ(a_first ? jb[0] : jb[1], a_first ? ab[0] : ab[1]);
        const auto bb = b.bounds();
        const auto &m = j->map;
        for (const T u : {bb[0], bb[1]})
            for (const T v : {bb[2], bb[3]})
            {
                const T x = m[0] * u + m[1] * v + m[2], y = m[3] * u + m[4] * v + m[5];
                EXPECT_GE(x, jb[0] - 4 * ulp);
                EXPECT_LE(x, jb[1] + 4 * ulp);
                EXPECT_GE(y, jb[2] - 4 * ulp);
                EXPECT_LE(y, jb[3] + 4 * ulp);
            }
        expect_on(*j, a, b, 1e-9);
    }
}

// Knots that are not clamped (along the side, or at the far end of a curve): nothing; degree
// elevation and the samples need clamped ends.
TEST(tests_bssjoin, not_clamped_refused)
{
    points_vector<T, 3> pa, pb;
    for (std::size_t j{}; j < 4; ++j)
        for (std::size_t i{}; i < 4; ++i)
        {
            pa.push_back({T(i), T(j), 0});
            pb.push_back({3 + T(i), T(j), 0});
        }
    const S a{pa, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 1, 2, 3, 4, 5, 6}, 3, 2}; // along v: 0, 1, ..., 6
    const S b{pb, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}, 3, 3};
    EXPECT_FALSE(join_surfaces(a, b, 1e-6));
    const BSCurve<T, 3> c{{{0, 0, 0}, {1, 0, 0}, {2, 1, 0}, {3, 1, 1}}, {0, 0, 0, 0, 1, 1, 1, 1}, 3};
    const BSCurve<T, 3> open_end{{{3, 1, 1}, {4, 1, 1}, {5, 0, 1}, {6, 0, 0}}, {1, 1, 1, 1, 2, 3, 4, 5}, 3};
    EXPECT_FALSE(join_curves(c, open_end, 1e-6));
}

// A side collapsed to a point (the apex of a dome) is no side to join along: the pieces of a dome
// cut along a meridian join along it, back to the dome, not through the apex.
TEST(tests_bssjoin, collapsed_side_skipped)
{
    const auto s = dome();
    const auto [a, b] = split_v(s, 2.0);
    const auto j = join_surfaces(a, b, 1e-9);
    ASSERT_TRUE(j);
    EXPECT_EQ(j->knots_removed, 2u); // (v = 2 a knot of multiplicity 1: back to it)
    EXPECT_EQ(j->side_a, SurfaceSide::v_max);
    EXPECT_EQ(j->side_b, SurfaceSide::v_min);
    EXPECT_TRUE(same_knots(j->surface.knotsFlatsV(), s.knotsFlatsV()));
    EXPECT_LE(gap(s.bounds(), [&](T u, T v) { return j->surface.value(u, v); }, [&](T u, T v) { return s.value(u, v); }), 1e-12);
    EXPECT_FALSE(join_surfaces(a, b, SurfaceSide::u_min, SurfaceSide::u_min, 1e-9));
}

// Two surfaces sharing two sides (the halves of a closed band), smooth along one and with a corner
// along the other, which comes first: the smooth side is chosen; the other one when asked for.
TEST(tests_bssjoin, smoothest_side_chosen)
{
    // the profile of a from P0 to P1, that of b from P1 back to P0, leaving P1 as a arrives there;
    // extruded along z
    const auto band = [](std::vector<std::array<T, 3>> profile) {
        points_vector<T, 3> poles;
        for (int j = 0; j < 2; ++j)
            for (auto p : profile)
            {
                p[2] = T(j);
                poles.push_back(p);
            }
        return S{poles, {0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}, 2, 1};
    };
    const auto a = band({{0, 0, 0}, {1, 1, 0}, {2, 0, 0}}), b = band({{2, 0, 0}, {3, -1, 0}, {0, 0, 0}});
    const auto j = join_surfaces(a, b, 1e-9);
    ASSERT_TRUE(j);
    EXPECT_EQ(j->knots_removed, 1u);
    EXPECT_EQ(j->side_a, SurfaceSide::u_max);
    EXPECT_EQ(j->side_b, SurfaceSide::u_min);
    expect_on(*j, a, b, 1e-12);
    const auto c = join_surfaces(a, b, SurfaceSide::u_min, SurfaceSide::u_max, 1e-9);
    ASSERT_TRUE(c);
    EXPECT_EQ(c->knots_removed, 0u);
    expect_on(*c, a, b, 1e-12);
    EXPECT_FALSE(join_surfaces(a, b, SurfaceSide::v_min, SurfaceSide::v_min, 1e-9)); // (not one curve)
}

// The deviation bounds the distance of the join from a and b at the same parameters, and stays
// within tol, whatever the poles of b beyond the junction (moved by up to 0.6 tol here).
TEST(tests_bssjoin, deviation_within_tol)
{
    std::mt19937 g{7};
    std::uniform_real_distribution<T> d{-1, 1};
    const T tol = 1e-3;
    const BSCurve<T, 3> c{{{0, 0, 0}, {1, 0.3, 0}, {2, -0.2, 0.4}, {3, 0.5, 0.1}, {4, 0, 0}, {5, 0.4, 0.2}}, {0, 0, 0, 0, 1, 2, 3, 3, 3, 3}, 3};
    std::size_t smooth{};
    for (int k = 0; k < 400; ++k)
    {
        auto a = c, b = c;
        a.trim(0, 1.5);
        b.trim(1.5, 3);
        auto pb = b.poles();
        for (std::size_t i = 1; i < pb.size(); ++i)
            for (auto &x : pb[i])
                x += 0.6 * tol * d(g);
        b = BSCurve<T, 3>{pb, b.knotsFlats(), 3};
        const auto j = join_curves(a, b, tol);
        ASSERT_TRUE(j);
        EXPECT_LE(j->deviation, tol);
        smooth += j->knots_removed > 0;
        T dist{};
        for (int i = 0; i <= 30; ++i)
        {
            const T t = 1.5 * T(i) / 30.;
            dist = std::max({dist, norm(j->curve.value(t) - a.value(t)), norm(j->curve.value(1.5 + t) - b.value(1.5 + t))});
        }
        EXPECT_LE(dist, j->deviation + 1e-15);
    }
    EXPECT_GT(smooth, 0u);
    const auto s = wavy();
    const auto [a, b0] = split_u(s, 1.5);
    for (int k = 0; k < 100; ++k)
    {
        auto poles = b0.poles();
        for (std::size_t i{}; i < poles.size(); ++i)
            if (i % b0.nPolesU() != 0) // (the common row kept)
                for (auto &x : poles[i])
                    x += 0.6 * tol * d(g);
        auto b = b0;
        b.copyPoles(poles);
        const auto j = join_surfaces(a, b, tol);
        ASSERT_TRUE(j);
        EXPECT_LE(j->deviation, tol);
        expect_on(*j, a, b, j->deviation + 1e-15);
    }
}

// Two curves end to end: the pieces of one curve joined back into it; b on other parameters and of
// another degree shifted to start where a ends, at its own speed; a kink kept; apart, nothing.
TEST(tests_bssjoin, curves)
{
    const BSCurve<T, 3> c{{{0, 0, 0}, {1, 0, 0}, {2, 1, 0}, {3, 1, 1}, {4, 0, 1}}, {0, 0, 0, 0, 1, 2, 2, 2, 2}, 3};
    auto a = c, b = c;
    a.trim(0, 0.7);
    b.trim(0.7, 2);
    const auto j = join_curves(a, b, 1e-9);
    ASSERT_TRUE(j);
    EXPECT_EQ(j->knots_removed, 3u);
    EXPECT_LE(j->deviation, 1e-12);
    EXPECT_TRUE(same_knots(j->curve.knotsFlats(), c.knotsFlats()));
    for (int i = 0; i <= 20; ++i)
        EXPECT_LE(norm(j->curve.value(0.1 * i) - c.value(0.1 * i)), 1e-12);
    // b elsewhere, a degree higher: still pieces of one polynomial, the four knots removed
    auto b2 = b;
    b2.changeBounds(10., 11.3);
    b2.increaseDegree();
    const auto j2 = join_curves(a, b2, 1e-9);
    ASSERT_TRUE(j2);
    EXPECT_EQ(j2->curve.degree(), 4u);
    EXPECT_EQ(j2->knots_removed, 4u);
    EXPECT_EQ(j2->curve.bounds()[1], 0.7 + (11.3 - 10.));
    for (int i = 0; i <= 13; ++i)
    {
        const T t = std::min(10. + 0.1 * i, 11.3);
        EXPECT_LE(norm(j2->curve.value(0.7 + (t - 10.)) - b2.value(t)), 1e-12);
    }
    // a corner: C0, b after a
    const BSCurve<T, 3> k1{{{0, 0, 0}, {1, 0, 0}}, {0, 0, 1, 1}, 1}, k2{{{1, 0, 0}, {1, 1, 0}}, {0, 0, 1, 1}, 1};
    const auto jk = join_curves(k1, k2, 1e-9);
    ASSERT_TRUE(jk);
    EXPECT_EQ(jk->knots_removed, 0u);
    EXPECT_EQ(jk->curve.bounds()[1], 2.);
    const BSCurve<T, 3> k3{{{1, 1e-3, 0}, {1, 1, 0}}, {0, 0, 1, 1}, 1};
    EXPECT_FALSE(join_curves(k1, k3, 1e-6));
    // gbs::join: b's knots shifted so that it ends at a's end + its own length, as the STEP reader
    // computes the bound of a curve across its closure
    auto a3 = a, b3 = b;
    a3.changeBounds(812.3, 845.1);
    b3.changeBounds(17.9, 64.2);
    EXPECT_EQ(join(a3, b3).bounds()[1], 845.1 + (64.2 - 17.9));
}
