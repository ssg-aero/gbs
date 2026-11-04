// Two B-spline surfaces joined along a side they share (gbs/bssjoin.h): the pieces of one surface
// joined back into it, whatever the way their parameters run; knots and degrees made the same
// along the side; a kink kept as a knot; refusals.
// Design: docs/sources/design/brep4_pr01ter_join_bsplines.md.
#include <doctest_gtest.hpp>
#include <gbs/bssjoin.h>
#include <gbs/bsstools.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
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
                const T u = box[0] + (box[1] - box[0]) * T(i) / T(n), v = box[2] + (box[3] - box[2]) * T(j) / T(n);
                d = std::max(d, norm(f(u, v) - g(u, v)));
            }
        return d;
    }

    auto same_knots(const std::vector<T> &a, const std::vector<T> &b) -> bool
    {
        return a.size() == b.size() && std::ranges::equal(a, b, [](T x, T y) { return std::abs(x - y) <= 1e-12; });
    }

    /// The point of the joined surface at the parameters of b, by the map.
    auto at_b(const SurfaceJoin<T, 3> &j, T u, T v) -> point<T, 3>
    {
        const auto &m = j.map;
        return j.surface.value(m[0] * u + m[1] * v + m[2], m[3] * u + m[4] * v + m[5]);
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
