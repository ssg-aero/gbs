#include <doctest_gtest.hpp>
#include <gbs-mesh/cdt.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <set>
#include <vector>

using namespace gbs;

namespace
{
    using T = double;
    using CDT = ConstrainedDelaunay<T>;
    using P = point<T, 2>;

    T area(const P &a, const P &b, const P &c) { return ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.; }

    /// Every triangle counter-clockwise, and the total area.
    T total_area(const CDT &cdt)
    {
        T s{};
        for (const auto &t : cdt.triangles())
        {
            const T a = area(cdt.point_at(t[0]), cdt.point_at(t[1]), cdt.point_at(t[2]));
            REQUIRE(a > 0.);
            s += a;
        }
        return s;
    }

    /// No point strictly inside the circumcircle of a triangle (relative margin).
    void expect_delaunay(const CDT &cdt)
    {
        const auto pts = cdt.points();
        for (const auto &t : cdt.triangles())
        {
            const auto &a = pts[t[0]], &b = pts[t[1]], &c = pts[t[2]];
            const T d = 2. * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
            const T ux = ((a[0] * a[0] + a[1] * a[1]) * (b[1] - c[1]) + (b[0] * b[0] + b[1] * b[1]) * (c[1] - a[1]) + (c[0] * c[0] + c[1] * c[1]) * (a[1] - b[1])) / d;
            const T uy = ((a[0] * a[0] + a[1] * a[1]) * (c[0] - b[0]) + (b[0] * b[0] + b[1] * b[1]) * (a[0] - c[0]) + (c[0] * c[0] + c[1] * c[1]) * (b[0] - a[0])) / d;
            const T r = std::hypot(a[0] - ux, a[1] - uy);
            for (std::size_t k = 0; k < pts.size(); ++k)
                if (k != t[0] && k != t[1] && k != t[2])
                    REQUIRE(std::hypot(pts[k][0] - ux, pts[k][1] - uy) >= r * (1. - 1e-9));
        }
    }

    /// The number of points on the convex hull (monotone chain, collinear points left out).
    std::size_t hull_size(std::vector<P> p)
    {
        std::ranges::sort(p);
        std::vector<P> h;
        for (int pass = 0; pass < 2; ++pass)
        {
            const auto start = h.size();
            for (const auto &x : p)
            {
                while (h.size() >= start + 2 && area(h[h.size() - 2], h.back(), x) <= 0.)
                    h.pop_back();
                h.push_back(x);
            }
            h.pop_back();
            std::ranges::reverse(p);
        }
        return h.size();
    }
}

TEST(tests_mesh_cdt, delaunay_of_random_points)
{
    std::mt19937 gen{42};
    std::uniform_real_distribution<T> u{0.01, 0.99};
    CDT cdt{{0., 1., 0., 1.}};
    std::vector<P> pts;
    for (int k = 0; k < 500; ++k)
    {
        pts.push_back({u(gen), u(gen)});
        ASSERT_EQ(cdt.insert(pts.back()), static_cast<CDT::index>(k));
    }
    ASSERT_TRUE(cdt.valid());
    expect_delaunay(cdt);
    // a point given again is the same point
    ASSERT_EQ(cdt.insert(pts[17]), 17u);
    ASSERT_EQ(cdt.size(), 500u);
    // in the square 0..1 (its sides constrained), n points with h = 4 on the hull: 2 n - 2 - h
    // triangles, still Delaunay (the sides are edges of the convex hull)
    std::vector<CDT::index> corners;
    for (const P &c : {P{0., 0.}, P{1., 0.}, P{1., 1.}, P{0., 1.}})
    {
        corners.push_back(cdt.insert(c));
        pts.push_back(c);
    }
    for (std::size_t k = 0; k < 4; ++k)
        cdt.constrain(corners[k], corners[(k + 1) % 4]);
    cdt.remove_outside();
    ASSERT_TRUE(cdt.valid());
    ASSERT_EQ(hull_size(pts), 4u);
    ASSERT_EQ(cdt.triangles().size(), 2 * pts.size() - 2 - 4);
    ASSERT_NEAR(total_area(cdt), 1., 1e-12);
    expect_delaunay(cdt);
}

TEST(tests_mesh_cdt, grid)
{
    // cocircular points everywhere: a valid triangulation of the 9 x 9 cells
    CDT cdt{{0., 9., 0., 9.}};
    for (int j = 0; j <= 9; ++j)
        for (int i = 0; i <= 9; ++i)
            cdt.insert({T(i), T(j)});
    ASSERT_EQ(cdt.triangles().size(), 2u * 81u);
    ASSERT_NEAR(total_area(cdt), 81., 1e-9);
}

TEST(tests_mesh_cdt, constraints)
{
    // random points and a long segment across them: recovered, then the triangulation is still
    // valid and Delaunay away from it
    std::mt19937 gen{7};
    std::uniform_real_distribution<T> u{0.05, 0.95};
    CDT cdt{{0., 1., 0., 1.}};
    const auto a = cdt.insert({0.01, 0.02}), b = cdt.insert({0.99, 0.97});
    for (int k = 0; k < 300; ++k)
        cdt.insert({u(gen), u(gen)});
    cdt.constrain(a, b);
    ASSERT_TRUE(cdt.constrained(a, b));
    ASSERT_TRUE(cdt.valid());
    total_area(cdt);
    // a segment through vertices is split at them: the diagonal of a grid
    CDT grid{{0., 4., 0., 4.}};
    std::vector<CDT::index> id;
    for (int j = 0; j <= 4; ++j)
        for (int i = 0; i <= 4; ++i)
            id.push_back(grid.insert({T(i), T(j)}));
    grid.constrain(id[0], id[24]);
    for (int k = 0; k < 4; ++k)
        ASSERT_TRUE(grid.constrained(id[k * 6], id[(k + 1) * 6]));
    // a point inserted on a constrained edge splits it
    const auto m = grid.insert({0.5, 0.5});
    ASSERT_TRUE(grid.constrained(id[0], m));
    ASSERT_TRUE(grid.constrained(m, id[6]));
    ASSERT_TRUE(grid.valid());
    // constraints crossing each other are refused: (2, 1) to (1, 2) crosses the diagonal at
    // (1.5, 1.5), not a vertex
    ASSERT_THROW(grid.constrain(id[7], id[11]), std::runtime_error);
}

TEST(tests_mesh_cdt, domains)
{
    // a concave L: three unit squares
    {
        const std::vector<P> l{{0., 0.}, {2., 0.}, {2., 1.}, {1., 1.}, {1., 2.}, {0., 2.}};
        CDT cdt{{0., 2., 0., 2.}};
        std::vector<CDT::index> id;
        for (const auto &p : l)
            id.push_back(cdt.insert(p));
        for (std::size_t k = 0; k < id.size(); ++k)
            cdt.constrain(id[k], id[(k + 1) % id.size()]);
        cdt.remove_outside();
        ASSERT_NEAR(total_area(cdt), 3., 1e-12);
        ASSERT_EQ(cdt.triangles().size(), 4u);
        ASSERT_THROW(cdt.insert({0.5, 0.5}), std::logic_error);
    }
    // a square 0..4 with a square hole 1..3, the boundaries refined, points inside: area 12
    {
        std::mt19937 gen{3};
        std::uniform_real_distribution<T> u{0., 4.};
        CDT cdt{{0., 4., 0., 4.}};
        auto loop = [&](T lo, T hi, int n) {
            std::vector<CDT::index> id;
            const std::array<P, 4> c{P{lo, lo}, P{hi, lo}, P{hi, hi}, P{lo, hi}};
            for (int s = 0; s < 4; ++s)
                for (int k = 0; k < n; ++k)
                {
                    const T f = T(k) / n;
                    id.push_back(cdt.insert({c[s][0] + f * (c[(s + 1) % 4][0] - c[s][0]), c[s][1] + f * (c[(s + 1) % 4][1] - c[s][1])}));
                }
            return id;
        };
        const auto outer = loop(0., 4., 8), hole = loop(1., 3., 5);
        for (int k = 0; k < 200; ++k)
        {
            const P p{u(gen), u(gen)};
            if (p[0] > 0.9 && p[0] < 3.1 && p[1] > 0.9 && p[1] < 3.1)
                continue;
            cdt.insert(p);
        }
        for (const auto *l : {&outer, &hole})
            for (std::size_t k = 0; k < l->size(); ++k)
                cdt.constrain((*l)[k], (*l)[(k + 1) % l->size()]);
        cdt.remove_outside();
        ASSERT_TRUE(cdt.valid());
        ASSERT_NEAR(total_area(cdt), 12., 1e-9);
        // no triangle inside the hole
        for (const auto &t : cdt.triangles())
        {
            const auto &a = cdt.point_at(t[0]), &b = cdt.point_at(t[1]), &c = cdt.point_at(t[2]);
            const P g{(a[0] + b[0] + c[0]) / 3., (a[1] + b[1] + c[1]) / 3.};
            ASSERT_FALSE(g[0] > 1. && g[0] < 3. && g[1] > 1. && g[1] < 3.);
        }
    }
}
