// n-sided surfaces from ribbons (stage 6, PR 2, docs/sources/design/brep6_pr02_nsided.md):
// fill_n_sided on holes of 3 to 6 sides cut in known surfaces (a plane, a sphere, a paraboloid,
// a saddle), the B-spline within the tolerances along the sides, near the surface the hole was
// cut from, the same under a rotation of the sides; refusals.
#include <doctest_gtest.hpp>
#include <gbs/nsided.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numbers>
#include <vector>

using namespace gbs;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;
    using P3 = std::array<T, 3>;

    /// A curve given by a function of t on [0, 1] and its derivative (by central differences).
    class FunctionCurve : public Curve<T, 3>
    {
    public:
        explicit FunctionCurve(std::function<P3(T)> f) : f_{std::move(f)} {}
        auto value(T u, size_t d = 0) const -> P3 override
        {
            if (d == 0)
                return f_(u);
            const T e = 1e-5;
            if (d == 1)
            {
                const auto a = f_(u + e), b = f_(u - e);
                return {(a[0] - b[0]) / (2 * e), (a[1] - b[1]) / (2 * e), (a[2] - b[2]) / (2 * e)};
            }
            const auto a = f_(u + e), b = f_(u), c = f_(u - e);
            return {(a[0] - 2 * b[0] + c[0]) / (e * e), (a[1] - 2 * b[1] + c[1]) / (e * e), (a[2] - 2 * b[2] + c[2]) / (e * e)};
        }
        auto bounds() const -> std::array<T, 2> override { return {0., 1.}; }

    private:
        std::function<P3(T)> f_;
    };

    /// A surface z = f(x, y) (or any map of the plane) with its normal by differences.
    struct Map
    {
        std::function<P3(T, T)> f;
        auto du(T x, T y) const
        {
            const T e = 1e-6;
            const auto a = f(x + e, y), b = f(x - e, y);
            return P3{(a[0] - b[0]) / (2 * e), (a[1] - b[1]) / (2 * e), (a[2] - b[2]) / (2 * e)};
        }
        auto dv(T x, T y) const
        {
            const T e = 1e-6;
            const auto a = f(x, y + e), b = f(x, y - e);
            return P3{(a[0] - b[0]) / (2 * e), (a[1] - b[1]) / (2 * e), (a[2] - b[2]) / (2 * e)};
        }
        auto normal(T x, T y) const
        {
            const auto n = cross(du(x, y), dv(x, y));
            return P3{n[0] / norm(n), n[1] / norm(n), n[2] / norm(n)};
        }
    };

    /**
     * The ribbons of the hole bounded, in the parameters of `s`, by the polygon `q` (counter-
     * clockwise): each side the image of a segment, its cross direction the image of the inner
     * normal of the segment.
     */
    std::vector<Ribbon<T>> hole(const Map &s, const std::vector<std::array<T, 2>> &q)
    {
        std::vector<Ribbon<T>> r;
        for (std::size_t i = 0; i < q.size(); ++i)
        {
            const auto a = q[i], b = q[(i + 1) % q.size()];
            const T nx = -(b[1] - a[1]), ny = b[0] - a[0];
            auto at = [a, b](T t) { return std::array<T, 2>{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])}; };
            Ribbon<T> rb;
            rb.curve = std::make_shared<FunctionCurve>([s, at](T t) { const auto p = at(t); return s.f(p[0], p[1]); });
            rb.cross = std::make_shared<FunctionCurve>([s, at, nx, ny](T t) {
                const auto p = at(t);
                const auto u = s.du(p[0], p[1]), v = s.dv(p[0], p[1]);
                return P3{nx * u[0] + ny * v[0], nx * u[1] + ny * v[1], nx * u[2] + ny * v[2]};
            });
            r.push_back(rb);
        }
        return r;
    }

    std::vector<std::array<T, 2>> polygon(std::size_t m, T radius, T turn = 0.)
    {
        std::vector<std::array<T, 2>> q;
        for (std::size_t k = 0; k < m; ++k)
            q.push_back({radius * std::cos(turn + 2 * pi * T(k) / T(m)), radius * std::sin(turn + 2 * pi * T(k) / T(m))});
        return q;
    }

    /// The largest distance from the trimmed surface (sampled inside its polygon) to the closest point of the map near the guess.
    T gap_inside(const NSidedSurface<T> &r, const std::function<T(const P3 &)> &dist)
    {
        T g{};
        for (int a = 0; a <= 40; ++a)
            for (int b = 0; b <= 40; ++b)
            {
                const std::array<T, 2> p{-1. + a / 20., -1. + b / 20.};
                const auto d = nsided_detail::distances(r.polygon, p);
                if (*std::ranges::min_element(d) <= 0.)
                    continue;
                g = std::max(g, dist(r.surface->value(p[0], p[1])));
            }
        return g;
    }

    void expect_tolerances(const NSidedSurface<T> &r, const NSidedOptions<T> &o)
    {
        INFO("spans " << r.spans << ", deviation " << r.deviation << ", angle " << r.angle);
        ASSERT_LE(r.deviation, o.approx_tol);
        ASSERT_LE(r.angle, o.angular_tol);
    }
}

TEST(tests_nsided, plane)
{
    // a hole in a plane: the surface is the plane
    const Map plane{[](T x, T y) { return P3{x, y, 0.5 * x - 0.25 * y + 1.}; }};
    for (std::size_t m : {3u, 4u, 5u, 6u})
    {
        INFO("sides " << m);
        const NSidedOptions<T> o;
        const auto r = fill_n_sided(hole(plane, polygon(m, 1.)), o);
        expect_tolerances(r, o);
        ASSERT_LE(gap_inside(r, [](const P3 &p) { return std::abs(p[2] - (0.5 * p[0] - 0.25 * p[1] + 1.)) / std::sqrt(1.3125); }), 1e-5);
    }
}

TEST(tests_nsided, paraboloid_and_saddle)
{
    // holes of 3 to 6 sides in z = (x^2 + y^2) / 2 and in z = (x^2 - y^2) / 2, bounded by the
    // images of regular polygons: within the tolerances along the sides, and the surface itself
    // inside (biharmonic, found again by the least bending energy)
    for (const T sign : {1., -1.})
    {
        const Map srf{[sign](T x, T y) { return P3{x, y, (x * x + sign * y * y) / 2.}; }};
        for (std::size_t m : {3u, 4u, 5u, 6u})
        {
            INFO("sign " << sign << ", sides " << m);
            const NSidedOptions<T> o;
            const auto r = fill_n_sided(hole(srf, polygon(m, 1., 0.3)), o);
            expect_tolerances(r, o);
            ASSERT_LE(gap_inside(r, [sign](const P3 &p) { return std::abs(p[2] - (p[0] * p[0] + sign * p[1] * p[1]) / 2.); }), 1e-4);
        }
    }
}

TEST(tests_nsided, sphere)
{
    // the octant of the unit sphere, bounded by three quarter great circles: within the
    // tolerances along them; inside, the surface of least bending energy is not the sphere: 8.3 %
    // of the radius off it at most (its right angles drawn from those of an equilateral triangle)
    const Map sphere{[](T u, T v) {
        // (u, v) in the triangle (0, 0), (1, 0), (0, 1): spherical coordinates of the octant
        const T w = std::max(1. - u - v, 0.);
        const P3 q{u, v, w};
        const T l = std::sqrt(u * u + v * v + w * w);
        return P3{q[0] / l, q[1] / l, q[2] / l};
    }};
    const NSidedOptions<T> o;
    const auto r = fill_n_sided(hole(sphere, {{0., 0.}, {1., 0.}, {0., 1.}}), o);
    expect_tolerances(r, o);
    ASSERT_LE(gap_inside(r, [](const P3 &p) { return std::abs(norm(p) - 1.); }), 0.085);
}

TEST(tests_nsided, rotation_of_the_sides)
{
    // the same hole, its sides numbered from another one: the same surface in space (the domain
    // turned by two sides), up to the grid of the B-spline, which does not turn with it
    const Map srf{[](T x, T y) { return P3{x, y, (x * x - 0.5 * y * y) / 2. + 0.1 * x * y * y}; }};
    const auto q = polygon(5, 1., 0.3);
    auto r1 = hole(srf, q);
    auto r2 = r1;
    std::rotate(r2.begin(), r2.begin() + 2, r2.end());
    const auto a = fill_n_sided(r1), b = fill_n_sided(r2);
    const T turn = 2. * 2. * pi / 5.;
    T gap{};
    for (int i = 0; i <= 10; ++i)
        for (int j = 0; j <= 10; ++j)
        {
            const std::array<T, 2> p{-0.6 + 0.12 * i, -0.6 + 0.12 * j};
            if (std::ranges::min(nsided_detail::distances(a.polygon, p)) <= 0.)
                continue;
            const std::array<T, 2> pr{std::cos(-turn) * p[0] - std::sin(-turn) * p[1], std::sin(-turn) * p[0] + std::cos(-turn) * p[1]};
            gap = std::max(gap, norm(a.surface->value(p[0], p[1]) - b.surface->value(pr[0], pr[1])));
        }
    ASSERT_LE(gap, 1e-4);
}

TEST(tests_nsided, refusals)
{
    const Map plane{[](T x, T y) { return P3{x, y, 0.}; }};
    auto r = hole(plane, polygon(4, 1.));
    // two ribbons
    ASSERT_THROW((void)fill_n_sided(std::vector<Ribbon<T>>{r[0], r[1]}), std::invalid_argument);
    // sides that do not follow one another
    auto s = r;
    std::swap(s[1], s[2]);
    ASSERT_THROW((void)fill_n_sided(s), std::invalid_argument);
    // a cross direction along its curve
    auto t = r;
    t[0].cross = t[0].curve->value(0.5, 1)[0] != 0. ? std::make_shared<FunctionCurve>([c = t[0].curve](T u) { return c->value(u, 1); }) : t[0].cross;
    ASSERT_THROW((void)fill_n_sided(t), std::invalid_argument);
    // ribbons whose normals disagree at a corner: a side tilted out of the plane
    auto u = r;
    u[1].cross = std::make_shared<FunctionCurve>([c = r[1].cross](T x) { const auto v = c->value(x); return P3{v[0], v[1], v[2] + 0.5}; });
    NSidedOptions<T> o;
    o.max_spans = 16;
    ASSERT_THROW((void)fill_n_sided(u, o), std::runtime_error);
}




