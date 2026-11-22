// n-sided surfaces from ribbons (stage 6, PR 2, docs/sources/design/brep6_pr02_nsided.md):
// fill_n_sided on holes of 3 to 6 sides cut in known surfaces (a plane, a sphere, a paraboloid,
// a saddle), the B-spline within the tolerances along the sides, near the surface the hole was
// cut from, the same under a rotation of the sides or a constant factor on the cross directions
// of a side; refusals, each with its message.
#include <doctest_gtest.hpp>
#include <gbs/nsided.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <string>
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

    /// The octant of the unit sphere over the triangle (0, 0), (1, 0), (0, 1): the point (u, v, 1 - u - v) on the sphere.
    Map octant()
    {
        return Map{[](T u, T v) {
            const T w = std::max(1. - u - v, 0.);
            const T l = std::sqrt(u * u + v * v + w * w);
            return P3{u / l, v / l, w / l};
        }};
    }

    /// The message of the exception E that fill_n_sided throws on the ribbons (empty if it does not throw).
    template <class E = std::invalid_argument>
    std::string refusal(const std::vector<Ribbon<T>> &r, const NSidedOptions<T> &o = {})
    {
        try
        {
            (void)fill_n_sided(r, o);
        }
        catch (const E &e)
        {
            return e.what();
        }
        return {};
    }

    bool says(const std::string &message, const char *what)
    {
        INFO("message: " << message);
        return message.find(what) != std::string::npos;
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
    const NSidedOptions<T> o;
    const auto r = fill_n_sided(hole(octant(), {{0., 0.}, {1., 0.}, {0., 1.}}), o);
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

TEST(tests_nsided, length_of_the_cross_directions)
{
    // a constant factor on the cross directions of a side does not count: the same surface with
    // them scaled by i + 1 on side i (only how their length varies along a side shapes the cross
    // derivatives)
    const Map srf{[](T x, T y) { return P3{x, y, (x * x - 0.5 * y * y) / 2. + 0.1 * x * y * y}; }};
    const auto r = hole(srf, polygon(5, 1., 0.3));
    auto s = r;
    for (std::size_t i = 0; i < s.size(); ++i)
        s[i].cross = std::make_shared<FunctionCurve>([c = r[i].cross, i](T t) {
            const auto v = c->value(t);
            const T k = T(i + 1);
            return P3{k * v[0], k * v[1], k * v[2]};
        });
    const auto a = fill_n_sided(r), b = fill_n_sided(s);
    ASSERT_EQ(a.spans, b.spans);
    // (inside the polygon; up to the rounding of a system of normal equations of large condition)
    T gap{};
    for (int i = 0; i <= 20; ++i)
        for (int j = 0; j <= 20; ++j)
            if (const std::array<T, 2> p{-1. + 0.1 * i, -1. + 0.1 * j}; std::ranges::min(nsided_detail::distances(a.polygon, p)) > 0.)
                gap = std::max(gap, norm(a.surface->value(p[0], p[1]) - b.surface->value(p[0], p[1])));
    INFO("gap " << gap);
    ASSERT_LE(gap, 1e-6);
}

TEST(tests_nsided, refusals)
{
    const Map plane{[](T x, T y) { return P3{x, y, 0.}; }};
    const auto r = hole(plane, polygon(4, 1.));
    // two ribbons
    ASSERT_TRUE(says(refusal({r[0], r[1]}), "fewer than three ribbons"));
    // options out of range, or not finite
    const T nan = std::numeric_limits<T>::quiet_NaN(), inf = std::numeric_limits<T>::infinity();
    for (const auto &set : std::vector<std::function<void(NSidedOptions<T> &)>>{
             [](auto &o) { o.approx_tol = 0.; }, [](auto &o) { o.angular_tol = -1e-5; }, [nan](auto &o) { o.approx_tol = nan; },
             [](auto &o) { o.side_weight = 0.; }, [inf](auto &o) { o.side_weight = inf; }, [](auto &o) { o.cross_weight = -1.; },
             [](auto &o) { o.shape_tol = -1e-3; }})
    {
        NSidedOptions<T> o;
        set(o);
        ASSERT_TRUE(says(refusal(r, o), "not positive"));
    }
    {
        NSidedOptions<T> o;
        o.max_spans = 4;
        ASSERT_TRUE(says(refusal(r, o), "fewer than 8 spans"));
    }
    // a ribbon without curve, without cross direction, or on an empty interval
    for (int k = 0; k < 3; ++k)
    {
        auto s = r;
        if (k == 0)
            s[2].curve = nullptr;
        else if (k == 1)
            s[2].cross = nullptr;
        else
            s[2].t1 = s[2].t0;
        ASSERT_TRUE(says(refusal(s), "without curve"));
    }
    // a curve whose derivative vanishes at its start (the side drawn as t^2)
    {
        auto s = r;
        const auto a = r[1].curve->value(0.), b = r[1].curve->value(1.);
        s[1].curve = std::make_shared<FunctionCurve>([a, b](T t) { return P3{a[0] + t * t * (b[0] - a[0]), a[1] + t * t * (b[1] - a[1]), 0.}; });
        ASSERT_TRUE(says(refusal(s), "derivative vanishes"));
    }
    // a cross direction along its curve
    {
        auto s = r;
        s[0].cross = std::make_shared<FunctionCurve>([c = r[0].curve](T u) { return c->value(u, 1); });
        ASSERT_TRUE(says(refusal(s), "along its curve"));
    }
    // sides that do not follow one another
    {
        auto s = r;
        std::swap(s[1], s[2]);
        ASSERT_TRUE(says(refusal(s), "do not follow"));
    }
    // a gap of 1e-4 at the corners of a side of a hole a thousand times as large: within 1e-6 of
    // the length of the sides, but above twice approx_tol
    {
        auto s = hole(plane, polygon(4, 1000.));
        s[1].curve = std::make_shared<FunctionCurve>([c = s[1].curve](T t) {
            auto p = c->value(t);
            p[2] += 1e-4;
            return p;
        });
        ASSERT_TRUE(says(refusal(s), "do not follow"));
    }
    // ribbons whose normals disagree at a corner: a side tilted out of the plane, or its cross
    // direction turned out of the hole (refused at once, not after the spans have been doubled
    // to max_spans)
    {
        auto s = r;
        s[1].cross = std::make_shared<FunctionCurve>([c = r[1].cross](T x) { const auto v = c->value(x); return P3{v[0], v[1], v[2] + 0.5}; });
        ASSERT_TRUE(says(refusal(s), "normals disagree"));
        s[1].cross = std::make_shared<FunctionCurve>([c = r[1].cross](T x) { const auto v = c->value(x); return P3{-v[0], -v[1], -v[2]}; });
        ASSERT_TRUE(says(refusal(s), "normals disagree"));
    }
    // a hole not convex at a corner: a dart, its last corner coming back inside
    ASSERT_TRUE(says(refusal(hole(plane, {{0., 0.}, {2., 1.}, {0., 2.}, {0.5, 1.}})), "not convex"));
}

TEST(tests_nsided, tolerances_not_reached)
{
    // the octant, 8 spans at most: the message gives the three criteria at the last doubling
    NSidedOptions<T> o;
    o.max_spans = 8;
    const auto m = refusal<std::runtime_error>(hole(octant(), {{0., 0.}, {1., 0.}, {0., 1.}}), o);
    ASSERT_TRUE(says(m, "with 8 spans"));
    ASSERT_TRUE(says(m, "in position for"));
    ASSERT_TRUE(says(m, "in angle for"));
    ASSERT_TRUE(says(m, "the shape moved by"));
}
