#include <doctest_gtest.hpp>
#include <gbs/offsets.h>
#include <gbs/curveoffset.h>
#include <gbs/intersect/surface_surface.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    /// z = (x² + y²) / 2 on [-a, a]², exact as a biquadratic Bézier patch on (u, v) = (x, y):
    /// curvature 1 at the origin, toward the normal S_u x S_v = (-x, -y, 1).
    auto paraboloid(T a = 1.) -> std::shared_ptr<BSSurface<T, 3>>
    {
        const std::array<T, 3> x{-a, 0., a}, c{a * a / 2, -a * a / 2, a * a / 2};
        std::vector<std::array<T, 3>> poles;
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                poles.push_back({x[i], x[j], c[i] + c[j]});
        const std::vector<T> k{-a, -a, -a, a, a, a};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 2, 2);
    }

    /// A wavy bicubic B-spline with interior knots.
    auto wavy() -> std::shared_ptr<BSSurface<T, 3>>
    {
        std::vector<std::array<T, 3>> poles;
        const int n = 6;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                poles.push_back({T(i), T(j), 0.3 * std::sin(1.3 * i) * std::cos(0.9 * j) + 0.05 * i * j});
        const std::vector<T> k{0., 0., 0., 0., 1., 2., 3., 3., 3., 3.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 3, 3);
    }

    /// Finite differences of the exact derivatives of order k - 1 against those of order k.
    void expect_consistent(const Surface<T, 3> &s, T u, T v, std::size_t k, T rel)
    {
        const auto b = s.bounds();
        const T hu = 1e-6 * (b[1] - b[0]), hv = 1e-6 * (b[3] - b[2]);
        for (std::size_t a = 0; a < k; ++a)
        {
            const std::size_t c = k - 1 - a;
            const auto du = (s.value(u + hu, v, a, c) - s.value(u - hu, v, a, c)) / (2 * hu);
            const auto dv = (s.value(u, v + hv, a, c) - s.value(u, v - hv, a, c)) / (2 * hv);
            const auto eu = s.value(u, v, a + 1, c), ev = s.value(u, v, a, c + 1);
            ASSERT_LT(norm(du - eu), rel * (1. + norm(eu))) << u << " " << v << " order " << k << " (" << a + 1 << ", " << c << ")";
            ASSERT_LT(norm(dv - ev), rel * (1. + norm(ev))) << "order " << k << " (" << a << ", " << c + 1 << ")";
        }
    }
}

TEST(tests_offsets, exact_derivatives_against_the_primitives)
{
    // the offset surface of a primitive against the exact offset primitive, on the same parameters:
    // every derivative up to the third
    const Frame3<T> f{ax2<T, 3>{point<T, 3>{0.1, -0.2, 0.3}, point<T, 3>{0.2, 0.3, 1.}, point<T, 3>{1., 0., 0.}}};
    std::vector<std::shared_ptr<ElementarySurface<T>>> prims{
        std::make_shared<CylindricalSurface<T>>(2., f, std::array<T, 4>{0., 2 * pi, -1., 1.}),
        std::make_shared<SphericalSurface<T>>(2., f, std::array<T, 4>{0., 2 * pi, -1.2, 1.2}),
        std::make_shared<ToroidalSurface<T>>(3., 1., f, std::array<T, 4>{0., 2 * pi, 0., 2 * pi}),
        std::make_shared<ConicalSurface<T>>(2., 0.4, f, std::array<T, 4>{0., 2 * pi, -1., 1.}),
        std::make_shared<Plane<T>>(f, std::array<T, 4>{-1., 1., -1., 1.})};
    for (const auto &p : prims)
        for (T d : {0.5, -0.4})
        {
            const OffsetSurface<T> o{p, d};
            const auto e = primitive_offset(*p, d, 1e-10);
            ASSERT_TRUE(e.has_value());
            ASSERT_TRUE(o.has_second_derivatives());
            for (T u : {0.3, 1.7, 4.})
                for (T v : {-0.5, 0.2, 0.9})
                    for (std::size_t a = 0; a <= 3; ++a)
                        for (std::size_t b = 0; a + b <= 3; ++b)
                        {
                            const auto x = o.value(u, v, a, b), y = e->surface->value(u, v + e->dv, a, b);
                            ASSERT_LT(norm(x - y), 1e-11 * (1. + norm(y))) << "(" << a << ", " << b << ") at " << u << ", " << v;
                        }
        }
}

TEST(tests_offsets, exact_derivatives_of_b_splines)
{
    // B-splines, polynomial and rational: each order the derivative of the one below
    auto torus = std::make_shared<ToroidalSurface<T>>(3., 1., Frame3<T>{}, std::array<T, 4>{0., 2 * pi, 0., 2 * pi});
    std::vector<std::shared_ptr<Surface<T, 3>>> basis{paraboloid(), wavy(), torus->to_nurbs()};
    for (const auto &s : basis)
    {
        const OffsetSurface<T> o{s, 0.2};
        const auto b = s->bounds();
        for (T fu : {0.21, 0.46, 0.83})
            for (T fv : {0.17, 0.46, 0.77})
            {
                const T u = b[0] + fu * (b[1] - b[0]), v = b[2] + fv * (b[3] - b[2]);
                for (std::size_t k = 1; k <= 3; ++k)
                    expect_consistent(o, u, v, k, 1e-6);
                // the one-pass derivatives are the values
                std::array<point<T, 3>, 9> D;
                o.derivatives(u, v, 2, 2, D.data());
                for (std::size_t a = 0; a <= 2; ++a)
                    for (std::size_t c = 0; c <= 2; ++c)
                        ASSERT_LT(norm(D[a * 3 + c] - o.value(u, v, a, c)), 1e-12 * (1. + norm(D[a * 3 + c])));
            }
    }
}

TEST(tests_offsets, nurbs_primitives)
{
    // the offset of a NURBS cylinder is a cylinder, its derivatives exact; the poles of a NURBS
    // sphere are limits
    const T r = 1., d = 0.5;
    auto cyl = std::make_shared<CylindricalSurface<T>>(r, Frame3<T>{}, std::array<T, 4>{0., 2 * pi, 0., 2.})->to_nurbs();
    const OffsetSurface<T> o{cyl, d};
    const auto b = cyl->bounds();
    for (int i = 0; i <= 20; ++i)
        for (int j = 0; j <= 4; ++j)
        {
            const T u = b[0] + (b[1] - b[0]) * i / 20., v = b[2] + (b[3] - b[2]) * j / 4.;
            const auto p = o.value(u, v);
            ASSERT_NEAR(std::hypot(p[0], p[1]), r + d, 1e-12);
            // the u derivative is tangent to the circle of radius r + d, of the length of the
            // derivative of the basis times (r + d) / r
            const auto pu = o.value(u, v, 1, 0), su = cyl->value(u, v, 1, 0);
            ASSERT_NEAR(norm(pu), norm(su) * (r + d) / r, 1e-10 * norm(pu));
        }
    auto sph = std::make_shared<SphericalSurface<T>>(1., Frame3<T>{}, std::array<T, 4>{0., 2 * pi, -pi / 2, pi / 2})->to_nurbs();
    const OffsetSurface<T> os{sph, 0.5};
    const auto c = sph->bounds();
    for (T v : {c[2], c[3]})
    {
        const auto p = os.value(0.3 * (c[0] + c[1]), v);
        ASSERT_NEAR(norm(p), 1.5, 1e-9);
        ASSERT_NEAR(std::abs(p[2]), 1.5, 1e-9);
        for (std::size_t a = 0; a <= 2; ++a)
            for (std::size_t e = 0; a + e <= 2; ++e)
            {
                const auto x = os.value(0.3 * (c[0] + c[1]), v, a, e);
                ASSERT_TRUE(std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]));
            }
    }
}

TEST(tests_offsets, fold)
{
    // the paraboloid: principal curvatures 1 at the origin, smaller away; its offset toward the
    // normal folds at d = 1
    auto s = paraboloid();
    const auto m = offset_margin(*s, 0.5, s->bounds());
    ASSERT_NEAR(m.margin, 0.5, 1e-9);
    ASSERT_NEAR(m.uv[0], 0., 1e-5);
    ASSERT_NEAR(m.uv[1], 0., 1e-5);
    // away from the normal: the smallest curvature, at a corner, (1 + 2a²)^(-3/2) along the diagonal
    const auto n = offset_margin(*s, -2., s->bounds());
    ASSERT_GT(n.margin, 1.);
    ASSERT_TRUE(offset_surface<T>(s, 0.9, 1e-7).has_value());
    ASSERT_FALSE(offset_surface<T>(s, 1., 1e-7).has_value());
    ASSERT_FALSE(offset_surface<T>(s, 1.5, 1e-7).has_value());
    ASSERT_TRUE(offset_surface<T>(s, -10., 1e-7).has_value());
    // on a box away from the origin the curvature is smaller: 1.2 does not fold there
    ASSERT_TRUE(offset_surface<T>(s, 1.2, 1e-7, std::array<T, 4>{0.5, 1., 0.5, 1.}).has_value());
    // a sphere: 1 + d / R outward
    auto sph = std::make_shared<SphericalSurface<T>>(2., Frame3<T>{}, std::array<T, 4>{0., 2 * pi, -1., 1.})->to_nurbs();
    ASSERT_NEAR(offset_margin(*sph, -1., sph->bounds()).margin, 0.5, 1e-9);
    ASSERT_FALSE(offset_surface<T>(sph, -2., 1e-7).has_value());
}

TEST(tests_offsets, exact_offsets)
{
    // a primitive gives a primitive, an offset the offset of its basis
    auto cyl = std::make_shared<CylindricalSurface<T>>(1., Frame3<T>{}, std::array<T, 4>{0., 2 * pi, 0., 1.});
    auto e = offset_surface<T>(cyl, 0.5, 1e-9);
    ASSERT_TRUE(e.has_value());
    ASSERT_NEAR(std::dynamic_pointer_cast<CylindricalSurface<T>>(e->surface)->radius(), 1.5, 1e-15);
    ASSERT_FALSE(offset_surface<T>(cyl, -1., 1e-9).has_value());
    auto s = wavy();
    auto o1 = offset_surface<T>(s, 0.1, 1e-9);
    ASSERT_TRUE(o1.has_value());
    auto o2 = offset_surface<T>(o1->surface, 0.15, 1e-9);
    ASSERT_TRUE(o2.has_value());
    const auto *o = dynamic_cast<const OffsetSurface<T> *>(o2->surface.get());
    ASSERT_NE(o, nullptr);
    ASSERT_EQ(o->basis(), s);
    ASSERT_NEAR(o->distance(), 0.25, 1e-15);
    ASSERT_LT(norm(o2->surface->value(1.3, 1.7) - OffsetSurface<T>{o1->surface, 0.15}.value(1.3, 1.7)), 1e-12);
}

TEST(tests_offsets, intersections)
{
    // an offset NURBS cylinder cut by a plane: the circle of radius r + d, by marching
    const T r = 1., d = 0.5;
    auto cyl = std::make_shared<CylindricalSurface<T>>(r, Frame3<T>{}, std::array<T, 4>{0., 2 * pi, 0., 2.})->to_nurbs();
    const OffsetSurface<T> o{cyl, d};
    Frame3<T> f;
    f.o = {0., 0., 0.7};
    const Plane<T> pl{f, {-3., 3., -3., 3.}};
    IntersectOptions<T> opt;
    opt.tol = 1e-7;
    const auto res = intersect(o, pl, opt);
    ASSERT_FALSE(res.curves.empty());
    T length{};
    for (const auto &c : res.curves)
    {
        for (int i = 0; i <= 50; ++i)
        {
            const T t = c.a + (c.b - c.a) * i / 50.;
            const auto p = c.curve->value(t);
            ASSERT_NEAR(std::hypot(p[0], p[1]), r + d, 1e-6);
            ASSERT_NEAR(p[2], 0.7, 1e-6);
            const auto q = c.pc1->value(t);
            ASSERT_LT(norm(o.value(q[0], q[1]) - p), 1e-6);
        }
        ASSERT_LT(c.deviation, 1e-6);
        T l{};
        for (int i = 0; i < 400; ++i)
            l += norm(c.curve->value(c.a + (c.b - c.a) * (i + 1) / 400.) - c.curve->value(c.a + (c.b - c.a) * i / 400.));
        length += l;
    }
    ASSERT_NEAR(length, 2 * pi * (r + d), 1e-4);

    // a plane tangent to the offset cylinder along a generatrix: a tangent line (the second
    // derivatives tell a tangency along a curve from an isolated contact; by finite differences,
    // two points were found)
    Frame3<T> g;
    g.o = {r + d, 0., 0.};
    g.z = {1., 0., 0.};
    g.x = {0., 1., 0.};
    g.y = {0., 0., 1.};
    const Plane<T> tp{g, {-1., 1., -0.5, 2.5}};
    const auto rt = intersect(o, tp, opt);
    ASSERT_EQ(rt.curves.size(), 1u);
    ASSERT_TRUE(rt.points.empty());
    ASSERT_TRUE(rt.curves[0].tangent);
    const auto &tc = rt.curves[0];
    ASSERT_LT(norm(tc.curve->value(tc.a) - point<T, 3>{r + d, 0., 0.}), 1e-7);
    ASSERT_NEAR(norm(tc.curve->value(tc.b) - tc.curve->value(tc.a)), 2., 1e-7);

    // two offset B-splines: the curve lies on both
    auto s1 = paraboloid(1.);
    auto s2 = wavy();
    const OffsetSurface<T> a{s1, -0.1};
    // the wavy surface brought over the paraboloid
    std::vector<std::array<T, 3>> poles = s2->poles();
    for (auto &p : poles)
        p = {(p[0] - 2.5) / 2., (p[1] - 2.5) / 2., p[2] + 0.25};
    auto s3 = std::make_shared<BSSurface<T, 3>>(poles, s2->knotsFlatsU(), s2->knotsFlatsV(), 3, 3);
    const OffsetSurface<T> b{s3, 0.1};
    const auto rb = intersect(a, b, opt);
    ASSERT_FALSE(rb.curves.empty());
    for (const auto &c : rb.curves)
        for (int i = 0; i <= 50; ++i)
        {
            const T t = c.a + (c.b - c.a) * i / 50.;
            const auto p = c.curve->value(t);
            const auto q1 = c.pc1->value(t), q2 = c.pc2->value(t);
            ASSERT_LT(norm(a.value(q1[0], q1[1]) - p), 1e-6);
            ASSERT_LT(norm(b.value(q2[0], q2[1]) - p), 1e-6);
        }
}

TEST(tests_offsets, contact_of_offset_spheres)
{
    // two NURBS spheres of radius 1 offset by 0.5, 3 apart: they touch at one point
    auto s1 = std::make_shared<SphericalSurface<T>>(1., Frame3<T>{}, std::array<T, 4>{0., 2 * pi, -pi / 2, pi / 2})->to_nurbs();
    Frame3<T> f;
    f.o = {3., 0., 0.};
    auto s2 = std::make_shared<SphericalSurface<T>>(1., f, std::array<T, 4>{0., 2 * pi, -pi / 2, pi / 2})->to_nurbs();
    const OffsetSurface<T> a{s1, 0.5}, b{s2, 0.5};
    ASSERT_NEAR(norm(a.value(0.4, 0.3)), 1.5, 1e-12);
    IntersectOptions<T> opt;
    opt.tol = 1e-7;
    const auto res = intersect(a, b, opt);
    ASSERT_TRUE(res.curves.empty());
    ASSERT_EQ(res.points.size(), 1u);
    ASSERT_LT(norm(res.points[0].p - point<T, 3>{1.5, 0., 0.}), 1e-6);
}

TEST(tests_offsets, offset_curves)
{
    // a circle offset inward in its plane, by D x C' with D its axis: the circle of radius r - f,
    // every derivative up to the third
    const T r = 2., f0 = 0.5;
    auto circle = std::make_shared<EllipseCurve<T>>(r, r, Frame3<T>{}, 0., 2 * pi);
    auto cst = [f0](T, std::size_t d) -> T { return d == 0 ? f0 : 0.; };
    const CurveOffset3D<T, decltype(cst)> o3{circle, cst, {0., 0., 1.}};
    const EllipseCurve<T> inner{r - f0, r - f0, Frame3<T>{}, 0., 2 * pi};
    for (T u : {0.3, 2.1, 5.})
        for (std::size_t d = 0; d <= 3; ++d)
            ASSERT_LT(norm(o3.value(u, d) - inner.value(u, d)), 1e-12) << d;

    // a variable offset of B-splines, in the plane and in space: each order the derivative of
    // the one below
    auto var = [](T u, std::size_t d) -> T { return d == 0 ? 0.2 + 0.1 * u * u : d == 1 ? 0.2 * u : d == 2 ? 0.2 : 0.; };
    const std::vector<T> k{0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
    auto c2 = std::make_shared<BSCurve<T, 2>>(std::vector<std::array<T, 2>>{{0., 0.}, {1., 2.}, {3., -1.}, {5., 1.}, {6., 3.}}, k, 3);
    auto c3 = std::make_shared<BSCurve<T, 3>>(std::vector<std::array<T, 3>>{{0., 0., 0.}, {1., 2., 1.}, {3., -1., 0.5}, {5., 1., 2.}, {6., 3., 1.}}, k, 3);
    const CurveOffset2D<T, decltype(var)> p2{c2, var};
    const CurveOffset3D<T, decltype(var)> p3{c3, var, {0., 0., 1.}};
    const T h = 1e-6;
    for (T u : {0.2, 0.7})
    {
        // the offset point: at f along the left normal (plane), along z x C' (space)
        const auto t2 = c2->value(u, 1);
        ASSERT_LT(norm(p2.value(u) - c2->value(u) - var(u, 0) * point<T, 2>{-t2[1], t2[0]} / norm(t2)), 1e-14);
        const auto m3 = cross(point<T, 3>{0., 0., 1.}, c3->value(u, 1));
        ASSERT_LT(norm(p3.value(u) - c3->value(u) - var(u, 0) * m3 / norm(m3)), 1e-14);
        for (std::size_t d = 1; d <= 3; ++d)
        {
            const auto e2 = p2.value(u, d);
            const auto e3 = p3.value(u, d);
            ASSERT_LT(norm((p2.value(u + h, d - 1) - p2.value(u - h, d - 1)) / (2 * h) - e2), 1e-6 * (1. + norm(e2))) << d;
            ASSERT_LT(norm((p3.value(u + h, d - 1) - p3.value(u - h, d - 1)) / (2 * h) - e3), 1e-6 * (1. + norm(e3))) << d;
        }
    }
}
