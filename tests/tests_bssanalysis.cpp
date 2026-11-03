#include <doctest_gtest.hpp>

#include <gbs/bssanalysis.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
using gbs::operator+;
using gbs::operator-;
using gbs::operator*;
using gbs::operator/;
const double tol = 1e-10;
TEST(tests_bssanalysis, discretize_basic)
{
    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;
    auto u = 2.5;
    auto v = 1.;


    std::vector<std::array<double,4> > poles_t = {
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,6,8,2}, {12,24,12,6}, {4,6,0,2},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1}
                                                    };

    //Pij avec i inner loop
    std::vector<std::array<double,4> > poles(poles_t.size());
    int ni = 5 , nj =8;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            poles[j + nj * i] = poles_t[i + ni * j];
        }
    }

    gbs::BSSurface<double,4> srf(poles,ku,kv,p,q);

    auto pts = gbs::discretize(srf,20,30);

}

namespace
{
    const double pi = std::numbers::pi;
    // a sphere of radius 2 about (1, -1, 0.5), turned, its axis along (1, 1, 1): poles at v = +-pi / 2
    auto sphere() -> gbs::SphericalSurface<double>
    {
        const gbs::Frame3<double> f{gbs::ax2<double, 3>{{{1., -1., 0.5}, {1., 1., 1.}, {1., -1., 0.}}}};
        return {2., f, {0., 2. * pi, -pi / 2., pi / 2.}};
    }
}

// S_u x S_v normed; at the poles, where S_u vanishes, the normal of the sphere all the same, for
// the primitive and for its NURBS (a collapsed row of poles).
TEST(tests_bssanalysis, normal_direction_on_a_sphere_and_at_its_poles)
{
    const auto s = sphere();
    const auto nurbs = s.to_nurbs();
    const gbs::point<double, 3> c{1., -1., 0.5};
    for (const gbs::Surface<double, 3> *srf : {static_cast<const gbs::Surface<double, 3> *>(&s), static_cast<const gbs::Surface<double, 3> *>(nurbs.get())})
    {
        const auto [u1, u2, v1, v2] = srf->bounds();
        // the side the normals are on: that of a regular point, then everywhere
        const double side = gbs::normal_direction(*srf, (u1 + u2) / 3, (v1 + v2) / 2) * (srf->value((u1 + u2) / 3, (v1 + v2) / 2) - c) > 0. ? 1. : -1.;
        for (double u : {u1, u1 + 0.3 * (u2 - u1), (u1 + u2) / 2, u2})
            for (double v : {v1, v1 + 1e-9 * (v2 - v1), v1 + 0.25 * (v2 - v1), (v1 + v2) / 2, v2})
            {
                const auto n = gbs::normal_direction(*srf, u, v);
                const auto radial = (srf->value(u, v) - c) / 2.;
                ASSERT_NEAR(gbs::norm(n), 1., 1e-14);
                ASSERT_LT(gbs::norm(n - side * radial), 1e-9);
            }
    }
}

TEST(tests_bssanalysis, normal_direction_of_a_plane_and_of_a_point)
{
    const gbs::BSSurface<double, 3> plane{{{0., 0., 1.}, {2., 0., 1.}, {0., 3., 1.}, {2., 3., 1.}}, {0., 0., 1., 1.}, {0., 0., 1., 1.}, 1, 1};
    for (double u : {0., 0.4, 1.})
        for (double v : {0., 0.7, 1.})
            ASSERT_EQ(gbs::normal_direction(plane, u, v), (gbs::point<double, 3>{0., 0., 1.}));
    // a surface reduced to a point has no normal anywhere
    const gbs::BSSurface<double, 3> point{{{1., 2., 3.}, {1., 2., 3.}, {1., 2., 3.}, {1., 2., 3.}}, {0., 0., 1., 1.}, {0., 0., 1., 1.}, 1, 1};
    ASSERT_THROW(gbs::normal_direction(point, 0.5, 0.5), std::domain_error);
    ASSERT_THROW(gbs::normal_direction(plane, std::nan(""), 0.5), std::domain_error);
}

// The apex of a cone has no normal: that of the side towards the centre of the domain, the limit
// along the generatrix through it.
TEST(tests_bssanalysis, normal_direction_at_the_apex_of_a_cone)
{
    const gbs::ConicalSurface<double> cone(0., pi / 6., gbs::Frame3<double>{}, {0., 2. * pi, 0., 1.});
    for (double u : {0., 1., 2.5})
    {
        const auto n = gbs::normal_direction(cone, u, 0.);
        // along the line to the centre (pi, 0.5), u stays u + t (pi - u): the normal there, near the apex
        const auto near = gbs::normal_direction(cone, u + 1e-7 * (pi - u), 1e-7 * 0.5);
        ASSERT_LT(gbs::norm(n - near), 1e-6);
    }
}

// The rule does not depend on the units nor on a scaling of the parameters: a cone 10 mm tall in
// metres and in millimetres, narrow cones at their apex, a cylinder of 2 km radius in millimetres
// and a strip 1e7 times longer than wide (where S_u and S_v differ by more than 1e6) have a normal.
TEST(tests_bssanalysis, normal_direction_whatever_the_units)
{
    const gbs::ConicalSurface<double> metres(0., pi / 6., gbs::Frame3<double>{}, {0., 2. * pi, 0., 0.01});
    const gbs::ConicalSurface<double> millimetres(0., pi / 6., gbs::Frame3<double>{}, {0., 2. * pi, 0., 10.});
    ASSERT_LT(gbs::norm(gbs::normal_direction(metres, 1., 0.) - gbs::normal_direction(millimetres, 1., 0.)), 1e-9);
    for (double deg : {1., 0.5, 0.1})
    {
        const gbs::ConicalSurface<double> cone(0., deg * pi / 180., gbs::Frame3<double>{}, {0., 2. * pi, 0., 1.});
        const auto n = gbs::normal_direction(cone, 1., 0.);
        // the generatrix at u = 1 + t (pi - 1), near the apex: its normal
        const gbs::point<double, 3> expected{-std::cos(1.) * std::cos(deg * pi / 180.), -std::sin(1.) * std::cos(deg * pi / 180.), std::sin(deg * pi / 180.)};
        ASSERT_LT(std::min(gbs::norm(n - expected), gbs::norm(n + expected)), 1e-6);
    }
    const gbs::CylindricalSurface<double> cylinder(2e6, gbs::Frame3<double>{}, {0., 2. * pi, 0., 1.});
    const auto n = gbs::normal_direction(cylinder, 1., 0.5);
    ASSERT_LT(std::min(gbs::norm(n - gbs::point<double, 3>{std::cos(1.), std::sin(1.), 0.}), gbs::norm(n + gbs::point<double, 3>{std::cos(1.), std::sin(1.), 0.})), 1e-12);
    const gbs::BSSurface<double, 3> strip{{{0., 0., 0.}, {1e4, 0., 0.}, {0., 1e-3, 0.}, {1e4, 1e-3, 0.}}, {0., 0., 1., 1.}, {0., 0., 1., 1.}, 1, 1};
    for (double u : {0., 0.5, 1.})
        ASSERT_EQ(gbs::normal_direction(strip, u, 0.5), (gbs::point<double, 3>{0., 0., 1.}));
}

// A double cone, its apex at v = 0 in the middle of the domain: the way to the centre runs along
// the apex, the way to a corner leaves it (the normal of the nappe it goes to).
TEST(tests_bssanalysis, normal_direction_at_the_apex_of_a_double_cone)
{
    const gbs::ConicalSurface<double> cone(0., pi / 6., gbs::Frame3<double>{}, {0., 2. * pi, -1., 1.});
    const auto n = gbs::normal_direction(cone, 1., 0.);
    // towards the corner (0, -1): u = 1 - t, v = -t
    ASSERT_LT(gbs::norm(n - gbs::normal_direction(cone, 1. - 1e-7, -1e-7)), 1e-6);
}

// A surface that throws at a node: its exception is rethrown out of discretize, in the parallel
// loop too (it terminated the process).
TEST(tests_bssanalysis, discretize_rethrows)
{
    struct Throwing : gbs::Surface<double, 3>
    {
        auto value(double u, double v, size_t du = 0, size_t dv = 0) const -> gbs::point<double, 3> override
        {
            if (u > 0.5 && v > 0.5)
                throw std::runtime_error("no value");
            return du + dv == 0 ? gbs::point<double, 3>{u, v, 0.} : gbs::point<double, 3>{du == 1 && dv == 0 ? 1. : 0., du == 0 && dv == 1 ? 1. : 0., 0.};
        }
        auto bounds() const -> std::array<double, 4> override { return {0., 1., 0., 1.}; }
    } throwing;
    ASSERT_THROW(gbs::discretize(throwing, 3, 3), std::runtime_error);
    ASSERT_THROW(gbs::discretize(throwing, 40, 40), std::runtime_error);
    ASSERT_THROW(gbs::discretize_with_normals(throwing, 40, 40), std::runtime_error);
}

// The same points as discretize, the normals of normal_direction at their parameters; past the
// size of the parallel loops too, where a surface without normal throws instead of terminating.
TEST(tests_bssanalysis, discretize_with_normals)
{
    const auto s = sphere();
    for (std::size_t n : {5u, 60u})
    {
        const auto [points, normals] = gbs::discretize_with_normals(s, n, n + 1);
        ASSERT_EQ(points, gbs::discretize(s, n, n + 1));
        ASSERT_EQ(normals.size(), n * (n + 1));
        const auto u = gbs::make_range(0., 2. * pi, n), v = gbs::make_range(-pi / 2., pi / 2., n + 1);
        for (std::size_t j = 0; j <= n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                ASSERT_EQ(normals[j * n + i], gbs::normal_direction(s, u[i], v[j]));
    }
    const gbs::BSSurface<double, 3> point{{{1., 2., 3.}, {1., 2., 3.}, {1., 2., 3.}, {1., 2., 3.}}, {0., 0., 1., 1.}, {0., 0., 1., 1.}, 1, 1};
    ASSERT_THROW(gbs::discretize_with_normals(point, 3, 3), std::domain_error);
    ASSERT_THROW(gbs::discretize_with_normals(point, 100, 100), std::domain_error);
}
