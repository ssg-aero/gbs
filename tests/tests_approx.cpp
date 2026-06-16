#include <doctest_gtest.hpp>
#include <gbs/bssapprox.h>
#include <gbs-render/vtkGbsRender.h>
#include <numbers>
#include <cmath>
namespace{
    const bool PLOT_ON=false;
}
using gbs::operator-;
TEST(tests_approx, surf_on_grid)
{
    size_t nu = 30;
    size_t nv = 25;
    std::vector<std::array<double,3> > points(nu*nv);
    std::vector<double> u(nu);
    std::vector<double> v(nv);
    auto x{0.};
    auto y{0.};
    auto z{0.};
    for (auto i = 0; i < nu; i++) u[i] = i / (nu -1.);
    for (auto j = 0; j < nv; j++) v[j] = j / (nv -1.);
    for (auto j = 0; j < nv; j++)
    {
        y = v[j] - 0.5;
        for (auto i = 0; i < nu; i++)
        {
            x = u[i] - 0.5;
            z = std::sin(std::numbers::pi * (std::sqrt(x*x+y*y)) );
            points[i+j*nu] = {x,y,z};
        }
    }

    size_t n_poles_u = 15;
    size_t n_poles_v = 10;
    size_t deg_u = 5;
    size_t deg_v = 4;
    auto ku = gbs::build_simple_mult_flat_knots(0.,1.,n_poles_u,deg_u);
    auto kv = gbs::build_simple_mult_flat_knots(0.,1.,n_poles_v,deg_v);
    auto poles = gbs::approx(points,ku,kv,u,v,deg_u,deg_v);
    gbs::BSSurface<double,3> srf(poles,ku,kv,deg_u,deg_v);

    auto d_avr {0.};
    for (auto j = 0; j < nv; j++)
    {
        for (auto i = 0; i < nu; i++)
        {
     
            d_avr+= gbs::norm( points[i+j*nu] - srf(u[i],v[j]) );
        }
    }
    d_avr/=(nu*nv);
    ASSERT_LT(d_avr,0.005);

    if (PLOT_ON)
        gbs::plot(
            points,
            srf
        );

}

// C3 (surface side): the grid approximation entry rejects ill-posed pole counts
// consistently with the curve entries: too few poles per direction (< deg+1) and
// more total poles than data points both throw.
TEST(tests_approx, surf_approx_rejects_illposed_pole_counts_C3)
{
    size_t deg_u = 3;
    size_t deg_v = 3;

    // a small valid grid of data points
    size_t nu = 6, nv = 6;
    std::vector<double> u(nu), v(nv);
    for (size_t i = 0; i < nu; ++i) u[i] = i / (nu - 1.);
    for (size_t j = 0; j < nv; ++j) v[j] = j / (nv - 1.);
    std::vector<std::array<double, 3>> points(nu * nv);
    for (size_t j = 0; j < nv; ++j)
        for (size_t i = 0; i < nu; ++i)
            points[i + j * nu] = {u[i], v[j], std::sin(u[i] * v[j])};

    // n_poles_u < deg_u + 1: a degenerate (too short) knot vector in u -> rejected.
    // {0,0,0,1,1,1,1} has 7 knots -> n_poles_u = 7 - (3+1) = 3 < 4.
    std::vector<double> ku_bad = {0., 0., 0., 1., 1., 1., 1.};
    auto kv_ok = gbs::build_simple_mult_flat_knots(0., 1., size_t{4}, deg_v);
    ASSERT_THROW(gbs::approx(points, ku_bad, kv_ok, u, v, deg_u, deg_v), std::length_error);

    // n_poles > n_pts: many poles, far fewer data points -> rejected.
    auto ku_big = gbs::build_simple_mult_flat_knots(0., 1., size_t{15}, deg_u);
    auto kv_big = gbs::build_simple_mult_flat_knots(0., 1., size_t{10}, deg_v); // 150 poles > 36 points
    ASSERT_THROW(gbs::approx(points, ku_big, kv_big, u, v, deg_u, deg_v), std::length_error);

    // a well-posed configuration still succeeds.
    auto ku_ok = gbs::build_simple_mult_flat_knots(0., 1., size_t{4}, deg_u);
    ASSERT_NO_THROW(gbs::approx(points, ku_ok, kv_ok, u, v, deg_u, deg_v));
}