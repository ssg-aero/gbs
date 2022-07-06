#include <gtest/gtest.h>
#include <gbs/bssapprox.h>
#include <gbs/render/vtkcurvesrender.h>
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