#include <doctest_gtest.hpp>
#include <gbs/curves>
#include <gbs-render/vtkGbsRender.h>
#include <gbs/bscbuild.h>
#include <gbs/bscshaping.h>
using namespace gbs;
using T = double;
const bool PLOT_ON = false;

TEST(curve_shaping, move_to_point)
{
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    T u = 2.3;

    auto crv = BSCurve<T,3>(poles,k,p);
    point<T,3> deltap{0.1,0.2,0.3}; 
    std::vector<std::array<T,3> > points{crv(u),crv(u)+deltap};

    auto new_crv = moved_to_point(crv, points[1], u);
    ASSERT_NEAR( distance(points[1], new_crv(u)), 0., 1e-9 );

    if(PLOT_ON)
    {
        auto seg = build_segment(points[0], points[1]);
        plot(
            gbs::crv_dsp<T, 3, false>{
                .c = &crv,
                .col_crv = {1., 0., 0.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },
            gbs::crv_dsp<T, 3, false>{
                .c = &new_crv,
                .col_crv = {0., 0., 1.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },

            gbs::crv_dsp<T, 3, false>{
                .c = &seg,
                .col_crv = {0., 1., 0.},
            },
            points);
    }
}

TEST(curve_shaping, singular_constraints_throw)
{
    // two constraints the same make the constraint system singular: it throws instead of
    // moving the poles by meaningless values
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T, 3>> poles{{0., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {1., 1., 1.}, {1., 1., 2.}, {3., 1., 1.}, {0., 4., 1.}};
    const auto crv = BSCurve<T, 3>(poles, k, 2);
    const point<T, 3> d{0.1, 0.2, 0.3};
    // fewer constraints than poles moved (minimal norm), as many (square)
    ASSERT_THROW(moved_to_constraints_delta(crv, std::vector<bsc_constraint<T, 3>>{{2.3, d, 0}, {2.3, d, 0}}, 1, 5), std::invalid_argument);
    ASSERT_THROW(moved_to_constraints_delta(crv, std::vector<bsc_constraint<T, 3>>{{2.3, d, 0}, {2.3, d, 0}, {3.1, d, 0}}, 2, 4), std::invalid_argument);
}
