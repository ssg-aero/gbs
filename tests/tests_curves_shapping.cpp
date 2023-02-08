#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-render/vtkcurvesrender.h>
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

TEST(curve_shaping, to_constraints)
{
    const size_t dim{3};
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T,dim> > poles =
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
    std::vector<std::array<T,dim> > points{crv(u),crv(u)+deltap};

    std::vector<bsc_constraint<T,dim>> constraints{
        {u,points[1],0},
        {u,crv(u,1),1},
    };

    auto new_crv = moved_to_constraints(crv, constraints);
    for(auto &cstr : constraints)
        ASSERT_NEAR( distance(std::get<1>(cstr), new_crv(std::get<0>(cstr),std::get<2>(cstr))), 0., 1e-9 );

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