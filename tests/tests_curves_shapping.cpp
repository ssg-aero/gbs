#include <gtest/gtest.h>
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

    std::vector<bsc_constraint<T,dim>> constraints1{
        {u,points[1],0},
        {u,crv(u,1),1},
    };

    auto new_crv1 = moved_to_constraints(crv, constraints1);
    for(auto &cstr : constraints1)
        ASSERT_NEAR( distance(std::get<1>(cstr), new_crv1(std::get<0>(cstr),std::get<2>(cstr))), 0., 1e-9 );

    std::vector<bsc_constraint<T,dim>> constraints2{
        {u,crv(u),0},
        {u,{0.3,0.,0.3},1}
    };
    auto new_crv2 = moved_to_constraints(crv, constraints2);
    for(auto &cstr : constraints2)
        ASSERT_NEAR( distance(std::get<1>(cstr), new_crv2(std::get<0>(cstr),std::get<2>(cstr))), 0., 1e-9 );


    std::vector<bsc_constraint<T,dim>> constraints3;
    int np = poles.size();
    for(int i = 0; i < np; i++)
    {
        constraints3.push_back(
            {i*5/(np-1.), {0.,i/(np-1.),2*i/(np-1.)}, 0}
        );
    }
    auto new_crv3 = moved_to_constraints(crv, constraints3);
    for(auto &cstr : constraints3)
        ASSERT_NEAR( distance(std::get<1>(cstr), new_crv3(std::get<0>(cstr),std::get<2>(cstr))), 0., 1e-9 );
    points_vector<T,3> pts_3(constraints3.size());
    std::transform(
        constraints3.begin(), constraints3.end(),
        pts_3.begin(),[](const auto &cstr){return std::get<1>(cstr);}
    );


    std::vector<bsc_constraint<T,dim>> constraints4;
    np = poles.size()*3;
    for(int i = 0; i < np; i++)
    {
        constraints4.push_back(
            {i*5/(np-1.), {0.,2*i/(np-1.),i/(np-1.)}, 0}
        );
    }
    auto new_crv4 = moved_to_constraints(crv, constraints4);
    for(auto &cstr : constraints4)
        ASSERT_NEAR( distance(std::get<1>(cstr), new_crv4(std::get<0>(cstr),std::get<2>(cstr))), 0., 1e-9 );
    points_vector<T,3> pts_4(constraints4.size());
    std::transform(
        constraints4.begin(), constraints4.end(),
        pts_4.begin(),[](const auto &cstr){return std::get<1>(cstr);}
    ); 

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
                .c = &new_crv1,
                .col_crv = {0., 0., 1.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },
            gbs::crv_dsp<T, 3, false>{
                .c = &new_crv2,
                .col_crv = {1., 0., 1.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },

            gbs::crv_dsp<T, 3, false>{
                .c = &new_crv3,
                .col_crv = {0., 1., 1.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },
            gbs::crv_dsp<T, 3, false>{
                .c = &new_crv4,
                .col_crv = {0., 1., 0.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },
            gbs::crv_dsp<T, 3, false>{
                .c = &seg,
                .col_crv = {0., 1., 0.},
            },
            points,
            pts_3,
            pts_4
            );
    }
}