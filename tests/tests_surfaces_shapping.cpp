#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs-render/vtkcurvesrender.h>
#include <gbs/bscbuild.h>
#include <gbs/bssshaping.h>
using namespace gbs;
using T = double;
const bool PLOT_ON = true;

TEST(surface_shaping, move_to_point)
{
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    gbs::points_vector<double,3> poles = 
    {                                       // ----U----
        {0,0,0},{1,0,0},                    // |
        {0,1,0},{1,1,1},                    // V
        {0,2,0},{1,2,0},                    // |
    };

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q);
    double u = 0.3;
    double v = 0.5;
    point<double,3> D{0, 0, 0.5};
    std::vector<std::array<T,3> > points{srf(u, v),srf(u, v)+D};
    auto srf1 = moved_to_point(srf, points[1] , u, v);
    ASSERT_NEAR( distance(points[1], srf1(u,v)), 0., 1e-9 );
    if(PLOT_ON)
    {
        auto seg = build_segment(points[0], points[1]);
        plot(
            srf,
            srf1,
            gbs::crv_dsp<T, 3, false>{
                .c = &seg,
                .col_crv = {0., 1., 0.},
            },
            points
        );
    }
}

TEST(surface_shaping, move_to_constraints)
{
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    gbs::points_vector<double,3> poles = 
    {                                       // ----U----
        {0,0,0},{1,0,0},                    // |
        {0,1,0},{1,1,1},                    // V
        {0,2,0},{1,2,0},                    // |
    };

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q);

    double u = 0.3;
    double v = 0.5;
    point<double,3> D{0, 0, 0.5};
    std::vector<std::array<T,3> > points{srf(u, v),srf(u, v)+D};

    std::vector<bss_constraint<double,3>> constraints1{
        {u,v,points[1],0,0},
        {u,v,srf(u,v,1,0),1,0},
        {u,v,srf(u,v,0,1),0,1}
    };

    auto srf1 = moved_to_constraints(srf, constraints1);

    // for(auto &cstr : constraints1)
        // ASSERT_NEAR( distance(std::get<2>(cstr), srf1(std::get<0>(cstr),std::get<1>(cstr),std::get<3>(cstr),std::get<4>(cstr))), 0., 1e-9 );

    srf.increaseDegreeU();
    srf.increaseDegreeU();
    srf.increaseDegreeV();
    srf.increaseDegreeV();
    srf.increaseDegreeU();
    srf.increaseDegreeV();
    srf.insertKnotU(0.5,4);
    srf.insertKnotV(0.5,4);

    auto srf2 =  moved_to_point(srf, points[1] , u, v);

    if(PLOT_ON)
    {
        auto seg = build_segment(points[0], points[1]);
        plot(
            srf,
            srf1,
            srf2,
            gbs::crv_dsp<T, 3, false>{
                .c = &seg,
                .col_crv = {0., 1., 0.},
            },
            points
        );
    }
}