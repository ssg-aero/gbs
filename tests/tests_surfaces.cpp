#include <gtest/gtest.h>

#include <gbs/surfaces>
#include <gbs/bscbuild.h>

#include <gbs-render/vtkcurvesrender.h>

TEST(tests_surfaces, surface_of_revolution)
{
    auto cir = gbs::build_circle<double,2>(0.2,{0.,1.});
    auto p_cir = std::make_shared<gbs::BSCurveRational<double,2>>(cir);
    gbs::SurfaceOfRevolution<double> sor1 {
        p_cir,
        {{
            {0.,0.,0.},
            {0.,0.,1.},
            {1.,0.,0.}
        }}
    };
    gbs::SurfaceOfRevolution<double> sor2 {
        p_cir,
        {{
            {0.,0.,0.},
            {1.,0.,0.},
            {0.,1.,0.}
        }}
    };
    gbs::SurfaceOfRevolution<double> sor3 {
        p_cir,
        {{
            {0.,0.,0.},
            {0.,1.,0.},
            {1.,0.,0.}
        }}
    };
    gbs::SurfaceOfRevolution<double> sor4 {
        p_cir,
        {{
            {0.,0.,2.},
            {0.,0.,1.},
            {1.,0.,0.}
        }}
        ,
        0.,
        0.5 * std::numbers::pi
    };


    auto Xaxis = gbs::build_segment<double,3>({0.,0.,0.},{1.,0.,0.});
    auto Yaxis = gbs::build_segment<double,3>({0.,0.,0.},{0.,1.,0.});
    auto Zaxis = gbs::build_segment<double,3>({0.,0.,0.},{0.,0.,1.});


    plot(
        sor1,
        sor2,
        sor3,
        sor4,
        gbs::crv_dsp<double,3,false>
        {
            .c = & Xaxis,
            .col_crv = {1.,0.,0.}
        },
        gbs::crv_dsp<double,3,false>
        {
            .c = & Yaxis,
            .col_crv = {0.,1.,0.}
        },
        gbs::crv_dsp<double,3,false>
        {
            .c = & Zaxis,
            .col_crv = {0.,0.,1.}
        }
        
    );
}