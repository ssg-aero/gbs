#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs/bscbuild.h>
#include <gbs/bscapprox.h>
#include <gbs/vecop.h>
#include <gbs/bscanalysis.h>
#include <gbs-render/vtkcurvesrender.h>

namespace
{
    const bool PLOT_ON = true;
}

using gbs::operator-;

TEST(tests_offset, curve2d_rational_offset)
{
    auto circle1 = gbs::build_circle<double, 2>(1.);
    auto f_offset = gbs::BSCfunction<double>(gbs::build_segment<double, 1>({1.}, {1.}));
    auto f_offset3 = gbs::BSCfunction<double>(gbs::build_segment<double, 1>({1.}, {2.}));
    auto f_offset4 = gbs::BSCfunction<double>(gbs::BSCurve<double, 1>(
        gbs::points_vector<double,1>{{0.},{0.},{1.},{0.},{0.}},
        {0.,0.,0.,0.1,0.5,1.,1.,1.},
        2
        )
    );
    auto p_circle1 = std::make_shared<gbs::BSCurveRational<double, 2>>(circle1);
    gbs::CurveOffset<double, 2> circle2{
        p_circle1,
        f_offset};
    gbs::CurveOffset<double, 2> circle3{
        p_circle1,
        f_offset3};
    gbs::CurveOffset<double, 2> circle4{
        p_circle1,
        f_offset4};
    auto u = gbs::deviation_based_params<double, 2>(circle1, 30, 0.01);
    for (auto u_ : u)
    {
        ASSERT_NEAR(gbs::norm(circle1(u_) - circle2(u_)), 1., 1e-6);
        ASSERT_NEAR(gbs::norm(circle1(u_) - circle3(u_)), f_offset3(u_), 1e-6);
        ASSERT_NEAR(gbs::norm(circle1(u_) - circle4(u_)), f_offset4(u_), 1e-6);
    }

    if (PLOT_ON)
        gbs::plot(
            gbs::crv_dsp<double, 2, true>{
                .c = &(circle1),
                .col_crv = {0., 0., 0.},
                .poles_on = true,
                .line_width = 3.,
            },
            gbs::make_actor(circle2, {1., 0., 0.}),
            gbs::make_actor(circle3, {0., 1., 0.}),
            gbs::make_actor(circle4, {0., 0., 1.}) );
}

TEST(tests_offset, curve2d_offset)
{
    std::string line;
    std::ifstream myfile("../tests/in/e1098.dat");
    // std::ifstream myfile("../tests/in/e817.dat");
    if (myfile.is_open())
    {
        std::vector<std::array<double, 2>> pts;
        getline(myfile, line);
        while (getline(myfile, line))
        {
            std::istringstream iss(line);
            std::string::size_type sz; // alias of size_t

            double x = std::stod(line, &sz);
            double y = std::stod(line.substr(sz));
            pts.push_back({x, y});
        }
        myfile.close();

        const auto d_offset {0.2};
        auto crv = gbs::approx(pts, 5, gbs::KnotsCalcMode::CHORD_LENGTH, true);
        auto f_offset = gbs::BSCfunction<double>(gbs::build_segment<double, 1>({d_offset}, {d_offset}));
        auto p_crv = std::make_shared<gbs::BSCurveRational<double, 2>>(crv);
        gbs::CurveOffset<double, 2> crv_offset{
            p_crv,
            f_offset};

        auto u = gbs::deviation_based_params<double,2>(crv,30,0.01);
        for(auto u_ :u)
        {
            ASSERT_NEAR(gbs::norm(crv(u_)-crv_offset(u_)),d_offset,1e-6);
        }

        if (PLOT_ON)
            gbs::plot(
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv),
            .col_crv = {0,0,1},
            // .poles_on = false,
            .poles_on = true,
            .line_width=3.,
            .show_curvature=false,
         }
                ,
                crv_offset);
    }
}