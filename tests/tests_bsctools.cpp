#include <gtest/gtest.h>
#include <gbs/bsctools.h>
#include <gbs/bscbuild.h>
#include <gbs/render/vtkcurvesrender.h>

using gbs::operator-;
using gbs::operator+;
TEST(tests_bsctools, trim)
{
    std::vector<double> k = {0., 0., 0., 1. / 4., 1. / 4., 1. / 2., 1. / 2., 3. / 4., 3. / 4., 1., 1., 1.};
    auto wi = sqrt(2.) / 2.;
    auto r = 1.23;
    std::vector<std::array<double, 4>> poles =
        {
            {r, 0., 0., 1.},
            {wi * r, wi * r, 0., wi},
            {0., r, 0., 1.},
            {-wi * r, wi * r, 0., wi},
            {-r, 0., 0., 1.},
            {-wi * r, -wi * r, 0., wi},
            {0., -r, 0., 1.},
            {wi * r, -wi * r, 0., wi},
            {r, 0., 0., 1.},
        };

    size_t p = 2;


    auto c1_3d_dp_ref = gbs::BSCurveRational<double, 3>(poles, k, p);

    gbs::trim(p, k, poles, 0., 3. / 8.);
    gbs::trim(p, k, poles, 0.1, 3. / 8.);

    auto c1_3d_dp = gbs::BSCurveRational<double, 3>(poles, k, p);

    ASSERT_TRUE(gbs::distance(c1_3d_dp.begin(),c1_3d_dp_ref(0.1)) < 1e-6);
    ASSERT_TRUE(gbs::distance(c1_3d_dp.end(),c1_3d_dp_ref(3. / 8.)) < 1e-6);

}

TEST(tests_bsctools, unify_curves_knots)
{

    std::vector<double> k1 = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<double> k2 = {0., 0., 0.1, 1.3, 1.8, 2.1, 2.6, 3., 3.};
    std::vector<std::array<double, 3>> poles =
        {
            {0., 0., 0.},
            {0., 1., 0.},
            {1., 1., 0.},
            {1., 1., 1.},
            {1., 1., 2.},
            {3., 1., 1.},
            {0., 4., 1.},
        };
    size_t p1 = 2;
    size_t p2 = 1;

    gbs::BSCurve3d_d c1(poles, k1, p1);
    gbs::BSCurve3d_d c2(poles, k2, p2);

    std::list<gbs::BSCurve3d_d> bsc_lst = {c1, c2};

    ASSERT_THROW(gbs::unify_curves_knots(bsc_lst), std::exception);

    gbs::unify_curves_degree(bsc_lst);

    gbs::unify_curves_knots(bsc_lst);

    auto k1_ = bsc_lst.front().knotsFlats();
    auto k2_ = bsc_lst.back().knotsFlats();

    ASSERT_TRUE(
        std::equal(
            k1_.begin(),
            k1_.end(),
            k2_.begin(),
            [](const auto &k1, const auto &k2) {
                return fabs(k1 - k2) < gbs::knot_eps;
            }));
}

TEST(tests_bsctools, join)
{

    std::vector<double> k1 = {0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1.};
    std::vector<double> k2 = {0., 0., 0., 0.33, 0.66, 1., 1., 1.};
    std::vector<std::array<double, 4>> poles1 =
        {
            {0., 1., 0., 1},
            {1., 1., 0., 1.5},
            {1., 2., 0., 1.},
            {2., 3., 0., 0.5},
            {3., 3., 0., 1.2},
            {4., 2., 0., 1.},
        };
    // gbs::scale_weights(poles1);
    // std::cerr << poles1.back().back() << std::endl;
    // gbs::scale_poles(poles1,1. / poles1.back().back());
    // std::cerr << poles1.back().back() << std::endl;

    // std::vector<std::array<double,3> > poles1 =
    // {
    //     {0.,1.,0.},
    //     {0.,1.,0.},
    //     {1.,2.,0.},
    //     {2.,3.,0.},
    //     {3.,3.,0.},
    //     {4.,2.,0.},
    // };

    std::vector<std::array<double, 3>> poles2 =
        {
            {4., 2., 0.},
            {5., 1., 0.},
            {5., 0., 0.},
            {4., -1., 0.},
            {3., -1., 0.},
        };
    size_t p1 = 2;
    size_t p2 = 2;

    gbs::BSCurveRational3d_d c1(poles1, k1, p1);
    gbs::BSCurve3d_d c2(poles2, k2, p2);
    auto c3 = gbs::join(&c1, &c2);
    auto c4 = gbs::join(c1, c2);
    // ASSERT_TRUE(std::is_sorted(k1.begin(), k1.end(), std::less_equal<double>()));

    ASSERT_LT(gbs::norm(c1.begin() - c3->begin()), 1e-6);
    ASSERT_LT(gbs::norm(c1.end() - c3->value(1.)), 1e-6);
    ASSERT_LT(gbs::norm(c2.begin() - c3->value(1.)), 1e-6);
    ASSERT_LT(gbs::norm(c2.end() - c3->end()), 1e-6);

    gbs::plot(
        c1, c2, *c3);
}

TEST(tests_bsctools, cn_connect)
{
    std::vector<double> k1 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<double> k2 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double, 3>> poles1 =
        {
            {0., 1., 0.},
            {1., 2., 0.},
            {2., 2., 0.},
            {3., 0., 0.},
        };

    std::vector<std::array<double, 3>> poles2 =
        {
            {0., 0.5, 0.},
            {1., 1., 0.},
            {2., 1.3, 0.},
            {3., 0., 0.},
        };
    size_t p1 = 3;
    size_t p2 = 3;

    gbs::BSCurve3d_d c1(poles1, k1, p1);
    gbs::BSCurve3d_d c2(poles2, k2, p2);
    c1.reverse();

    auto c3 = gbs::c2_connect(c1, c2);
    auto c4 = gbs::c3_connect(c1, c2);

    gbs::plot(
        c1, c2, c4);

    for (int i = 0; i < 3; i++)
    {
        ASSERT_LT(gbs::norm(c1.end(i) - c3.begin(i)), 1e-6);
        ASSERT_LT(gbs::norm(c2.begin(i) - c3.end(i)), 1e-6);
        ASSERT_LT(gbs::norm(c1.end(i) - c3.begin(i)), 1e-6);
        ASSERT_LT(gbs::norm(c2.begin(i) - c3.end(i)), 1e-6);
    }

    for (int i = 0; i < 4; i++)
    {
        ASSERT_LT(gbs::norm(c1.end(i) - c4.begin(i)), 1e-6);
        ASSERT_LT(gbs::norm(c2.begin(i) - c4.end(i)), 1e-6);
        ASSERT_LT(gbs::norm(c1.end(i) - c4.begin(i)), 1e-6);
        ASSERT_LT(gbs::norm(c2.begin(i) - c4.end(i)), 1e-6);
    }
}

TEST(tests_bsctools, c2_connect_2d)
{
    // auto crv1 = gbs::build_segment<double, 2>({1., 1.0}, {0., 1.});
    // auto crv2 = gbs::build_segment<double, 2>({0., -1.}, {1., -1.});
    auto crv1 = gbs::BSCurve<double,2>(
        {
            {2., 1.5},
            // {1.2, 1.5},
            {1., 1.5}, 
            {0., 1.},
        },
        {0.,0.,0.,1.,1.,1.},
        2
    );
    auto crv2 = gbs::BSCurve<double,2>(
        {
            {0.2, -1.},
            {1., -1}, 
            {2., -1.5},
        },
        {0.,0.,0.,1.,1.,1.},
        2
    );

    auto c1 = gbs::c2_connect(crv1, crv2, 1.);
    auto c2 = gbs::c2_connect(crv1, crv2, 2.);
    auto c3 = gbs::c2_connect(crv1,crv2, 3.);

    auto c2_1 = gbs::c2_connect(crv1, crv2,crv1.bounds()[0]+0.1,crv2.bounds()[0],true,false, 2.);

    gbs::plot(
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv1),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
            .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv2),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
            .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c1),
            .col_crv = {1,0,0},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c2),
            .col_crv = {0,1,0},
            // .col_crv = {1,0,0},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
            // .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c2_1),
            .col_crv = {0,1,1},
            // .col_crv = {1,0,0},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
            // .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c3),
            .col_crv = {0,0,1},
            // .poles_on = false,
            .poles_on = true,
            .line_width=3.,
            .show_curvature=true,
         }
    );

}

TEST(tests_bsctools, c3_connect_2d)
{
    // auto crv1 = gbs::build_segment<double, 2>({1., 1.5}, {0., 1.});
    // auto crv2 = gbs::build_segment<double, 2>({0.2, -1.}, {1., -1.});
    // auto crv1 = gbs::BSCurve<double,2>(
    //     {
    //         {2., 1.8},
    //         {1.2, 1.8},
    //         {1., 1.5}, 
    //         {0., 1.},
    //     },
    //     {0.,0.,0.,0.,1.,1.,1.,1.},
    //     3
    // );
    // auto crv2 = gbs::BSCurve<double,2>(
    //     {
    //         {0.2, -1.},
    //         // {1., -1.5}, 
    //         {1.5, -1.5}, 
    //         {2., -1.5},
    //     },
    //     {0.,0.,0.,1.,1.,1.},
    //     2
    // );

        auto crv1 = gbs::BSCurve<double,2>(
        {
            {2., 1.5},
            // {1.2, 1.5},
            {1., 1.5}, 
            {0., 1.},
        },
        {0.,0.,0.,1.,1.,1.},
        2
    );
    auto crv2 = gbs::BSCurve<double,2>(
        {
            {0.2, -1.},
            {1., -1}, 
            {2., -1.5},
        },
        {0.,0.,0.,1.,1.,1.},
        2
    );

    auto c1 = gbs::c3_connect(crv1, crv2, 1.);
    auto c2 = gbs::c3_connect(crv1, crv2, 2.);
    auto c3 = gbs::c3_connect(crv1,crv2, 3.);

    gbs::plot(
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv1),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
            .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv2),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
            .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c1),
            .col_crv = {1,0,0},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c2),
            .col_crv = {0,1,0},
            // .col_crv = {1,0,0},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
            // .show_curvature=true,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c3),
            .col_crv = {0,0,1},
            // .poles_on = false,
            .poles_on = true,
            .line_width=3.,
            .show_curvature=true,
         }
    );

}

TEST(tests_bsctools, extend_to_point)
{
    using T = double;
    // using T = float;
    std::vector<T> k = {-0.5,-0.5, -0.5, 0.25, 0.5, 0.75, 1., 1., 1.};
    std::vector<std::array<T, 3>> poles =
        {
            {0., 1., 1},
            {1., .9, 1.5},
            {1., 2., 1.},
            {2., 3., 0.5},
            {3., 3., 1.2},
            {4., 2., 2.},
        };
    size_t p = 2;
    gbs::BSCurve<T,3> c1(poles, k, p);

    gbs::point<T,3> pt2 {4.,3.,1.};
    gbs::point<T,3> pt3 {0.,1.3,0.};
    bool natural_end = true;
    auto c2 = extended_to_point(c1,pt2,natural_end);
    auto c3 = extended_to_point(c1,pt3,false,natural_end);
    auto c4 = extended(c1,1.,true,false,natural_end);
    size_t max_cont = 1;
    auto c5 = extended(c1,0.2,false,true,natural_end,max_cont);

    ASSERT_LT(gbs::norm(c2.end()-pt2),1e-8);
    ASSERT_LT(gbs::norm(c2(c1.bounds()[1])-c1.end()),1e-8);
    ASSERT_LT(gbs::norm(c2(c1.bounds()[1],1)-c1.end(1)),1e-8);
    // ASSERT_LT(gbs::norm(c2(c1.bounds()[1],2)-c1.end(2)),1e-8); // p = 2 -> C1 at knots
    ASSERT_LT(gbs::norm(c3.begin()-pt3),1e-8);
    ASSERT_LT(gbs::norm(c3(c1.bounds()[0])-c1.begin()),1e-8);
    ASSERT_LT(gbs::norm(c3(c1.bounds()[0],1)-c1.begin(1)),1e-8);
    ASSERT_LT(gbs::norm(c1.end()-c4.end())-1.0,1e-8);
    ASSERT_LT(gbs::norm(c4(c1.bounds()[0])-c1.begin()),1e-8);
    ASSERT_LT(gbs::norm(c4(c1.bounds()[0],1)-c1.begin(1)),1e-8);
    ASSERT_LT(gbs::norm(c4(c1.bounds()[1])-c1.end()),1e-8);
    ASSERT_LT(gbs::norm(c4(c1.bounds()[1],1)-c1.end(1)),1e-8);
    ASSERT_LT(gbs::norm(c5(c1.bounds()[0])-c1.begin()),1e-8);
    ASSERT_LT(gbs::norm(c5(c1.bounds()[0],1)-c1.begin(1)),1e-8);

    gbs::plot(
        // gbs::crv_dsp<T, 3, false>{
        //     .c = &(c1),
        //     .col_crv = {0,0,0},
        //     // .poles_on = true,
        //     // .line_width=3.,
        //     .show_curvature=true,
        //  },
        //  gbs::crv_dsp<T, 3, false>{
        //     .c = &(c2),
        //     .col_crv = {1,0,0},
        //     // .poles_on = true,
        //     // .line_width=3.,
        //     .show_curvature=true,
        //  },
        c1,
        c2,
        c3,
        c4,
        c5,
        gbs::points_vector<T,3>{pt2,c1.end(),pt3,c1.begin(),c4.end(),c5.begin()}
    );

}

TEST(tests_bsctools, extend_to_point_rational)
{
    using T = double;
    // using T = float;
    std::vector<T> k = {-0.5,-0.5, -0.5,-0.5, 0.25,  0.75, 1., 1.,1., 1.};
    std::vector<std::array<T, 3>> poles =
        {
            {0., 1., 1},
            {1., .9, 1.5},
            {1., 2., 1.},
            {2., 3., 0.5},
            {3., 3., 1.2},
            {4., 2., 2.},
        };
    size_t p = 3;
    gbs::BSCurveRational<T,2> c1(poles, k, p);

    gbs::point<T,2> pt2 {4.4,3.};
    gbs::point<T,2> pt3 {0.,0.};
    bool natural_end = true;
    auto c2 = extended_to_point(c1,pt2,natural_end);
    size_t max_cont = 1;
    auto c3 = extended_to_point(c1,pt3,false,natural_end,max_cont);

    ASSERT_LT(gbs::norm(c2.end()-pt2),1e-8);
    ASSERT_LT(gbs::norm(c2(c1.bounds()[1])-c1.end()),1e-8);
    ASSERT_LT(gbs::norm(c2(c1.bounds()[1],1)-c1.end(1)),1e-8);

    ASSERT_LT(gbs::norm(c3.begin()-pt3),1e-8);
    ASSERT_LT(gbs::norm(c3(c1.bounds()[0])-c1.begin()),1e-8);
    ASSERT_LT(gbs::norm(c3(c1.bounds()[0],1)-c1.begin(1)),1e-8);

    gbs::plot(
        c2,
        c3,
        gbs::points_vector<T,2>{c1.end(),c2.end(),c1.begin(),c3.begin()}
    );
}
