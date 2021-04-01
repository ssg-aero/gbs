#include <gtest/gtest.h>
#include <gbs/bsctools.h>
#include <gbs/bscbuild.h>
#include <gbs-render/vtkcurvesrender.h>

#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/export.h>

using gbs::operator-;
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

    // gbs::trim(p,k,poles,0.1,3./8.);
    gbs::trim(p, k, poles, 0., 3. / 8.);
    gbs::trim(p, k, poles, 0.1, 3. / 8.);

    auto c1_3d_dp = gbs::BSCurveRational<double, 3>(poles, k, p);

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(occt_utils::NURBSplineCurve(c1_3d_dp));

    occt_utils::to_iges(crv_lst, "trim.igs");
}

TEST(tests_bsctools, unify_knots)
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

    ASSERT_THROW(gbs::unify_knots(bsc_lst), std::exception);

    gbs::unify_degree(bsc_lst);

    gbs::unify_knots(bsc_lst);

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

    std::vector<Handle(Geom_Geometry)> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(c1));
    crv_lst.push_back(occt_utils::BSplineCurve(c2));
    crv_lst.push_back(occt_utils::BSplineCurve(bsc_lst.front()));
    crv_lst.push_back(occt_utils::BSplineCurve(bsc_lst.back()));

    occt_utils::to_iges(crv_lst, "unify_knots.igs");
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
    ASSERT_TRUE(gbs::check_curve(c1.poles(), c1.knotsFlats(), c1.degree()));
    ASSERT_TRUE(gbs::check_curve(c3->poles(), c3->knotsFlats(), c3->degree()));

    ASSERT_LT(gbs::norm(c1.begin() - c3->begin()), 1e-6);
    ASSERT_LT(gbs::norm(c1.end() - c3->value(1.)), 1e-6);
    ASSERT_LT(gbs::norm(c2.begin() - c3->value(1.)), 1e-6);
    ASSERT_LT(gbs::norm(c2.end() - c3->end()), 1e-6);

    std::vector<Handle(Geom_Geometry)> crv_lst;
    crv_lst.push_back(occt_utils::NURBSplineCurve(c1));
    crv_lst.push_back(occt_utils::BSplineCurve(c2));
    crv_lst.push_back(occt_utils::NURBSplineCurve(*c3));

    occt_utils::to_iges(crv_lst, "join.igs");

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
    std::vector<Handle(Geom_Geometry)> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(c1));
    crv_lst.push_back(occt_utils::BSplineCurve(c2));
    crv_lst.push_back(occt_utils::BSplineCurve(c3));
    crv_lst.push_back(occt_utils::BSplineCurve(c4));

    occt_utils::to_iges(crv_lst, "cn_connect.igs");
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

    gbs::plot(
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv1),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv2),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
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
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c3),
            .col_crv = {0,0,1},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         }
    );

    std::vector<Handle(Geom_Geometry)> occt_crv_lst;
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c1,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c2,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c3,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(crv1,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(crv2,0.)));

    occt_utils::to_iges(occt_crv_lst, "c2_connect_2d.igs");
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
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(crv2),
            .col_crv = {0,0,0},
            .poles_on = false,
            .line_width=3.,
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
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         },
        gbs::crv_dsp<double, 2, false>{
            .c = &(c3),
            .col_crv = {0,0,1},
            // .poles_on = false,
            .poles_on = true,
            .line_width=1.,
         }
    );

    std::vector<Handle(Geom_Geometry)> occt_crv_lst;
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c1,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c2,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(c3,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(crv1,0.)));
    occt_crv_lst.push_back(occt_utils::BSplineCurve(gbs::add_dimension(crv2,0.)));

    occt_utils::to_iges(occt_crv_lst, "c3_connect_2d.igs");
}