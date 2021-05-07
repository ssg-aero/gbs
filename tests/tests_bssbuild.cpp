#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bssbuild.h>
#include <gbs-render/vtkcurvesrender.h>

const double tol = 1e-6;

using gbs::operator-;

TEST(tests_bssbuild, has_nurbs)
{
    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_3d_d poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {3.,1.,0.},
        {0.,4.,0.},
    };
    gbs::points_vector<double,4> poles2 =
    {
        {0.,0.,0.,1.1},
        {0.,1.,0.,1.1},
        {1.,1.,0.,1.2},
        {1.,1.,0.,1.3},
        {1.,1.,0.,1.},
        {3.,1.,0.,1.},
        {0.,4.,0.,1.},
    };
    gbs::BSCurve3d_d         c1(poles1,k,p);
    gbs::BSCurveRational3d_d c2(poles2,k,p);

    std::list<gbs::Curve<double, 3> *> bs_lst = {&c1, &c2};

    ASSERT_TRUE(gbs::has_nurbs(bs_lst));
    bs_lst.pop_back();
    ASSERT_FALSE(gbs::has_nurbs(bs_lst));
}

TEST(tests_bssbuild, loft)
{
    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_3d_d poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {3.,1.,0.},
        {0.,4.,0.},
    };
    gbs::points_vector_3d_d poles2 =
    {
        {0.,0.,1.},
        {0.,1.,1.},
        {1.,1.,1.},
        {1.,1.,1.},
        {1.,1.,1.},
        {3.,1.,1.},
        {2.,4.,1.},
    };
    gbs::points_vector_3d_d poles3 =
    {
        {0.,0.,2.},
        {0.,1.,2.},
        {1.,1.,2.},
        {1.,1.,2.},
        {1.,1.,2.},
        {3.,1.,2.},
        {1.,4.,2.},
    };
    gbs::BSCurve3d_d c1(poles1,k,p);
    gbs::BSCurve3d_d c2(poles2,k,p);
    gbs::BSCurve3d_d c3(poles3,k,p);

    std::list<gbs::BSCurveGeneral<double,3,false>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst );
    std::list<gbs::BSCurve3d_d> bs_lst2 = {c1,c2,c3};
    auto s2 = gbs::loft( bs_lst2 );
    gbs::plot(s,c1,c2,c3);
}

TEST(tests_bssbuild, loft_rational)
{
    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_4d_d poles1 =
    {
        {0.,0.,0.,1.},
        {0.,1.,0.,1.1},
        {1.,1.,0.,1.},
        {1.,1.,0.0,1.},
        {2.,1.,0.0,1.},
        {3.,1.,0.,1.1},
        {0.,4.,0.,1.},
    };
    gbs::points_vector_4d_d poles2 =
    {
        {0.,0.,1.,1.},
        {0.,1.,1.3,1.},
        {1.,1.,1.2,1.},
        {1.,2.,1.4,1.},
        {2.,1.,1.3,1.2},
        {3.,1.,1.1,1.},
        {2.,4.,1.,1.},
    };
    gbs::points_vector_4d_d poles3 =
    {
        {0.,0.,2.,1.},
        {0.,1.,2.2,1.},
        {1.,1.,2.2,1.2},
        {1.,2.,2.3,1.},
        {2.,1.,2.1,1.},
        {3.,1.,2.,1.},
        {1.,4.,2.,1.},
    };
    gbs::BSCurveRational3d_d c1(poles1,k,p);
    gbs::BSCurveRational3d_d c2(poles2,k,p);
    gbs::BSCurveRational3d_d c3(poles3,k,p);

    std::list<gbs::BSCurveGeneral<double,3,true>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst );
    gbs::plot(s,c1,c2,c3);
}

TEST(tests_bssbuild, loft_with_spine)
{
    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_3d_d poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {0.7,1.,0.},
        {1.,1.3,0.},
        {1.8,1.,0.3},
        {3.,1.,0.},
        {0.,4.,0.},
    };
    gbs::points_vector_3d_d poles2 =
    {
        {0.,0.,1.},
        {0.,1.,1.},
        {1.,1.,1.},
        {1.3,0.4,1.},
        {1.5,0.5,1.},
        {3.,1.,1.5},
        {2.,4.,1.},
    };
    gbs::points_vector_3d_d poles3 =
    {
        {0.,0.,2.},
        {0.,1.,2.},
        {0.5,1.,2.},
        {1.,1.,2.5},
        {1.5,1.,2.5},
        {3.,1.,2.},
        {1.,4.,2.},
    };
    gbs::BSCurve3d_d c1(poles1,k,p);
    gbs::BSCurve3d_d c2(poles2,k,p);
    gbs::BSCurve3d_d c3(poles3,k,p);
    gbs::BSCurve3d_d sp({{1.5,1.5,0.},{1.5,1.5,3.}},{0.,0.,3.,3.},1);

    std::list<gbs::BSCurveGeneral<double,3,false>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst, sp );
    gbs::plot(s,c1,c2,c3,sp);
    // gbs::plot(s);
}

// TEST(tests_bssbuild, loft_rational_with_spine)
// {
//     size_t p = 2;
//     std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
//     gbs::points_vector_4d_d poles1 =
//     {
//         {0.,0.,0.,1.},
//         {0.,1.,0.,1.1},
//         {1.,1.,0.,1.},
//         {1.,1.,0.0,1.},
//         {2.,1.,0.0,1.},
//         {3.,1.,0.,1.1},
//         {0.,4.,0.,1.},
//     };
//     gbs::points_vector_4d_d poles2 =
//     {
//         {0.,0.,1.,1.},
//         {0.,1.,1.3,1.},
//         {1.,1.,1.2,1.},
//         {1.,2.,1.4,1.},
//         {2.,1.,1.3,1.2},
//         {3.,1.,1.1,1.},
//         {2.,4.,1.,1.},
//     };
//     gbs::points_vector_4d_d poles3 =
//     {
//         {0.,0.,2.,1.},
//         {0.,1.,2.2,1.},
//         {1.,1.,2.2,1.2},
//         {1.,2.,2.3,1.},
//         {2.,1.,2.1,1.},
//         {3.,1.,2.,1.},
//         {1.,4.,2.,1.},
//     };
//     gbs::BSCurveRational3d_d c1(poles1,k,p);
//     gbs::BSCurveRational3d_d c2(poles2,k,p);
//     gbs::BSCurveRational3d_d c3(poles3,k,p);
//     gbs::BSCurve3d_d sp({{1.5,1.5,0.},{1.5,1.5,3.}},{0.,0.,1.,1.},1);

//     std::list<gbs::BSCurveGeneral<double,3,true>*> bs_lst = {&c1,&c2,&c3};
//     auto s = gbs::loft( bs_lst , sp);
//     gbs::plot(s,c1,c2,c3);
// }