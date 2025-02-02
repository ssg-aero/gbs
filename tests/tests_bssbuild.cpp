#include <gbs/bscurve.h>
#include <gbs/bsctools.h>
#include <gbs/bssbuild.h>
#include <gbs-io/fromtext.h>
#include <gbs-io/fromjson.h>
#include <gbs-io/fromjson2.h>
#include <gbs-render/vtkGbsRender.h>
#include <gtest/gtest.h>
#include "tests_helpers.h"


const double tol = 1e-6;

using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

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

TEST(tests_bssbuild, loft2d_sharp)
{
    gbs::BSCurve2d_d crv1{
        {{0.,0.},{0.5,0.5},{1.,0.}},
        {0.,0.5,1.},
        {2,1,2},
        1
    };

    gbs::BSCurve2d_d crv2{
        {{0.,2.},{1.,2.}},
        {0.,1.},
        {2,2},
        1
    };

    auto s = gbs::loft( std::list<gbs::BSCurve2d_d>{crv1, crv2}, 1 );
    if(PLOT_ON)
        gbs::plot(s,crv1, crv2, s.isoU(0.5), s.isoV(0.5));
}

TEST(tests_bssbuild, loft2d_sharp_partial)
{
    gbs::BSCurve2d_d crv1{
        {{0.3611771432346219,0.10511112605663975},{0.5,0.135},{0.53,0.15},{0.6179529455155457,0.19125644692155786}},
        {0.0,0.10653709580151404,0.13170091743633459,0.20458556663833828},
        {2,1,1,2},
        1
    };

    gbs::BSCurve2d_d crv2{
        {{0.4,0.25},{0.4,0.25},{0.6000000000000001,0.29307266043245908},{0.6000000000000001,0.29307266043245908}},
        {0.0,0.20458556663833828},
        {4,4},
        3
    };

    gbs::BSCurve2d_d crv3{
        {{0.4388228567653781,0.39488887394336028},{0.43882285676537816,0.39488887394336028},{0.5820470544844545,0.39488887394336028},{0.5820470544844545,0.39488887394336028}},
        {0.0,0.20458556663833828},
        {4,4},
        3
    };

    auto s = gbs::loft( std::list<gbs::BSCurve<double,2>>{ crv1, crv2, crv3}, std::vector<double>{0.,0.5,1.}, 2 );
    if(PLOT_ON)
        gbs::plot(s,crv1, crv2, s.isoU(0.1), s.isoV(0.5));
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

    auto s = gbs::loft_generic<double,3>( std::list<gbs::BSCurve3d_d>{ c1, c2, c3} , 3);
    std::list<gbs::BSCurve3d_d> bs_lst2 = {c1,c2,c3};
    auto s2 = gbs::loft( bs_lst2 );
    if(PLOT_ON)
        gbs::plot(s,c1,c2,c3);
}

TEST(tests_bssbuild, loft_u)
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
    k[3] = 1.2;
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
    k[5] = 2.8;
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
    auto c1 {std::make_shared<gbs::BSCurve3d_d>(poles1,k,p)};
    auto c2 {std::make_shared<gbs::BSCurve3d_d>(poles2,k,p)};
    auto c3 {std::make_shared<gbs::BSCurve3d_d>(poles3,k,p)};

    std::vector<std::shared_ptr<gbs::BSCurveGeneral<double,3,false>>> bs_lst = {c1,c2,c3};
    auto s = gbs::loft_generic<double,3>( bs_lst, {0.,0.5,1.}, 2 );
    // auto s = gbs::loft<double,3,false>( bs_lst.begin(),bs_lst.end(), {0.,0.5,1.}, 2 );
    if(PLOT_ON)
        gbs::plot(s,c1,c2,c3);
}

TEST(tests_bssbuild, loft1d)
{
    size_t p = 1;
    std::vector<double> k = {0., 0., 1., 1.};
    std::vector<std::array<double,1>> poles1 =
    {
        {1.0},
        {1.2},
    };
    std::vector<std::array<double,1>> poles2 =
    {
        {2.0},
        {1.8},
    };
    gbs::BSCurve<double,1> c1{poles1,k,p};
    gbs::BSCurve<double,1> c2{poles2,k,p};
    auto s = gbs::loft( {c1, c2}, 3 );

    ASSERT_DOUBLE_EQ(s(0., 0.)[0], poles1[0][0]);
    ASSERT_DOUBLE_EQ(s(1., 0.)[0], poles1[1][0]);
    ASSERT_DOUBLE_EQ(s(0., 1.)[0], poles2[0][0]);
    ASSERT_DOUBLE_EQ(s(1., 1.)[0], poles2[1][0]);
    ASSERT_DOUBLE_EQ(s(0.5, 0.5)[0], 0.25 * ( poles2[1][0] + poles2[0][0] + poles1[1][0] + poles1[0][0]) );
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

    std::list<gbs::BSCurveRational<double,3>> bs_lst = {c1,c2,c3};
    auto s = gbs::loft( bs_lst, {0., 0.5, 1.0}, 3 );
    if(PLOT_ON)
        gbs::plot(s,c1,c2,c3);
}

// TODO: Fix shape
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
    gbs::BSCurve3d_d sp({{1.5,1.5,0.},{1.5,1.0,1.5},{1.5,1.5,3.}},{0.,0.,0.,1.,1.,1.},2);

    std::list<gbs::BSCurveGeneral<double,3,false>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst, sp );
    if(PLOT_ON)
        gbs::plot(s,c1,c2,c3,sp);

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

TEST(tests_bssbuild, gordon_bss)
{
    using namespace gbs;
    using T = double;
    const size_t dim = 3;
    size_t p = 3;
    size_t q = 2;

    std::vector<T> u{0.,0.33,0.66,1.};
    std::vector<T> v{0.,0.5,1.};
    auto ku = build_simple_mult_flat_knots<T>(u, p);
    auto kv = build_simple_mult_flat_knots<T>(v, q);
    points_vector<T,dim> Q{
        {0.,0.,0.}, {0.33,0.,0.}, {0.66,0.,0.}, {1.,0.,0.},
        {0.,0.5,0.}, {0.33,0.5,0.1}, {0.66,0.5,0.2}, {1.,0.5,0.},
        {0.,1.,0.}, {0.33,1.,0.}, {0.66,1.,0.}, {1.,1.,0.},
    };
    auto poles_T = build_poles(Q, ku, kv, u, v, p, q);
    BSSurface<T,dim> Tuv{
        poles_T,
        ku,kv,
        p,q
    };

    std::vector<gbs::constrType<T,dim,1> > Q_u_crv1(4),Q_u_crv2(4),Q_u_crv3(4);
    std::transform(Q.begin(),std::next(Q.begin(),4),Q_u_crv1.begin(),[](const auto Q_){return gbs::constrType<T,dim,1>{Q_};});
    std::transform(std::next(Q.begin(),4),std::next(Q.begin(),8),Q_u_crv2.begin(),[](const auto Q_){return gbs::constrType<T,dim,1>{Q_};});
    std::transform(std::next(Q.begin(),8),std::next(Q.begin(),12),Q_u_crv3.begin(),[](const auto Q_){return gbs::constrType<T,dim,1>{Q_};});
    auto u_crv1 = BSCurve<T,dim>(build_poles(Q_u_crv1, ku, u, p), ku, p);
    auto u_crv2 = BSCurve<T,dim>(build_poles(Q_u_crv2, ku, u, p), ku, p);
    auto u_crv3 = BSCurve<T,dim>(build_poles(Q_u_crv3, ku, u, p), ku, p);

    std::vector<gbs::constrType<T,dim,1> > Q_v_crv1(3),Q_v_crv2(3),Q_v_crv3(3),Q_v_crv4(3);
    Q_v_crv1[0] = gbs::constrType<T,dim,1>{Q[0]};
    Q_v_crv1[1] = gbs::constrType<T,dim,1>{Q[4]};
    Q_v_crv1[2] = gbs::constrType<T,dim,1>{Q[8]};

    Q_v_crv2[0] = gbs::constrType<T,dim,1>{Q[1]};
    Q_v_crv2[1] = gbs::constrType<T,dim,1>{Q[5]};
    Q_v_crv2[2] = gbs::constrType<T,dim,1>{Q[9]};

    Q_v_crv3[0] = gbs::constrType<T,dim,1>{Q[2]};
    Q_v_crv3[1] = gbs::constrType<T,dim,1>{Q[6]};
    Q_v_crv3[2] = gbs::constrType<T,dim,1>{Q[10]};

    Q_v_crv4[0] = gbs::constrType<T,dim,1>{Q[3]};
    Q_v_crv4[1] = gbs::constrType<T,dim,1>{Q[7]};
    Q_v_crv4[2] = gbs::constrType<T,dim,1>{Q[11]};

    auto v_crv1 = BSCurve<T,dim>(build_poles(Q_v_crv1, kv, v, q), kv, q);
    auto v_crv2 = BSCurve<T,dim>(build_poles(Q_v_crv2, kv, v, q), kv, q);
    auto v_crv3 = BSCurve<T,dim>(build_poles(Q_v_crv3, kv, v, q), kv, q);
    auto v_crv4 = BSCurve<T,dim>(build_poles(Q_v_crv4, kv, v, q), kv, q);

    std::vector<BSCurve<T,dim>> u_crv_lst{ u_crv1, u_crv2, u_crv3};
    std::vector<BSCurve<T,dim>> v_crv_lst{ v_crv1, v_crv2, v_crv3, v_crv4};

    auto G2 =gbs::gordon<T,dim>(u_crv_lst.begin(), u_crv_lst.end(), v_crv_lst.begin(),v_crv_lst.end());

    for( auto u_ : gbs::make_range(u.front(), u.back(), 100))
    {
        ASSERT_LT( gbs::distance( G2(u_,v[0]), u_crv1(u_) ), 1e-6);
        ASSERT_LT( gbs::distance( G2(u_,v[1]), u_crv2(u_) ), 1e-6);
        ASSERT_LT( gbs::distance( G2(u_,v[2]), u_crv3(u_) ), 1e-6);
    }
    for( auto v_ : gbs::make_range(v.front(), v.back(), 100))
    {
        ASSERT_LT( gbs::distance( G2(u[0],v_), v_crv1(v_) ), 1e-6);
        ASSERT_LT( gbs::distance( G2(u[1],v_), v_crv2(v_) ), 1e-6);
        ASSERT_LT( gbs::distance( G2(u[2],v_), v_crv3(v_) ), 1e-6);
        ASSERT_LT( gbs::distance( G2(u[3],v_), v_crv4(v_) ), 1e-6);
    }

    if(PLOT_ON)
        gbs::plot(
            G2,
            u_crv1, u_crv2, u_crv3, 
            v_crv1, v_crv2, v_crv3, v_crv4
        );
}

TEST(tests_bssbuild,gordon_foils)
{
    GTEST_SKIP() << "Skipping this test for now.";
    using namespace gbs;
    using T = double;

    
    std::string dir = get_directory(__FILE__);

    auto crv1_2d = bscurve_approx_from_points<T,2>(dir +"/in/e1098.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv2_2d = bscurve_approx_from_points<T,2>(dir +"/in/e817.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv3_2d = bscurve_approx_from_points<T,2>(dir +"/in/e186.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    translate(crv2_2d,{-0.5,0.});
    rotate(crv2_2d,std::numbers::pi/8.);
    translate(crv2_2d,{0.5,0.});

    translate(crv3_2d,{-0.5,0.});
    rotate(crv3_2d,std::numbers::pi/6.);
    translate(crv3_2d,{0.5,0.});

    auto crv1 = add_dimension(crv1_2d,0.0);
    auto crv2 = add_dimension(crv2_2d,1.0);
    auto crv3 = add_dimension(crv3_2d,2.0);

    std::vector<BSCurve<T,3>> u_crv_lst = {crv1, crv2, crv3};

    // auto Lu = loft<T,3,false>(u_crv_lst.begin(),u_crv_lst.end(),{0., 0.5, 1.},2);
    auto Lu = loft(u_crv_lst,{0., 0.5, 1.},2);

    
    auto g1 = interpolate(points_vector<T,3>{crv1.begin(),crv2.begin(), crv3.begin()},{0.,0.5,1.},2);
    auto ule1 = max_curvature_pos(crv1, 1e-6)[0];
    auto ule2 = max_curvature_pos(crv2, 1e-6)[0];
    auto ule3 = max_curvature_pos(crv3, 1e-6)[0];
    auto g2 = interpolate(points_vector<T,3>{crv1(0.5),crv2(0.5), crv3(0.5)},{0.,0.5,1.},2);
    auto g3 = interpolate(points_vector<T,3>{crv1.end(),crv2.end(), crv3.end()},{0.,0.5,1.},2);

    if(PLOT_ON)
        gbs::plot( crv1, crv2, crv3, g1, g2, g3, Lu );

}

TEST(tests_bssbuild, gordon_prop)
{

    using namespace gbs;
    using T = double;
    constexpr size_t dim = 3;

    rapidjson::Document document;
    std::string dir = get_directory(__FILE__);
    gbs::parse_file((dir + "/in/prop3.json").c_str(), document);

    auto s1 = gbs::make_surface<T, dim>(document["s1"].GetObject());
    auto s2 = gbs::make_surface<T, dim>(document["s2"].GetObject());
    auto le = gbs::make_surface<T, dim>(document["le"].GetObject());
    auto te = gbs::make_surface<T, dim>(document["te"].GetObject());

    auto spans = make_vec<T>(document["spans"]);

    auto s1_bs = std::dynamic_pointer_cast<gbs::BSSurface<T, dim>>(s1);
    auto s2_bs = std::dynamic_pointer_cast<gbs::BSSurface<T, dim>>(s2);
    auto le_bs = std::dynamic_pointer_cast<gbs::BSSurface<T, dim>>(le);
    auto te_bs = std::dynamic_pointer_cast<gbs::BSSurface<T, dim>>(te);

    s2_bs->changeUBounds(s1_bs->boundsU()[0], s1_bs->boundsU()[1]);
    te_bs->changeUBounds(le_bs->boundsU()[0], le_bs->boundsU()[1]);

    auto le1 = s1_bs->isoU(0.);
    auto le2 = s2_bs->isoU(0.);
    auto te1 = s1_bs->isoU(s1_bs->boundsU()[1]);
    auto te2 = s2_bs->isoU(s2_bs->boundsU()[1]);


    // auto j1 = c2_connect(le1, te1.reversed());
    // auto j2 = c2_connect(le2, te2.reversed());
    auto t1 = cn_connect(le1, te1.reversed(), 1., 3, 2);
    auto t2 = cn_connect(le2, te2.reversed(), 1., 3, 2);

    auto j1 = join(join(le1, t1), te1.reversed());
    auto j2 = join(join(le2, t2), te2.reversed());

    // std::vector<gbs::BSCurve<T, dim>> v_crv{s1_bs->isoU(0.), s2_bs->isoU(0.)}, u_crv{le_bs->isoV(0.), le_bs->isoV(1.)};
    // {
    //     auto G = gbs::gordon<double, 3>(u_crv.begin(), u_crv.end(), v_crv.begin(), v_crv.end(), 1e-6);

    //     T u_ = le->boundsU()[1];
    //     T v_ = le->boundsV()[1];
    //     ASSERT_LE(distance(G(0., 0.5), le->value(0., 0.5)), 1e-6);
    //     ASSERT_LE(distance(G(u_, 0.5), le->value(u_, 0.5)), 1e-6);

    //     if (PLOT_ON)
    //         gbs::plot(s1, s2, G, te, v_crv, u_crv);
    // }

    {
        auto le_root = le_bs->isoV(0.);
        auto te_root = te_bs->isoV(0.);
        std::vector<gbs::BSCurve<T, dim>> 
            v_crv{j1, j2},
            u_crv{le_root, te_root};
        

        auto G = gbs::gordon<double, 3>(u_crv.begin(), u_crv.end(), v_crv.begin(), v_crv.end(), 1e-6);

        gbs::plot(le1, le2, j1, j2, te1, te2, G);
    }

    // {
    //     std::vector<gbs::BSCurve<T, dim>> u_crv(2*spans.size());
    //     std::transform(
    //         spans.begin(), spans.end(),
    //         u_crv.begin(),
    //         [&le_bs](auto v_){
    //             return le_bs->isoV(v_); 
    //         }
    //     );
    //     std::transform(
    //         spans.rbegin(), spans.rend(),
    //         std::next(u_crv.begin(), spans.size() ),
    //         [&te_bs](auto v_){
    //             return te_bs->isoV(v_); 
    //         }
    //     );
    //     std::vector<gbs::BSCurve<T, dim>> v_crv{j1, j2};
    //     auto G = gbs::gordon<double, 3>(u_crv.begin(), u_crv.end(), v_crv.begin(), v_crv.end(), 1e-6);
    //     gbs::plot(le1, le2, j1, j2, te1, te2, G);
    // }
}

TEST(tests_bssbuild, loft_algo)
{
    using T = double;
    const size_t d = 3;

    size_t p = 2;
    std::vector<T> ku = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T, d>> poles1 =
        {
            {0., 0., 0.},
            {0., 1., 0.},
            {1., 1., 0.},
            {1., 1., 0.},
            {1., 1., 0.},
            {3., 1., 0.},
            {0., 4., 0.},
        };
    gbs::BSCurve<T, d> c1(poles1, ku, p);

    ku[3] = 1.2;
    std::vector<std::array<T, d>> poles2 =
        {
            {0., 0., 1.},
            {0., 1., 1.},
            {1., 1., 1.},
            {1., 1., 1.},
            {1., 1., 1.},
            {3., 1., 1.},
            {2., 4., 1.},
        };
    gbs::BSCurve<T, d> c2(poles2, ku, p);
    ku[5] = 2.8;

    std::vector<std::array<T, d>> poles3 =
        {
            {0., 0., 2.},
            {0., 1., 2.},
            {1., 1., 2.},
            {1., 1., 2.},
            {1., 1., 2.},
            {3., 1., 2.},
            {1., 4., 2.},
        };
    gbs::BSCurve<T, d> c3(poles3, ku, p);

    std::vector<T> v{0., 0.5, 1.0};
    size_t q = 2;
    auto flat_v = gbs::build_simple_mult_flat_knots(v, q);

    ///////////////////////
    // Loft construction //
    ///////////////////////
    std::vector<gbs::BSCurve<T, d>> bsc_lst_original{c1, c2, c3};
    auto bsc_lst = gbs::unified_curves(bsc_lst_original);
    // store unified values
    auto ncr= std::distance(bsc_lst.begin(), bsc_lst.end());
    auto nv = ncr;
    auto nu = bsc_lst.front().poles().size();
    std::vector<std::vector<std::array<T, 3>>> poles_curves(ncr);
    std::transform(
        bsc_lst.begin(), bsc_lst.end(),
        poles_curves.begin(),
        [](const auto &curve){return curve.poles();}
    );
    auto flat_u = bsc_lst.front().knotsFlats();
    // Build surface
    auto poles = gbs::build_loft_surface_poles(poles_curves, flat_v, v, flat_u, p, q);
    gbs::BSSurface<T, d> srf(poles, flat_u, flat_v, p, q);
    // Check Loft definition
    auto [u1, u2] = c1.bounds();
    for (size_t j{}; j < ncr; j++)
    {
        auto v_ = v[j];
        const auto &crv = bsc_lst_original[j]; 
        for (T u : gbs::make_range(u1, u2, 100))
        {
            ASSERT_LT(gbs::distance(srf(u, v_), crv(u)), std::numeric_limits<T>::epsilon()*10);
        }
    }

    if(PLOT_ON)
        gbs::plot(gbs::make_actor(srf, { 51./255.,  161./255.,  201./255.}, 100, 100), c1, c2, c3);
}


TEST(tests_bssbuild, loft2)
{
    using namespace gbs;
    using T = double;
    const size_t d = 3;


    std::string dir = get_directory(__FILE__);
    auto crv1_2d = bscurve_approx_from_points<T,2>(dir + "/in/e1098.dat",3,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv2_2d = bscurve_approx_from_points<T,2>(dir + "/in/e817.dat",4,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv3_2d = bscurve_approx_from_points<T,2>(dir + "/in/e186.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    translate(crv2_2d,{-0.5,0.});
    rotate(crv2_2d,std::numbers::pi/8.);
    translate(crv2_2d,{0.5,0.});

    translate(crv3_2d,{-0.5,0.});
    rotate(crv3_2d,std::numbers::pi/6.);
    translate(crv3_2d,{0.5,0.});

    auto c1 = add_dimension(crv1_2d,0.0);
    auto c2 = add_dimension(crv2_2d,1.0);
    auto c3 = add_dimension(crv3_2d,2.0);

    std::vector<T> v{0., 0.5, 1.0};


    ///////////////////////
    // Loft construction //
    ///////////////////////
    size_t q = 2;
    auto flat_v = gbs::build_simple_mult_flat_knots(v, q);
    std::vector<gbs::BSCurve<T, d>> bsc_lst_original{c1, c2, c3};

    auto curves_info = gbs::get_bs_curves_info<T,d>(bsc_lst_original.begin(), bsc_lst_original.end());

    gbs::unify_degree(curves_info);
    auto degree = std::get<2>(curves_info.front());
    for(auto [poles, knots, p] : curves_info)
    {
        ASSERT_EQ(p, degree);
    }
    gbs::unify_knots(curves_info);
    auto &[poles0, k0, degree0] = curves_info.front();
    std::vector<T> delta(k0.size());
    for(auto &curve_info : curves_info)
    {
        const auto &k = std::get<1>(curve_info);
        ASSERT_EQ(k.size(), k0.size());
        for(size_t i{}; i < k0.size(); i++)
            ASSERT_DOUBLE_EQ(k[i], k0[i]);   
    }

    auto p = std::get<2>(curves_info.front());
    auto ncr= std::distance(bsc_lst_original.begin(), bsc_lst_original.end());
    auto nv = ncr;
    auto nu = std::get<0>(curves_info.front()).size();
    auto flat_u = std::get<1>(curves_info.front());
    std::vector<std::vector<std::array<T, 3>>> poles_curves(ncr);
    std::transform(
        curves_info.begin(), curves_info.end(),
        poles_curves.begin(),
        [](const auto &curve_info){return std::get<0>(curve_info);}
    );

    // Build surface
    auto poles = build_loft_surface_poles(poles_curves, flat_v, v, flat_u, p, q);
    gbs::BSSurface<T, d> srf(poles, flat_u, flat_v, p, q);
    // Check Loft definition
    auto [u1, u2] = c1.bounds();
    auto bsc_lst = gbs::unified_curves(bsc_lst_original);
    for (size_t j{}; j < ncr; j++)
    {
        auto v_ = v[j];
        const auto &crv = bsc_lst[j]; 
        for (T u : gbs::make_range(u1, u2, 100))
        {
            ASSERT_LT(gbs::distance(srf(u, v_), crv(u)), 1e-7);
        }
    }

    auto [poles_, flat_u_, p_]  = loft(curves_info, v, flat_v, q);

    gbs::BSSurface<T, d> srf_(poles_, flat_u_, flat_v, p_, q);
        for (size_t j{}; j < ncr; j++)
    {
        auto v_ = v[j];
        const auto &crv = bsc_lst[j]; 
        for (T u : gbs::make_range(u1, u2, 100))
        {
            ASSERT_LT(gbs::distance(srf_(u, v_), crv(u)), 1e-7);
        }
    }

    if(PLOT_ON)
        gbs::plot(gbs::make_actor(srf_, { 51./255.,  161./255.,  201./255.}, 500, 500), c1, c2, c3);

}