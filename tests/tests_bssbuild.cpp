#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bssbuild.h>
#include <gbs-render/vtkcurvesrender.h>
#include <gbs-io/fromtext.h>

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
    // auto s = gbs::loft( bs_lst, {0.,0.5,1.}, 2 );
    auto s = gbs::loft<double,3,false>( bs_lst.begin(),bs_lst.end(), {0.,0.5,1.}, 2 );
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
    auto s = gbs::loft<double,1>( {c1, c2} );

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

    std::list<gbs::BSCurveGeneral<double,3,true>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst );
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
    auto ku = build_mult_flat_knots<T>({0.,1.}, p, 2);
    auto kv = build_mult_flat_knots<T>({0.,1.}, q, 2);
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

    std::vector<std::shared_ptr<gbs::BSCurveGeneral<T,dim,false>>> u_crv_lst = {
        std::make_shared<BSCurve<T,dim>>(u_crv1),
        std::make_shared<BSCurve<T,dim>>(u_crv2),
        std::make_shared<BSCurve<T,dim>>(u_crv3),
    };
    std::vector<std::shared_ptr<gbs::BSCurveGeneral<T,dim,false>>> v_crv_lst = {
        std::make_shared<BSCurve<T,dim>>(v_crv1),
        std::make_shared<BSCurve<T,dim>>(v_crv2),
        std::make_shared<BSCurve<T,dim>>(v_crv3),
        std::make_shared<BSCurve<T,dim>>(v_crv4),
    };

    auto Lu = loft<T,dim,false>(u_crv_lst.begin(),u_crv_lst.end(),v,q);
    auto Lv = loft<T,dim,false>(v_crv_lst.begin(),v_crv_lst.end(),u,p);
    Lv.invertUV();

    auto poles_gordon = Lu.poles();

    std::transform(
        Lv.poles().begin(),Lv.poles().end(),
        poles_gordon.begin(),
        poles_gordon.begin(),
        [](const auto &v_, const auto &g_){return g_ + v_;}
    );
    std::transform(
        Tuv.poles().begin(),Tuv.poles().end(),
        poles_gordon.begin(),
        poles_gordon.begin(),
        [](const auto &uv_, const auto &g_){return g_ - uv_;}
    );

    gbs::BSSurface<T,dim> G {
        poles_gordon,
        ku, kv,
        p, q
    };

    auto G2 =gbs::gordon<T,dim>(u_crv_lst.begin(), u_crv_lst.end(), v_crv_lst.begin(),v_crv_lst.end());

    gbs::plot(
        // Lu, 
        // Lv, 
        // Tuv, 
        // G,
        G2,
        u_crv1, u_crv2, u_crv3, 
        v_crv1, v_crv2, v_crv3, v_crv4,
        G.poles()
    );
}

TEST(tests_bssbuild,gordon_foils)
{
    using namespace gbs;
    using T = double;

    auto crv1_2d = bscurve_approx_from_points<T,2>("../tests/in/e1098.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv2_2d = bscurve_approx_from_points<T,2>("../tests/in/e817.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    auto crv3_2d = bscurve_approx_from_points<T,2>("../tests/in/e186.dat",5,KnotsCalcMode::CHORD_LENGTH,1);
    translate(crv2_2d,{-0.5,0.});
    rotate(crv2_2d,std::numbers::pi/8.);
    translate(crv2_2d,{0.5,0.});

    translate(crv3_2d,{-0.5,0.});
    rotate(crv3_2d,std::numbers::pi/6.);
    translate(crv3_2d,{0.5,0.});

    auto crv1 = add_dimension(crv1_2d,0.0);
    auto crv2 = add_dimension(crv2_2d,1.0);
    auto crv3 = add_dimension(crv3_2d,2.0);

    std::vector<std::shared_ptr<BSCurveGeneral<T,3,false>>> u_crv_lst = {
        std::make_shared<BSCurve<T,3>>(crv1),
        std::make_shared<BSCurve<T,3>>(crv2),
        std::make_shared<BSCurve<T,3>>(crv3),
    };

    auto Lu = loft<T,3,false>(u_crv_lst.begin(),u_crv_lst.end(),{0., 0.5, 1.},2);

    
    auto g1 = interpolate(points_vector<T,3>{crv1.begin(),crv2.begin(), crv3.begin()},{0.,0.5,1.},2);
    auto ule1 = max_curvature_pos(crv1, 1e-6)[0];
    auto ule2 = max_curvature_pos(crv2, 1e-6)[0];
    auto ule3 = max_curvature_pos(crv3, 1e-6)[0];
    auto g2 = interpolate(points_vector<T,3>{crv1(0.5),crv2(0.5), crv3(0.5)},{0.,0.5,1.},2);
    auto g3 = interpolate(points_vector<T,3>{crv1.end(),crv2.end(), crv3.end()},{0.,0.5,1.},2);

    plot( crv1, crv2, crv3, g1, g2, g3, Lu );

}