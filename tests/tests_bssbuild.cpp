#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bsctools.h>
#include <gbs/bssbuild.h>
#include <gbs-render/vtkGbsRender.h>
#include <gbs-io/fromtext.h>

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

    auto s = gbs::loft( std::list<gbs::BSCurve2d_d>{ crv1, crv2, crv3}, std::vector<double>{0.,0.5,1.}, 2 );
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

    std::list<gbs::BSCurveGeneral<double,3,false>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst );
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
    // auto s = gbs::loft( bs_lst, {0.,0.5,1.}, 2 );
    auto s = gbs::loft<double,3,false>( bs_lst.begin(),bs_lst.end(), {0.,0.5,1.}, 2 );
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
    auto s = gbs::loft<double,1>( {c1, c2} );

    ASSERT_DOUBLE_EQ(s(0., 0.)[0], poles1[0][0]);
    ASSERT_DOUBLE_EQ(s(1., 0.)[0], poles1[1][0]);
    ASSERT_DOUBLE_EQ(s(0., 1.)[0], poles2[0][0]);
    ASSERT_DOUBLE_EQ(s(1., 1.)[0], poles2[1][0]);
    ASSERT_DOUBLE_EQ(s(0.5, 0.5)[0], 0.25 * ( poles2[1][0] + poles2[0][0] + poles1[1][0] + poles1[0][0]) );
}

TEST(tests_bssbuild, loft_rational)
{
    GTEST_SKIP();
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

    if(PLOT_ON)
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
    GTEST_SKIP() << "Skipping this test for now.";
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

    if(PLOT_ON)
        plot( crv1, crv2, crv3, g1, g2, g3, Lu );

}
/**
 * Transposes a matrix of poles.
 * 
 * @param poles A 2D vector representing a matrix of poles.
 * @return A 2D vector representing the transposed matrix of poles.
 */
template <typename T, size_t dim>
std::vector<std::vector<std::array<T, 3>>> transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles) {
    if (poles.empty()) return {};

    size_t nu = poles.size();
    size_t nv = poles[0].size();

    std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

    for (size_t i = 0; i < nu; ++i) {
        for (size_t j = 0; j < nv; ++j) {
            poles_t[j][i] = poles[i][j];
        }
    }

    return poles_t;
}
/**
 * Transposes a block of a matrix of poles.
 * 
 * @param input The input matrix of poles to transpose.
 * @param output The output matrix where the transposed block will be stored.
 * @param startRow The starting row index for the block.
 * @param startCol The starting column index for the block.
 * @param blockSize The size of the block to be transposed.
 */
template <typename T, size_t dim>
// Function to transpose a small block of the matrix
void transpose_block(const std::vector<std::vector<std::array<T, dim>>>& input,
                    std::vector<std::vector<std::array<T, dim>>>& output,
                    size_t startRow, size_t startCol, size_t blockSize) {
    size_t endRow = std::min(startRow + blockSize, input.size());
    size_t endCol = std::min(startCol + blockSize, input[0].size());

    for (size_t i = startRow; i < endRow; ++i) {
        for (size_t j = startCol; j < endCol; ++j) {
            output[j][i] = input[i][j];
        }
    }
}
/**
 * Transposes a matrix of poles using block-based transposition for efficiency.
 * 
 * @param poles The matrix of poles to be transposed.
 * @param blockSize The size of each block used in the transposition.
 * @return The transposed matrix of poles.
 */
template <typename T, size_t dim>
auto transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles, size_t blockSize) {

    size_t nu = poles.size();
    size_t nv = poles[0].size();

    std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

    for (size_t i = 0; i < nu; i += blockSize) {
        for (size_t j = 0; j < nv; j += blockSize) {
            transpose_block(poles, poles_t, i, j, blockSize);
        }
    }

    return poles_t;
}
/**
 * Flattens a 2D vector of poles into a 1D vector.
 * 
 * @param poles_curves The 2D vector of poles to be flattened.
 * @return A 1D vector containing all the poles.
 */
template <typename T, size_t dim>
auto flatten_poles(const std::vector<std::vector<std::array<T, dim>>>& poles_curves) {
    std::vector<std::array<T, dim>> flattened;

    // Calculate total size for memory reservation
    size_t totalSize = std::accumulate(poles_curves.begin(), poles_curves.end(), size_t(0),
        [](size_t sum, const auto& row) { return sum + row.size(); });
    flattened.reserve(totalSize);

    // Flatten the vector
    for (const auto& row : poles_curves) {
        std::ranges::copy(row, std::back_inserter(flattened));
    }

    return flattened;
}
/**
 * Builds the poles for a loft surface from a set of curve poles, using NURBS mathematics.
 * 
 * @param poles_curves The poles of the curves used to build the loft surface.
 * @param flat_v The flattened knot vector in the v direction.
 * @param v The knot vector in the v direction.
 * @param flat_u The flattened knot vector in the u direction.
 * @param p The degree of the curve in the u direction.
 * @param q The degree of the curve in the v direction.
 * @return A flattened vector of poles representing the loft surface.
 */
template <typename T, size_t dim>
auto build_loft_surface_poles(const std::vector<std::vector<std::array<T, dim>>> &poles_curves,
                           const std::vector<T> &flat_v, const std::vector<T> &v, const std::vector<T> &flat_u, size_t p,
                           size_t q)
{
    auto poles_curves_t = transpose_poles(poles_curves);
    size_t ncr = poles_curves.size();
    size_t nu = poles_curves[0].size();
    size_t nv = poles_curves.size();

    gbs::MatrixX<T> N(ncr, ncr);
    gbs::build_poles_matrix<T, 1>(flat_v, v, q, ncr, N);
    auto N_inv = N.partialPivLu();

    std::vector<std::vector<std::array<T, dim>>> poles_t(nu, std::vector<std::array<T, dim>>(nv));

    for (size_t i = 0; i < nu; ++i)
    {
        for (size_t d = 0; d < dim; ++d)
        {
            gbs::VectorX<T> b(nv);
            std::transform(poles_curves_t[i].begin(), poles_curves_t[i].end(), b.begin(),
                           [d](const auto &arr)
                           { return arr[d]; });

            auto x = N_inv.solve(b);

            for (size_t j = 0; j < nv; ++j)
            {
                poles_t[i][j][d] = x(j);
            }
        }
    }

    return flatten_poles(transpose_poles(poles_t));
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
    // store unifed values
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
    auto poles = build_loft_surface_poles(poles_curves, flat_v, v, flat_u, p, q);
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
}           