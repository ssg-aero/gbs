#include <doctest_gtest.hpp>
#include <gbs/bscinterp.h>
#include <gbs/bscapprox.h>
#ifdef GBS_USE_MODULES
    import knots_functions;
#else
    #include <gbs/knotsfunctions.ixx>
#endif
#include <gbs/bscurve.h>
#include <gbs/extrema.h>
#include <gbs/bscanalysis.h>
#include <gbs-render/vtkGbsRender.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <numbers>

const double tol = 1e-10;

using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

template <typename T, size_t dim, size_t nc>
inline void test_crv(const std::vector<gbs::constrType<T, dim, nc>> Q, const gbs::BSCurve<T, dim> &c, const std::vector<T> &u)
{
    for (int i = 0; i < Q.size(); i++)
    {
        for (int j = 0; j < Q[i].size(); j++)
        {
            auto p1 = c.value(u[i], j);
            auto p2 = Q[i][j];
            ASSERT_LT(gbs::distance(p1, p2), tol);
        }
    }
}

template <typename T, size_t dim, size_t nc>
inline bool test_crv(const std::vector<gbs::constrType<T, dim, nc>> Q, const gbs::BSCurve<T, dim> &c, gbs::KnotsCalcMode mode, bool adim = false)
{
    std::vector<std::array<T, dim>> pts(Q.size());
    std::transform(Q.begin(), Q.end(), pts.begin(), [](const auto &q_) { return q_.front(); });
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, adim);
    for (int i = 0; i < Q.size(); i++)
    {
        for (int j = 0; j < Q[i].size(); j++)
        {
            auto p1 = c.value(u[i], j);
            auto p2 = Q[i][j];
            if (gbs::distance(p1, p2) > tol)
                return false;
        }
    }
    return true;
}

TEST(tests_bscurve, interp1)
{
    std::vector<gbs::constrType<double, 3, 1>> Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
        };

    std::vector<std::array<double, 3>> pts(Q.size());
    std::transform(Q.begin(), Q.end(), pts.begin(), [](const auto &q_) { return q_.front(); });

    // interp bezier

    auto p = pts.size() - 1;
    std::vector<size_t> m = {p + 1, p + 1};
    auto k_flat = gbs::flat_knots(std::vector<double>{0., 1.}, m);
    std::vector<double> u = {0., 0.5, 1.};

    auto c_bez = gbs::BSCurve<double,3>(gbs::build_poles(Q, k_flat, u, p), k_flat, p);

    test_crv(Q, c_bez, u);

    Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
            {{1., 1., 2.}},
            {{0., 1., 1.}},
            {{0., -1., 4.}},
        };

    for (size_t p = 1; p < Q.size(); p++)
    {
        auto c = gbs::interpolate(Q, p, gbs::KnotsCalcMode::CHORD_LENGTH);
        ASSERT_TRUE(test_crv(Q, c, gbs::KnotsCalcMode::CHORD_LENGTH, true));
    }
    // p = 2 produit les courbes les + smmooths
}

TEST(tests_bscurve, interp_C1)
{
    std::vector<gbs::constrType<double, 3, 2>> Q =
        {
            {{{0., 0., 0.}, {0., 1., 0.}}},
            {{{1., 0., 0.}, {1., 0., 0.}}},
            {{{1., 1., 0.}, {0., 1., 0.}}},
        };

    auto mode = gbs::KnotsCalcMode::CHORD_LENGTH;
    auto c1 = gbs::interpolate(Q, mode);

    ASSERT_TRUE(test_crv(Q, c1, mode));
}

TEST(tests_bscurve, interp_C1_to_C2)
{
    std::vector<gbs::constrType<double, 3, 2>> Q =
        {
            {{{0., 0., 0.}, {0., 1., 0.}}},
            {{{1., 0., 0.}, {1., 0., 0.}}},
            {{{1., 1., 0.}, {0., 1., 0.}}},
        };

    auto mode = gbs::KnotsCalcMode::CHORD_LENGTH;
    auto c1 = gbs::interpolate(Q, mode,true);

    ASSERT_TRUE(test_crv(Q, c1, mode));
    ASSERT_TRUE(c1.degree() == 5);
    // visual stuff
    auto c2 = gbs::interpolate(Q, mode);
    if(PLOT_ON)
        gbs::plot(
            // c1, c2,
            gbs::crv_dsp<double,3,false>{
                .c =&c1,
                .col_crv = {1.,0.,0.},
                .poles_on = true,
                .col_poles = {0.,1.,0.},
                // .col_ctrl = {0.,0.,0.},
                .show_curvature=true,
                } // c++20
            , gbs::crv_dsp<double,3,false>{
                .c =&c2,
                .col_crv = {0.,1.,0.},
                .poles_on = true,
                .col_poles = {0.,1.,1.},
                // .col_ctrl = {0.,0.,0.},
                .show_curvature=true,
                } // c++20
        );
}

TEST(tests_bscurve, interp_C2)
{

    // https://stackoverflow.com/questions/49689336/initialize-struct-using-initializer-list-fails-for-stdarray
    std::vector<gbs::constrType<double, 3, 3>> Q =
        {
            {{{0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}}},
            {{{1., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}},
            {{{1., 0.3, 1.}, {1., 0., 0.}, {0., 1., 0.}}},
            {{{1., 1., 0.}, {0., 1., 0.}, {0., 0., 1.}}},
        };

    gbs::KnotsCalcMode mode = gbs::KnotsCalcMode::CHORD_LENGTH;
    ASSERT_TRUE(test_crv(Q, gbs::interpolate(Q, mode), mode));

    std::vector<gbs::constrType<double, 2, 3>> Q2d =
        {
            {{{0., 0.}, {0., 1.}, {1., 0.}}},
            {{{1., 0.}, {1., 0.}, {0., 1.}}},
            {{{1., 0.3}, {1., 0.}, {0., 1.}}},
            {{{1., 1.}, {0., 1.}, {0.2, 0.}}},
        };

    gbs::KnotsCalcMode mode2d = gbs::KnotsCalcMode::CHORD_LENGTH;
    ASSERT_TRUE(test_crv(Q, gbs::interpolate(Q, mode2d), mode2d));

    std::vector<gbs::constrType<double, 1, 3>> Q1d =
        {
            {{{0.}, {0.}, {1.}}},
            {{{1.}, {1.}, {0.}}},
            {{{1.}, {1.}, {0.}}},
            {{{1.}, {0.}, {0.2}}},
        };

    gbs::KnotsCalcMode mode1d = gbs::KnotsCalcMode::CHORD_LENGTH;
    ASSERT_TRUE(test_crv(Q, gbs::interpolate(Q, mode1d), mode1d));

    std::vector<gbs::constrType<double, 4, 3>> Q4d =
        {
            {{{0., 0., 0., 1}, {0., 1., 0., 0}, {1., 0., 0., 0}}},
            {{{1., 0., 0., 2}, {1., 0., 0., 0}, {0., 1., 0., 0}}},
            {{{1., 0.3, 1., 3}, {1., 0., 0., 0}, {0., 1., 0., 0}}},
            {{{1., 1., 0., 4}, {0., 1., 0., 0}, {0., 0., 1., 0}}},
        };
    gbs::KnotsCalcMode mode4d = gbs::KnotsCalcMode::CHORD_LENGTH;
    ASSERT_TRUE(test_crv(Q, gbs::interpolate(Q, mode4d), mode4d));
}

TEST(tests_bscurve, build_3pt_tg)
{
    std::vector<std::array<double, 3>> pts =
        {
            {0., 0., 0.},
            {0., 1., 0.},
            {1., 1., 0.},
            // {1.,1.,1.},
            // {1.,1.,2.},
            // {3.,1.,1.},
            // {0.,4.,1.},
            {3., 1., 0.},
            {4., 0., 0.},
        };
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH);

    auto D = gbs::build_3pt_tg_vec(pts, u);
    auto T = gbs::build_3pt_tg_dir(pts, u);

    std::vector<gbs::constrType<double, 3, 2>> Q(pts.size());
    std::transform(
        pts.begin(), pts.end(),
        // T.begin(),
        D.begin(),
        Q.begin(),
        [](const auto &pt_, const auto &tg_) {
            return gbs::constrType<double, 3, 2>({pt_, tg_});
        });

    gbs::KnotsCalcMode mode = gbs::KnotsCalcMode::CHORD_LENGTH;
    auto crv1 = gbs::interpolate(Q, mode);
    ASSERT_TRUE(test_crv(Q, crv1, mode));
}

TEST(tests_bscurve, multi_constrained)
{
    std::vector<gbs::constrPoint<double, 2>> cstr{
        {{0., 0. }, 0. , 0},
        {{.5, 0.3}, 0.5, 0},
        {{1., 0. }, 1. , 0},
        {{0., 0.2}, 0. , 1},
        {{0.,-0.1}, 1. , 1},
    };
    size_t p = 2;
    auto k_flat = gbs::build_simple_mult_flat_knots(0.,1.,cstr.size(),p);
    auto poles = gbs::build_poles(cstr,k_flat,p);
    auto crv = gbs::BSCurve2d_d(poles,k_flat,p);

    std::for_each(
        cstr.begin(),
        cstr.end(),
        [&crv](const auto &c_)
        {
            ASSERT_LT(gbs::norm( crv.value(c_.u,c_.d) - c_.v ),1e-6);
        }
    );
}


TEST(tests_bscurve, multi_constrained_general)
{

    // "Natural" cubic bspline interp
    gbs::bsc_bound<double,2> pt_begin{0.,{0.,0.}};
    gbs::bsc_bound<double,2> pt_end{1.,{1.,0.}};
    gbs::bsc_constraint<double,2> pt_int{0.5,{0.3,0.2},0};
    gbs::bsc_constraint<double,2> cr_begin{0.0,{0.0,0.0},2};
    gbs::bsc_constraint<double,2> cr_end{1.0,{0.0,0.0},2};

    std::vector<gbs::bsc_constraint<double,2>> cstr_lst = {
        pt_int,
        cr_begin,   
        cr_end
    };

    auto crv = gbs::interpolate(pt_begin,pt_end,cstr_lst,3);

    ASSERT_LT(gbs::norm( crv.value(std::get<0>(pt_begin)) - std::get<1>(pt_begin) ),1e-6);
    ASSERT_LT(gbs::norm( crv.value(std::get<0>(pt_end)) - std::get<1>(pt_end) ),1e-6);
    ASSERT_LT(gbs::norm( crv.value(std::get<0>(pt_int)) - std::get<1>(pt_int) ),1e-6);
    ASSERT_LT(gbs::norm( crv.value(std::get<0>(pt_begin),2) ),1e-6);
    ASSERT_LT(gbs::norm( crv.value(std::get<0>(pt_end),2) ),1e-6);
    if(PLOT_ON)
        gbs::plot(
            crv
        );
}

TEST(tests_bscurve, perf_build_poles_unif_constr)
{

    using T = double;
    size_t n = 360;
    std::vector<gbs::constrType<T,2,3> > Q(n);
    auto th = gbs::make_range<T>(0., std::numbers::pi,n);
    std::transform(
        th.begin(), th.end(), Q.begin(),
        [](T th){
            return gbs::constrType<T,2,3>{
                std::array<T,2>{ std::cos(th), std::sin(th)},
                std::array<T,2>{-std::sin(th), std::cos(th)},
                std::array<T,2>{-std::cos(th),-std::sin(th)}
            };
            
        }
    );

    size_t p = 5;
    auto k_flat = gbs::build_mult_flat_knots<T>(th, p, 3);
    auto poles =  gbs::build_poles(Q,k_flat,th,p);

    std::for_each(
        th.begin(), th.end(),
        [&](T th_)
        {
            ASSERT_LT(
                gbs::distance(
                    gbs::eval_value_decasteljau(th_,k_flat, poles, p, 0),
                    std::array<T,2>{ std::cos(th_), std::sin(th_)}
                ),
                1e-6
            );
            ASSERT_LT(
                gbs::distance(
                    gbs::eval_value_decasteljau(th_,k_flat, poles, p, 1),
                    std::array<T,2>{-std::sin(th_), std::cos(th_)}
                ),
                1e-6
            );
            ASSERT_LT(
                gbs::distance(
                    gbs::eval_value_decasteljau(th_,k_flat, poles, p, 2),
                    std::array<T,2>{-std::cos(th_),-std::sin(th_)}
                ),
                1e-6
            );
        }
    );
    
}
// Large-count curve interpolation goes through the banded sparse solve
// (issue #34, plan #3). Verify the curve interpolates all 201 points to
// tolerance — i.e. the SparseLU solve matches the data, at a size where the
// old dense O(n^3) factorization dominated.
TEST(tests_bsinterp, large_point_interpolation_banded)
{
    using T = double;
    const size_t n = 201, p = 3;
    gbs::points_vector<T, 3> pts(n);
    for (size_t i = 0; i < n; ++i)
    {
        T t = T(i) / T(n - 1);
        pts[i] = {std::cos(6.0 * t), std::sin(6.0 * t), 2.0 * t};
    }
    auto crv = gbs::interpolate(pts, p, gbs::KnotsCalcMode::CHORD_LENGTH);
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    ASSERT_EQ(u.size(), n);
    for (size_t i = 0; i < n; ++i)
        ASSERT_LT(gbs::distance(crv.value(u[i]), pts[i]), 1e-9) << "i=" << i;
}

// Large scattered constraint set exercises the sparse branch of
// build_poles(constrPoint) (issue #34): above the size threshold it solves with
// a sparse LU. Verify the curve interpolates all constraints to tolerance.
TEST(tests_bsinterp, scattered_constr_points_sparse_path)
{
    using T = double;
    const size_t n = 201, p = 3;
    std::vector<T> u(n);
    std::vector<gbs::constrPoint<T, 3>> Q(n);
    for (size_t i = 0; i < n; ++i)
    {
        T t = T(i) / T(n - 1);
        u[i] = t;
        Q[i] = gbs::constrPoint<T, 3>{ {std::cos(6.0 * t), std::sin(6.0 * t), 2.0 * t}, t, 0 };
    }
    auto k = gbs::build_simple_mult_flat_knots<T>(u, p);
    auto poles = gbs::build_poles(Q, k, p); // n >= threshold -> sparse path
    ASSERT_EQ(poles.size(), n);
    gbs::BSCurve<T, 3> crv(poles, k, p);
    for (size_t i = 0; i < n; ++i)
        ASSERT_LT(gbs::distance(crv.value(u[i]), Q[i].v), 1e-9) << "i=" << i;
}

// ===========================================================================
// Coverage edge cases (issue #50): the interpolation argument guards. These
// reject malformed inputs before the solve, and were uncovered by the baseline.
// ===========================================================================
TEST(tests_bsinterp, interpolation_argument_guards)
{
    using T = double;
    using gbs::point;

    // build_3pt_tg_dir needs at least 3 points (Bessel 3-point tangents).
    {
        std::vector<std::array<T, 2>> pts = {{0., 0.}, {1., 0.}};
        std::vector<T> u = {0., 1.};
        ASSERT_THROW(gbs::build_3pt_tg_dir(pts, u), std::invalid_argument);
    }

    // interpolate(Q, p, mode): degree must be strictly below the point count.
    {
        std::vector<gbs::constrType<T, 2, 1>> Q = {
            {point<T, 2>{0., 0.}}, {point<T, 2>{1., 1.}}, {point<T, 2>{2., 0.}}};
        ASSERT_THROW((gbs::interpolate<T, 2>(Q, size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH)),
                     std::domain_error);
    }

    // interpolate(begin, end, constraints, p): the number of constraints (here
    // none, so n = 2 endpoints) must be >= p+1, else the system is rank-deficient.
    {
        gbs::bsc_bound<T, 2> pt_begin{0., point<T, 2>{0., 0.}};
        gbs::bsc_bound<T, 2> pt_end{1., point<T, 2>{1., 1.}};
        std::vector<gbs::bsc_constraint<T, 2>> no_cstr;
        ASSERT_THROW(gbs::interpolate(pt_begin, pt_end, no_cstr, size_t{3}),
                     std::invalid_argument);
    }
}
