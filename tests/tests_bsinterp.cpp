#include <gtest/gtest.h>
#include <gbs/bscinterp.h>
#include <gbs/bscapprox.h>
#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>
#include <gbs/extrema.h>
#include <gbs/bscanalysis.h>
#include <gbs-occt/export.h>
#include <gbs-occt/curvesbuild.h>
#include <gbs-render/vtkcurvesrender.h>

#include <algorithm>
#include <iostream>
#include <fstream>

#include <GeomTools.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

const double tol = 1e-10;

using gbs::operator-;

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

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(c_bez));
    for (size_t p = 1; p < Q.size(); p++)
    {
        auto c = gbs::interpolate(Q, p, gbs::KnotsCalcMode::CHORD_LENGTH);
        ASSERT_TRUE(test_crv(Q, c, gbs::KnotsCalcMode::CHORD_LENGTH, true));
        crv_lst.push_back(occt_utils::BSplineCurve(c));
    }
    // p = 2 produit les courbes les + smmooths
    occt_utils::to_iges(crv_lst, "tests/out/interp1.igs");
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

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(c1));
    occt_utils::to_iges(crv_lst, "tests/out/interpC1.igs");
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

    // gbs::interpolate(Q,p,gbs::KnotsCalcMode::CHORD_LENGTH);

    std::vector<Handle_Geom_Geometry> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(crv1));
    occt_utils::to_iges(crv_lst, "tests/out/build_3pt_tg.igs");
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
    gbs::extr<double,2> pt_begin{0.,{0.,0.}};
    gbs::extr<double,2> pt_end{1.,{1.,0.}};
    gbs::cstr<double,2> pt_int{0.5,{0.3,0.2},0};
    gbs::cstr<double,2> cr_begin{0.0,{0.0,0.0},2};
    gbs::cstr<double,2> cr_end{1.0,{0.0,0.0},2};

    std::vector<gbs::cstr<double,2>> cstr_lst = {
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

    gbs::plot(
        crv
    );
}