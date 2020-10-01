#include <gtest/gtest.h>
#include <gbslib/bscinterp.h>
#include <gbslib/bscapprox.h>
#include <gbslib/knotsfunctions.h>
#include <gbslib/bscurve.h>
#include <gbslib/extrema.h>
#include <gbslib/bscanalysis.h>
#include <occt-utils/export.h>
#include <occt-utils/curvesbuild.h>

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

TEST(tests_bscurve, approx_simple)
{

    std::string line;
    std::ifstream myfile("../tests/in/e1098.dat");
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

        // auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
        // auto crv = gbs::approx(pts, 5, 10, u);
        auto crv = gbs::approx(pts, 5, 10, gbs::KnotsCalcMode::CHORD_LENGTH);

        std::vector<Handle_Geom2d_Curve> crv_lst;
        crv_lst.push_back(occt_utils::BSplineCurve(crv));
        occt_utils::to_iges(crv_lst, "approx_simple.igs");
        auto res = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << res.d_avg << ", d_max:" << res.d_max << ", u_max:" << res.u_max << std::endl;
    }
    else
        std::cout << "Unable to open file";
}

TEST(tests_bscurve, approx_refined)
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

        auto crv = gbs::approx(pts, 5, gbs::KnotsCalcMode::CHORD_LENGTH,true);
        std::cout << "n poles: " << crv.poles().size() << ", n_flat: " << crv.knotsFlats().size() <<std::endl;

        std::vector<Handle_Geom2d_Curve> crv_lst;
        crv_lst.push_back(occt_utils::BSplineCurve(crv));
        // GeomTools::Dump( occt_utils::BSplineCurve(crv), std::cout );
        occt_utils::to_iges(crv_lst, "approx_refined.igs");

        auto res = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << res.d_avg << ", d_max:" << res.d_max << ", u_max:" << res.u_max << std::endl;

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_PC(crv, pnt, u0, 1e-6);
                u0 = res.u;
                // std::cout << gbs::norm(crv.value(u0)-pnt) << " , " << crv_lst.back()->Value(u0).Distance(occt_utils::point(pnt)) << std::endl;
            });

        std::vector<TopoDS_Shape> pt_lst(pts.size());
        std::transform(
            pts.begin(), pts.end(),
            pt_lst.begin(),
            [](const auto &p_) {
                return BRepBuilderAPI_MakeVertex (gp_Pnt(p_[0],p_[1],0.));
            });
        occt_utils::to_iges(pt_lst, "foilpoints.igs");
    }
    else
        std::cout << "Unable to open file";
}

TEST(tests_bscurve, approx_simple_nurbs)
{

    std::string line;
    std::ifstream myfile("../tests/in/e1098.dat");
    if (myfile.is_open())
    {
        std::vector<std::array<double, 3>> pts;
        getline(myfile, line);
        while (getline(myfile, line))
        {
            std::istringstream iss(line);
            std::string::size_type sz; // alias of size_t

            double x = std::stod(line, &sz);
            double y = std::stod(line.substr(sz));
            pts.push_back({x, y,1});
        }
        myfile.close();

        // auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
        // auto crv = gbs::approx(pts, 5, 10, u);
        auto crv = gbs::approx(pts, 5, 10, gbs::KnotsCalcMode::CHORD_LENGTH);

        std::vector<Handle_Geom_Curve> crv_lst;
        crv_lst.push_back(occt_utils::BSplineCurve(crv));
        occt_utils::to_iges(crv_lst, "approx_simple_nurbs.igs");
        auto res = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << res.d_avg << ", d_max:" << res.d_max << ", u_max:" << res.u_max << std::endl;
    }
    else
        std::cout << "Unable to open file";
}

TEST(tests_bscurve, approx_refined_nurbs)
{

    std::string line;
    std::ifstream myfile("../tests/in/e1098.dat");
    // std::ifstream myfile("../tests/in/e817.dat");
    if (myfile.is_open())
    {
        std::vector<std::array<double, 3>> pts;
        getline(myfile, line);
        while (getline(myfile, line))
        {
            std::istringstream iss(line);
            std::string::size_type sz; // alias of size_t

            double x = std::stod(line, &sz);
            double y = std::stod(line.substr(sz));
            pts.push_back({x, y,1});
        }
        myfile.close();

        auto crv = gbs::approx(pts, 5, gbs::KnotsCalcMode::CHORD_LENGTH,true);
        std::cout << "n poles: " << crv.poles().size() << ", n_flat: " << crv.knotsFlats().size() <<std::endl;

        std::vector<Handle_Geom_Curve> crv_lst;
        crv_lst.push_back(occt_utils::BSplineCurve(crv));
        // GeomTools::Dump( occt_utils::BSplineCurve(crv), std::cout );
        occt_utils::to_iges(crv_lst, "approx_refined_nurbs.igs");

        auto res = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << res.d_avg << ", d_max:" << res.d_max << ", u_max:" << res.u_max << std::endl;

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_PC(crv, pnt, u0, 1e-6);
                u0 = res.u;
                // std::cout << gbs::norm(crv.value(u0)-pnt) << " , " << crv_lst.back()->Value(u0).Distance(occt_utils::point(pnt)) << std::endl;
            });

        std::vector<TopoDS_Shape> pt_lst(pts.size());
        std::transform(
            pts.begin(), pts.end(),
            pt_lst.begin(),
            [](const auto &p_) {
                return BRepBuilderAPI_MakeVertex (gp_Pnt(p_[0],p_[1],0.));
            });
        occt_utils::to_iges(pt_lst, "foilpoints_nurbs.igs");
    }
    else
        std::cout << "Unable to open file";
}