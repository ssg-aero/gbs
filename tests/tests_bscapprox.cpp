#include <gtest/gtest.h>
#include <gbs/bscapprox.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>
#include <gbs/vecop.h>
#include <iostream>
#include <fstream>
#include <gbs-occt/export.h>
#include <gbs-occt/curvesbuild.h>
#include <gbs-render/vtkcurvesrender.h>
#include <Geom2d_Curve.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

TEST(tests_bscapprox, approx_simple)
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
        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
    }
    else
        std::cout << "Unable to open file";
}

TEST(tests_bscapprox, approx_refined)
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

        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res[0];
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

TEST(tests_bscapprox, approx_simple_nurbs)
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
        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
    }
    else
        std::cout << "Unable to open file";
}

TEST(tests_bscapprox, approx_refined_nurbs)
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

        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res[0];
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

TEST(tests_bscapprox, approx_curve)
{
    auto circle = gbs::build_circle<double,2>(1.);
    size_t p = 5;
    size_t n_poles = 36;
    double dev = 0.05;
    auto circle_approx1 = gbs::approx(circle,dev,n_poles,p,gbs::KnotsCalcMode::CENTRIPETAL);
    auto u = gbs::make_range(circle_approx1.bounds() ,30);
    for(auto u_ : u)
    {
        auto pt = circle_approx1(u_);
        ASSERT_LT(gbs::extrema_curve_point(circle,pt,1e-6)[1],1e-4);
    }

    auto circle_approx2 = gbs::approx(circle,dev,p,gbs::KnotsCalcMode::CHORD_LENGTH);
    u = gbs::make_range(circle_approx2.bounds() ,30);
    for(auto u_ : u)
    {
        auto pt = circle_approx2(u_);
        ASSERT_LT(gbs::extrema_curve_point(circle,pt,1e-6)[1],1e-3);
    }

    std::vector<double> k = {1., 1., 1., 1.5, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    p = 2;

    auto c1_3d_dp        = gbs::BSCurve3d_d(poles, k, p);
    auto c1_3d_dp_approx = gbs::approx(c1_3d_dp,0.01,1.e-6,p);
    ASSERT_DOUBLE_EQ(c1_3d_dp.bounds()[0],c1_3d_dp_approx.bounds()[0]);
    ASSERT_DOUBLE_EQ(c1_3d_dp.bounds()[1],c1_3d_dp_approx.bounds()[1]);
    auto u_lst = gbs::deviation_based_params(c1_3d_dp,30,0.01);
    using gbs::operator-;
    for(auto u_ : u_lst)
    {
        ASSERT_LT(gbs::norm(c1_3d_dp(u_)-c1_3d_dp_approx(u_)),1e-5);
    }

}

