#include <doctest_gtest.hpp>
#include <gbs/bscapprox.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>

#ifdef GBS_USE_MODULES
    import vecop;
 #else
    #include <gbs/vecop.ixx>
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <gbs-render/vtkGbsRender.h>

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

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

        size_t p = 5; // degree
        size_t n = 10;// n poles
        auto crv = gbs::approx(pts, p, n, gbs::KnotsCalcMode::CHORD_LENGTH);

        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
        ASSERT_LT(d_avg, 1.e-4);
        ASSERT_LT(d_max, 1.e-3);
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

        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
        ASSERT_TRUE(d_avg<1.e-6);
        ASSERT_TRUE(d_max<1.e-5);

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res[0];
                // std::cout << gbs::norm(crv.value(u0)-pnt) << " , " << crv_lst.back()->Value(u0).Distance(occt_utils::point(pnt)) << std::endl;
            });
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
        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
        ASSERT_LT(d_avg, 1.e-4);
        ASSERT_LT(d_max, 1.e-3);
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

        auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts, crv);
        std::cout << "d_avg: " << d_avg << ", d_max:" << d_max << ", u_max:" << u_max << std::endl;
        ASSERT_TRUE(d_avg<1.e-6);
        ASSERT_TRUE(d_max<1.e-5);

        auto u0 = crv.knotsFlats().front();
        std::for_each(
            pts.begin(),
            pts.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res[0];
                // std::cout << gbs::norm(crv.value(u0)-pnt) << " , " << crv_lst.back()->Value(u0).Distance(occt_utils::point(pnt)) << std::endl;
            });
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
    double deviation{0.01};
    double tol{1e-6};
    size_t np{30};

    auto c1_3d_dp        = gbs::BSCurve3d_d(poles, k, p);
    auto c1_3d_dp_approx = gbs::approx(c1_3d_dp, deviation, tol,p, np);
    ASSERT_DOUBLE_EQ(c1_3d_dp.bounds()[0],c1_3d_dp_approx.bounds()[0]);
    ASSERT_DOUBLE_EQ(c1_3d_dp.bounds()[1],c1_3d_dp_approx.bounds()[1]);
    auto u_lst = gbs::deviation_based_params(c1_3d_dp,np, deviation);
    using gbs::operator-;
    for(auto u_ : u_lst)
    {
        ASSERT_LT(gbs::norm(c1_3d_dp(u_)-c1_3d_dp_approx(u_)),1.e-4);
    }
    if(PLOT_ON)
    {
        gbs::plot(
            c1_3d_dp,
            c1_3d_dp_approx
        );
    }
}

namespace
{
    // A smooth, multi-scale 2D point set: two sinusoids of different frequency. It is
    // not fit by the initial handful of poles, so refine_approx must iterate; yet it
    // is well-behaved enough that both the average and the maximum deviation shrink
    // monotonically with added knots, so a moderate tolerance is reachable well below
    // the pole cap. Worst-fit points migrate across the domain as knots are inserted,
    // exercising the "error near a knot" path (C2).
    std::vector<std::array<double, 2>> wavy_point_set(size_t n = 120)
    {
        std::vector<std::array<double, 2>> pts;
        pts.reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
            double t = double(i) / (n - 1);
            double y = 0.3 * std::sin(5.0 * t) + 0.15 * std::sin(13.0 * t);
            pts.push_back({t, y});
        }
        return pts;
    }

    // A 2D point set with a localized feature: a moderate Gaussian bump on a gentle
    // slope. refine_approx concentrates inserted knots around the bump, so the worst
    // residuals end up sitting NEXT TO those knots -- exactly the situation the C2
    // bug mishandled (it updated the max only for points well inside a span).
    std::vector<std::array<double, 2>> localized_point_set(size_t n = 120)
    {
        std::vector<std::array<double, 2>> pts;
        pts.reserve(n);
        for (size_t i = 0; i < n; ++i)
        {
            double t = double(i) / (n - 1);
            double y = 0.1 * t + 0.6 * std::exp(-std::pow((t - 0.5) / 0.08, 2));
            pts.push_back({t, y});
        }
        return pts;
    }

    // True deviation over ALL points, using the SAME metric as refine_approx
    // (parametric distance at the fitting parameter u[j]), returns {d_max, d_avg}.
    std::array<double, 2> true_deviation(
        const gbs::BSCurve<double, 2> &crv,
        const std::vector<std::array<double, 2>> &pts,
        const std::vector<double> &u)
    {
        using gbs::operator-;
        double d_max = 0., d_avg = 0.;
        for (size_t j = 0; j < pts.size(); ++j)
        {
            double d = gbs::norm(crv(u[j]) - pts[j]);
            d_avg += d;
            d_max = std::max(d_max, d);
        }
        d_avg /= pts.size();
        return {d_max, d_avg};
    }
}

// C1: the average-deviation stop criterion must reflect the TRUE mean over all
// points. Before the fix, `d_avg_` accumulated only record-breaking distances, so a
// large true mean was masked and refine stopped early. Here d_avg is the binding
// criterion (d_max is left loose): refine MUST keep inserting until the true mean is
// actually met.
TEST(tests_bscapprox, refine_approx_honors_average_tolerance_C1)
{
    auto pts = wavy_point_set();
    size_t p = 3;
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);

    const double req_d_max = 1.e-1; // loose -> the average criterion binds
    const double req_d_avg = 1.e-3;
    auto crv = gbs::approx(pts, u, p, true, req_d_max, req_d_avg, 200);

    auto [d_max, d_avg] = true_deviation(crv, pts, u);
    std::cout << "[C1] n_poles=" << crv.poles().size()
              << " d_max=" << d_max << " d_avg=" << d_avg << std::endl;

    // It did not stop at the initial pole count: it actually refined.
    ASSERT_GT(crv.poles().size(), p * 2);
    // It did not declare convergence before the true mean met the tolerance.
    ASSERT_LT(d_avg, req_d_avg);
}

// C2: the maximum-deviation stop criterion must reflect the TRUE max over all
// points, including those near a knot. Before the fix, `d_max_`/`u_max` were updated
// only for points "well inside" a span, so a large error next to a knot was ignored
// and convergence was falsely declared. Here d_max is the binding criterion.
TEST(tests_bscapprox, refine_approx_honors_max_tolerance_near_knots_C2)
{
    auto pts = localized_point_set();
    size_t p = 3;
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);

    const double req_d_max = 2.e-3;
    const double req_d_avg = 1.e-1; // loose -> the maximum criterion binds
    auto crv = gbs::approx(pts, u, p, true, req_d_max, req_d_avg, 200);

    auto [d_max, d_avg] = true_deviation(crv, pts, u);
    std::cout << "[C2] n_poles=" << crv.poles().size()
              << " d_max=" << d_max << " d_avg=" << d_avg << std::endl;

    // For this dataset the worst residuals sit next to inserted knots. The buggy
    // criterion tracked the max only for span-interior points, so it declared
    // convergence at ~0.0026 (true) while believing ~0.0012; the fix tracks the true
    // max over ALL points and refines until it is genuinely below req_d_max.
    ASSERT_GT(crv.poles().size(), p * 2);
    ASSERT_LT(d_max, req_d_max);
}

// C3: every point-set approximation entry rejects ill-posed pole counts
// consistently: n_poles < p+1 and n_poles > n_pts both throw.
TEST(tests_bscapprox, approx_rejects_illposed_pole_counts_C3)
{
    auto pts = wavy_point_set(16);
    size_t p = 3;
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    auto k_flat = gbs::build_simple_mult_flat_knots(u.front(), u.back(), 6, p);

    // n_poles < p + 1 -> rejected on every entry point
    ASSERT_THROW(gbs::approx(pts, p, p, u, true), std::length_error);            // dispatcher (fix_bound)
    ASSERT_THROW(gbs::approx(pts, p, p, u, false), std::length_error);           // dispatcher (free)
    ASSERT_THROW(gbs::approx(pts, p, p, gbs::KnotsCalcMode::CHORD_LENGTH), std::length_error); // public
    ASSERT_THROW(gbs::approx_bound_fixed(pts, p, p, u, k_flat), std::length_error);

    // n_poles > n_pts -> rejected on every entry point
    size_t too_many = pts.size() + 1;
    ASSERT_THROW(gbs::approx(pts, p, too_many, u, true), std::length_error);
    ASSERT_THROW(gbs::approx(pts, p, too_many, u, false), std::length_error);

    // a well-posed count still succeeds
    ASSERT_NO_THROW(gbs::approx(pts, p, size_t{6}, u, true));
}

// A1 (assembly): approx_bound_fixed now assembles the collocation rows with the
// banded one-pass fill_basis_row instead of the recursive basis_function. The
// produced poles must match the original dense/recursive assembly to ~1e-9.
TEST(tests_bscapprox, approx_bound_fixed_banded_matches_recursive_A1)
{
    using gbs::operator-;
    auto pts = wavy_point_set(20);
    size_t p = 3, n_poles = 8;
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    auto k_flat = gbs::build_simple_mult_flat_knots(u.front(), u.back(), n_poles, p);

    // Reference: the original DENSE assembly via the recursive basis_function,
    // same endpoint-fixed least-squares solve.
    int n_params = int(u.size());
    gbs::MatrixX<double> N(n_params - 2, int(n_poles) - 2);
    for (int i = 0; i < n_params - 2; ++i)
        for (int j = 0; j < int(n_poles) - 2; ++j)
            N(i, j) = gbs::basis_function(u[i + 1], size_t(j + 1), p, size_t{0}, k_flat);
    std::vector<double> Nbegin(n_params - 2), Nend(n_params - 2);
    for (int i = 0; i < n_params - 2; ++i)
    {
        Nbegin[i] = gbs::basis_function(u[i + 1], size_t{0}, p, size_t{0}, k_flat);
        Nend[i]   = gbs::basis_function(u[i + 1], n_poles - 1, p, size_t{0}, k_flat);
    }
    auto qr = N.colPivHouseholderQr();
    int n_pt = int(pts.size());
    std::vector<std::array<double, 2>> poles_ref(n_poles);
    gbs::VectorX<double> b(n_pt - 2);
    for (int d = 0; d < 2; ++d)
    {
        for (int i = 0; i < n_pt - 2; ++i)
            b(i) = pts[i + 1][d] - pts.front()[d] * Nbegin[i] - pts.back()[d] * Nend[i];
        auto x = qr.solve(b);
        for (int i = 1; i < int(n_poles) - 1; ++i)
            poles_ref[i][d] = x(i - 1);
    }
    poles_ref.front() = pts.front();
    poles_ref.back()  = pts.back();

    auto crv = gbs::approx_bound_fixed(pts, p, n_poles, u, k_flat);
    ASSERT_EQ(crv.poles().size(), poles_ref.size());
    for (size_t i = 0; i < poles_ref.size(); ++i)
        ASSERT_LT(gbs::norm(crv.poles()[i] - poles_ref[i]), 1e-9);
}

// A3 (structured solve): the plain curve approx now solves the banded normal
// equations (NᵀN, SPD) via Cholesky instead of a dense colPivHouseholderQr on the
// rectangular system. With 150 points it crosses approx_sparse_threshold, so the
// SPARSE SimplicialLDLT branch runs. The poles must match the dense QR LS minimizer
// to ~1e-9.
TEST(tests_bscapprox, approx_plain_structured_matches_dense_qr_A3)
{
    using gbs::operator-;
    auto pts = wavy_point_set(150); // > approx_sparse_threshold -> sparse Cholesky path
    size_t p = 3, n_poles = 30;
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    auto k_flat = gbs::build_simple_mult_flat_knots(u.front(), u.back(), n_poles, p);

    // Dense QR reference on the rectangular collocation system.
    gbs::MatrixX<double> N(Eigen::Index(u.size()), Eigen::Index(n_poles));
    gbs::build_poles_matrix<double, 1>(k_flat, u, p, n_poles, N);
    auto qr = N.colPivHouseholderQr();
    std::vector<std::array<double, 2>> poles_ref(n_poles);
    gbs::VectorX<double> b(Eigen::Index(pts.size()));
    for (int d = 0; d < 2; ++d)
    {
        for (size_t i = 0; i < pts.size(); ++i) b(Eigen::Index(i)) = pts[i][d];
        auto x = qr.solve(b);
        for (size_t i = 0; i < n_poles; ++i) poles_ref[i][d] = x(Eigen::Index(i));
    }

    auto crv = gbs::approx(pts, p, n_poles, u, k_flat);
    ASSERT_EQ(crv.poles().size(), poles_ref.size());
    for (size_t i = 0; i < n_poles; ++i)
        ASSERT_LT(gbs::norm(crv.poles()[i] - poles_ref[i]), 1e-9);
}

// PR4 dedup: the deviation-driven Curve overloads now share sample_by_deviation /
// refine_from_points. These overloads had no in-repo callers; this exercises each so
// the refactor (and their default n_poles paths) is actually covered. Each must
// preserve the source curve's bounds and approximate it within tolerance.
TEST(tests_bscapprox, curve_approx_overloads_smoke_PR4)
{
    using gbs::operator-;
    std::vector<double> k = {0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
    std::vector<std::array<double, 3>> poles =
        {{0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 1, 2}, {4, 0, 1}};
    size_t p = 3;
    gbs::BSCurve<double, 3> crv(poles, k, p);

    const double deviation = 0.01, tol = 1.e-5;
    const size_t np = 60;

    auto check = [&](const gbs::BSCurve<double, 3> &a, double max_dev)
    {
        ASSERT_DOUBLE_EQ(crv.bounds()[0], a.bounds()[0]);
        ASSERT_DOUBLE_EQ(crv.bounds()[1], a.bounds()[1]);
        auto us = gbs::deviation_based_params(crv, np, deviation);
        for (auto u_ : us)
            ASSERT_LT(gbs::norm(crv(u_) - a(u_)), max_dev);
    };

    // 443: explicit n_poles, refined to tol.
    check(gbs::approx(crv, deviation, tol, p, size_t{8}, np), 1.e-3);
    // 466: (u1,u2) overload (whole curve sampled), default n_poles, refined to tol.
    check(gbs::approx(crv, crv.bounds()[0], crv.bounds()[1], deviation, tol, p, np), 1.e-3);
    // 477: per-point transform (identity here), default n_poles, refined to tol.
    std::function<std::array<double, 3>(const std::array<double, 3> &)> identity =
        [](const std::array<double, 3> &x) { return x; };
    check(gbs::approx(crv, identity, deviation, tol, p, np), 1.e-3);
    // 512: uniformly sampled, no refine (coarser) -> bounds preserved, loose deviation.
    check(gbs::approx(crv, p, np), 1.e-1);
}

