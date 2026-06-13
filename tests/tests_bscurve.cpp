#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/curvereparametrized.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>

#include <numbers>
#include <chrono>
#include <limits>

#ifdef GBS_USE_MODULES
    import knots_functions;
    import vecop;
 #else
    #include <gbs/vecop.ixx>
    #include <gbs/knotsfunctions.ixx>
#endif

const double tol = 1e-10;
const double tol_confusion = 1e-7;


using gbs::operator-;

TEST(tests_bscurve, perf_long_double)
{
    using T = long double;
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    double u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<T,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = static_cast<T>(count) / static_cast<T>(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<T, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs took " << ms.count() << " ms\n";

}

TEST(tests_bscurve, perf_double)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
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
    size_t p = 2;
    double u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = double(count) / double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs took " << ms.count() << " ms\n";

}

TEST(tests_bscurve, perf_float)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    float u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = double(count) / double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs float took " << ms.count() << " ms\n";

}
TEST(tests_bscurve, curve_values)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    float u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);

    size_t n{1000};
    auto u_lst = gbs::make_range(k.front(), k.back(), n);
    const auto t1 = std::chrono::high_resolution_clock::now();
    auto pts = c1_3d_dp.values(u_lst);
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs float took " << ms.count() << " ms\n";
                  
    for (size_t i{}; i < n; i++)
    {
        ASSERT_LT( gbs::distance( c1_3d_dp(u_lst[i]), pts[i] ), tol);
    }

}
TEST(tests_bscurve, curve_parametrization)
{
    std::vector<std::array<double,3> > pt =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    auto k1 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::EQUALLY_SPACED);
    std::for_each(k1.begin(), k1.end(), [](auto k_) { printf("k=%f\n", k_); });

    auto k2 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::CHORD_LENGTH);
    std::for_each(k2.begin(), k2.end(), [](auto k_) { printf("k=%f\n", k_); });

    auto k3 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::CENTRIPETAL);
    std::for_each(k3.begin(), k3.end(), [](auto k_) { printf("k=%f\n", k_); });
}

TEST(tests_bscurve, curve_increase_degree)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);
    auto c2_3d_dp = gbs::BSCurve<float,3>(poles,k,p);
    c1_3d_dp.increaseDegree();
    auto u = gbs::make_range(k.front(),k.back(),100);
    for(auto u_ : u) 
    {
        ASSERT_LT(gbs::norm(c1_3d_dp.value(u_) - c2_3d_dp.value(u_)), 1e-6);
    }
}

TEST(tests_bscurve, curve_reparametrized)
{
    auto circle = gbs::build_circle<float,2>(1.);
    
    auto f_u = gbs::BSCfunction{ gbs::abs_curv<float,2,10>(circle,36*30) };

    gbs::CurveReparametrized<float,2> reparam(std::make_shared<gbs::BSCurveRational<float,2>>(circle),f_u);

    // Reparametrization goes through a float arc-length integration, so the
    // error budget is many epsilons of accumulated quadrature/interpolation
    // error, not single rounding. Express it relative to float epsilon to make
    // the precision basis explicit and portable: 1e-6 (~8*eps) was below the
    // achievable precision and failed on MSVC (passed on gcc/clang by luck of
    // rounding).
    const float tol = 1.e3f * std::numeric_limits<float>::epsilon(); // ~1.2e-4

    ASSERT_NEAR(reparam(std::numbers::pi/2.)[1], 1.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi/2.)[0], 0.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi   )[1], 0.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi   )[0],-1.,tol);
}

TEST(tests_bscurve, eval_except)
{
    auto crv =  gbs::build_segment<double,2>({0.,0.},{1.,0.});
    ASSERT_THROW(crv.value(2.), gbs::OutOfBoundsCurveEval<double>);
}

TEST(tests_bscurve, d_dm_vectorized)
{
    using T = double;
    std::vector<T> k = {0., 0., 0., 1., 2., 3., 4., 5., 5., 5.};
    std::vector<std::array<T,3>> poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    gbs::BSCurve<T,3> crv{poles, k, p};

    auto u_lst = gbs::make_range(k.front(), k.back(), 50);

    // The vectorized derivatives w.r.t. the curvilinear abscissa must match the
    // scalar overloads parameter by parameter.
    auto d1 = crv.d_dms(u_lst);
    auto d2 = crv.d_dm2s(u_lst);

    ASSERT_EQ(d1.size(), u_lst.size());
    ASSERT_EQ(d2.size(), u_lst.size());

    for (size_t i{}; i < u_lst.size(); ++i)
    {
        ASSERT_LT(gbs::norm(d1[i] - crv.d_dm(u_lst[i])),  tol);
        ASSERT_LT(gbs::norm(d2[i] - crv.d_dm2(u_lst[i])), tol);
    }
}

// Span reduction (summing only over the p+1 non-zero basis functions) must
// give the same result as the full-span sum for DERIVATIVE orders too, even
// when sampled exactly on interior knots of multiplicity > 1 (C^0 / reduced
// continuity). This requires find_span to use the half-open convention
// (k[s] <= u < k[s+1]) consistent with basis_function; a plain lower_bound
// returns the left span at multiple knots and yields the wrong one-sided
// derivative. Guards the curve eval path behind d_dm / d_dm2 and the
// equivalent surface path.
TEST(tests_bscurve, deriv_span_reduction_matches_full_span)
{
    using T = double;
    // Degree 3 with interior knot u=4 of multiplicity 3 (== p) => C^0 joint,
    // discontinuous derivatives there, like a grafted curve extension.
    std::vector<T> k = {0., 0., 0., 0., 1., 2., 4., 4., 4., 5., 6., 7., 7., 7., 7.};
    std::vector<std::array<T,3>> poles =
    {
        {0.,0.,0.}, {0.,1.,0.}, {1.,1.,0.}, {1.,2.,1.}, {2.,2.,2.},
        {3.,1.,1.}, {3.,4.,1.}, {4.,4.,2.}, {5.,3.,0.}, {6.,5.,1.}, {7.,4.,2.},
    };
    size_t p = 3;
    ASSERT_EQ(poles.size(), k.size() - p - 1);

    // Sample the grid plus every distinct knot value (hits the multiple knot).
    auto u_lst = gbs::make_range(k.front(), k.back(), 50);
    for (auto kv : k) u_lst.push_back(kv);

    for (size_t d : {size_t{0}, size_t{1}, size_t{2}, size_t{3}})
    {
        for (auto u : u_lst)
        {
            auto reduced = gbs::eval_value_decasteljau(u, k, poles, p, d, /*use_span_reduction=*/true);
            auto full    = gbs::eval_value_decasteljau(u, k, poles, p, d, /*use_span_reduction=*/false);
            ASSERT_LT(gbs::norm(reduced - full), tol) << "d=" << d << " u=" << u;
        }
    }
}

// Ground-truth guard for the allocation-free A2.3 evaluator: its value and every
// derivative order must match the textbook recursive Cox-de Boor basis summed
// over ALL poles, to machine precision, including at knots of multiplicity > 1.
// This is the regression net for issue #30 (the self-comparison above only checks
// that use_span_reduction is a no-op; both branches now share the same kernel).
TEST(tests_bscurve, eval_matches_recursive_basis_ground_truth)
{
    using T = double;
    for (size_t p : {size_t{1}, size_t{2}, size_t{3}, size_t{5}, size_t{7}})
    {
        // Clamped knot vector with an interior knot of multiplicity p (a C^0
        // joint, where one-sided derivatives are discontinuous).
        std::vector<T> k;
        for (size_t i = 0; i <= p; ++i) k.push_back(0.);
        k.push_back(1.); k.push_back(2.);
        for (size_t i = 0; i < p; ++i) k.push_back(3.); // interior mult == p
        k.push_back(4.); k.push_back(5.);
        for (size_t i = 0; i <= p; ++i) k.push_back(6.);

        const size_t n_poles = k.size() - p - 1;
        std::vector<std::array<T,3>> poles(n_poles);
        for (size_t i = 0; i < n_poles; ++i)
            poles[i] = { T(i), T((i * 7) % 11), T((i * 3) % 5) };

        auto u_lst = gbs::make_range(k.front(), k.back(), 73);
        for (auto kv : k) u_lst.push_back(kv); // include every distinct knot

        for (size_t d = 0; d <= p + 1; ++d) // include d > p (must be exact zero)
        {
            for (auto u : u_lst)
            {
                std::array<T,3> expected{0.,0.,0.};
                for (size_t i = 0; i < n_poles; ++i)
                {
                    T N = gbs::basis_function(u, i, p, d, k);
                    for (size_t c = 0; c < 3; ++c) expected[c] += N * poles[i][c];
                }
                auto got = gbs::eval_value_decasteljau(u, k, poles, p, d);
                ASSERT_LT(gbs::norm(got - expected), tol)
                    << "p=" << p << " d=" << d << " u=" << u;
            }
        }
    }
}