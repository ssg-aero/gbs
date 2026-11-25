#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs-render/vtkfunctionrender.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iterator>
#include <algorithm>

#ifdef GBS_USE_MODULES
    import math;
    import basis_functions;
#endif


const double tol = gbs::knot_eps<double>;
using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

TEST(tests_basis_functions, eval_basis)
{
    std::vector<double> k1 = {0.,0.,0.,1.,1.,1.};
    auto u1_test = gbs::make_range<double>(0.,1.,10);
    size_t p = 1;
    size_t i = 0;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                0.,
                tol
                );
    }  
    );

    p = 1;
    i = 1;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                1-u,
                tol
                );
    }  
    );

    i = 2;
    p = 1;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                u,
                tol
                );
    }  
    );

    i = 3;
    p = 1;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                0.,
                tol
                );
    }  
    );

    i = 0;
    p = 2;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                (1-u)*(1-u),
                tol
                );
    }  
    );

    i = 0;
    p = 2;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                (1-u)*(1-u),
                tol
                );
    }  
    );

    i = 1;
    p = 2;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                2*u*(1-u),
                tol
                );
    }  
    );

    i = 2;
    p = 2;

    std::for_each( u1_test.begin(), u1_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k1.begin(), i), p, k1.end()), 
                u*u,
                tol
                );
    }  
    );

    auto k2 = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    auto u2_test = gbs::make_range(0.,5.,30);

    i = 2;
    p = 0;

    std::for_each( u2_test.begin(), u2_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k2.begin(), i), p, k2.end()), 
                0<=u && u<1 ? 1. : 0.,
                tol
                );
    }  
    );  

    i = 4;
    p = 0;

    std::for_each( u2_test.begin(), u2_test.end(), [&](double u)
    {
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k2.begin(), i), p, k2.end()), 
                2.<=u && u<3. ? 1. : 0.,
                tol
                );
    }  
    );     


    i = 4;
    p = 1;

    std::for_each( u2_test.begin(), u2_test.end(), [&](double u)
    {
        double val = 0.;
        if(2.<=u && u<3.) val =u-2;
        if(3.<=u && u<4.) val =4-u;
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k2.begin(), i), p, k2.end()), 
                val,
                tol
                );
    }  
    );    


    i = 3;
    p = 2;

    std::for_each( u2_test.begin(), u2_test.end(), [&](double u)
    {
        double val = 0.;
        if(1.<=u && u<2.) val =0.5*(u-1.)*(u-1.);
        if(2.<=u && u<3.) val =-11./2+5.*u-u*u;
        if(3.<=u && u<4.) val =0.5*(4-u)*(4-u);
        ASSERT_NEAR( 
                gbs::basis_function(u, std::next(k2.begin(), i), p, k2.end()), 
                val,
                tol
                );
    }  
    );  
}

// Pins find_span to Piegl & Tiller A2.1: the returned span is always within
// [p, n-1] (n = number of poles), so the downstream min(i, k.size()-p-2) clamp
// is no longer needed (#42) and evaluation never runs past the last pole at
// u == knots.back() (#46). Covers the last knot and an interior knot of
// multiplicity > 1.
TEST(tests_basis_functions, eval_span)
{
    // Distinct interior knots; n = 7 poles -> valid spans [p, n-1] = [2, 6].
    {
        std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
        size_t p = 2;
        size_t n = k.size() - p - 1;
        auto span = [&](double u) { return size_t(gbs::find_span(n, p, u, k) - k.begin()); };

        ASSERT_EQ(span(0.5), 2u);
        ASSERT_EQ(span(1.5), 3u);
        ASSERT_EQ(span(2.5), 4u);
        ASSERT_EQ(span(3.5), 5u);
        ASSERT_EQ(span(4.5), 6u);
        // At interior knots the half-open convention lands on the right span.
        ASSERT_EQ(span(1.0), 3u);
        ASSERT_EQ(span(4.0), 6u);
        // At the last knot u == k.back(): span must be the LAST valid one (n-1),
        // not n. Returning n here was the off-by-one behind the #46 OOB.
        ASSERT_EQ(span(5.0), n - 1);
        // Out-of-domain above the last knot clamps to the last valid span.
        ASSERT_EQ(span(6.5), n - 1);
        // Whole domain: span stays within [p, n-1] (poles[span-p .. span] valid).
        for (auto u : gbs::make_range<double>(k.front(), k.back(), 100))
        {
            size_t s = span(u);
            CAPTURE(u);
            CAPTURE(s);
            ASSERT_GE(s, p);
            ASSERT_LE(s, n - 1);
        }
    }

    // Interior knot of multiplicity > 1 (C0 at u = 2). n = 7 poles, spans [3, 6].
    {
        size_t p = 3;
        std::vector<double> k = {0., 0., 0., 0., 2., 2., 6., 8., 8., 8., 8.};
        size_t n = k.size() - p - 1;
        auto span = [&](double u) { return size_t(gbs::find_span(n, p, u, k) - k.begin()); };

        ASSERT_EQ(span(1.5), 3u);
        // Half-open: at the double knot u = 2 we skip the empty span [2,2) and
        // land on [k[5], k[6]) = [2, 6) -> span 5.
        ASSERT_EQ(span(2.0), 5u);
        ASSERT_EQ(span(5.0), 5u);
        ASSERT_EQ(span(6.0), 6u);
        // Last knot: span n-1, no OOB.
        ASSERT_EQ(span(8.0), n - 1);

        // The non-zero basis functions of u are exactly N[span-p .. span]; every
        // other index must evaluate to 0 (consistency of span with the support).
        for (auto u_ : gbs::make_range<double>(k.front(), k.back(), 100))
        {
            size_t s = span(u_);
            ASSERT_GE(u_, k[s]);
            if (std::abs(u_ - k.back()) > gbs::knot_eps<double>)
                ASSERT_LT(u_, k[s + 1]);
            for (size_t i = 0; i < n; i++)
            {
                auto N = gbs::basis_function(u_, i, p, 0, k);
                CAPTURE(u_);
                CAPTURE(i);
                CAPTURE(s);
                if (i + p < s || i > s)
                    ASSERT_NEAR(N, 0., 1e-10);
            }
        }
    }
}

TEST(tests_basis_functions, eval_span_perf)
{

    std::vector<double> k1 = {0.,0.,0.,1,2,3,4,5.,5.,5.};   
    auto p = 2;
    auto n = k1.size() - p - 1;
	auto count = 100000;
	const auto t1_ref = std::chrono::high_resolution_clock::now();
    double u;
	while (count)
	{
		u = (rand() % 1000) / 999.;
		gbs::find_span(n,p,u,k1);
		count--;
	}
	const auto t2_ref = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> ms_ref = t2_ref - t1_ref;
	std::cout << std::fixed
		<< " took " << ms_ref.count() << " ms\n";

}

TEST(tests_basis_functions, eval_surf)
{
    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;
    auto u = 2.5;
    auto v = 1.;

    //Pij avec j inner loop
    const std::vector<std::array<double,4> > poles_t = {
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,6,8,2}, {12,24,12,6}, {4,6,0,2},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1}
                                                    };
    //Pij avec i inner loop
    std::vector<std::array<double,4> > poles(poles_t.size());
    int ni = 5 , nj =8;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            poles[j + nj * i] = poles_t[i + ni * j];
        }
    }
    auto pt2 = gbs::eval_value_decasteljau(u,v,ku,kv,poles,p,q);
    ASSERT_LT(gbs::norm(pt2 - std::array<double,4>({54/8.,98/8.,68/8.,27/8.})),tol);

    //Evalutaion rationnelle
    std::array<double, 3> r;
    std::transform(pt2.begin(), std::next(pt2.end(), - 1), r.begin(), [&pt2](const auto &p_) { return p_ / pt2.back(); });
    ASSERT_LT(gbs::norm(r-std::array<double,3>({2.,98/27.,68./27.})),tol);
}

TEST(tests_basis_functions, scale_weights)
{
    std::vector<double> k1 = {0., 0., 0., 0.25, 0.5 , 0.75, 1., 1., 1.};
    std::vector<std::array<double,4> > poles1 =
    {
        {0.,1.,0.,1},
        {0.,1.,0.,1.5},
        {1.,2.,0.,1.},
        {2.,3.,0.,0.5},
        {3.,3.,0.,1.2},
        {4.,2.,0.,1.},
    };
    size_t p1 = 2;
    gbs::BSCurveRational3d_d c1(poles1,k1,p1);
    gbs::scale_weights(poles1);
    gbs::BSCurveRational3d_d c2(poles1,k1,p1);

    ASSERT_LT(gbs::norm(c1.value(0.5)-c2.value(0.5)),tol);

}

TEST(tests_basis_functions, basis_funcs)
{
    using T = double;
    using namespace gbs;
    size_t p = 3;
    std::vector<std::vector<T>> N_plot;
    std::vector<T> k = {0., 0., 0., 0., 1., 5. , 6., 8., 8., 8., 8.};
    auto n_poles =k.size() - p - 1;
    size_t np = 1000;
    auto u = make_range<T>(0.,8.,np);
    for( auto u_ : u)
    {
        std::vector<T> N(p+1);
        size_t i = find_span(n_poles, p, u_, k) - k.begin();
        basis_funcs( i, p, u_, k, N );
        std::vector<T> N_full(n_poles,static_cast<T>(0.));
        for( size_t n{} ; n <=p ; n++ )
        {
            N_full[i-p] = N[n];
            i++;
        }
        N_plot.push_back( N_full );
        for( size_t n{}; n < n_poles; n++)
        {
            ASSERT_NEAR( basis_function( u_, n, p, 0, k), N_full[n], 1e-6 );
        }

    }

    if(PLOT_ON)
        plot_basis_funcs(p, u, N_plot);
}

TEST(tests_basis_functions, basis_funcs_perf)
{
    // Re-enabled after #42/#46: this benchmark sweeps u over [0, 8] inclusive
    // and used to abort at u = 8.0 (last knot) because find_span returned a span
    // one too high, indexing `poles` out of bounds under _GLIBCXX_ASSERTIONS.
    // With find_span faithful to A2.1 the boundary span is in range.
    using T = double;
    using namespace gbs;

    size_t p = 3;
    std::vector<T> k = {0., 0., 0., 0., 1., 5. , 6., 8., 8., 8., 8.};
    points_vector<T,3> poles =
    {
        {0.,1.,0.},
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,3.,0.},
        {3.,3.,0.},
        {4.,2.,0.},
    };
    size_t np = 1000000;
    auto u = make_range<T>(0.,8.,np);
    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for( auto u_ : u)
        {
            eval_value_deboor_cox(u_, k, poles, p);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms_ref = t2 - t1;
        std::cout << std::fixed
            << "eval_value_deboor_cox took " << ms_ref.count() << " ms\n";
    }
    
    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for( auto u_ : u)
        {
            eval_value_decasteljau(u_, k, poles, p);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms_ref = t2 - t1;
        std::cout << std::fixed
            << "eval_value_decasteljau took " << ms_ref.count() << " ms\n";
    }

}

// Cross-check the two evaluators agree over the whole domain, INCLUDING the
// last knot u = 8.0 where both used to index `poles` out of bounds (#46).
TEST(tests_basis_functions, eval_value_deboor_cox)
{
    using T = double;
    using namespace gbs;

    size_t p = 3;
    std::vector<T> k = {0., 0., 0., 0., 1., 5. , 6., 8., 8., 8., 8.};
    points_vector<T,3> poles =
    {
        {0.,1.,0.},
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,3.,0.},
        {3.,3.,0.},
        {4.,2.,0.},
    };
    size_t np = 100;
    auto u = make_range<T>(0.,8.,np);
    for( auto u_ : u)
    {
       auto p1 = eval_value_decasteljau(u_, k, poles, p);
       auto p2 = eval_value_deboor_cox(u_, k, poles, p);
       ASSERT_NEAR( norm(p1-p2), 0., 1e-6);
    }
}

// Non-regression for #42/#46: evaluate EXACTLY at u == knots.back() (and the
// surface corner u==ku.back(), v==kv.back()), value AND derivatives. Before the
// find_span A2.1 fix this indexed `poles` out of bounds and aborted under
// _GLIBCXX_ASSERTIONS. The fast evaluators are checked against the independent
// recursive basis_function reference (which uses no find_span / span reduction),
// and the clamped-endpoint identity C(u_max) == last pole is pinned directly.
TEST(tests_basis_functions, eval_at_last_knot)
{
    using T = double;
    using namespace gbs;

    // ---- Curve: clamped cubic, well-formed (6 poles, 10 knots) -------------
    {
        size_t p = 3;
        std::vector<T> k = {0., 0., 0., 0., 1., 2., 3., 3., 3., 3.};
        points_vector<T, 3> poles = {
            {0., 0., 0.}, {1., 2., 0.}, {2., 3., 1.}, {3., 1., 2.}, {4., 2., 1.}, {5., 0., 0.}};
        const size_t n_poles = poles.size();
        const T u_max = k.back();

        // Independent reference: full recursive sum over every pole.
        auto ref = [&](T u, size_t d) {
            std::array<T, 3> r{0., 0., 0.};
            for (size_t i = 0; i < n_poles; ++i)
            {
                const T N = basis_function(u, i, p, d, k);
                for (size_t c = 0; c < 3; ++c)
                    r[c] += N * poles[i][c];
            }
            return r;
        };

        // Endpoint identity: a clamped curve interpolates its last pole.
        ASSERT_LT(norm(eval_value_decasteljau(u_max, k, poles, p) - poles.back()), 1e-12);

        // Value + each derivative order at the last knot, vs. the reference.
        for (size_t d = 0; d <= p; ++d)
        {
            CAPTURE(d);
            ASSERT_LT(norm(eval_value_decasteljau(u_max, k, poles, p, d) - ref(u_max, d)), 1e-9);
        }

        // One-pass all-orders evaluator must agree order-by-order too.
        std::array<std::array<T, 3>, 4> CK;
        eval_ders_decasteljau<T, 3>(u_max, k, poles, p, p, CK.data());
        for (size_t d = 0; d <= p; ++d)
        {
            CAPTURE(d);
            ASSERT_LT(norm(CK[d] - ref(u_max, d)), 1e-9);
        }
    }

    // ---- Surface: clamped biquadratic, well-formed (4x4 poles) -------------
    {
        size_t p = 2, q = 2;
        std::vector<T> ku = {0., 0., 0., 1., 2., 2., 2.};
        std::vector<T> kv = {0., 0., 0., 1., 2., 2., 2.};
        const size_t nU = ku.size() - p - 1; // 4
        const size_t nV = kv.size() - q - 1; // 4
        std::vector<std::array<T, 3>> poles(nU * nV);
        for (size_t j = 0; j < nV; ++j)
            for (size_t i = 0; i < nU; ++i)
                poles[i + nU * j] = {T(i), T(j), T(i * j)};
        const T u_max = ku.back();
        const T v_max = kv.back();

        auto ref = [&](T u, T v, size_t du, size_t dv) {
            std::array<T, 3> r{0., 0., 0.};
            for (size_t j = 0; j < nV; ++j)
                for (size_t i = 0; i < nU; ++i)
                {
                    const T N = basis_function(u, i, p, du, ku) * basis_function(v, j, q, dv, kv);
                    for (size_t c = 0; c < 3; ++c)
                        r[c] += N * poles[i + nU * j][c];
                }
            return r;
        };

        // Corner identity: clamped surface interpolates its last corner pole.
        ASSERT_LT(norm(eval_value_decasteljau(u_max, v_max, ku, kv, poles, p, q) - poles.back()), 1e-12);

        // Value + mixed derivatives at the corner, vs. the reference.
        for (size_t du = 0; du <= p; ++du)
            for (size_t dv = 0; dv <= q; ++dv)
            {
                CAPTURE(du);
                CAPTURE(dv);
                ASSERT_LT(norm(eval_value_decasteljau(u_max, v_max, ku, kv, poles, p, q, du, dv) - ref(u_max, v_max, du, dv)), 1e-9);
            }

        // One-pass mixed-derivative evaluator must agree on every (du,dv) too.
        std::array<std::array<T, 3>, 9> SKL; // (p+1)*(q+1)
        eval_ders_decasteljau<T, 3>(u_max, v_max, ku, kv, poles, p, q, p, q, SKL.data());
        for (size_t du = 0; du <= p; ++du)
            for (size_t dv = 0; dv <= q; ++dv)
            {
                CAPTURE(du);
                CAPTURE(dv);
                ASSERT_LT(norm(SKL[du * (q + 1) + dv] - ref(u_max, v_max, du, dv)), 1e-9);
            }
    }
}

// ===========================================================================
// Coverage edge cases (issue #50): low-degree evaluation, degenerate weight
// and out-of-domain find_span. These exercise branches the baseline left
// uncovered.
// ===========================================================================

// Degree 0 and degree 1 were never evaluated through the fast path: p == 0 is
// a piecewise-constant (value == pole of the half-open span, derivative zero),
// p == 1 is the control polyline (interpolates poles at knots, segment slope as
// first derivative). Exercises basis_funcs/ders_basis_funs at their trivial
// loop bounds.
TEST(tests_basis_functions, eval_degree_0_and_1)
{
    using T = double;
    using namespace gbs;

    // ---- degree 0: 3 poles, knots {0,1,2,3}; C(u) == pole of its span -------
    {
        size_t p = 0;
        std::vector<T> k = {0., 1., 2., 3.};
        points_vector<T, 2> poles = {{0., 0.}, {1., 10.}, {2., 20.}};
        ASSERT_LT(norm(eval_value_decasteljau<T, 2>(0.5, k, poles, p) - poles[0]), 1e-12);
        ASSERT_LT(norm(eval_value_decasteljau<T, 2>(1.5, k, poles, p) - poles[1]), 1e-12);
        ASSERT_LT(norm(eval_value_decasteljau<T, 2>(2.5, k, poles, p) - poles[2]), 1e-12);
        // a derivative of a piecewise-constant (d > p) is identically zero
        ASSERT_LT(norm(eval_value_decasteljau<T, 2>(1.5, k, poles, p, 1)), 1e-12);
    }

    // ---- degree 1: clamped linear, 3 poles; the control polyline -----------
    {
        size_t p = 1;
        std::vector<T> k = {0., 0., 1., 2., 2.};
        points_vector<T, 2> poles = {{0., 0.}, {1., 2.}, {2., 0.}};
        BSCurve<T, 2> crv(poles, k, p);
        ASSERT_LT(norm(crv.value(0.) - poles.front()), 1e-12);
        ASSERT_LT(norm(crv.value(1.) - poles[1]), 1e-12);
        ASSERT_LT(norm(crv.value(2.) - poles.back()), 1e-12);
        // midpoint of the first segment is the average of its endpoints
        ASSERT_LT(norm(crv.value(0.5) - point<T, 2>{0.5, 1.}), 1e-12);
        // first derivative on [0,1) is the segment slope poles[1]-poles[0]
        ASSERT_LT(norm(crv.value(0.5, 1) - point<T, 2>{1., 2.}), 1e-12);
    }
}

// A homogeneous point with zero weight projects to the origin (degenerate
// guard in weight_projection); a non-zero weight divides through.
TEST(tests_basis_functions, weight_projection_zero_weight_is_origin)
{
    using T = double;
    std::array<T, 4> hp0 = {3., 4., 5., 0.};
    auto r0 = gbs::weight_projection(hp0); // -> std::array<T,3>
    ASSERT_EQ(r0[0], 0.);
    ASSERT_EQ(r0[1], 0.);
    ASSERT_EQ(r0[2], 0.);

    std::array<T, 4> hp = {3., 4., 5., 2.};
    auto r = gbs::weight_projection(hp);
    ASSERT_NEAR(r[0], 1.5, 1e-12);
    ASSERT_NEAR(r[1], 2.0, 1e-12);
    ASSERT_NEAR(r[2], 2.5, 1e-12);
}

// find_span on an out-of-domain parameter below the first knot must not loop
// forever: the "no progress" break returns the lowest valid span p (the
// closest one), instead of running off the knot vector.
TEST(tests_basis_functions, find_span_out_of_domain_below)
{
    using T = double;
    std::vector<T> k = {0., 0., 0., 1., 2., 3., 3., 3.};
    size_t p = 2;
    size_t n = k.size() - p - 1; // 5 poles, valid spans [2, 4]
    size_t s = size_t(gbs::find_span(n, p, T(-1.0), k) - k.begin());
    ASSERT_EQ(s, p);
}
// ===========================================================================
// Degrees above bspline_stack_max_degree (24): the evaluators run the same
// A2.2/A2.3 on heap buffers, O(p^2) per evaluation. They fell back on the
// recursive basis_function, whose cost doubles with each degree (0.4 s per
// evaluation at degree 26, about two hours at degree 40).
// ===========================================================================

namespace
{
    // Independent reference: de Boor's algorithm (affine combinations of the
    // poles, no basis function) on the curve of the d-th derivative, whose
    // poles come from d hodographs: P'_i = p (P_{i+1} - P_i) / (k_{i+p+1} - k_{i+1}),
    // knots without the first and the last.
    template <size_t dim>
    auto deboor_reference(double u, std::vector<double> k, gbs::points_vector<double, dim> P, size_t p, size_t d) -> std::array<double, dim>
    {
        for (size_t r = 0; r < d; ++r)
        {
            if (p == 0)
                return std::array<double, dim>{};
            gbs::points_vector<double, dim> Q(P.size() - 1);
            for (size_t i = 0; i + 1 < P.size(); ++i)
            {
                const double den = k[i + p + 1] - k[i + 1];
                for (size_t c = 0; c < dim; ++c)
                    Q[i][c] = den > 0. ? double(p) * (P[i + 1][c] - P[i][c]) / den : 0.;
            }
            k = std::vector<double>(k.begin() + 1, k.end() - 1);
            P = Q;
            --p;
        }
        const size_t n = P.size();
        size_t s = p; // k[s] <= u < k[s+1], the last span at u == k.back()
        while (s + 1 < n && u >= k[s + 1])
            ++s;
        gbs::points_vector<double, dim> D(P.begin() + std::ptrdiff_t(s - p), P.begin() + std::ptrdiff_t(s + 1));
        for (size_t r = 1; r <= p; ++r)
            for (size_t j = p; j >= r; --j)
            {
                const double a = (u - k[j + s - p]) / (k[j + 1 + s - r] - k[j + s - p]);
                for (size_t c = 0; c < dim; ++c)
                    D[j][c] = (1. - a) * D[j - 1][c] + a * D[j][c];
            }
        return D[p];
    }

    // Clamped knots on [0, 2] with interior knots 0.5, 0.9 (double) and 1.4.
    auto high_degree_knots(size_t p) -> std::vector<double>
    {
        std::vector<double> k(p + 1, 0.);
        for (double x : {0.5, 0.9, 0.9, 1.4})
            k.push_back(x);
        k.insert(k.end(), p + 1, 2.);
        return k;
    }
}

// Degree 26 against the recursive basis_function, the former fallback: values,
// derivatives and the basis functions themselves, on a curve and a surface.
// One Bezier span, where the recursion costs about 0.4 s per evaluation: its
// basis vectors are computed once and shared.
TEST(tests_basis_functions, high_degree_matches_recursive_basis)
{
    using T = double;
    using namespace gbs;
    const size_t p = 26;
    static_assert(26 > bspline_stack_max_degree);
    std::vector<T> k(p + 1, 0.);
    k.insert(k.end(), p + 1, 1.);
    points_vector<T, 3> P(p + 1);
    for (size_t i = 0; i <= p; ++i)
        P[i] = {T(i) / T(p), std::sin(0.7 * T(i)), std::cos(0.3 * T(i))};

    const T u = 0.3;
    std::vector<std::vector<T>> Nref(3, std::vector<T>(p + 1));
    for (size_t d = 0; d <= 2; ++d)
        for (size_t i = 0; i <= p; ++i)
            Nref[d][i] = basis_function(u, i, p, d, k);
    auto rel = [](const std::array<T, 3> &a, const std::array<T, 3> &b) { return norm(a - b) / std::max(T(1), norm(b)); };

    // the basis functions, function by function (the span of u is p)
    std::vector<T> N(p + 1), D(3 * (p + 1));
    basis_funcs(p, p, u, k, N.data());
    ders_basis_funs(p, p, 2, u, k, D.data());
    for (size_t d = 0; d <= 2; ++d)
    {
        T scale = 0;
        for (T x : Nref[d])
            scale = std::max(scale, std::abs(x));
        for (size_t i = 0; i <= p; ++i)
        {
            CAPTURE(d);
            CAPTURE(i);
            if (d == 0)
                ASSERT_LT(std::abs(N[i] - Nref[0][i]), 1e-12 * scale);
            ASSERT_LT(std::abs(D[d * (p + 1) + i] - Nref[d][i]), 1e-12 * scale);
        }
    }

    // curve: value, each derivative, all orders in one pass
    std::array<std::array<T, 3>, 3> CK;
    eval_ders_decasteljau<T, 3>(u, k, P, p, 2, CK.data());
    for (size_t d = 0; d <= 2; ++d)
    {
        CAPTURE(d);
        std::array<T, 3> ref{0., 0., 0.};
        for (size_t i = 0; i <= p; ++i)
            for (size_t c = 0; c < 3; ++c)
                ref[c] += Nref[d][i] * P[i][c];
        ASSERT_LT(rel(eval_value_decasteljau<T, 3>(u, k, P, p, d), ref), 1e-12);
        ASSERT_LT(rel(CK[d], ref), 1e-12);
    }

    // surface: degree 26 in u (the same knots), 2 in v
    const size_t q = 2;
    const std::vector<T> kv = {0., 0., 0., 0.5, 1., 1., 1.};
    const size_t nU = p + 1, nV = kv.size() - q - 1;
    std::vector<std::array<T, 3>> S(nU * nV);
    for (size_t j = 0; j < nV; ++j)
        for (size_t i = 0; i < nU; ++i)
            S[i + nU * j] = {T(i) / T(p), T(j), std::sin(0.5 * T(i)) * std::cos(0.9 * T(j))};
    const T v = 0.62;
    std::array<std::array<T, 3>, 9> SKL;
    eval_ders_decasteljau<T, 3>(u, v, k, kv, S, p, q, 2, 2, SKL.data());
    for (size_t du = 0; du <= 2; ++du)
        for (size_t dv = 0; dv <= 2; ++dv)
        {
            CAPTURE(du);
            CAPTURE(dv);
            std::array<T, 3> ref{0., 0., 0.};
            for (size_t j = 0; j < nV; ++j)
            {
                const T Nv = basis_function(v, j, q, dv, kv);
                for (size_t i = 0; i < nU; ++i)
                    for (size_t c = 0; c < 3; ++c)
                        ref[c] += Nref[du][i] * Nv * S[i + nU * j][c];
            }
            ASSERT_LT(rel(eval_value_decasteljau<T, 3>(u, v, k, kv, S, p, q, du, dv), ref), 1e-12);
            ASSERT_LT(rel(SKL[du * 3 + dv], ref), 1e-12);
        }
}

// Degrees 30 and 40 against de Boor's algorithm, at 101 points across four
// spans (knots and ends included): values and derivatives up to order 3 by
// each evaluator, and the partition of unity. Milliseconds in all; the
// recursion took minutes at degree 30, hours at degree 40.
TEST(tests_basis_functions, high_degree_curves_against_de_boor)
{
    using T = double;
    using namespace gbs;
    for (size_t p : {30, 40})
    {
        CAPTURE(p);
        const auto k = high_degree_knots(p);
        const size_t np = k.size() - p - 1;
        points_vector<T, 3> P(np);
        for (size_t i = 0; i < np; ++i)
            P[i] = {std::cos(0.4 * T(i)), std::sin(0.9 * T(i)), 0.05 * T(i)};
        const BSCurve<T, 3> crv(P, k, p);

        std::vector<T> us;
        for (int s = 0; s <= 100; ++s)
            us.push_back(2. * s / 100.);
        std::array<std::vector<std::array<T, 3>>, 4> ref;
        std::array<T, 4> scale{};
        for (size_t d = 0; d <= 3; ++d)
            for (T u : us)
            {
                ref[d].push_back(deboor_reference<3>(u, k, P, p, d));
                scale[d] = std::max(scale[d], norm(ref[d].back()));
            }

        for (size_t s = 0; s < us.size(); ++s)
        {
            const T u = us[s];
            CAPTURE(u);
            std::array<std::array<T, 3>, 4> CK;
            crv.derivatives(u, 3, CK.data());
            for (size_t d = 0; d <= 3; ++d)
            {
                CAPTURE(d);
                ASSERT_LT(norm(crv.value(u, d) - ref[d][s]), 1e-12 * scale[d]);
                ASSERT_LT(norm(CK[d] - ref[d][s]), 1e-12 * scale[d]);
            }
            std::vector<T> N(p + 1);
            basis_funcs(size_t(find_span(np, p, u, k) - k.begin()), p, u, k, N.data());
            T sum = 0;
            for (T x : N)
                sum += x;
            ASSERT_NEAR(sum, 1., 1e-13);
        }
    }
}

// Surfaces of degree 30 in u (3 in v) and of degree 30 in v (2 in u): each
// evaluator against de Boor along the high degree, the recursive basis (cheap
// at low degree) across it.
TEST(tests_basis_functions, high_degree_surfaces_against_de_boor)
{
    using T = double;
    using namespace gbs;
    const std::vector<T> k_low = {0., 0., 0., 0., 1., 2., 2., 2., 2.}; // degree 3
    for (bool high_in_u : {true, false})
    {
        CAPTURE(high_in_u);
        const size_t p = high_in_u ? 30 : 3;
        const size_t q = high_in_u ? 3 : 30;
        const auto ku = high_in_u ? high_degree_knots(30) : k_low;
        const auto kv = high_in_u ? k_low : high_degree_knots(30);
        const size_t nU = ku.size() - p - 1, nV = kv.size() - q - 1;
        std::vector<std::array<T, 3>> S(nU * nV);
        for (size_t j = 0; j < nV; ++j)
            for (size_t i = 0; i < nU; ++i)
                S[i + nU * j] = {T(i) / T(nU), T(j) / T(nV), std::sin(0.5 * T(i)) * std::cos(0.7 * T(j))};
        const BSSurface<T, 3> srf(S, ku, kv, p, q);

        auto reference = [&](T u, T v, size_t du, size_t dv) {
            std::array<T, 3> r{0., 0., 0.};
            if (high_in_u)
                for (size_t j = 0; j < nV; ++j)
                {
                    points_vector<T, 3> row(S.begin() + std::ptrdiff_t(nU * j), S.begin() + std::ptrdiff_t(nU * (j + 1)));
                    const auto C = deboor_reference<3>(u, ku, row, p, du);
                    const T Nv = basis_function(v, j, q, dv, kv);
                    for (size_t c = 0; c < 3; ++c)
                        r[c] += Nv * C[c];
                }
            else
                for (size_t i = 0; i < nU; ++i)
                {
                    points_vector<T, 3> col(nV);
                    for (size_t j = 0; j < nV; ++j)
                        col[j] = S[i + nU * j];
                    const auto C = deboor_reference<3>(v, kv, col, q, dv);
                    const T Nu = basis_function(u, i, p, du, ku);
                    for (size_t c = 0; c < 3; ++c)
                        r[c] += Nu * C[c];
                }
            return r;
        };

        for (int a = 0; a <= 8; ++a)
            for (int b = 0; b <= 8; ++b)
            {
                const T u = 2. * a / 8., v = 2. * b / 8.;
                CAPTURE(u);
                CAPTURE(v);
                std::array<std::array<T, 3>, 9> SKL;
                srf.derivatives(u, v, 2, 2, SKL.data());
                for (size_t du = 0; du <= 2; ++du)
                    for (size_t dv = 0; dv <= 2; ++dv)
                    {
                        CAPTURE(du);
                        CAPTURE(dv);
                        const auto ref = reference(u, v, du, dv);
                        const T tol = 1e-12 * std::max(T(1), norm(ref));
                        ASSERT_LT(norm(srf.value(u, v, du, dv) - ref), tol);
                        ASSERT_LT(norm(SKL[du * 3 + dv] - ref), tol);
                    }
            }
    }
}

// A rational curve of degree 30 (the homogeneous evaluation of degree 30 in 4
// dimensions), and the derivatives of a rational curve to an order above
// bspline_stack_max_degree: one pass on a heap table, where it was the
// per-order recursion (2^n evaluations, so far out of reach at n = 30).
TEST(tests_basis_functions, high_degree_rational)
{
    using T = double;
    using namespace gbs;
    {
        const size_t p = 30;
        const auto k = high_degree_knots(p);
        const size_t np = k.size() - p - 1;
        points_vector<T, 3> P(np);
        std::vector<T> w(np);
        points_vector<T, 4> Pw(np);
        for (size_t i = 0; i < np; ++i)
        {
            P[i] = {std::cos(0.4 * T(i)), std::sin(0.9 * T(i)), 0.05 * T(i)};
            w[i] = 1. + 0.5 * std::sin(1.3 * T(i));
            Pw[i] = {P[i][0] * w[i], P[i][1] * w[i], P[i][2] * w[i], w[i]};
        }
        const BSCurveRational<T, 3> crv(P, k, w, p);
        for (int s = 0; s <= 40; ++s)
        {
            const T u = 2. * s / 40.;
            CAPTURE(u);
            // C = A / w, C' = (A' - w' C) / w
            const auto A0 = deboor_reference<4>(u, k, Pw, p, 0);
            const auto A1 = deboor_reference<4>(u, k, Pw, p, 1);
            std::array<T, 3> C0, C1;
            for (size_t c = 0; c < 3; ++c)
            {
                C0[c] = A0[c] / A0[3];
                C1[c] = (A1[c] - A1[3] * C0[c]) / A0[3];
            }
            std::array<std::array<T, 3>, 2> CK;
            crv.derivatives(u, 1, CK.data());
            ASSERT_LT(norm(crv.value(u) - C0), 1e-12 * std::max(T(1), norm(C0)));
            ASSERT_LT(norm(CK[0] - C0), 1e-12 * std::max(T(1), norm(C0)));
            ASSERT_LT(norm(crv.value(u, 1) - C1), 1e-12 * std::max(T(1), norm(C1)));
            ASSERT_LT(norm(CK[1] - C1), 1e-12 * std::max(T(1), norm(C1)));
        }
    }
    {
        // a rational cubic, derivatives to order 30 (> 24)
        const std::vector<T> k = {0., 0., 0., 0., 0.4, 1., 1., 1., 1.};
        const points_vector<T, 2> P = {{0., 0.}, {1., 2.}, {2., -1.}, {3., 1.}, {4., 0.}};
        const std::vector<T> w = {1., 0.6, 1.5, 0.8, 1.};
        const BSCurveRational<T, 2> crv(P, k, w, 3);
        const size_t n = 30;
        static_assert(30 > bspline_stack_max_degree);
        for (T u : {0., 0.2, 0.4, 0.75, 1.})
        {
            CAPTURE(u);
            std::array<std::array<T, 2>, n + 1> CK;
            crv.derivatives(u, n, CK.data());
            // the orders the stack table also computes are the same, to the bit
            std::array<std::array<T, 2>, 7> CK6;
            crv.derivatives(u, 6, CK6.data());
            for (size_t d = 0; d <= 6; ++d)
            {
                CAPTURE(d);
                ASSERT_EQ(CK[d][0], CK6[d][0]);
                ASSERT_EQ(CK[d][1], CK6[d][1]);
            }
            // the higher ones as the per-order quotient rule (2^d evaluations)
            for (size_t d : {8, 10})
            {
                CAPTURE(d);
                const auto ref = eval_rational_value_simple<T, 2>(u, k, crv.poles(), 3, d);
                ASSERT_LT(norm(CK[d] - ref), 1e-12 * std::max(T(1), norm(ref)));
            }
            // value(u, d) takes order d of the same pass (it was that recursion)
            for (size_t d : {1, 2, 10, 30})
            {
                CAPTURE(d);
                const auto x = crv.value(u, d);
                ASSERT_EQ(x[0], CK[d][0]);
                ASSERT_EQ(x[1], CK[d][1]);
            }
            for (size_t d = 0; d <= n; ++d)
                ASSERT_TRUE(std::isfinite(CK[d][0]) && std::isfinite(CK[d][1]));
        }
    }
}
