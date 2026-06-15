#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>
#include <gbs-render/vtkfunctionrender.h>
#include <chrono>
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