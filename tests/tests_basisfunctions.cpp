#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/basisfunctions.h>
#include <gbs/maths.h>
#include <gbs-render/vtkfunctionrender.h>
#include <chrono>
#include <algorithm>
#include <iterator>

const double tol = 1e-10;
using gbs::operator-;

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

TEST(tests_basis_functions, eval_span)
{
    {
        std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
        auto p = 2;
        auto n = k.size() - p - 1;
        auto u = 0.5;
        auto it = std::next(k.begin(),2);
        for (int i = 0; i < n; i++)
        {
            auto it1 = gbs::find_span(n, p, u, k);
            auto it2 = gbs::find_span2(n, p, u, k);
            if (i == n - 1)
            {
                ASSERT_EQ(it1, std::next(it, -1));
                ASSERT_EQ(it2, std::next(it, -1));
            }
            else
            {
                ASSERT_EQ(it1, it);
                ASSERT_EQ(it2, it);
            }
            std::cout << "index: " << it1 - k.begin() << " " << *it1 << " " << *(++it1) << std::endl;
            std::cout << "index: " << it2 - k.begin() << " " << *it2 << " " << *(++it2) << std::endl;
            u += 1;
            it = std::next(it);
        }
    }

    {
        size_t p = 3;
        std::vector<double> k = {0., 0., 0., 0., 2., 2. , 6., 8., 8., 8., 8.};
        auto n = k.size() - p - 1;
        auto u = gbs::make_range<double>(k.front(), k.back(), 10);
        for( auto u_ : u){
            auto first= k.cbegin();
            auto it1  = gbs::find_span(n, p, u_, k);
            auto it2  = gbs::find_span2(n, p, u_, k);
            auto span1 =  std::distance( first, it1 );
            auto span2 =  std::distance( first, it2 );
                 span1 =  std::min<size_t>( span1, k.size() - p - 2);
                 span2 =  std::min<size_t>( span2, k.size() - p - 2);
            ASSERT_EQ(span1, span2);
            ASSERT_GE(u_, k[span1]);
            if(std::abs(u_-k.back()) > gbs::knot_eps<double>)
                ASSERT_LT(u_, k[span1+1]);
            for (int i = 0; i < n; i++)
            {
                auto N = gbs::basis_function( u_, i, p, 0, k);
                if(i<span1-p || i > span1)
                    ASSERT_NEAR( N, 0., 1e-10);
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

    const auto t1 = std::chrono::high_resolution_clock::now();
    count = 100000;
    while(count)
    {
        u = (rand() % 1000) / 999.;
        gbs::find_span2(n,p,u,k1);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << " took " << ms.count() << " ms\n";
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
        auto i =  std::min<size_t>( find_span2(n_poles, p, u_, k) - k.begin(), k.size() - p - 2);
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


    plot_basis_funcs(p, u, N_plot);
}

TEST(tests_basis_functions, basis_funcs_perf)
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