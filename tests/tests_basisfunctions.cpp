#include <gtest/gtest.h>

#include <gbs/basisfunctions.h>
#include <gbs/math.h>
#include <gbs-occt/geomprim.h>
#include <gbs-occt/curvesbuild.h>

#include <chrono>

const double tol = 1e-10;
// using namespace gbs;
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

TEST(tests_basis_functions, eval_curve)
{
    std::vector<double> k1 = {0.,0.,0.,1.,1.,1.};
    auto u = 0.3;
    auto p = 2;

    auto crv = occt_utils::BSplineCurve(
        {gp_Pnt{0., 0., 0.},
         gp_Pnt{1., 0., 0.},
         gp_Pnt{1., 1., 1.}},
        {0., 1.},
        {3, 3},
        2);
    auto pt1 = crv->Value(u);
    const std::vector<std::array<double,3> > poles = {{0., 0., 0.}, {1., 0., 0.}, {1., 1., 1.}};
    auto pt2 = occt_utils::point(gbs::eval_value_simple(u, k1, poles , 2));
    auto v1 = crv->DN(u,1);
    auto v2 = occt_utils::vector(gbs::eval_value_simple(u, k1, poles , 2,1));

    ASSERT_LT(pt1.Distance(pt2),tol);
    ASSERT_LT((v1-v2).Magnitude(),tol);


}

TEST(tests_basis_functions, eval_curve_perf)
{

    std::vector<double> k1 = {0.,0.,0.,1.,1.,1.};
    const std::vector<std::array<double,3> > poles = {{0., 0., 0.}, {1., 0., 0.}, {1., 1., 1.}};
    auto u = 0.3;
    auto p = 2;

    auto crv = occt_utils::BSplineCurve(
        {gp_Pnt{0., 0., 0.},
         gp_Pnt{1., 0., 0.},
         gp_Pnt{1., 1., 1.}},
        {0., 1.},
        {3, 3},
        2);

	auto count = 100000;
	const auto t1_ref = std::chrono::high_resolution_clock::now();
	while (count)
	{
		u = (rand() % 1000) / 999.;
		crv->Value(u);
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
        gbs::eval_value_simple(u, k1, poles , 2);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << " took " << ms.count() << " ms\n";
}

TEST(tests_basis_functions, eval_span)
{
    std::vector<double> k1 = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    auto p = 2;
    auto n = k1.size() - p - 1;
    auto u = 0.5;
    auto it = std::next(k1.begin(),2);
    for (int i = 0; i < n; i++)
    {
        auto it1 = gbs::find_span(n, p, u, k1);
        auto it2 = gbs::find_span2(n, p, u, k1);
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
        std::cout << "index: " << it1 - k1.begin() << " " << *it1 << " " << *(++it1) << std::endl;
        std::cout << "index: " << it2 - k1.begin() << " " << *it2 << " " << *(++it2) << std::endl;
        u += 1;
        it = std::next(it);
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
    auto pt2 = gbs::eval_value_simple(u,v,ku,kv,poles,p,q);
    ASSERT_LT(gbs::norm(pt2 - std::array<double,4>({54/8.,98/8.,68/8.,27/8.})),tol);

    //Evalutaion rationnelle
    std::array<double, 3> r;
    std::transform(pt2.begin(), std::next(pt2.end(), - 1), r.begin(), [&pt2](const auto &p_) { return p_ / pt2.back(); });
    ASSERT_LT(gbs::norm(r-std::array<double,3>({2.,98/27.,68./27.})),tol);
}