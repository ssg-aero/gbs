#include <gtest/gtest.h>
#include <gbs/bssurf.h>
#include <gbs/vecop.h>
#include <gbs/bssinterp.h>
#include <gbs/bssapprox.h>
#include <gbs/render/vtkcurvesrender.h>

const double tol = 1e-10;
using gbs::operator-;
TEST(tests_bssurf, ctor)
{
    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;
    auto u = 2.5;
    auto v = 1.;


    std::vector<std::array<double,4> > poles_t = {
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

    gbs::BSSurface<double,4> srf(poles,ku,kv,p,q);
    gbs::BSSurfaceRational<double,3> srfNURBS(poles,ku,kv,p,q);

    auto pt1 = srf.value(u,v);
    ASSERT_LT(gbs::norm(pt1 - std::array<double,4>({54/8.,98/8.,68/8.,27/8.})),tol); // NURBS's Book



    auto ptr1 = srfNURBS.value(u,v);
    ASSERT_LT(gbs::norm(ptr1-std::array<double,3>({2.,98/27.,68./27.})),tol); // NURBS's Book

    auto polesU = srfNURBS.polesU(0);
    auto polesV = srfNURBS.polesV(0);

    // gbs::NURBSSurface<double,3> srf_nurbs(poles,ku,kv,p,q);
    // auto ptr1 = srf_nurbs.value2(u,v);
    // ASSERT_LT(gbs::norm(ptr1-std::array<double,3>({2.,98/27.,68./27.})),tol);


    // occt_utils::
}

TEST(tests_bssurf, increaseDegree)
{
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    gbs::points_vector<double,3> poles = 
    {                                       // ----U----
        {0,0,0},{1,0,0},                    // |
        {0,1,0},{1,1,1},                    // V
        {0,2,0},{1,2,0},                    // |
    };

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q);

    srf.increaseDegreeU();
    ASSERT_EQ(srf.degreeU(),p+1);
    srf.increaseDegreeU();
    ASSERT_EQ(srf.degreeU(),p+2);

    srf.increaseDegreeV();
    ASSERT_EQ(srf.degreeV(),q+1);

    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[0]-srf(0,0)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[1]-srf(1,0)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[4]-srf(0,1)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[5]-srf(1,1)), 0.
    );

    // gbs::plot(
    //     srf,
    //     srf.poles()
    // );

}

TEST(tests_bssurf, invertUV)
{
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    gbs::points_vector<double,3> poles = 
    {                                       // ----U----
        {0,0,0},{1,0,0},                    // |
        {0,1,0},{1,1,1},                    // V
        {0,2,0},{1,2,0},                    // |
    };

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q);
    // gbs::BSSurface srf_orig(poles,ku,kv,p,q);

    srf.invertUV();

    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[0]-srf(0,0)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[1]-srf(0,1)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[4]-srf(1,0)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[5]-srf(1,1)), 0.
    );

    // gbs::plot(
    //     srf,
    //     srf_orig,
    //     srf.poles()
    // );

}

TEST(tests_bssurf, reverseU)
{
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    gbs::points_vector<double,3> poles = 
    {                                       // ----U----
        {0,0,0},{1,0,0},                    // |
        {0,1,0},{1,1,1},                    // V
        {0,2,0},{1,2,0},                    // |
    };

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q);

    srf.reverseU();
    srf.reverseV();


    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[0]-srf(1,1)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[1]-srf(0,1)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[4]-srf(1,0)), 0.
    );
    ASSERT_DOUBLE_EQ
    (
        gbs::norm(poles[5]-srf(0,0)), 0.
    );

    // gbs::plot(
    //     srf,
    //     srf.poles()
    // );

}

TEST(tests_bssurf,interp1)
{
    //Pij avec j inner loop
    // ---U--
    // |
    // V
    // |
    const std::vector<std::array<double,3> > points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };
    std::vector<double> ku = {0.,0.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,0.5,1.,1.,1.};
    size_t p = 1;
    size_t q = 2;
    std::vector<double> u = {0.,1.};
    std::vector<double> v = {0.,0.33,0.66,1.};
    auto poles = gbs::build_poles(points,ku,kv,u,v,p,q);

    gbs::BSSurface<double,3> srf(poles,ku,kv,p,q) ;

    for (int j = 0; j < v.size(); j++)
    {
        for (int i = 0; i < u.size(); i++)
        {
            auto pt = srf.value(u[i],v[j]);
            ASSERT_LT(gbs::norm(points[i+u.size()*j] - pt ), tol);
        }
    }

}



TEST(tests_bssurf,interp_general_cstr)
{
    // ---U--
    // |
    // V
    // |
    using T = double;
    const size_t dim =3;

    auto cstr = [](T u,T v, const gbs::point<T,dim> &x, size_t du, size_t dv)
    {
        return std::make_tuple(u, v, x, du, dv);
    };

    gbs::points_vector<T,dim> pts{
        {0,0,0},{1,0,0},{2,0,1},{3,0,0},
        {0,1,0},{1,1,1},{2,1,2},{3,1,0},
    };
    gbs::points_vector<T,dim> tgs{
        {0,1,0},{0,1,0},{0,1,0},{0,1,0},
        {0,1,0},{0,1,0},{0,1,0},{0,1,0},
    };
    std::vector<T> u{0.1,0.33,0.66,0.9};
    std::vector<T> v{0.1,0.9};

    const std::vector<gbs::bss_constrain<double,3>> Q =
    {
        cstr(u[0],v[0],pts[0],0,0), cstr(u[1],v[0],pts[1],0,0),cstr(u[2],v[0],pts[2],0,0), cstr(u[3],v[0],pts[3],0,0),
        cstr(u[0],v[1],pts[4],0,0), cstr(u[1],v[1],pts[5],0,0),cstr(u[2],v[1],pts[6],0,0), cstr(u[3],v[1],pts[7],0,0),
        cstr(u[0],v[0],tgs[0],0,1), cstr(u[1],v[0],tgs[1],0,1),cstr(u[2],v[0],tgs[2],0,1), cstr(u[3],v[0],tgs[3],0,1),
        cstr(u[0],v[1],tgs[4],0,1), cstr(u[1],v[1],tgs[5],0,1),cstr(u[2],v[1],tgs[6],0,1), cstr(u[3],v[1],tgs[7],0,1),
    };

    auto [u_min, u_max, v_min, v_max] = gbs::get_constrains_bounds(Q);
    ASSERT_NEAR(u_min,u.front(),1e-6);
    ASSERT_NEAR(u_max,u.back(),1e-6);
    ASSERT_NEAR(v_min,v.front(),1e-6);
    ASSERT_NEAR(v_max,v.back(),1e-6);

    size_t p = 3;
    size_t q = 2;
    size_t n_polesv = 4;
    auto srf = gbs::interpolate(Q,n_polesv,p,q) ;

    for (int j = 0; j < v.size(); j++)
    {
        for (int i = 0; i < u.size(); i++)
        {
            auto pt = srf.value(u[i],v[j]);
            auto tg = srf.value(u[i],v[j],0,1);
            ASSERT_LT(gbs::norm(pts[i+u.size()*j] - pt ), tol);
            ASSERT_LT(gbs::norm(tgs[i+u.size()*j] - tg ), tol);
        }
    }

    gbs::plot(srf,pts);

}

TEST(tests_bssurf,extract_row_col)
{
    const std::vector<std::array<double,3> > points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };

    for (auto i = 0; i < 4; i++)
    {
        for (auto j = 0; j < 2; j++)
        {
            auto point_U = gbs::extract_U(i, points, 2);
            ASSERT_LT(gbs::norm(points[j + 2 * i] - point_U[j]), tol);
        }
    }
}
TEST(cpp_algo, reduce)
{
    const std::vector<double> v(10'000'007, 0.5);

    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double result = std::accumulate(v.cbegin(), v.cend(), 0.0);
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << std::fixed << "std::accumulate result " << result
                  << " took " << ms.count() << " ms/n";
    }

    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double result = std::reduce(std::execution::par, v.cbegin(), v.cend());
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "std::reduce result "
                  << result << " took " << ms.count() << " ms/n";
    }
}

TEST(cpp_algo, par_vs_seq)
{
    {
        std::vector<double> v1(10);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::fill(
            std::execution::seq,
            v1.begin(), v1.end(), 1.);
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "fill took " << ms.count() << " ms/n in seq for a vec of size:" << v1.size() << std::endl;
    }
    {
        std::vector<double> v1(10);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::fill(
            std::execution::par,
            v1.begin(), v1.end(), 1.);
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "fill took " << ms.count() << " ms/n in // for a vec of size:" << v1.size() << std::endl;
    }
    {
        std::vector<double> v1(10);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::transform(
            std::execution::seq,
            v1.begin(), v1.end(),
            v1.begin(),
            [](const auto v_){return 2*v_;});
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "transform took " << ms.count() << " ms/n in seq for a vec of size:" << v1.size() << std::endl;
    }
    {
        std::vector<double> v1(10);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::transform(
            std::execution::par,
            v1.begin(), v1.end(),
            v1.begin(),
            [](const auto v_){return 2*v_;});
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "transform took " << ms.count() << " ms/n in // for a vec of size:" << v1.size() << std::endl;
    }
    size_t n = 100000000;
    {
        std::vector<double> v1(n);
        std::vector<double> v2(n);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::transform(
            std::execution::seq,
            v1.begin(), v1.end(),
            v2.begin(),
            [](const auto v_){return 2*v_;});
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "transform took " << ms.count() << " ms/n in seq for a vec of size:" << v1.size() << std::endl;
    }
    {
        std::vector<double> v1(n);
        std::vector<double> v2(n);
        const auto t1 = std::chrono::high_resolution_clock::now();
        std::transform(
            std::execution::par,
            v1.begin(), v1.end(),
            v2.begin(),
            [](const auto v_){return 2*v_;});
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms = t2 - t1;
        std::cout << "transform took " << ms.count() << " ms/n in // for a vec of size:" << v1.size() << std::endl;
    }
}

TEST(tests_bssurf, eval_except)
{
    gbs::BSSurface<double,2> srf {
        {
            {0., 0.},
            {1., 0.},
            {0., 1.},
            {1., 1.},
        },
        {0., 0., 1., 1.},
        {0., 0., 1., 1.},
        1,
        1
    };
    ASSERT_THROW(srf.value(2.0,0.5), gbs::OutOfBoundsSurfaceUEval<double>);
    ASSERT_THROW(srf.value(0.5,2.0), gbs::OutOfBoundsSurfaceVEval<double>);
}

TEST(tests_bssurf, approx_constrained)
{
    // ---U--
    // |
    // V
    // |
    using T = double;
    const size_t dim =3;

    auto cstr = [](T u,T v, const gbs::point<T,dim> &x, size_t du, size_t dv)
    {
        return std::make_tuple(u, v, x, du, dv);
    };

    gbs::points_vector<T,dim> pts{
        {0,0,0},{1,0,0},{2,0,1},{3,0,0},
        {0,1,0},{1,1,1},{2,1,2},{3,1,0},
    };
    gbs::points_vector<T,dim> tgs{
        {0,1,0},{0,1,0},{0,1,0},{0,1,0},
        {0,1,0},{0,1,0},{0,1,0},{0,1,0},
    };
    std::vector<T> u{0.1,0.33,0.66,0.9};
    std::vector<T> v{0.1,0.9};

    const std::vector<gbs::bss_constrain<double,3>> Q =
    {
        cstr(u[0],v[0],pts[0],0,0), cstr(u[1],v[0],pts[1],0,0),cstr(u[2],v[0],pts[2],0,0), cstr(u[3],v[0],pts[3],0,0),
        cstr(u[0],v[1],pts[4],0,0), cstr(u[1],v[1],pts[5],0,0),cstr(u[2],v[1],pts[6],0,0), cstr(u[3],v[1],pts[7],0,0),
        cstr(u[0],v[0],tgs[0],0,1), cstr(u[1],v[0],tgs[1],0,1),cstr(u[2],v[0],tgs[2],0,1), cstr(u[3],v[0],tgs[3],0,1),
        cstr(u[0],v[1],tgs[4],0,1), cstr(u[1],v[1],tgs[5],0,1),cstr(u[2],v[1],tgs[6],0,1), cstr(u[3],v[1],tgs[7],0,1),
    };

    size_t p = 3;
    size_t q = 3;

    std::vector<T> ku{0.,0.,0.,0.,1.,1.,1.,1.};
    std::vector<T> kv{0.,0.,0.,0.,1.,1.,1.,1.};

    auto srf = gbs::approx(Q,ku,kv,p,q);
    gbs::plot(srf,pts);
}
