#include <gtest/gtest.h>
#include <gbslib/bssurf.h>
#include <gbslib/vecop.h>
#include <gbslib/bssinterp.h>
#include <occt-utils/surfacesbuild.h>
#include <occt-utils/export.h>

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

    gbs::BSSurface srf(poles,ku,kv,p,q);

    auto pt1 = srf.value(u,v);
    ASSERT_LT(gbs::norm(pt1 - std::array<double,4>({54/8.,98/8.,68/8.,27/8.})),tol); // NURBS's Book

    auto ptr1 = srf.valueRational(u,v);
    ASSERT_LT(gbs::norm(ptr1-std::array<double,3>({2.,98/27.,68./27.})),tol); // NURBS's Book

    // gbs::NURBSSurface<double,3> srf_nurbs(poles,ku,kv,p,q);
    // auto ptr1 = srf_nurbs.value2(u,v);
    // ASSERT_LT(gbs::norm(ptr1-std::array<double,3>({2.,98/27.,68./27.})),tol);


    // occt_utils::
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

    gbs::BSSurface srf(poles,ku,kv,p,q) ;

    for (int j = 0; j < v.size(); j++)
    {
        for (int i = 0; i < u.size(); i++)
        {
            auto pt = srf.value(u[i],v[j]);
            ASSERT_LT(gbs::norm(points[i+u.size()*j] - pt ), tol);
        }
    }

    std::list<Handle(Geom_Geometry)> lst;
    lst.push_back(occt_utils::BSplineSurface(srf));

    occt_utils::to_iges(lst,"C:/Users/sebastien/workspace2/gbslib/tests/out/srf_interp1.igs");

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