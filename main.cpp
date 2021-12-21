#include <vector>
#include <array>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/maths.h>

int main()
{
    //using T = double;
    //size_t p = 5;
    //auto u = gbs::make_range<T>(0., 1.,10);
    //auto k = gbs::build_simple_mult_flat_knots(u, p);
    //auto n = k.size() - p - 1;
    //gbs::points_vector_3d_d poles(n);
    //auto dx = 1. / (n - 1);
    ////T u_ = 0.3;
    //gbs::BSCurve<T, 3> c1_3d_dp{ poles, k, p };
    ////c1_3d_dp.value(u_);

    //const auto count_max = 10;

    //auto count = count_max;
    //while (count)
    //{
    //    c1_3d_dp.value(double(count) / double(count_max));
    //    count--;
    //}

    using T = double;
    using namespace gbs;

    size_t p = 3;
    std::vector<T> k = { 0., 0., 0., 0., 1., 5. , 6., 8., 8., 8., 8. };
    points_vector<T, 3> poles =
    {
        {0.,1.,0.},
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,3.,0.},
        {3.,3.,0.},
        {4.,2.,0.},
    };
    size_t np = 1000000;
    auto u = make_range<T>(0., 8., np);
    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for (auto u_ : u)
        {
            eval_value_deboor_cox(u_, k, poles, p);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms_ref = t2 - t1;
        std::cout << std::fixed
            << " took " << ms_ref.count() << " ms\n";
    }

    {
        const auto t1 = std::chrono::high_resolution_clock::now();
        for (auto u_ : u)
        {
            eval_value_decasteljau(u_, k, poles, p);
        }
        const auto t2 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double, std::milli> ms_ref = t2 - t1;
        std::cout << std::fixed
            << " took " << ms_ref.count() << " ms\n";
    }

	return 0;
}