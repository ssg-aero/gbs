#include <gtest/gtest.h>
#include <gbs-cuda/basisfunctions.cuh>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <gbs/maths.h>
#include <gbs/basisfunctions.h>

namespace gbs
{
    struct f_gbs_cuda_basis_function
    {
        Id p_;
        device_vector<Real>::iterator knots_begin_;
        device_vector<Real>::iterator knots_end_;
        f_gbs_cuda_basis_function(
            Id p, 
            const device_vector<Real>::iterator &knots_begin,
            const device_vector<Real>::iterator &knots_end
        ) : p_{p}, knots_begin_{knots_begin}, knots_end_{knots_end} {}
        __device__
        auto operator()(Real u)
        {
            return basis_function(
                p_, 
                u, 
                knots_begin_, 
                knots_end_
            );
        }
    };
}


TEST(gbs_cuda, basis_function)
{
    using namespace gbs;

    const size_t p = 3;
    const size_t n = 100;
    std::vector<Real> knots_h
    {
        0., 0., 0.,
        0.3,
        0.7,
        1., 1., 1.
    };

    auto u_h = make_range<Real>(0., 1., n);
    std::vector<Real> eval_h(n);
    std::transform(
        u_h.begin(), u_h.end(),
        eval_h.begin(),
        [&](Real u){
            return basis_function(
                u, 
                std::begin(knots_h), 
                p, 
                std::end(knots_h)
            );
        }
    );


    device_vector<Real> u_d{ u_h };
    device_vector<Real> knots_d{knots_h};
    device_vector<Real> eval_d(n);

    thrust::transform(
        u_d.begin(), u_d.end(),
        eval_d.begin(),
        f_gbs_cuda_basis_function(p,knots_d.begin(), knots_d.end())
    );

    for(int i{}; i < n ; i++)
    {
        ASSERT_NEAR( eval_d[i] , eval_h[i], 1e-6 );
    }

}