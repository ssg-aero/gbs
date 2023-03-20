#pragma once
#include <thrust/device_vector.h>
#include <numeric>
#define KNOT_EPS 1e-6
namespace gbs
{

    template <typename T>
    using device_vector = thrust::device_vector<T>;
    using Real = float;
    using Id = unsigned int;
    // const Real KNOT_EPS = std::numeric_limits<Real>::min();

    /**
     * @brief Basis function used to compute BSpline
     * 
     * @tparam _InIt 
     * @param u  Parameter on BSpline object
     * @param it Flat knots iterator
     * @param p  BSpline object degree
     * @param _Last Container's end
     * @return T 
     */
    template <typename _InIt>
    __host__ __device__
    Real basis_function(Id p, Real u, const _InIt &first, const _InIt &last)
    {
        auto u_last = *(last - 1);

        auto ui = *first;
        auto ui1 = *(first + 1);
        if (p == 0)
        {
            return ((ui <= u) && (u < ui1)) || (fabs(ui1 - u_last) < KNOT_EPS && fabs(u - u_last) < KNOT_EPS)
                       ? Real(1.)
                       : Real(0.);
        }
        else
        {
            auto uip  = *(first + p);
            auto ui1p = *(first + p + 1);
            auto C1   = (uip - ui);
            auto C2   = (ui1p - ui1);
            if (C1 > KNOT_EPS)
            {
                C1 = (u - ui) / C1;
                C1 *= basis_function( p - 1, u, first, last);
            }
            if (C2 > KNOT_EPS)
            {
                C2 = (ui1p - u) / C2;
                C2 *= basis_function( p - 1, u, first + 1, last);
            }
            return C1 + C2;
        }
    }

    /**
     * @brief Basis function used to compute BSpline's derivatives
     * 
     * @tparam _InIt 
     * @tparam T 
     * @param u   Parameter on BSpline object
     * @param it  Flat knots iterator
     * @param p   BSpline object degree
     * @param d   Derivative order
     * @param _Last 
     * @return T 
     */
    template <typename _InIt>
    __host__ __device__
    Real basis_function(Id p, Real u, Id d, const _InIt &it, const _InIt &last)
    {
        if (d == 0)
        {
            return basis_function(p, u, it, last);
        }
        else if (d > p)
        {
            return 0.;
        }
        else
        {
            auto ui   = *it;
            auto ui1  = *(it + 1);
            auto uip  = *(it + p);
            auto ui1p = *(it + p + 1);
            auto C1   = (uip - ui);
            auto C2   = (ui1p - ui1);
            if (C1 > KNOT_EPS)
            {
                C1 = basis_function( p - 1, u, d - 1, it, last) / C1;
            }
            if (C2 > KNOT_EPS)
            {
                C2 = basis_function(p, u, d - 1, it, last) / C2;
            }
            return p * (C1 - C2);
        }
    }
    
} // namespace gbs
