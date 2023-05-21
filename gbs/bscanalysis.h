#pragma once

#include <list>
#include <nlopt.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "extrema.h"
#include "bscurve.h"
#include "bscinterp.h"
#include "bscapprox.h"

// GSL_INTEG_GAUSS15

// GSL_INTEG_GAUSS21

// GSL_INTEG_GAUSS31

// GSL_INTEG_GAUSS41

// GSL_INTEG_GAUSS51

// GSL_INTEG_GAUSS61

namespace{
    constexpr size_t N_gauss_pt = 31;
}
namespace gbs
{
/**
 * @brief Computes the deviation from a set of points
 * 
 * @tparam T 
 * @tparam dim 
 * @param points 
 * @param crv 
 * @return (pos of max, max distance, average distance) 
 */
    template <typename T, size_t dim>
    auto dev_from_points(const std::vector<std::array<T, dim>> &points, const BSCurve<T, dim> &crv)
    {
        auto d_avg = 0., d_max = 0., u_max = 0.;
        auto u0 = crv.knotsFlats().front();
        std::for_each(
            points.begin(),
            points.end(),
            [&](const auto &pnt) {
                auto [res_u,res_d] = extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res_u;
                if (res_d > d_max)
                {
                    d_max = res_d;
                    u_max = res_u;
                }
                d_avg += res_d;
            });
        d_avg /= points.size();
        return std::make_tuple(u_max, d_max, d_avg);
    }
    /**
     * \brief Compute the length of a curve using Gauss-Kronrod quadrature or Gauss-Legendre quadrature.
     * The optimal number of points for the Gauss-Kronrod quadrature depends on the specific function being integrated and the desired level of accuracy. It is not possible to determine a single set of optimal numbers that will work for all integrands.
     * 
     *  However, here are some common choices for the degree (number of points) of the Gauss-Kronrod rule, based on the desired level of accuracy:
     * 
     *    Degree 3: This rule uses 7 points (3 Gauss points and 4 Kronrod points) and provides an accuracy of about 2 decimal places.
     * 
     *    Degree 5: This rule uses 11 points (5 Gauss points and 6 Kronrod points) and provides an accuracy of about 4 decimal places.
     * 
     *    Degree 7: This rule uses 15 points (7 Gauss points and 8 Kronrod points) and provides an accuracy of about 6 decimal places.
     * 
     *    Degree 9: This rule uses 21 points (9 Gauss points and 12 Kronrod points) and provides an accuracy of about 8 decimal places.
     * 
     *    Degree 11: This rule uses 27 points (11 Gauss points and 16 Kronrod points) and provides an accuracy of about 10 decimal places.
     *
     * \tparam T The type of the curve's parameter values.
     * \tparam dim The dimension of the curve.
     * \tparam N The number of Gauss-Kronrod or Gauss-Legendre integration points to use. Defaults to 31.
     * Internally class gauss_kronrod has pre-computed tables of abscissa and weights for 15, 31, 41, 51 and 61 Gauss-Kronrod points at up to 100-decimal digit precision. That means that using for example, gauss_kronrod<double, 31>::integrate incurs absolutely zero set-up overhead from computing the abscissa/weight pairs. When using multiprecision types with less than 100 digits of precision, then there is a small initial one time cost, while the abscissa/weight pairs are constructed from strings.
     * \param crv The curve.
     * \param u1 The starting parameter value for the integration.
     * \param u2 The ending parameter value for the integration.
     * \param d The derivative to compute. 0 for the length, 1 for the tangent vector at the end point, 2 for the curvature at the end point.
     * \param adaptive Whether to use Gauss-Kronrod quadrature (true) or Gauss-Legendre quadrature (false).
     *
     * \returns The length of the curve or the tangent vector or the curvature, depending on the value of the `d` parameter.
     */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto length(const Curve<T,dim> &crv,T u1 , T u2, size_t d = 0, bool adaptive=true) -> T
    {
        using namespace boost::math::quadrature;

        auto f = [&](const auto& u) { return norm(crv.value(u,1)); };
        switch (d)
        {
        case 0:
            return  adaptive ? gauss_kronrod<T, N>::integrate(f, u1, u2) : gauss<T, N>::integrate(f, u1, u2);
            break;
        case 1:
            return f(u2);
            break;
        case 2:
            return (crv.value(u2,1) * crv.value(u2,2)) / f(u2);
            break;
        default:
            {
                throw std::runtime_error("Not implemented.");
                return  0.;
            }
            break;
        }
    }
/**
 * @brief Compute full curve length
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param d  : derivate order max 2 respective to u2
 * @return T 
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto length(const Curve<T,dim> &crv, size_t d = 0, bool adaptive=true) -> T
    {
        auto [u1,u2] = crv.bounds();
        return length<T,dim,N>(crv,u1,u2,d, adaptive);
    }
/**
 * Calculate the absolute curvature of a curve.
 * @tparam T The floating-point type.
 * @tparam dim The dimension of the curve.
 * @tparam N The number of divisions for length calculation.
 * @param crv The curve object.
 * @param u1 The starting parameter value.
 * @param u2 The ending parameter value.
 * @param n The number of divisions to calculate the curvature (default: 30).
 * @param p The degree of the B-spline curve (default: 3).
 * @return A BSCfunction object representing the absolute curvature.
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30, size_t p = 3) -> BSCfunction<T>
    {
        points_vector<T, 1> u;

        u = make_range(point<T, 1>{u1}, point<T, 1>{u2}, n);

        std::vector<T> dm(u.size() - 1);
        std::transform(
            std::execution::par,
            u.begin(), std::next(u.end(), -1), std::next(u.begin(), 1), dm.begin(),
            [&crv](const auto &u1_,const auto &u2_) {
                return length<T, dim, N>(crv, u1_[0], u2_[0]);
            });

        std::vector<T> m(u.size());
        m.front() = 0.;
        std::transform(dm.begin(), dm.end(), m.begin(), std::next(m.begin()),
                       [](T dm_, T sum_) {
                           return sum_ + dm_;
                       });


        return  BSCfunction<T>{ interpolate<T, 1>(u, m, fmin(p, u.size()-1)) };
    }

/**
 * @brief Computes the adaptive discretization of a curve based on its absolute curvature.
 *
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the adaptive discretization will be calculated
 * @param u1 The lower bound of the curve parameter
 * @param u2 The upper bound of the curve parameter
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return std::pair<points_vector<T, 1>, std::vector<T>> A pair of vectors, where the first vector contains the parameter values (u) and the second vector contains the accumulated arc length (m) at each parameter value.
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv_adaptive_discretzation(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30, size_t p = 3, T tolerance = 0.001)
    {
        std::vector<T> u;

        u = make_range(u1, u2, n);

        // Adaptive refinement based on local curvature
        bool inserted = true;
        while (inserted)
        {
            inserted = false;
            auto it_u1 = u.begin();
            auto it_u3 = std::next(it_u1);
            while (it_u3 != u.end())
            {
                T mid_u = 0.5 * ((*it_u1) + (*it_u3));
                T l12 = length<T, dim, N>(crv, (*it_u1), mid_u);
                T l23 = length<T, dim, N>(crv, mid_u, (*it_u3));
                T l13 = length<T, dim, N>(crv, (*it_u1), (*it_u3));
                T curv_diff = std::abs( l12 + l23 - l13) / l13;
                if (curv_diff > tolerance)
                {
                    inserted = true;
                    u.insert(it_u3, mid_u);
                }
                it_u1 = it_u3;
                it_u3 = std::next(it_u1);
            }
        }

        std::vector<T> dm(u.size() - 1);
        std::transform(
            std::execution::par,
            u.begin(), std::next(u.end(), -1), std::next(u.begin(), 1), dm.begin(),
            [&crv](const auto &u1_, const auto &u2_) {
                return length<T, dim, N>(crv, u1_, u2_);
            });

        std::vector<T> m(u.size());
        m.front() = 0.;
        std::transform(dm.begin(), dm.end(), m.begin(), std::next(m.begin()),
                    [](T dm_, T sum_) {
                        return sum_ + dm_;
                    });
        return std::make_pair(u,m);
    }
/**
 * @brief Computes the curvilinear abscissa of a curve using an adaptive refinement strategy.
 * 
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the curvilinear abscissa will be calculated
 * @param u1 The lower bound of the curve parameter
 * @param u2 The upper bound of the curve parameter
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return BSCfunction<T> The curvilinear abscissa as a BSCfunction of curve parametrization
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto abs_curv_adaptive(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30, size_t p = 3, T tolerance = 0.001) -> BSCfunction<T>
    {
        auto [u, m] = abs_curv_adaptive_discretzation<T, dim, N>(crv, u1, u2, n, p, tolerance);
        return interpolate<T>(u, m, std::min(p, u.size() - 1));
    }
/**
 * @brief Computes the curve length function using an adaptive discretization strategy based on the curve's absolute curvature.
 *
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the length function will be calculated
 * @param u1 The lower bound of the curve parameter
 * @param u2 The upper bound of the curve parameter
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return BSCfunction<T> The curve length function as a BSCfunction
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto curve_length_function(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30, size_t p = 3, T tolerance = 0.001) -> BSCfunction<T>
    {
        auto [u, m] = abs_curv_adaptive_discretzation<T, dim, N>(crv, u1, u2, n, p, tolerance);
        return interpolate( m, u, std::min(p, u.size() - 1));
    }
/**
 * @brief Computes the curve length function using an adaptive discretization strategy based on the curve's absolute curvature.
 *
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the length function will be calculated
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return BSCfunction<T> The curve length function as a BSCfunction
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto curve_length_function(const Curve<T, dim> &crv, size_t n = 30, size_t p = 3, T tolerance = 0.001) -> BSCfunction<T>
    {
        auto [u1, u2] = crv.bounds();
        auto [u, m] = abs_curv_adaptive_discretzation<T, dim, N>(crv, u1, u2, n, p, tolerance);
        return interpolate( m, u, std::min(p, u.size() - 1));
    }
/**
 * @brief Computes the curvilinear abscissa of a curve using an adaptive refinement strategy.
 * 
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the curvilinear abscissa will be calculated
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return BSCfunction<T> The curvilinear abscissa as a BSCfunction of curve parametrization
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv_adaptive(const Curve<T, dim> &crv, size_t n = 30, size_t p = 3) -> BSCfunction<T>
    {
        auto [u1, u2] = crv.bounds();
        return abs_curv<T,dim,N>(crv,u1,u2,n,p);
    }

/**
 * @brief Computes the curvilinear abscissa of a curve.
 * 
 * @tparam T Floating-point type for calculations
 * @tparam dim Dimension of the curve
 * @tparam N Number of segments for the length computation
 * @param crv The curve for which the curvilinear abscissa will be calculated
 * @param n Initial number of sampling points
 * @param p Degree for the BSCfunction
 * @param tolerance Tolerance for adaptive refinement based on local curvature
 * @return BSCfunction<T> The curvilinear abscissa as a BSCfunction of curve parametrization
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv(const Curve<T, dim> &crv, size_t n = 30, size_t p = 3) -> BSCfunction<T>
    {
        auto [u1, u2] = crv.bounds();
        return abs_curv<T,dim,N>(crv,u1,u2,n,p);
    }
/**
 * @brief Generate a list of uniformly distributed parameters along the curve based on its arc length
 * 
 * @tparam T Numeric type
 * @tparam dim Dimension of the curve
 * @tparam N Number of Gauss quadrature points for arc length computation in abs_curv function
 * @param crv The curve object
 * @param u1 Starting parameter value
 * @param u2 Ending parameter value
 * @param n Number of parameters to generate
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> A list of uniformly distributed parameters along the curve
 */
    template <typename T, size_t dim, size_t N = 10>
    auto uniform_distrib_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, size_t n_law = 30) -> std::list<T>
    {
        // Initialize the list of parameters and set the first and last elements
        std::list<T> u_lst(n);
        u_lst.front() = u1;
        u_lst.back() = u2;

        // Generate the arc length function
        auto f_u = abs_curv<T, dim, N>(crv, u1, u2, n_law);

        // Calculate the uniform step in arc length
        auto dm = f_u.bounds()[1] / (n - T(1));

        // Generate the remaining parameters based on the arc length function
        std::generate(
            std::next(u_lst.begin(), 1),
            std::prev(u_lst.end()),
            [&, m_ = T(0.)]() mutable {
                m_ += dm;
                return f_u(m_);
            });

        // Verify if the parameters are sorted, otherwise, the arc length function construction failed
        if (!std::is_sorted(u_lst.begin(), u_lst.end()))
            throw std::length_error("Building abs curve fails, please refine n_law");

        return u_lst;
    }
/**
 * @brief Generate a list of uniformly distributed parameters along the curve based on its arc length
 * 
 * @tparam T Numeric type
 * @tparam dim Dimension of the curve
 * @tparam N Number of Gauss quadrature points for arc length computation in abs_curv function
 * @param crv The curve object
 * @param n Number of parameters to generate
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> A list of uniformly distributed parameters along the curve
 */
    template <typename T, size_t dim, size_t N = 10>
    auto uniform_distrib_params(const Curve<T, dim> &crv, size_t n, size_t n_law = 30) -> std::list<T>
    {
        auto [u1, u2] = crv.bounds();
        return uniform_distrib_params(crv,u1,u2,n,n_law);
    }
/**
 * @brief Refines the range of curve parameters recursively based on the deviation.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param crv The curve to discretize
 * @param dev_max The maximum allowed deviation
 * @param it_u1 Iterator pointing to the start of the range in the parameter list
 * @param it_u3 Iterator pointing to the end of the range in the parameter list
 * @param u_lst The list of curve parameters
 */
    template <typename T, size_t dim>
    void refine_params_recursive(const Curve<T, dim> &crv, T dev_max, typename std::list<T>::iterator it_u1, typename std::list<T>::iterator it_u3, std::list<T> &u_lst)
    {
        auto u_ = 0.5 * (*it_u1 + *it_u3);
        auto v1 = crv(u_) - crv(*it_u1);
        auto v2 = crv(*it_u3) - crv(*it_u1);
        auto dev_ = norm(cross(v1, v2)) / (norm(v1) * norm(v2));

        if (dev_ > dev_max)
        {
            auto it_u2 = u_lst.insert(it_u3, u_);
            refine_params_recursive(crv, dev_max, it_u1, it_u2, u_lst);
            refine_params_recursive(crv, dev_max, it_u2, it_u3, u_lst);
        }
    }

/**
 * @brief Generates a list of curve parameters based on deviation using a recursive approach.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param crv The curve to discretize
 * @param u1 Start value of the parameter range
 * @param u2 End value of the parameter range
 * @param n Initial number of points to create
 * @param dev_max The maximum allowed deviation
 * @param n_max_pts Maximum number of points allowed
 * @return A list of curve parameters with the refined deviation
 */
    template <typename T, size_t dim>
    auto deviation_based_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, T dev_max, size_t n_max_pts = 5000) -> std::list<T>
    {
        // Create initial list of parameters
        std::list<T> u_lst;
        for (size_t i{}; i < n; i++)
            u_lst.push_back(u1 + (u2 - u1) * i / (n - 1.0));

        // Refine the list recursively if the maximum number of points is not reached
        if (u_lst.size() < n_max_pts)
        {
            auto it_u1 = u_lst.begin();
            auto it_u3 = std::next(it_u1);
            while (it_u3 != u_lst.end())
            {
                refine_params_recursive(crv, dev_max, it_u1, it_u3, u_lst);

                it_u1 = it_u3;
                it_u3 = std::next(it_u1);
            }
        }
        return u_lst;
    }

/**
 * @brief 
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @param dev_max 
 * @param n_max_pts 
 * @return std::list<T> 
 */
    template <typename T, size_t dim>
    auto deviation_based_params(const Curve<T, dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000) -> std::list<T>
    {

        auto [u1, u2] = crv.bounds();
        return deviation_based_params(crv, u1, u2, n, dev_max, n_max_pts );
    }
/**
 * @brief Create a vector of points at curve's parameters positions
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam _Container 
 * @param crv Curve
 * @param u_lst Poisitions on curve
 * @return gbs::points_vector<T,dim> 
 */
    template <typename T, size_t dim, typename _Container>
    auto make_points(const Curve<T,dim> &crv,const _Container &u_lst, size_t d = 0) -> gbs::points_vector<T,dim>
    {
        points_vector<T,dim> points(u_lst.size());

        std::transform(
            std::execution::par,
            u_lst.begin(),u_lst.end(),points.begin(),
            [d,&crv](T u_){return crv(u_,d);}
        );

        return points;
    }
/**
 * @brief Uniformly spaced discretization
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @return PointArray<T,dim> 
 */
    template <typename T, size_t dim>
    auto discretize(const Curve<T,dim> &crv, size_t n) -> gbs::points_vector<T,dim>
    {
        return make_points(crv,uniform_distrib_params<T,dim>(crv,n));
    }
    /**
     * @brief Build curve's discretization with a nim/max number of points and refined to respect the maximal deviation
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv curve
     * @param n minimum number of points
     * @param dev_max maximal deviation
     * @param n_max_pts maximal number of points
     * @return gbs::points_vector<T,dim> 
     */
    template <typename T, size_t dim>
    auto discretize(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000) -> gbs::points_vector<T,dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, n,dev_max,n_max_pts);
        // build points
        return make_points(crv,u_lst);

    }
/**
 * @brief Build curve's discretization with a nim/max number of points and refined to respect the maximal deviation, both points and corresponding paramets are returned
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @param dev_max 
 * @param n_max_pts 
 * @return (gbs::points_vector<T,dim>, std::vector<T>) 
 */
    template <typename T, size_t dim>
    auto discretize_with_params(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000)
    {
        auto u_lst = deviation_based_params<T, dim>(crv, n,dev_max,n_max_pts);
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        // build points
        return std::make_tuple( make_points(crv,u_lst) , u );
    }
/**
 * @brief Compute normalized tangential direction of the curve
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return point<T,dim> 
 */
    template <typename T, size_t dim>
    auto tangential_direction(const Curve<T,dim> &crv,T u) -> point<T,dim>
    {
        auto tg = crv(u, 1);
        return tg / norm(tg);
    }
/**
 * @brief Compute curve direction information at given parameter
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return ax1<T,dim> 
 */
    template <typename T,size_t dim>
    auto tangential_line(const Curve<T,dim> &crv,T u) -> ax1<T,dim>
    {
        return {crv(u),tangential_direction(crv,u)};
    }
/**
 * @brief Compute normalized normal direction of the curve
 * 
 * @tparam T 
 * @param crv 
 * @param u 
 * @return point<T,2> 
 */
    template <typename T>
    auto normal_direction(const Curve<T,2> &crv,T u) -> point<T,2>
    {
        auto tg = crv(u, 1);
        point<T, 2> n;
        n[0] = -tg[1];
        n[1] =  tg[0];
        return n / norm(n);
    }
/**
 * @brief Compute normalized normal direction of the curve using curvature, thus is the later is null result is +/-infinity
 * 
 * @tparam T 
 * @param crv 
 * @param u 
 * @return point<T,3> 
 */
    template <typename T>
    auto normal_direction(const Curve<T,3> &crv,T u) -> point<T,3>
    {
        auto tg = crv(u, 1);
        auto cu = crv(u, 2);
        auto n_pln = tg^cu; // normal osculating plane
        n_pln = n_pln / norm(n_pln);
        return tg^n_pln / norm(tg);
    }
/**
 * @brief Compute curve's normal direction information at given parameter
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return ax1<T,dim> 
 */
    template <typename T,size_t dim>
    auto normal_line(const Curve<T,dim> &crv,T u) -> ax1<T,dim>
    {
        return {crv(u),normal_direction(crv,u)};
    }
/**
 * @brief find the closest curve to point
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt : the point
 * @param curve_begin 
 * @param curve_end 
 * @param tol 
 * @return auto : curve's iterator
 */
    template <typename T, size_t dim>
    auto closest_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = 1e-6)
    {
        return std::min_element(
            std::execution::par,
            curve_begin, curve_end,
            [&pt, tol](const auto &c1, const auto &c2)
            {
                return extrema_curve_point(c1, pt, c1.midPoint(), tol)[1] < extrema_curve_point(c2, pt, c2.midPoint(), tol)[1];
            });
    }
/**
 * @brief  find the closest curve's pointer to point
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt : the point
 * @param curve_begin 
 * @param curve_end 
 * @param tol 
 * @return auto : curve pointer's iterator
 */
    template <typename T, size_t dim>
    auto closest_p_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = 1e-6)
    {
        return std::min_element(
            std::execution::par,
            curve_begin, curve_end,
            [&pt, tol](const auto &c1, const auto &c2)
            {
                return extrema_curve_point(*c1, pt, c1->midPoint(), tol)[1] < extrema_curve_point(*c2, pt, c2->midPoint(), tol)[1];
            });
    }

    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T u1, T u2, T tol_x, nlopt::algorithm solver=nlopt::LN_COBYLA)
    {
        auto f = [&crv](const std::vector<double> &x)
        {
            return std::vector<double>{1./gbs::sq_norm(crv(x[0],2))};
        };
        auto gf = [&crv](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto f2 = crv(x[0],2);
            auto f3 = crv(x[0],3);
            return std::vector<double>{ -2 * ( f2 * f3 ) / ( ( f2 * f2 ) * ( f2 * f2 ) )};
        };
        std::vector<T> x{ 0.5 * (u1 + u2)};
        std::vector<T> lb{u1};
        std::vector<T> hb{u2};
        auto minf = gbs::solve_D_nlop(
            f,gf,
            x, lb, hb,
            tol_x*tol_x,
            solver
        );

        return std::array<T,2> { x[0], 1 / sqrt(minf)};
    }

    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T tol_x, nlopt::algorithm solver=nlopt::LN_COBYLA)
    {
        auto [ u1, u2 ] = crv.bounds();
        return max_curvature_pos(crv, u1, u2, tol_x, solver);

    }

    template <typename T>
    auto compute_rolling_ball(const gbs::Curve<T, 2> &crv1, const gbs::Curve<T, 2> &crv2, T u1)
    {
        using gbs::operator-;
        using gbs::operator+;
        using gbs::operator*;
        T tol_x = 1.e-6;
        T tol_a = 1.e-3;

        auto bisec = [&crv1,&crv2](T u1_, T u2_)
        {
            return
            std::make_tuple(
                gbs::Line<T, 2>{gbs::normal_line(crv1, u1_)},
                gbs::Line<T, 2>{gbs::normal_line(crv2, u2_)}
            );
        };
        // Checking if curves are // at this point
        auto u20 =  extrema_curve_point(crv2,crv1(u1), u1, tol_x)[0];
        auto [L1, L2] = bisec(u1,u20);
        if(fabs(fabs(L1.getAx()[1]*L2.getAx()[1])-1.)<tol_a)
            return std::make_tuple(
                u20,
                0.5 * (crv1(u1) + crv2(u20)),
                true
            );
        // Normal case
        auto f_ = [u1, tol_x, &crv1, &crv2,&bisec](const std::vector<double> &x) { // nlopt works only with doubles
            T u2 = x[0];
            auto [L1, L2] = bisec(u1,u2);
            auto [u1_, u2_] = extrema_curve_curve<T>(L1, L2);
            T res = fabs(gbs::sq_norm(L1(u1_) - L1(0.)) - gbs::sq_norm(L2(u2_) - L2(0)));
            // std::cout << u2 << " " << u1_ << " " << u2_ << " " << res << std::endl;
            return res;
        };

        auto [u1_2, u2_2] = crv2.bounds();
        std::vector<T> x{u20};
        std::vector<T> lb{u1_2}, hb{u2_2};
        gbs::solve_N_nlop(
            f_,
            x,lb,hb,
            tol_x,
            nlopt::LN_COBYLA // <- better if use gbs::extrema_curve_point(crv2,crv1(u1), u1, tol_x).u to init
        );

        auto u2 = x[0];
        auto [L1_, L2_] = bisec(u1,u2);
        auto [u1__, u2__] = extrema_curve_curve(L1_, L2_);

        return std::make_tuple(
            x[0],
            T(0.5) * (L1_(u1__) + L2_(u2__)),
            false
        );
    }

} // namespace gbs