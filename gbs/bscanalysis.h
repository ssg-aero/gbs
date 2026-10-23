#pragma once

#include <concepts>
#include <iterator>
#include <list>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "extrema.h"
#include <gbs/intersect/curve_pieces.h> // directly: extrema.h may be in progress (extrema.h -> curves -> curveoffset.h -> here)
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
                auto [res_u,res_d] = extrema_curve_point(crv, pnt, u0, extrema_tol<T>);
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
     * @brief Hodograph: the k-th derivative of a B-spline curve as a B-spline curve
     * (NURBS Book A3.3 CurveDerivCpts, Eq. 3.8).
     *
     * The result has degree p-k, the knot vector of crv without its k first and
     * k last knots, and poles given by the recurrence
     * P^(j)_i = (p-j+1) (P^(j-1)_{i+1} - P^(j-1)_i) / (u_{i+p+1} - u_{i+j}).
     * A zero denominator (knot of multiplicity > p-j+1) only multiplies a basis
     * function that vanishes identically, so the pole is set to 0 (0/0 := 0).
     * For k > p the derivative vanishes: a degree-0 zero curve is returned.
     *
     * Non-rational curves only: the derivative of a NURBS is not a NURBS with the
     * same weights.
     *
     * @param crv The curve
     * @param k   Derivative order (default 1)
     * @return BSCurve<T, dim> such that derivative_curve(crv, k).value(u) == crv.value(u, k)
     */
    template <typename T, size_t dim>
    auto derivative_curve(const BSCurve<T, dim> &crv, size_t k = 1) -> BSCurve<T, dim>
    {
        const auto &U = crv.knotsFlats();
        const size_t p = crv.degree();
        const size_t kk = std::min(k, p); // orders beyond p: same layout, zero poles
        points_vector<T, dim> PK{crv.poles()};
        for (size_t j = 1; j <= kk; j++)
        {
            const T tmp = static_cast<T>(p - j + 1);
            for (size_t i = 0; i + 1 < PK.size(); i++)
            {
                const T den = U[i + p + 1] - U[i + j];
                for (size_t c = 0; c < dim; c++)
                    PK[i][c] = den == T(0) ? T(0) : tmp * (PK[i + 1][c] - PK[i][c]) / den;
            }
            PK.pop_back();
        }
        if (k > p)
            std::fill(PK.begin(), PK.end(), std::array<T, dim>{});
        std::vector<T> UK(std::next(U.begin(), kk), std::prev(U.end(), kk));
        return BSCurve<T, dim>(PK, UK, p - kk);
    }
    /// Rational curves are out of scope (the derivative of a NURBS is not a NURBS
    /// with the same weights).
    template <typename T, size_t dim>
    auto derivative_curve(const BSCurveRational<T, dim> &crv, size_t k = 1) = delete;
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
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto abs_curv(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30, size_t p = 3) -> BSCfunction<T>
    {
        points_vector<T, 1> u;

        u = make_range(point<T, 1>{u1}, point<T, 1>{u2}, n);

        std::vector<T> dm(u.size() - 1);
        std::transform(
            GBS_PAR_EXEC
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
    template <typename T, size_t dim, size_t N = N_gauss_pt>
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
            GBS_PAR_EXEC
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
    auto abs_curv_adaptive(const Curve<T, dim> &crv, size_t n = 30, size_t p = 3, T tolerance = 0.001) -> BSCfunction<T>
    {
        auto [u1, u2] = crv.bounds();
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
 * @brief Curve parameters at prescribed normalized arc lengths
 *
 * Each value s of `s_lst` is a fraction of the arc length between u1 and u2:
 * the returned parameter u satisfies length(u1, u) = s * length(u1, u2).
 * s = 0 and s = 1 give exactly u1 and u2.
 *
 * @tparam T Numeric type
 * @tparam dim Dimension of the curve
 * @tparam N Number of Gauss quadrature points for arc length computation in abs_curv function
 * @param crv The curve object
 * @param u1 Starting parameter value
 * @param u2 Ending parameter value
 * @param s_lst Normalized arc lengths, in [0, 1] and non decreasing
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> The curve parameters, one per value of s_lst
 */
    template <typename T, size_t dim, size_t N = 10>
    auto arc_length_distrib_params(const Curve<T, dim> &crv, T u1, T u2, const std::vector<T> &s_lst, size_t n_law = 30) -> std::list<T>
    {
        if (!std::is_sorted(s_lst.begin(), s_lst.end()) ||
            (!s_lst.empty() && (s_lst.front() < T(0) || s_lst.back() > T(1))))
            throw std::invalid_argument("Normalized arc lengths shall be sorted and within [0, 1]");

        // Arc length function: arc length -> parameter
        auto f_u = abs_curv<T, dim, N>(crv, u1, u2, n_law);
        auto l = f_u.bounds()[1];

        std::list<T> u_lst;
        std::transform(
            s_lst.begin(), s_lst.end(), std::back_inserter(u_lst),
            [&](T s_) {
                if (s_ <= T(0)) return u1;
                if (s_ >= T(1)) return u2;
                return f_u(s_ * l);
            });

        // Verify if the parameters are sorted, otherwise, the arc length function construction failed
        if (!std::is_sorted(u_lst.begin(), u_lst.end()))
            throw std::length_error("Building abs curve fails, please refine n_law");

        return u_lst;
    }
/**
 * @brief Curve parameters at prescribed normalized arc lengths, between the curve bounds
 *
 * @param crv The curve object
 * @param s_lst Normalized arc lengths, in [0, 1] and non decreasing
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> The curve parameters, one per value of s_lst
 */
    template <typename T, size_t dim, size_t N = 10>
    auto arc_length_distrib_params(const Curve<T, dim> &crv, const std::vector<T> &s_lst, size_t n_law = 30) -> std::list<T>
    {
        auto [u1, u2] = crv.bounds();
        return arc_length_distrib_params<T, dim, N>(crv, u1, u2, s_lst, n_law);
    }
/**
 * @brief n curve parameters distributed in arc length following a law
 *
 * The law maps the node index fraction ξ = i / (n - 1) in [0, 1] to the
 * normalized arc length s(ξ), with s(0) = 0, s(1) = 1 and s non decreasing
 * (e.g. s(ξ) = ξ: uniform distribution, s(ξ) = ξ²: nodes clustered at u1).
 *
 * @tparam T Numeric type
 * @tparam dim Dimension of the curve
 * @tparam N Number of Gauss quadrature points for arc length computation in abs_curv function
 * @param crv The curve object
 * @param u1 Starting parameter value
 * @param u2 Ending parameter value
 * @param n Number of parameters to generate (at least 2)
 * @param law Callable T(T), normalized arc length s as a function of ξ
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> The n curve parameters
 */
    template <typename T, size_t dim, size_t N = 10, typename Law>
        requires std::invocable<Law, T>
    auto arc_length_distrib_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, Law &&law, size_t n_law = 30) -> std::list<T>
    {
        if (n < 2)
            throw std::invalid_argument("At least 2 parameters are required");
        std::vector<T> s_lst(n);
        for (size_t i{}; i < n; i++)
            s_lst[i] = static_cast<T>(law(T(i) / T(n - 1)));
        return arc_length_distrib_params<T, dim, N>(crv, u1, u2, s_lst, n_law);
    }
/**
 * @brief n curve parameters distributed in arc length following a law, between the curve bounds
 *
 * @param crv The curve object
 * @param n Number of parameters to generate (at least 2)
 * @param law Callable T(T), normalized arc length s as a function of ξ = i / (n - 1)
 * @param n_law Number of points to sample for the construction of the arc length function
 * @return std::list<T> The n curve parameters
 */
    template <typename T, size_t dim, size_t N = 10, typename Law>
        requires std::invocable<Law, T>
    auto arc_length_distrib_params(const Curve<T, dim> &crv, size_t n, Law &&law, size_t n_law = 30) -> std::list<T>
    {
        auto [u1, u2] = crv.bounds();
        return arc_length_distrib_params<T, dim, N>(crv, u1, u2, n, std::forward<Law>(law), n_law);
    }
/**
 * @brief Refines the range of curve parameters recursively based on the deviation.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param crv The curve to discretize
 * @param dev_max The maximum allowed deviation
 * @param it_u1 Iterator pointing to the start of the range in the parameter list
 * @param it_u2 Iterator pointing to the end of the range in the parameter list
 * @param u_lst The list of curve parameters
 */
    template <typename T, size_t dim>
    void refine_params_recursive(const Curve<T, dim> &crv, T dev_max, typename std::list<T>::iterator it_u1, typename std::list<T>::iterator it_u2, std::list<T> &u_lst)
    {
        auto u1 = *it_u1;
        auto u2 = *it_u2;

        if (std::abs(u2 - u1) > knot_eps<T>)
        {
            auto um = static_cast<T>(0.5) * (u1 + u2);
            auto v1 = crv(um) - crv(u1);
            auto v2 = crv(u2) - crv(u1);
            auto dev = norm(cross(v1, v2)) / (norm(v1) * norm(v2));

            if (dev > dev_max)
            {
                auto it_um = u_lst.insert(it_u2, um);
                refine_params_recursive(crv, dev_max, it_u1, it_um, u_lst);
                refine_params_recursive(crv, dev_max, it_um, it_u2, u_lst);
            }
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
    auto deviation_based_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points) -> std::list<T>
    {

        // Create initial list of parameters
        std::list<T> u_lst;
        for (size_t i{}; i < n; i++)
            u_lst.push_back(u1 + (u2 - u1) * i / (n - 1.0));

        // Try to cast to a B-spline curve
        auto *bsc = dynamic_cast<const BSCurve<T, dim>*>(&crv);
        auto *bscr= dynamic_cast<const BSCurveRational<T, dim>*>(&crv);
        // Add potentially hard points to the list
        std::vector<T> knots;

        if(bsc  != nullptr) knots = bsc->knots();
        if(bscr != nullptr) knots = bscr->knots();

        for (auto u : knots) // Non flat knots
        {
            if( u > u1 && u < u2)
            {
                u_lst.push_back(u);
            }
        }


        // Sort the list
        u_lst.sort();

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
    auto deviation_based_params(const Curve<T, dim> &crv, size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points) -> std::list<T>
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
 * @return points_vector<T,dim> 
 */
    template <typename T, size_t dim, typename _Container>
    auto make_points(const Curve<T,dim> &crv,const _Container &u_lst, size_t d = 0) -> points_vector<T,dim>
    {
        points_vector<T,dim> points(u_lst.size());

        std::transform(
            GBS_PAR_EXEC
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
    auto discretize(const Curve<T,dim> &crv, size_t n) -> points_vector<T,dim>
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
     * @return points_vector<T,dim> 
     */
    template <typename T, size_t dim>
    auto discretize(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points) -> points_vector<T,dim>
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
 * @return (points_vector<T,dim>, std::vector<T>) 
 */
    template <typename T, size_t dim>
    auto discretize_with_params(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points)
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
        auto n_pln = cross(tg, cu); // normal osculating plane
        n_pln = n_pln / norm(n_pln);
        return cross(tg, n_pln) / norm(tg);
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
    namespace analysis_detail
    {
        /// Samples of [u1, u2]: 8 per knot span of a B-spline (its pieces), so that a feature on a
        /// short span (a fillet among many spans) is not stepped over, at least 64 in all (64 on
        /// another curve), then u2.
        template <typename T, size_t dim>
        auto span_samples(const Curve<T, dim> &crv, T u1, T u2) -> std::vector<T>
        {
            std::vector<T> u;
            const auto pieces = CurvePieces<T, dim>{crv}.pieces(u1, u2);
            const std::size_t per_piece = std::max<std::size_t>(8, (64 + pieces.size() - 1) / pieces.size());
            for (const auto &pc : pieces)
                for (std::size_t k{}; k < per_piece; ++k)
                    u.push_back(pc.a + (pc.b - pc.a) * T(k) / T(per_piece));
            u.push_back(u2);
            return u;
        }

        /// Iterator on the curve closest to pt in [first, last), `curve(*it)` giving the curve
        /// (global projection on each curve, in parallel; `last` for an empty range).
        template <typename T, size_t dim, typename It, typename F>
        auto closest(const point<T, dim> &pt, It first, It last, T tol, const F &curve) -> It
        {
            std::vector<T> d(static_cast<std::size_t>(std::ranges::distance(first, last)));
            transform_rethrow(first, last, d.begin(), [&](const auto &c) { return extrema_curve_point(curve(c), pt, tol)[1]; });
            return std::ranges::next(first, std::ranges::distance(d.begin(), std::ranges::min_element(d)));
        }
    }
/**
 * @brief Finds the curve closest to a point (global projection on each curve, in parallel).
 *
 * @param pt : the point
 * @param curve_begin, curve_end : range of curves
 * @param tol : tolerance on the parameters of the projections
 * @return iterator on the closest curve (curve_end for an empty range)
 */
    template <typename T, size_t dim>
    auto closest_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = extrema_tol<T>)
    {
        return analysis_detail::closest(pt, curve_begin, curve_end, tol, [](const Curve<T, dim> &c) -> const Curve<T, dim> & { return c; });
    }
/**
 * @brief Finds the curve closest to a point, in a range of pointers to curves.
 *
 * @param pt : the point
 * @param curve_begin, curve_end : range of pointers to curves
 * @param tol : tolerance on the parameters of the projections
 * @return iterator on the pointer to the closest curve (curve_end for an empty range)
 */
    template <typename T, size_t dim>
    auto closest_p_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = extrema_tol<T>)
    {
        return analysis_detail::closest(pt, curve_begin, curve_end, tol, [](const auto &c) -> const Curve<T, dim> & { return *c; });
    }

    /**
     * @brief Curvature of a curve: |C' x C''| / |C'|^3 (0 where C' vanishes). The norm of the
     * cross product is exact in 2D (a determinant) and 3D; in other dimensions it comes from
     * Lagrange's identity sqrt(|C'|^2 |C''|^2 - (C'.C'')^2), which cancels where the curvature
     * is small (C' and C'' nearly parallel: about 1e-8 relative at best).
     */
    template <typename T, size_t dim>
    auto curvature(const Curve<T, dim> &crv, T u) -> T
    {
        const auto d1 = crv(u, 1), d2 = crv(u, 2);
        const T s2 = d1 * d1;
        if (!(s2 > T(0)))
            return T(0);
        T c;
        if constexpr (dim == 2)
            c = std::abs(d1[0] * d2[1] - d1[1] * d2[0]);
        else if constexpr (dim == 3)
            c = norm(cross(d1, d2));
        else
            c = std::sqrt(std::max(T(0), s2 * (d2 * d2) - (d1 * d2) * (d1 * d2)));
        return c / (s2 * std::sqrt(s2));
    }

    /**
     * @brief Position of the maximum of the curvature on [u1, u2]: the largest of samples taken
     * on every knot span of a B-spline (8 per span, at least 64 in all; 64 on another curve),
     * then Brent (golden section and parabolic interpolation) between its neighbours.
     *
     * @param crv    : the curve (second derivatives required)
     * @param u1, u2 : the range, u1 < u2 (std::invalid_argument otherwise)
     * @param tol_x  : tolerance on the parameter (limited by the precision of Brent, about the
     *                 square root of the machine epsilon relative to the range)
     * @return {u, curvature at u}
     */
    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T u1, T u2, T tol_x)
    {
        if (!(u1 < u2) || !std::isfinite(u2 - u1))
            throw std::invalid_argument("max_curvature_pos: the range must be finite and u1 < u2");
        const auto u = analysis_detail::span_samples(crv, u1, u2);
        const auto c = u | std::views::transform([&crv](T u_) { return curvature(crv, u_); }) | std::ranges::to<std::vector>();
        const auto k = static_cast<std::size_t>(std::ranges::distance(c.begin(), std::ranges::max_element(c)));
        const T a = u[k == 0 ? 0 : k - 1], b = u[std::min(k + 1, u.size() - 1)];
        // bits of precision: those of tol_x relative to the range, at most half of the mantissa
        const int bits = std::clamp(int(-std::log2(std::clamp(tol_x / (u2 - u1), std::numeric_limits<T>::min(), T(1)))), 8,
                                    std::numeric_limits<T>::digits / 2);
        std::uintmax_t it = 100;
        const auto [ub, minus_c] = boost::math::tools::brent_find_minima([&crv](T u_) { return -curvature(crv, u_); }, a, b, bits, it);
        // a maximum at a knot, where C'' jumps, may be better sampled than refined
        return -minus_c >= c[k] ? std::array<T, 2>{ub, -minus_c} : std::array<T, 2>{u[k], c[k]};
    }

    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T tol_x)
    {
        auto [ u1, u2 ] = crv.bounds();
        return max_curvature_pos(crv, u1, u2, tol_x);

    }

    /**
     * @brief Rolling ball between two 2D curves: the largest circle tangent to crv1 at u1 that
     * holds no point of crv2, on a side of crv1 where it holds no point of crv1 either.
     *
     * The circle tangent to crv1 at p = crv1(u1), of centre p + r n (n a unit normal of crv1),
     * passes through a point q of crv2 for r = rho(q) = |q - p|² / (2 (q - p).n),
     * and holds q for any larger r: the ball has the radius min rho over crv2. A minimum inside
     * crv2 is a tangency (the radii along both normals are equal); a minimum at an end of crv2 is
     * a contact on that end, where the ball rests on a corner. rho is sampled on every knot span
     * of crv2 (8 per span, at least 64; 64 on another curve), the smallest sample refined by
     * Brent, and the ball certified by the global projection of its centre on crv2: closer than
     * r, the minimum is refined again from there (crv2 passing close to p, between samples).
     * The projection is exact on a B-spline; on another curve it samples too, and a feature
     * narrower than its samples may be stepped over. The side of the projection of crv1(u1) on
     * crv2 is tried first, then the other one.
     *
     * @throws std::invalid_argument if crv2 is unbounded (a Line: pass a trimmed curve)
     * @throws std::runtime_error if there is no such circle: on each side of crv1, crv2 is behind
     * its tangent at u1, or the ball holds points of crv1 (crv1 bends more than the ball, or comes
     * back across it)
     * @return {u2, centre of the ball, true if the normals of crv1 at u1 and of crv2 at the
     * projection of crv1(u1) are parallel (within 1e-3 on the cosine)}; u2 is the contact on crv2,
     * an end of crv2 for a contact on a corner
     */
    template <typename T>
    auto compute_rolling_ball(const Curve<T, 2> &crv1, const Curve<T, 2> &crv2, T u1)
    {
        T tol_x = 1.e-6;
        T tol_a = 1.e-3;
        // the search runs over the range of crv2: an unbounded curve (Line) cannot be searched
        check_bracket_bounds(crv2.bounds());

        const auto p = crv1(u1);
        // the normals at p and at its projection on crv2: their parallelism is reported, the ball
        // is computed the same way (min rho is the common normal then)
        auto u20 = extrema_curve_point(crv2, p, tol_x)[0];
        const Line<T, 2> L1{normal_line(crv1, u1)}, L2{normal_line(crv2, u20)};
        const bool parallel = fabs(fabs(L1.getAx()[1] * L2.getAx()[1]) - 1.) < tol_a;

        const auto [a, b] = crv2.bounds();
        const auto u = analysis_detail::span_samples(crv2, a, b);
        const auto [b1a, b1b] = crv1.bounds();
        const bool crv1_bounded = std::isfinite(b1a) && std::isfinite(b1b);

        // the ball on the side n of crv1: {u2, centre}, empty when crv2 is behind the tangent of
        // crv1 or when the ball holds points of crv1
        auto ball = [&](const point<T, 2> &n) -> std::optional<std::pair<T, point<T, 2>>> {
            // radius of the circle tangent to crv1 at p through crv2(u2), +inf for a point behind p
            auto rho = [&](T u2) {
                const auto d = crv2(u2) - p;
                const T s = d * n;
                return s > T(0) ? (d * d) / (T(2) * s) : std::numeric_limits<T>::infinity();
            };
            // Brent on rho capped: two infinite values among its points would make its parabolic
            // fit NaN; and its iterations are bounded
            constexpr T cap = std::numeric_limits<T>::max() / T(4);
            const auto rhos = u | std::views::transform(rho) | std::ranges::to<std::vector>();
            // the smallest rho near u[k]: Brent between the neighbours of u[k], or the sample
            auto refine = [&](std::size_t k) -> std::pair<T, T> {
                const T lo = u[k == 0 ? 0 : k - 1], hi = u[std::min(k + 1, u.size() - 1)];
                if (!(hi > lo) || !std::isfinite(rhos[k]))
                    return {u[k], rhos[k]};
                std::uintmax_t it = 100;
                const auto [ub, rb] = boost::math::tools::brent_find_minima([&](T v) { return std::min(rho(v), cap); }, lo, hi,
                                                                            std::numeric_limits<T>::digits / 2, it);
                return rb < rhos[k] ? std::pair{ub, rb} : std::pair{u[k], rhos[k]};
            };
            auto [u2, r] = refine(static_cast<std::size_t>(std::ranges::distance(rhos.begin(), std::ranges::min_element(rhos))));
            if (!std::isfinite(r))
                return std::nullopt;
            // certificate: crv2 must keep out of the ball. A minimum of rho between two samples,
            // sharper than the one they point to (crv2 passing close to p), leaves points of crv2
            // inside: the centre projects on crv2 closer than r (subdivision, exact on a B-spline),
            // and the minimum is refined again from the sample next to the projection
            for (int it{}; it < 4; ++it)
            {
                const auto [uc, dc] = extrema_curve_point(crv2, p + r * n, tol_x);
                if (!(dc < r * (T(1) - T(1e-9))))
                    break;
                const auto kc = static_cast<std::size_t>(std::ranges::distance(u.begin(), std::ranges::lower_bound(u, uc)));
                auto [u2c, rc] = refine(std::min(kc, u.size() - 1));
                if (const T rp = rho(uc); rp < rc)
                {
                    u2c = uc;
                    rc = rp;
                }
                if (!(rc < r))
                    break;
                u2 = u2c;
                r = rc;
            }
            const auto centre = p + r * n;
            // the ball must not hold points of crv1 either (crv1 bending more than the ball on that side)
            if (crv1_bounded && extrema_curve_point(crv1, centre, tol_x)[1] < r * (T(1) - T(1e-9)))
                return std::nullopt;
            // the contact, where the radius is normal to crv2: Brent locates the minimum of rho to
            // the square root of the machine epsilon only (the radius is exact), the root of
            // (crv2 - centre).crv2' next to it is solved by TOMS 748; none on a corner (an end of crv2)
            {
                auto g = [&](T v) { return (crv2(v) - centre) * crv2(v, 1); };
                const T h = T(1e-6) * (b - a);
                const T lo = std::max(a, u2 - h), hi = std::min(b, u2 + h);
                const T glo = g(lo), ghi = g(hi);
                if (lo < hi && (glo > T(0)) != (ghi > T(0)) && glo != T(0) && ghi != T(0))
                {
                    std::uintmax_t it = 100;
                    const auto [r_lo, r_hi] = boost::math::tools::toms748_solve(g, lo, hi, glo, ghi, boost::math::tools::eps_tolerance<T>(std::numeric_limits<T>::digits - 3), it);
                    u2 = (r_lo + r_hi) / T(2);
                }
            }
            return std::pair{u2, centre};
        };

        // the side of crv2 first (the projection of p), then the other one
        auto n = L1.getAx()[1];
        if ((crv2(u20) - p) * n < T(0))
            n = T(-1) * n;
        auto found = ball(n);
        if (!found)
            found = ball(T(-1) * n);
        if (!found)
            throw std::runtime_error("compute_rolling_ball: no circle tangent to both curves (on each side of crv1, crv2 is behind its tangent or the ball holds points of crv1)");
        return std::make_tuple(found->first, found->second, parallel);
    }

    /**
     * @brief One piece of a curvature comb: its teeth as one polyline (p0, tip0, p0, p1, tip1, p1, ...) and their
     * envelope (tip0, tip1, ...).
     */
    template <typename T>
    struct CurvatureCombPiece
    {
        points_vector<T, 2> teeth;
        points_vector<T, 2> envelope;
    };

    /**
     * @brief The curvature comb of a chain of planar curves: along each piece, teeth on the normal, of length
     * proportional to the square root of the signed curvature (so that gentle pieces stay visible beside tight
     * ones), the longest `length` long over all the pieces (one scale: the pieces compare across their junctions).
     *
     * `teeth` teeth are shared out by length (8 at least per piece when there are several), then refined by
     * deviation_based_params (`dev_max`): more where the curve turns faster. A chain run counterclockwise has its
     * teeth outwards where it is convex: a tooth inwards is a curvature inversion. A piece of zero length (a point)
     * and the points where a piece has no tangent get no tooth.
     *
     * @param pieces the curves, in their order along the chain
     * @param length the longest tooth
     * @param teeth the teeth before refinement
     * @param dev_max the deviation (sine of the angle between successive chords) above which the teeth are refined
     * @return one CurvatureCombPiece per curve
     */
    template <typename T>
    auto curvature_comb(const std::vector<const Curve<T, 2> *> &pieces, T length, size_t teeth = 200,
                        T dev_max = T(0.05)) -> std::vector<CurvatureCombPiece<T>>
    {
        struct Sample
        {
            point<T, 2> p;
            point<T, 2> normal; // left of the run
            T k;
        };
        std::vector<T> lengths;
        T total{};
        for (const auto *c : pieces)
            total += lengths.emplace_back(gbs::length(*c));
        std::vector<std::vector<Sample>> samples;
        T largest{};
        for (size_t j = 0; j < pieces.size(); ++j)
        {
            const auto *c = pieces[j];
            auto &out = samples.emplace_back();
            if (!(lengths[j] > T(0)))
                continue;
            const size_t n = pieces.size() == 1
                                 ? teeth + 1
                                 : static_cast<size_t>(std::max(8L, std::lround(T(teeth) * lengths[j] / total))) + 1;
            for (T u : deviation_based_params(*c, n, dev_max))
            {
                const auto d1 = c->value(u, 1), d2 = c->value(u, 2);
                const T speed = std::hypot(d1[0], d1[1]);
                if (!(speed > T(1e-12)))
                    continue;
                const T k = (d1[0] * d2[1] - d1[1] * d2[0]) / (speed * speed * speed);
                out.push_back({c->value(u), {-d1[1] / speed, d1[0] / speed}, k});
                largest = std::max(largest, std::abs(k));
            }
        }
        std::vector<CurvatureCombPiece<T>> comb;
        const T scale = largest > T(0) ? length / std::sqrt(largest) : T(0);
        for (const auto &piece : samples)
        {
            auto &out = comb.emplace_back();
            for (const auto &[p, n, k] : piece)
            {
                // outwards where convex: against the left normal for k > 0
                const T t = scale * std::copysign(std::sqrt(std::abs(k)), k);
                const point<T, 2> tip{p[0] - t * n[0], p[1] - t * n[1]};
                out.teeth.insert(out.teeth.end(), {p, tip, p});
                out.envelope.push_back(tip);
            }
        }
        return comb;
    }

} // namespace gbs
