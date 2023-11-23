#pragma once
#include "loftGeneric.h"

namespace gbs{
   /**
     * Overload of the loft function for a vector of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst Vector of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurve<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, v, q);
    }

    /**
     * Overload of the loft function for a list of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::list<BSCurve<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, v, q);
    }

    /**
     * Overload of the loft function for a initializer_list of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::initializer_list<BSCurve<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(std::vector<BSCurve<T, dim>>{bs_lst}, v, q);
    }

    /**
     * Overload of the loft function for a list of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::list<BSCurve<T, dim>> &bs_lst, size_t q_max=3) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

    /**
     * Overload of the loft function for a vector of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurve<T, dim>> &bs_lst, size_t q_max=3) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

    /**
     * Overload of the loft function for a initializer_list of B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSSurface object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::initializer_list<BSCurve<T, dim>> &bs_lst, size_t q_max=3) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(std::vector<BSCurve<T, dim>>{bs_lst}, q_max);
    }

    /**
     * Overload of the loft function for a vector of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst Vector of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurveRational<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, v, q);
    }

    /**
     * Overload of the loft function for a list of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::list<BSCurveRational<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, v, q);
    }

    /**
     * Overload of the loft function for a initializer_list of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::initializer_list<BSCurveRational<T, dim>> &bs_lst, const std::vector<T> &v, size_t q) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(std::vector<BSCurveRational<T, dim>>{bs_lst}, v, q);
    }

    /**
     * Overload of the loft function for a list of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::list<BSCurveRational<T, dim>> &bs_lst, size_t q_max=3) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

    /**
     * Overload of the loft function for a vector of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurveRational<T, dim>> &bs_lst, size_t q_max=3) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

    /**
     * Overload of the loft function for a initializer_list of rational B-spline curves.
     *
     * @tparam T Floating point type.
     * @tparam dim Dimensionality.
     * @param bs_lst List of rational B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSCurveRational object representing the lofted surface.
     */
    template <typename T, size_t dim>
    auto loft(const std::initializer_list<BSCurveRational<T, dim>> &bs_lst, size_t q_max=3) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(std::vector<BSCurveRational<T, dim>>{bs_lst}, q_max);
    }

    /**
     * Constructs a lofted surface from a list of curves.
     * 
     * This function takes a vector of shared pointers to Curves and constructs a lofted 
     * surface (B-spline surface) by converting each Curve to a B-spline curve and then 
     * lofting them together. The conversion of Curves to B-spline curves is done using 
     * the 'to_bs_curve' function.
     * 
     * @param crv_lst A vector of shared pointers to Curves to be lofted.
     * @param v_degree_max The maximum degree for the V-direction of the lofted surface (default: 3).
     * @param dev The maximum deviation for curve approximation (default: 0.01).
     * @param np The number of points to be used in the approximation (default: 100).
     * @param deg_approx The degree of the B-spline curve to be used in the approximation (default: 5).
     * @return A lofted surface created from the input curves.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, size_t v_degree_max = 3,T dev = 0.01, size_t np=100, size_t deg_approx=5)
    {
        std::vector<std::shared_ptr < BSCurve<T, dim> >> bs_lst(crv_lst.size());
        std::transform(
            crv_lst.begin(),
            crv_lst.end(),
            bs_lst.begin(),
            [dev, np, deg_approx](const auto &crv)
            {return to_bs_curve(crv, dev, np, deg_approx);}
        );
        return loft_generic<T,dim>(bs_lst, v_degree_max);
    }

    /**
     * Constructs a lofted surface from a list of curves with specified knot vector and degree in V-direction.
     * 
     * This overloaded function takes a vector of shared pointers to Curves, a knot vector, and a degree 
     * for the V-direction, and constructs a lofted surface. Each Curve is converted to a B-spline curve 
     * using the 'to_bs_curve' function, and then they are lofted together using the specified knot vector 
     * and degree in V-direction.
     * 
     * @param crv_lst A vector of shared pointers to Curves to be lofted.
     * @param v A vector representing the knot vector for the V-direction of the lofted surface.
     * @param q The degree for the V-direction of the lofted surface.
     * @param dev The maximum deviation for curve approximation (default: 0.01).
     * @param np The number of points to be used in the approximation (default: 100).
     * @param deg_approx The degree of the B-spline curve to be used in the approximation (default: 5).
     * @return A lofted surface created from the input curves with specified knot vector and degree.
     */
    template <typename T, size_t dim>
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, const std::vector<T> &v, size_t q,T dev = 0.01, size_t np=100, size_t deg_approx=5)
    {
        std::vector<std::shared_ptr < BSCurve<T, dim> >> bs_lst(crv_lst.size());
        std::transform(
            crv_lst.begin(),
            crv_lst.end(),
            bs_lst.begin(),
            [dev, np, deg_approx](const auto &crv)
            {return to_bs_curve(crv, dev, np, deg_approx);}
        );
        return loft_generic<T,dim>(bs_lst ,v, q);
    }
}