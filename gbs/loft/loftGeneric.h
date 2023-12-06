#pragma once
#include <gbs/bsctools.h>
#include "loftBase.h"

namespace gbs{


    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This template function generates a lofted B-spline surface or B-spline curve based on a sequence of input B-spline curves.
     * The lofting process essentially creates a new surface or curve that passes through or approximates the given sequence of curves.
     *
     * @param first An iterator pointing to the beginning of the sequence of B-spline curves.
     * @param last An iterator pointing to the end of the sequence of B-spline curves.
     * @param v A vector of parameter values in the V-direction for the lofting process.
     * @param flat_v A vector of flattened knot values in the V-direction.
     * @param q The degree of the B-spline surface or curve in the V-direction.
     * @return A lofted B-spline surface or curve, depending on the type of input curves.
     *
     * Note: The iterators `first` and `last` should point to a sequence of B-spline curves with compatible dimensions
     * and types. The function relies on the type of the curve pointed to by the iterator for determining the return type.
     */
    template <typename T, size_t dim, typename InputIt>
    auto loft_generic(InputIt first, InputIt last, const std::vector<T> &v, const std::vector<T> &flat_v, size_t q)
    {
        auto curves_info = get_bs_curves_info<T, dim>(first, last);
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);
        using CurveType = typename std::iterator_traits<InputIt>::value_type;
        if constexpr (std::is_same_v<CurveType, BSCurveRational<T, dim>>) {
            return BSCurveRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
        }
    }
    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This function creates a lofted surface by connecting a set of B-spline curves contained within a generic container.
     * It works with any container type (e.g., std::vector, std::list) that holds the B-spline curves. Based on the 
     * information extracted from these curves, it computes the necessary parameters and delegates the construction
     * of the lofted surface to the 'loft' function. The function also handles different curve types (rational or non-rational),
     * creating the appropriate surface type based on the curve type in the container.
     *
     * @tparam T Floating point type used for curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @tparam Container The container type holding the B-spline curves.
     * @param bs_lst Container of B-spline curves.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param flat_v A flattened knot vector in the 'v' direction for the lofted surface.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A lofted surface represented in a suitable format. The type of the returned object depends on whether
     *         the input curves are rational or non-rational. If the curves are rational, a `BSCurveRational` object
     *         is returned; otherwise, a `BSSurface` object is returned.
     */
    template <typename T, size_t dim, typename Container>
    auto loft_generic(const Container &bs_lst, const std::vector<T> &v, const std::vector<T> &flat_v, size_t q)
    {
        return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), v, flat_v, q);
    }

    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This template function generates a lofted B-spline surface or B-spline curve based on a sequence of input B-spline curves.
     * The lofting process essentially creates a new surface or curve that passes through or approximates the given sequence of curves.
     *
     * @param first An iterator pointing to the beginning of the sequence of B-spline curves.
     * @param last An iterator pointing to the end of the sequence of B-spline curves.
     * @param v A vector of parameter values in the V-direction for the lofting process.
     * @param q The degree of the B-spline surface or curve in the V-direction.
     * @return A lofted B-spline surface or curve, depending on the type of input curves.
     *
     * Note: The iterators `first` and `last` should point to a sequence of B-spline curves with compatible dimensions
     * and types. The function relies on the type of the curve pointed to by the iterator for determining the return type.
     */
    template <typename T, size_t dim, typename InputIt>
    auto loft_generic(InputIt first, InputIt last, const std::vector<T> &v, size_t q)
    {
        auto curves_info = get_bs_curves_info<T, dim>(first, last);
        auto [poles, flat_u, flat_v, p] = loft(curves_info, v, q);

        using CurveType = typename std::iterator_traits<InputIt>::value_type;
        if constexpr (std::is_same_v<CurveType, BSCurveRational<T, dim>>) {
            return BSSurfaceRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
        }
    }
    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This function is a generic implementation that creates a lofted surface by connecting a set of B-spline curves.
     * It works with any container type (e.g., std::vector, std::list) that holds the B-spline curves. The function 
     * computes the necessary information from the curves, and based on the type of curves (rational or non-rational),
     * it constructs the appropriate lofted surface.
     *
     * @tparam T Floating point type used for curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @tparam Container The container type holding the B-spline curves.
     * @param bs_lst Container of B-spline curves.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A lofted surface represented in a suitable format.
     *         The type of the returned object depends on whether the input curves are rational or non-rational.
     */
    template <typename T, size_t dim, typename Container>
    auto loft_generic(const Container &bs_lst, const std::vector<T> &v, size_t q)
    {
        return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), v, q);
    }

    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This template function generates a lofted B-spline surface or B-spline curve based on a sequence of input B-spline curves.
     * The lofting process essentially creates a new surface or curve that passes through or approximates the given sequence of curves.
     *
     * @param first An iterator pointing to the beginning of the sequence of B-spline curves.
     * @param last An iterator pointing to the end of the sequence of B-spline curves.
     * @param q_max The maximal degree of the lofted surface in the 'v' direction.
     * @return A lofted B-spline surface or curve, depending on the type of input curves.
     *
     * Note: The iterators `first` and `last` should point to a sequence of B-spline curves with compatible dimensions
     * and types. The function relies on the type of the curve pointed to by the iterator for determining the return type.
     */
    template <typename T, size_t dim, typename InputIt>
    auto loft_generic(InputIt first, InputIt last, size_t q_max)
    {
        auto curves_info = get_bs_curves_info<T, dim>(first, last);
        auto [poles, flat_u, flat_v, p, q] = loft(curves_info, q_max);

        using CurveType = typename std::iterator_traits<InputIt>::value_type;
        if constexpr (std::is_same_v<CurveType, BSCurveRational<T, dim>>) {
            return BSSurfaceRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
        }
    }
    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This function is a generic implementation that creates a lofted surface by connecting a set of B-spline curves.
     * It works with any container type (e.g., std::vector, std::list) that holds the B-spline curves. The function 
     * computes the necessary information from the curves, and based on the type of curves (rational or non-rational),
     * it constructs the appropriate lofted surface.
     *
     * @tparam T Floating point type used for curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @tparam Container The container type holding the B-spline curves.
     * @param bs_lst Container of B-spline curves.
     * @param q_max The maximal degree of the lofted surface in the 'v' direction.
     * @return A lofted surface represented in a suitable format.
     *         The type of the returned object depends on whether the input curves are rational or non-rational.
     */
    template <typename T, size_t dim, typename Container>
    auto loft_generic(const Container &bs_lst, size_t q_max)
    {
        return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), q_max);
    }
}