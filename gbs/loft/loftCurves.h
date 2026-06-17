#pragma once
#include "loftGeneric.h"
#include <gbs/gbsconstants.h>
#include <ranges>
#include <initializer_list>

namespace gbs{

    namespace detail {
        // Extract (T, dim) from a loftable curve type. Specialized only for the
        // curve kinds loft accepts, so it doubles as the membership test below.
        template <typename C> struct loft_curve_traits;
        template <typename T, size_t dim> struct loft_curve_traits<BSCurve<T, dim>> {
            using value_type = T; static constexpr size_t dimension = dim;
        };
        template <typename T, size_t dim> struct loft_curve_traits<BSCurveRational<T, dim>> {
            using value_type = T; static constexpr size_t dimension = dim;
        };

        template <typename> inline constexpr bool is_initializer_list_v = false;
        template <typename U> inline constexpr bool is_initializer_list_v<std::initializer_list<U>> = true;

        // A non-rational or rational B-spline curve loft knows how to consume.
        template <typename C>
        concept LoftableCurve = requires { typename loft_curve_traits<std::remove_cvref_t<C>>::value_type; };

        // Any range of such curves (vector, list, …) — initializer_list is handled
        // by its own overloads because a braced-init-list can't deduce a range.
        template <typename R>
        concept LoftableCurveRange =
            !is_initializer_list_v<std::remove_cvref_t<R>> &&
            std::ranges::input_range<R> &&
            LoftableCurve<std::ranges::range_value_t<R>>;
    }

    /**
     * Loft a range of B-spline curves through given v-parameters and v-degree.
     *
     * Handles every curve container (vector, list, …) and both rational and
     * non-rational curves: the curve type drives degree/dimension deduction and the
     * surface type (BSSurface vs BSSurfaceRational) via loft_generic's `if
     * constexpr`. Collapses the former dozen near-identical overloads (issue #40
     * PR3 / #44).
     *
     * @param bs_lst Range of B-spline curves.
     * @param v Parameter values in the 'v' direction.
     * @param q Degree of the lofted surface in the 'v' direction.
     * @return A BSSurface / BSSurfaceRational representing the lofted surface.
     */
    template <typename R>
        requires detail::LoftableCurveRange<R>
    auto loft(const R &bs_lst,
              const std::vector<typename detail::loft_curve_traits<std::ranges::range_value_t<R>>::value_type> &v,
              size_t q)
    {
        using Tr = detail::loft_curve_traits<std::ranges::range_value_t<R>>;
        return loft_generic<typename Tr::value_type, Tr::dimension>(bs_lst.begin(), bs_lst.end(), v, q);
    }

    /**
     * Loft a range of B-spline curves, choosing the v-parametrization automatically
     * and capping the v-degree at q_max. See the (range, v, q) overload above.
     *
     * @param bs_lst Range of B-spline curves.
     * @param q_max Maximal degree of the lofted surface in the 'v' direction. Default 3.
     * @return A BSSurface / BSSurfaceRational representing the lofted surface.
     */
    template <typename R>
        requires detail::LoftableCurveRange<R>
    auto loft(const R &bs_lst, size_t q_max = 3)
    {
        using Tr = detail::loft_curve_traits<std::ranges::range_value_t<R>>;
        return loft_generic<typename Tr::value_type, Tr::dimension>(bs_lst.begin(), bs_lst.end(), q_max);
    }

    /**
     * initializer_list overload of loft (v, q) — a braced-init-list cannot deduce a
     * generic range, so it gets its own entry point that materializes a vector.
     */
    template <detail::LoftableCurve C>
    auto loft(std::initializer_list<C> bs_lst,
              const std::vector<typename detail::loft_curve_traits<C>::value_type> &v, size_t q)
    {
        using Tr = detail::loft_curve_traits<C>;
        return loft_generic<typename Tr::value_type, Tr::dimension>(std::vector<C>{bs_lst}, v, q);
    }

    /**
     * initializer_list overload of loft (q_max) — see the (init_list, v, q) overload.
     */
    template <detail::LoftableCurve C>
    auto loft(std::initializer_list<C> bs_lst, size_t q_max = 3)
    {
        using Tr = detail::loft_curve_traits<C>;
        return loft_generic<typename Tr::value_type, Tr::dimension>(std::vector<C>{bs_lst}, q_max);
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
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, size_t v_degree_max = loft_v_degree_max,T dev = bs_conversion_dev<T>, size_t np = bs_conversion_np, size_t deg_approx = bs_conversion_degree)
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
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, const std::vector<T> &v, size_t q,T dev = bs_conversion_dev<T>, size_t np = bs_conversion_np, size_t deg_approx = bs_conversion_degree)
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