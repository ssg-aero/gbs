#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
#include <gbs/bscapprox.h>
#include <gbs/bscanalysis.h>
#include <boost/range/combine.hpp>
// #include <boost/foreach.hpp>
namespace gbs
{
    ///
    ///@brief Search in container if at least one curve pointer has a rational definition
    ///
    ///@tparam Ptr_Container 
    ///@param p_crv_lst 
    ///@return auto 
    ///
    template <typename T, size_t  dim>
    auto has_nurbs(const std::list<Curve<T, dim>*> &p_crv_lst)
    {
        return 
        std::find_if(
            p_crv_lst.begin(),
            p_crv_lst.end(),
            [](const auto p_c){
                auto p_nurbs = dynamic_cast<BSCurveRational<T,dim>*>(p_c);
                return  p_nurbs != nullptr;
                }
        )
        != p_crv_lst.end();
    }

    /**
     * @brief Get the BSCurves cpy object
     * 
     * @tparam T 
     * @tparam dim 
     * @param bs_lst 
     * @return std::list<BSCurveGeneral<T, dim,rational>>
    **/
    template <typename T, size_t  dim, bool rational>
    auto get_BSCurves_ptr_cpy(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst)
    {
        using bs_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        std::list<bs_type> bs_lst_cpy(bs_lst.size());
        std::transform(
            std::execution::par,
            bs_lst.begin(),
            bs_lst.end(),
            bs_lst_cpy.begin(),
            [](const auto bs) { return bs_type(*bs); });
        return bs_lst_cpy;
    }

    template <typename T, size_t dim,bool rational,typename Container>
    auto get_v_aligned_poles(size_t n_poles_u,size_t n_poles_v,const Container &bs_lst_cpy)
    {
        // get v aligned poles
        std::vector< points_vector<T,dim + rational> > poles_v_lst{ n_poles_u };
        std::vector< std::vector<T> > poles_v_params{ n_poles_u };
        for(auto i = 0 ; i < n_poles_u; i++)
        {
            points_vector<T,dim + rational> poles_v{ n_poles_v };
            std::transform(
                std::execution::par,
                bs_lst_cpy.begin(),
                bs_lst_cpy.end(),
                poles_v.begin(),
                [&i](const auto &crv)
                {
                    return crv.poles()[i];
                }
            );
            poles_v_params[i] = curve_parametrization(poles_v, KnotsCalcMode::CHORD_LENGTH, true);
            std::swap( poles_v_lst[i] , poles_v );
        }
        // build parametrization
        std::vector<T> v(n_poles_v,0.);
        for(auto j = 0 ; j < n_poles_v; j++)
        {
            for(auto i = 0 ; i < n_poles_u; i++)
            {
                v[j] += poles_v_params[i][j];
            }
            v[j] /=  n_poles_u;
        }
        return std::make_tuple(poles_v_lst,v);
    }

    template <typename ForwardIt>
    auto get_max_degree(ForwardIt first, ForwardIt last) ->size_t
    {
        auto max_deg_curve = *std::max_element(
            first, last,
            [](const auto &bsc1, const auto &bsc2){return bsc1->degree() < bsc2->degree();}
        );
        return max_deg_curve->degree();
    }

    /**
     * Transposes a matrix of poles.
     * 
     * @param poles A 2D vector representing a matrix of poles.
     * @return A 2D vector representing the transposed matrix of poles.
     */
    template <typename T, size_t dim>
    std::vector<std::vector<std::array<T, dim>>> transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles) {
        if (poles.empty()) return {};

        size_t nu = poles.size();
        size_t nv = poles[0].size();

        std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

        for (size_t i = 0; i < nu; ++i) {
            for (size_t j = 0; j < nv; ++j) {
                poles_t[j][i] = poles[i][j];
            }
        }

        return poles_t;
    }

    /**
     * Transposes a block of a matrix of poles.
     * 
     * @param input The input matrix of poles to transpose.
     * @param output The output matrix where the transposed block will be stored.
     * @param startRow The starting row index for the block.
     * @param startCol The starting column index for the block.
     * @param blockSize The size of the block to be transposed.
     */
    template <typename T, size_t dim>
    // Function to transpose a small block of the matrix
    void transpose_block(const std::vector<std::vector<std::array<T, dim>>>& input,
                        std::vector<std::vector<std::array<T, dim>>>& output,
                        size_t startRow, size_t startCol, size_t blockSize) {
        size_t endRow = std::min(startRow + blockSize, input.size());
        size_t endCol = std::min(startCol + blockSize, input[0].size());

        for (size_t i = startRow; i < endRow; ++i) {
            for (size_t j = startCol; j < endCol; ++j) {
                output[j][i] = input[i][j];
            }
        }
    }

    /**
     * Transposes a matrix of poles using block-based transposition for efficiency.
     * 
     * @param poles The matrix of poles to be transposed.
     * @param blockSize The size of each block used in the transposition.
     * @return The transposed matrix of poles.
     */
    template <typename T, size_t dim>
    auto transpose_poles(const std::vector<std::vector<std::array<T, dim>>>& poles, size_t blockSize) {

        size_t nu = poles.size();
        size_t nv = poles[0].size();

        std::vector<std::vector<std::array<T, dim>>> poles_t(nv, std::vector<std::array<T, dim>>(nu));

        for (size_t i = 0; i < nu; i += blockSize) {
            for (size_t j = 0; j < nv; j += blockSize) {
                transpose_block(poles, poles_t, i, j, blockSize);
            }
        }

        return poles_t;
    }

    /**
     * Flattens a 2D vector of poles into a 1D vector.
     * 
     * @param poles_curves The 2D vector of poles to be flattened.
     * @return A 1D vector containing all the poles.
     */
    template <typename T, size_t dim>
    auto flatten_poles(const std::vector<std::vector<std::array<T, dim>>>& poles_curves) {
        std::vector<std::array<T, dim>> flattened;

        // Calculate total size for memory reservation
        size_t totalSize = std::accumulate(poles_curves.begin(), poles_curves.end(), size_t(0),
            [](size_t sum, const auto& row) { return sum + row.size(); });
        flattened.reserve(totalSize);

        // Flatten the vector
        for (const auto& row : poles_curves) {
            std::ranges::copy(row, std::back_inserter(flattened));
        }

        return flattened;
    }

    /**
     * Builds the poles for a loft surface from a set of curve poles, using NURBS mathematics.
     * 
     * @param poles_curves The poles of the curves used to build the loft surface.
     * @param flat_v The flattened knot vector in the v direction.
     * @param v The knot vector in the v direction.
     * @param flat_u The flattened knot vector in the u direction.
     * @param p The degree of the curve in the u direction.
     * @param q The degree of the curve in the v direction.
     * @return A flattened vector of poles representing the loft surface.
     */
    template <typename T, size_t dim>
    auto build_loft_surface_poles(const std::vector<std::vector<std::array<T, dim>>> &poles_curves,
                            const std::vector<T> &flat_v, const std::vector<T> &v, const std::vector<T> &flat_u, size_t p,
                            size_t q)
    {
        auto poles_curves_t = transpose_poles(poles_curves);
        size_t ncr = poles_curves.size();
        size_t nu = poles_curves[0].size();
        size_t nv = poles_curves.size();

        gbs::MatrixX<T> N(ncr, ncr);
        gbs::build_poles_matrix<T, 1>(flat_v, v, q, ncr, N);
        auto N_inv = N.partialPivLu();

        std::vector<std::vector<std::array<T, dim>>> poles_t(nu, std::vector<std::array<T, dim>>(nv));

        for (size_t i = 0; i < nu; ++i)
        {
            for (size_t d = 0; d < dim; ++d)
            {
                gbs::VectorX<T> b(nv);
                std::transform(poles_curves_t[i].begin(), poles_curves_t[i].end(), b.begin(),
                            [d](const auto &arr)
                            { return arr[d]; });

                auto x = N_inv.solve(b);

                for (size_t j = 0; j < nv; ++j)
                {
                    poles_t[i][j][d] = x(j);
                }
            }
        }

        return flatten_poles(transpose_poles(poles_t));
    }

    /**
     * Creates a lofted surface from a series of B-spline curves.
     *
     * The function performs a loft operation, which is a common technique in computer-aided design
     * and computer graphics for creating a surface that smoothly interpolates between a series of curves.
     * This is achieved by first unifying the degrees and knot vectors of the input B-spline curves.
     * The control points for the lofted surface are then computed based on these unified curves.
     *
     * @tparam T The floating point type used for curve coordinates and knots.
     * @tparam dim The dimensionality of the control points.
     * @param curves_info A vector of gbs::BSCurveInfo structures, each representing a B-spline curve
     *                    with control points, a knot vector, and a degree.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param flat_v A flattened knot vector in the 'v' direction for the lofted surface.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A tuple containing the control points of the lofted surface, the knot vector in the 'u' 
     *         direction, and the degree of the surface in the 'u' direction.
     *
     * Example Usage:
     *     auto surface = loft(curves, v_values, flat_v, degree_v);
     */
    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, const std::vector<T> &v, const std::vector<T> &flat_v, size_t q)
    {
        // Unify the degrees and knot vectors of the input curves
        unify_degree(curves_info);
        unify_knots(curves_info);

        // Store the unified degree and the size of the knot vector
        auto p = std::get<2>(curves_info.front());
        auto ncr = std::distance(curves_info.begin(), curves_info.end());
        auto nu = std::get<0>(curves_info.front()).size();
        auto flat_u = std::get<1>(curves_info.front());

        // Extract the control points from each curve
        std::vector<std::vector<std::array<T, dim>>> poles_curves(ncr);
        std::transform(
            curves_info.begin(), curves_info.end(),
            poles_curves.begin(),
            [](const auto &curve_info) { return std::get<0>(curve_info); }
        );

        // Build and return the lofted surface
        return std::make_tuple(build_loft_surface_poles(poles_curves, flat_v, v, flat_u, p, q), flat_u, p);
    }

    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This function creates a lofted surface by connecting a set of B-spline curves specified by the `curves_info`.
     * It first computes a simple multiple flat knot vector in the 'v' direction based on the provided parameters. 
     * The function then uses this knot vector, along with the curve information, to create the lofted surface. 
     * The result includes the control points of the surface, the knot vectors in both the 'u' and 'v' directions, 
     * and the degree in the 'u' direction.
     *
     * @tparam T Floating point type used for curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @param curves_info A vector of BSCurveInfo<T, dim>, each containing essential information of a B-spline curve.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A tuple containing the following elements:
     *         - The control points of the lofted surface.
     *         - The knot vector in the 'u' direction.
     *         - The knot vector in the 'v' direction.
     *         - The degree of the surface in the 'u' direction.
     *
     * Example Usage:
     *     std::vector<BSCurveInfo<float, 3>> curves_info = ...; // Curve information
     *     std::vector<float> v = ...; // Parameter values in the 'v' direction
     *     auto [poles, flat_u, flat_v, p] = loft(curves_info, v, degree_v); // Lofted surface
     */
    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, const std::vector<T> &v, size_t q)
    {
        // Build a simple multiple flat knot vector in the 'v' direction
        auto flat_v = gbs::build_simple_mult_flat_knots(v, q);

        // Delegate to the main loft function to create the lofted surface
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);
        return std::make_tuple(poles, flat_u, flat_v, p);
    }

    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, size_t q_max)
    {
        // TODO: Compute an optimal parametrization
        auto n_curves = curves_info.size();
        auto v = make_range(T{}, T{1}, n_curves);
        auto q = std::min(q_max, n_curves-1);
        // Build a simple multiple flat knot vector in the 'v' direction
        auto flat_v = gbs::build_simple_mult_flat_knots(v, q);

        // Delegate to the main loft function to create the lofted surface
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);
        return std::make_tuple(poles, flat_u, flat_v, p, q);
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
        auto curves_info = get_bs_curves_info<T, dim>(bs_lst.begin(), bs_lst.end());
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);

        if constexpr (std::is_same_v<typename Container::value_type, BSCurveRational<T, dim>>) {
            return gbs::BSCurveRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return gbs::BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
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
        auto curves_info = get_bs_curves_info<T, dim>(bs_lst.begin(), bs_lst.end());
        auto [poles, flat_u, flat_v, p] = loft(curves_info, v, q);

        if constexpr (std::is_same_v<typename Container::value_type, BSCurveRational<T, dim>>) {
            return gbs::BSSurfaceRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return gbs::BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
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
        auto curves_info = get_bs_curves_info<T, dim>(bs_lst.begin(), bs_lst.end());
        auto [poles, flat_u, flat_v, p, q] = loft(curves_info, q_max);

        if constexpr (std::is_same_v<typename Container::value_type, BSCurveRational<T, dim>>) {
            return gbs::BSSurfaceRational<T, dim>(poles, flat_u, flat_v, p, q);
        } else {
            return gbs::BSSurface<T, dim>(poles, flat_u, flat_v, p, q);
        }
    }

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

    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurve<T, dim>> &bs_lst, size_t q_max=3) -> BSSurface<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

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

    template <typename T, size_t dim>
    auto loft(const std::vector<BSCurveRational<T, dim>> &bs_lst, size_t q_max=3) -> BSSurfaceRational<T, dim>
    {
        return loft_generic<T, dim>(bs_lst, q_max);
    }

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

    template <typename T, size_t dim, bool rational,typename Container>
    auto get_tg_from_spine(const std::vector<points_vector<T, dim + rational>> &poles_v_lst, const BSCurve<T, dim> &spine,const Container &bs_lst_cpy)
    {
        std::vector< points_vector<T,dim + rational> > tg_v_lst{ poles_v_lst.size() };
        // compute positions on spine
        std::vector<T> v_spine(bs_lst_cpy.size());
        std::transform(
                std::execution::par,
                bs_lst_cpy.begin(),
                bs_lst_cpy.end(),
                v_spine.begin(),
                [&spine](const Curve<T, dim> &profile_)
                {
                    return extrema_curve_curve<T,dim>(spine,profile_,1e-6)[0];
                }
        );
        // compute spine distances
        std::vector<T> d(v_spine.size() - 1);
        std::transform(
            std::execution::par,
            std::next(v_spine.begin()),
            v_spine.end(),
            v_spine.begin(),
            d.begin(),
            [&](const auto &v2_, const auto &v1_) {
                return norm(spine.value(v2_) - spine.value(v1_));
            });
        // compute spine tg
        std::vector<std::array<T,dim+rational>> D(v_spine.size());
        std::transform(
            std::execution::par,
            v_spine.begin(),
            v_spine.end(),
            D.begin(),
            [&](const auto &v_) {
                return spine.value(v_,1);
            });
        // create weighted tangents
        std::transform(
            poles_v_lst.begin(),
            poles_v_lst.end(),
            tg_v_lst.begin(),
            [&](const auto &pts)
            {
                auto n_poles_v = pts.size();
                points_vector<T,dim + rational> tg_i(n_poles_v);
                tg_i[0] = norm(pts[1]-pts[0])*D[0]/d[0];
                for(auto k = 1 ; k < n_poles_v - 1; k++)
                {
                    tg_i[k] = ( norm(pts[k+1]-pts[k]) +  norm(pts[k]-pts[k-1]) ) * D[k] / (d[k]+d[k-1]);
                }
                tg_i[n_poles_v-1] = norm(pts[n_poles_v-1]-pts[n_poles_v-2])*D[n_poles_v-1]/d[n_poles_v-2];
                return tg_i;
            });
            return tg_v_lst;
    }

    template <typename T, size_t dim,bool rational>
    auto loft(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst, const BSCurve<T, dim> &spine, size_t v_degree_max = 3)
    {
        if(bs_lst.size()<2)
        {
            throw std::length_error("loft needs at least 2 curves.");
        }
        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        std::list<bsc_type> bs_lst_cpy = get_BSCurves_ptr_cpy(bs_lst);

        // static auto dim_poles = dim + rational; 

        unify_curves_degree(bs_lst_cpy);
        unify_curves_knots(bs_lst_cpy);

        auto n_poles_v = bs_lst_cpy.size(); 
        auto n_poles_u = bs_lst_cpy.front().poles().size();
        auto ku_flat = bs_lst_cpy.front().knotsFlats();
        auto p = bs_lst_cpy.front().degree();
        auto [ poles_v_lst, v] = get_v_aligned_poles<T, dim, rational>(n_poles_u, n_poles_v, bs_lst_cpy);

        auto tg_i = get_tg_from_spine<T, dim, rational>(poles_v_lst,spine,bs_lst_cpy);

        // interpolate poles
        size_t q = fmax(fmin(v_degree_max,n_poles_v-1),1); // degree V // diff
        auto kv_flat = build_simple_mult_flat_knots<T>(v.front(),v.back(),n_poles_v * 2,q); //dif
        points_vector<T,dim + rational> poles;
        std::vector<constrType<T,dim + rational,2> > Q(n_poles_v); //diff

        auto pt_and_tg_i = boost::combine(poles_v_lst,tg_i);
        for (auto pt_and_tg : pt_and_tg_i)
        {
            auto [pt, tg] = pt_and_tg;
            std::transform(
                // std::execution::par, // not working with boost::combine
                pt.begin(),
                pt.end(),
                tg.begin(),
                Q.begin(),
                [&](const auto &pt_,const auto &tg_) {
                    return constrType<T, dim + rational, 2>{pt_, tg_};
                });
            auto poles_v = build_poles(Q, kv_flat, v, q);
            poles.insert(poles.end(), poles_v.begin(), poles_v.end());
        }

        // build surf
        using bss_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;
        return bss_type( inverted_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
    }

    template <typename T, size_t dim>
    auto loft(const std::list<BSCurve<T, dim>> &bs_lst, const BSCurve<T, dim> &spine, size_t v_degree_max = 3)
    {
        std::list<gbs::BSCurveGeneral<T, dim, false> *> p_curve_lst(bs_lst.size());
        std::transform(
            bs_lst.begin(),
            bs_lst.end(),
            p_curve_lst.begin(),
            [](auto &bs) 
                { 
                    auto p_bs = static_cast<const gbs::BSCurveGeneral<T, dim, false> *>(&bs);
                    return const_cast<gbs::BSCurveGeneral<T, dim, false> *>(p_bs); 
                }
            ); // pointer to temporary address works within the scope

        return gbs::loft<T, dim, false>(p_curve_lst, spine, v_degree_max);
    }

/*
    template <typename T, size_t dim, typename ForwardIt>
    auto gordon(ForwardIt first_u,ForwardIt last_u,ForwardIt first_v,ForwardIt last_v, T tol = 1e-6) -> BSSurface<T,dim>
    {
        // gather informations
        auto nu = std::distance(first_v, last_v);
        auto nv = std::distance(first_u, last_u);
        auto p  = get_max_degree( first_u, last_u);
        auto q  = get_max_degree( first_v, last_v);
        // compute and check intersections
        std::vector<T> u(nu);
        std::transform(
            std::execution::par,
            first_v, last_v,
            u.begin(),
            [first_u,first_v,tol](const auto &crv_v){
                auto [u_, v_, d]  =extrema_curve_curve(**first_u, *crv_v, tol);
                if( d > tol )
                {
                    throw std::invalid_argument("V curve is not touching first U curve");
                }
                return u_;
            }
        );

        std::vector<T> v(nv);
        std::transform(
            std::execution::par,
            first_u, last_u,
            v.begin(),
            [first_u,first_v,tol](const auto &crv_u){
                auto [u_, v_, d]  =extrema_curve_curve(*crv_u, **first_v, tol);
                if( d > tol )
                {
                    throw std::invalid_argument("U curve is not touching first V curve");
                }
                return v_;
            }
        );

        points_vector<T,dim> Q(nu*nv);
        auto it_u = first_u;
        for(size_t j = 0; j < nv; j++)
        {
            auto v_ = v[j];
            auto it_v = first_v;
            for(size_t i = 0; i < nu; i++)
            {
                auto u_ = u[i];
                auto Q_v = (*it_v)->value(v_);
                auto Q_u = (*it_u)->value(u_);
                if(distance(Q_v, Q_v)>tol)
                {
                    throw std::invalid_argument("U V curve are not intersecting.");
                }
                Q[i + j * nu] = Q_v;
                it_v = std::next(it_v);
            }
            it_u = std::next(it_u);
        }

        // build intermediates surfaces
        auto Lu = loft<T,dim,false>(first_u, last_u,v,q);
        auto Lv = loft<T,dim,false>(first_v, last_v,u,p);
        Lv.invertUV();
        // auto ku = Lu.knotsFlatsU();
        // auto kv = Lu.knotsFlatsV();
        auto ku = build_simple_mult_flat_knots(u.front(), u.back(), nu, p);
        auto kv = build_simple_mult_flat_knots(v.front(), v.back(), nv, q);

        auto poles_T = build_poles(Q, ku, kv, u, v, p, q);
        BSSurface<T,dim> Tuv{
            poles_T,
            ku,kv,
            p,q
        };

        // // uniformize knots
        // auto [ku_Lu, mu_Lu] = unflat_knots(Lu.knotsFlatsU());
        // auto [kv_Lu, mv_Lu] = unflat_knots(Lu.knotsFlatsV());
        // auto [ku_Lv, mu_Lv] = unflat_knots(Lv.knotsFlatsU());
        // auto [kv_Lv, mv_Lv] = unflat_knots(Lv.knotsFlatsV());



        auto poles_gordon = Lu.poles();

        std::transform(
            Lv.poles().begin(),Lv.poles().end(),
            poles_gordon.begin(),
            poles_gordon.begin(),
            [](const auto &v_, const auto &g_){return g_ + v_;}
        );
        std::transform(
            Tuv.poles().begin(),Tuv.poles().end(),
            poles_gordon.begin(),
            poles_gordon.begin(),
            [](const auto &uv_, const auto &g_){return g_ - uv_;}
        );

        return gbs::BSSurface<T,dim>{
            poles_gordon,
            ku, kv,
            p, q
        };
    }
*/
} // namespace gbs