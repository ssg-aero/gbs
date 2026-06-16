#pragma once
#include <cassert>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
#include <gbs/bscapprox.h>
#include <gbs/bscanalysis.h>
#include <boost/range/combine.hpp>
#include "loft/loftBase.h"
#include "loft/loftGeneric.h"
#include "loft/loftCurves.h"
#include <gbs/bsstools.h>
namespace gbs
{

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
            GBS_PAR_EXEC
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
                GBS_PAR_EXEC
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

    /**
     * @brief Finds the maximum degree among a range of curves.
     * 
     * This function iterates over a range of curves, applying conditional_dereference to handle
     * both pointer and non-pointer types seamlessly. It calculates the degree of each curve in the range
     * and returns the maximum degree found.
     * 
     * @tparam InputIt The type of the iterators defining the range, deduced automatically.
     * @param first Iterator to the first element in the range.
     * @param last Iterator to the element following the last element in the range.
     * @return The maximum degree among the curves in the given range.
     */
    template <typename InputIt>
    auto get_max_degree(InputIt first, InputIt last) -> size_t
    {
        auto get_degree = [](const auto &crv){
            return conditional_dereference(crv).degree();
        };

        auto max_deg_curve = *std::max_element(
            first, last,
            [&get_degree](const auto &bsc1, const auto &bsc2){return get_degree(bsc1) < get_degree(bsc2);}
        );
        return get_degree(max_deg_curve);
    }

    /**
     * @brief Lift a model-space spine tangent into pole (homogeneous when rational) space.
     *
     * The spine is non-rational, so `spine.value(v,1)` is a model-space tangent
     * `V` (dim components). For a rational loft the surface control net lives in
     * homogeneous space, so the v-tangent the Hermite solve interpolates must be
     * homogeneous too.
     *
     * Moving a pole along the spine direction is a CONSTANT-weight model-space
     * displacement: with the homogeneous pole `P^w = (w*P, w)`,
     *     P^w(t) = ( w*(P + t*V), w )  =>  dP^w/dt = ( w*V, 0 ).
     * Hence the spatial part is scaled by the pole weight `w` and the weight
     * component is null (the spine carries no weight variation) — the 4th
     * component is zero by derivation, not to silence the compiler.
     *
     * For the non-rational case this is just `s * V` (no weight component), which
     * reproduces the previous behaviour bit-for-bit.
     *
     * @param s         scalar magnitude (chord-length / spine-distance scaling).
     * @param tg_model  model-space spine tangent `V` (dim components).
     * @param pole      the homogeneous (or plain) pole the tangent is attached to.
     */
    template <typename T, size_t dim, bool rational>
    inline auto spine_tangent_to_pole(T s, const std::array<T, dim> &tg_model, const point<T, dim + rational> &pole)
        -> point<T, dim + rational>
    {
        point<T, dim + rational> r;
        if constexpr (rational)
        {
            const T w = pole.back();
            for (size_t i = 0; i < dim; ++i)
                r[i] = s * w * tg_model[i];
            r[dim] = T{0}; // constant weight along the spine direction
        }
        else
        {
            for (size_t i = 0; i < dim; ++i)
                r[i] = s * tg_model[i];
        }
        return r;
    }

    template <typename T, size_t dim, bool rational,typename Container>
    auto get_tg_from_spine(const std::vector<points_vector<T, dim + rational>> &poles_v_lst, const BSCurve<T, dim> &spine,const Container &bs_lst_cpy)
    {
        std::vector< points_vector<T,dim + rational> > tg_v_lst{ poles_v_lst.size() };
        // compute positions on spine
        std::vector<T> v_spine(bs_lst_cpy.size());
        std::transform(
                GBS_PAR_EXEC
                bs_lst_cpy.begin(),
                bs_lst_cpy.end(),
                v_spine.begin(),
                [&spine](const Curve<T, dim> &profile_)
                {
                    return extrema_curve_curve<T,dim>(spine,profile_,1e-6)[0];
                }
        );
        // compute spine distances (model space — the spine is non-rational)
        std::vector<T> d(v_spine.size() - 1);
        std::transform(
            GBS_PAR_EXEC
            std::next(v_spine.begin()),
            v_spine.end(),
            v_spine.begin(),
            d.begin(),
            [&](const auto &v2_, const auto &v1_) {
                return norm(spine.value(v2_) - spine.value(v1_));
            });
        // compute spine tangents in MODEL space (dim components); they are lifted
        // into homogeneous pole space per-pole below (weight depends on the pole).
        std::vector<std::array<T,dim>> D(v_spine.size());
        std::transform(
            GBS_PAR_EXEC
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
                // Project the column's poles to model space so the chord-length
                // scaling is taken consistently with the model-space spine
                // distances `d` (for the rational case `pts` are homogeneous, so a
                // raw norm would mix in the weights — exactly the #58 inconsistency).
                points_vector<T, dim> pm(n_poles_v);
                std::transform(pts.begin(), pts.end(), pm.begin(),
                               [](const auto &p) { return weight_projected_pole<T, dim, rational>(p); });

                points_vector<T,dim + rational> tg_i(n_poles_v);
                tg_i[0] = spine_tangent_to_pole<T, dim, rational>(
                    norm(pm[1] - pm[0]) / d[0], D[0], pts[0]);
                for(auto k = 1 ; k < n_poles_v - 1; k++)
                {
                    tg_i[k] = spine_tangent_to_pole<T, dim, rational>(
                        (norm(pm[k+1] - pm[k]) + norm(pm[k] - pm[k-1])) / (d[k] + d[k-1]), D[k], pts[k]);
                }
                tg_i[n_poles_v-1] = spine_tangent_to_pole<T, dim, rational>(
                    norm(pm[n_poles_v-1] - pm[n_poles_v-2]) / d[n_poles_v-2], D[n_poles_v-1], pts[n_poles_v-1]);
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
        constexpr size_t dim_t = dim + rational;
        size_t q = std::max(std::min(v_degree_max, n_poles_v - 1), size_t{1}); // degree V, clamped to a valid range
        auto kv_flat = build_simple_mult_flat_knots<T>(v.front(),v.back(),n_poles_v * 2,q); //dif

        // One constraint column per u-pole: at each section, position + spine
        // tangent (nc == 2). Assemble them straight into the shared loft core's RHS
        // (rows = 2*n_poles_v v-DOFs, cols = n_poles_u*dim_t) and let it factorize
        // the v-system ONCE, batch-solve every column, and emit poles already in
        // the surface layout — the same core the plain loft uses, so spine and
        // non-spine no longer carry parallel pole builders (issue #40 PR2/PR3 / #44).
        MatrixX<T> B(Eigen::Index(2 * n_poles_v), Eigen::Index(n_poles_u * dim_t));
        for (size_t i = 0; i < n_poles_u; ++i)
            for (size_t j = 0; j < n_poles_v; ++j)
                for (size_t d = 0; d < dim_t; ++d)
                {
                    B(Eigen::Index(2 * j + 0), Eigen::Index(i * dim_t + d)) = poles_v_lst[i][j][d];
                    B(Eigen::Index(2 * j + 1), Eigen::Index(i * dim_t + d)) = tg_i[i][j][d];
                }
        auto poles = build_loft_surface_poles_core<T, dim_t, 2>(B, kv_flat, v, q);

        // build surf
        using bss_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;
        return bss_type( poles,  ku_flat, kv_flat, p , q );
    }

    template <typename T, size_t dim>
    auto loft(const std::list<BSCurve<T, dim>> &bs_lst, const BSCurve<T, dim> &spine, size_t v_degree_max = 3)
    {
        std::list<BSCurveGeneral<T, dim, false> *> p_curve_lst(bs_lst.size());
        std::transform(
            bs_lst.begin(),
            bs_lst.end(),
            p_curve_lst.begin(),
            [](auto &bs) 
                { 
                    auto p_bs = static_cast<const BSCurveGeneral<T, dim, false> *>(&bs);
                    return const_cast<BSCurveGeneral<T, dim, false> *>(p_bs); 
                }
            ); // pointer to temporary address works within the scope

        return loft<T, dim, false>(p_curve_lst, spine, v_degree_max);
    }

    /**
     * @brief Builds a vector of intersection points between a given curve and a range of other curves.
     *
     * This function calculates intersection points between a provided curve 'crv_u' and a range of curves
     * defined by iterators 'first_curve' and 'last_curve'. If the distance is greater than a specified tolerance, 
     * it throws an exception.
     * 
     * @tparam T The numeric type used in the curves (e.g., double, float).
     * @tparam dim The dimensionality of the curve (e.g., 2D, 3D).
     * @tparam InputIt The type of the iterators defining the range of curves.
     * @param first_curve Iterator to the first curve in the range.
     * @param last_curve Iterator to the element following the last curve in the range.
     * @param crv The curve to intersect with each curve in the range.
     * @param tol The tolerance value used to determine if two curves are touching.
     * @return std::vector<T> A vector containing parameters (of type T) at each intersection point with 'crv_u'.
     * 
     * @throws std::invalid_argument if any curve in the range does not touch 'crv' within the specified tolerance.
     */
    template <typename T, size_t dim, typename InputIt>
    auto build_intersections(InputIt first_curve, InputIt last_curve, const Curve<T, dim> &crv, T tol = 1e-6) -> std::vector<T>
    {
        auto nu = std::distance(first_curve, last_curve);
        std::vector<T> u(nu);

        // Compute intersections
        auto f_intersection = [tol, &crv](const auto &crv_) {
            auto [u_, v_, d] = extrema_curve_curve( crv, conditional_dereference(crv_), tol);
            if (d > tol) {
                throw std::invalid_argument("A curve from set is not touching the curve. Distance is: "+std::to_string(d));
            }
            return u_;
        };

        std::transform(
            GBS_PAR_EXEC
            first_curve, last_curve,
            u.begin(),
            f_intersection
        );

        return u;
    }

    /**
     * @brief Builds intersections between two sets of curves.
     * 
     * This function calculates the intersection points between every curve in the first set (from first_u to last_u)
     * with every curve in the second set (from first_v to last_v). It uses the build_intersections function to find
     * intersection parameters, then calculates the actual points of intersection, checking if they are within
     * the specified tolerance.
     * 
     * @tparam T The numeric type used in the curves (e.g., double, float).
     * @tparam dim The dimensionality of the curve (e.g., 2D, 3D).
     * @tparam InputIt The type of the iterators defining the range of curves.
     * @param first_u Iterator to the first curve in the first set.
     * @param last_u Iterator to the element following the last curve in the first set.
     * @param first_v Iterator to the first curve in the second set.
     * @param last_v Iterator to the element following the last curve in the second set.
     * @param tol The tolerance value used to determine if two curves are intersecting.
     * @return A tuple containing the intersection parameters and points.
     * 
     * @throws std::invalid_argument if any pair of curves from the two sets do not intersect within the specified tolerance.
     */
    template <typename T, size_t dim, typename InputIt>
    auto build_intersections(InputIt first_u, InputIt last_u, InputIt first_v, InputIt last_v, T tol = 1e-6)
    {
        auto nu = std::distance(first_u, last_u);
        auto nv = std::distance(first_v, last_v);

        // Precompute intersections
        auto intersections_v = build_intersections<T, dim>(first_u, last_u, conditional_dereference(*first_v), tol);
        auto intersections_u = build_intersections<T, dim>(first_v, last_v, conditional_dereference(*first_u), tol);

        points_vector<T,dim> Q(nu * nv);

        for (size_t i = 0; i < nu; ++i)
        {
            const auto &curve_u = conditional_dereference(*(std::next(first_u, i)));
            for (size_t j = 0; j < nv; ++j)
            {
                const auto &curve_v = conditional_dereference(*(std::next(first_v, j)));
                auto Q_u = curve_u.value(intersections_u[j]);
                auto Q_v = curve_v.value(intersections_v[i]);

                auto d = distance(Q_u, Q_v);
                if (d > tol)
                {
                    throw std::invalid_argument("U V curve are not intersecting. Distance is: " + std::to_string(d));
                }

                Q[i * nv + j] = 0.5 * (Q_u + Q_v);
            }
        }

        return std::make_tuple(intersections_u, intersections_v, Q);
    }

    template <typename T, size_t dim, typename ForwardIt>
    auto gordon(ForwardIt first_u,ForwardIt last_u,ForwardIt first_v,ForwardIt last_v, T tol = 1e-6) -> BSSurface<T,dim>
    {

        // gather informations
        auto p  = get_max_degree( first_u, last_u);
        auto q  = get_max_degree( first_v, last_v);
        // compute and check intersections
        auto [u, v, Q] = build_intersections<T, dim>(first_u, last_u, first_v, last_v, tol);
        // Cap max degree according to curves numbers
        auto p_ = std::min(p, u.size() - 1);
        auto q_ = std::min(q, v.size() - 1);

        // build intermediates surfaces
        // Lu
        auto Lu = loft_generic<T, dim>(first_u, last_u, v, q_);
        for(auto j{Lu.degreeV()}; j<q; j++)
        {
            Lu.increaseDegreeV();
        }
        // Tu
        auto Lv = loft_generic<T, dim>(first_v, last_v, u, p_);
        for( auto i{Lv.degreeV()}; i<p; i++) // Lv is UV inverted
        {
            Lv.increaseDegreeV();
        }
        Lv.invertUV();

        assert(Lv.degreeU() == Lu.degreeU());
        assert(Lv.degreeV() == Lu.degreeV());

        // Tuv
        auto ku = build_simple_mult_flat_knots(u, p_);
        auto kv = build_simple_mult_flat_knots(v, q_);
        auto poles_T = build_poles(Q, ku, kv, u, v, p_, q_);
        BSSurface<T,dim> Tuv{
            poles_T,
            ku,kv,
            p_,q_
        };
        for(auto i{p_}; i<p; i++)
        {
            Tuv.increaseDegreeU();
        }
        for(auto j{q_}; j<q; j++)
        {
            Tuv.increaseDegreeV();
        }

        // uniformize knots
        auto Lu_info = Lu.info();
        auto Lv_info = Lv.info();
        auto Tuv_info = Tuv.info();
        std::vector<BSSurfaceInfo<T, dim >> surfaces_info{Lu_info, Lv_info, Tuv_info};
        
        unify_knots(surfaces_info);

        auto poles_u = flatten_poles(std::get<0>(surfaces_info[0]));
        auto poles_v = flatten_poles(std::get<0>(surfaces_info[1]));
        auto poles_uv= flatten_poles(std::get<0>(surfaces_info[2]));

        auto ku_gordon = std::get<1>(surfaces_info[0]);
        auto kv_gordon = std::get<2>(surfaces_info[0]);

        auto poles_gordon = poles_u;

        std::ranges::transform(
            poles_gordon,
            poles_v,
            poles_gordon.begin(),
            [](const auto &g_, const auto &v_){return g_ + v_;}
        );
        std::ranges::transform(
            poles_gordon,
            poles_uv,
            poles_gordon.begin(),
            [](const auto &g_, const auto &uv_){return g_ - uv_;}
        );

        return BSSurface<T,dim>{
            poles_gordon,
            ku_gordon, kv_gordon,
            p, q
        };
    }
// */
} // namespace gbs