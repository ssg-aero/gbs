#pragma once
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
            std::execution::par,
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

        for(size_t j = 0; j < nv; ++j) {
            const auto &curve_v = conditional_dereference(*(std::next(first_v, j)));
            for(size_t i = 0; i < nu; ++i) {
                const auto &curve_u = conditional_dereference(*(std::next(first_u, i)));
                auto Q_u = curve_u.value(intersections_u[j]);
                auto Q_v = curve_v.value(intersections_v[i]);

                auto d = distance(Q_u, Q_v);
                if( d > tol) {
                    throw std::invalid_argument("U V curve are not intersecting. Distance is: "+std::to_string(d));
                }

                // Q[i + j * nu] = 0.5 * ( Q_u + Q_v );
                Q[i * nv + j ] = 0.5 * ( Q_u + Q_v );
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

        // build intermediates surfaces
        auto Lu = loft_generic<T,dim>(first_u, last_u,v,q);
        auto Lv = loft_generic<T,dim>(first_v, last_v,u,p);
        Lv.invertUV();

        auto ku = build_simple_mult_flat_knots(u, p);
        auto kv = build_simple_mult_flat_knots(v, q);

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

        return BSSurface<T,dim>{
            poles_gordon,
            ku, kv,
            p, q
        };
    }
// */
} // namespace gbs