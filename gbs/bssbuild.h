#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
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
     * @return auto 
    **/
    template <typename T, size_t  dim, bool rational>
    auto get_BSCurves_cpy(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst)
    {
        typedef std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type bs_type;
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
        return std::make_tuple(poles_v_params,poles_v_lst,v);
    }

    // TODO try to make template container work
    // template <typename T, size_t dim,bool rational, template<typename> typename Container >
    // auto loft(const Container<BSCurveGeneral<T, dim,rational>*> &bs_lst,size_t v_degree_max = 3)
    template <typename T, size_t dim,bool rational>
    auto loft(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst,size_t v_degree_max = 3)
    {
        if(bs_lst.size()<2)
        {
            throw std::length_error("loft needs at least 2 curves.");
        }
        auto bs_lst_cpy = get_BSCurves_cpy(bs_lst);
        // static auto dim_poles = dim + rational; 

        unify_degree(bs_lst_cpy);
        unify_knots(bs_lst_cpy);

        auto n_poles_v = bs_lst_cpy.size();
        auto n_poles_u = bs_lst_cpy.front().poles().size();
        auto ku_flat = bs_lst_cpy.front().knotsFlats();
        auto p = bs_lst_cpy.front().degree();
        
        auto [poles_v_params,poles_v_lst,v] = get_v_aligned_poles<T,dim,rational>(n_poles_u,n_poles_v,bs_lst_cpy);

        // interpolate poles
        size_t q = fmax(fmin(v_degree_max,n_poles_v-1),1); // degree V
        auto kv_flat = build_simple_mult_flat_knots<T>(v,n_poles_v,q);
        points_vector<T,dim + rational> poles;
        std::vector<constrType<T,dim + rational,1> > Q(n_poles_v);
        for(auto pts : poles_v_lst)
        {
            std::transform(
                std::execution::par,
                pts.begin(),
                pts.end(),
                Q.begin(),
                [](const auto &p){constrType<T,dim + rational,1> c{p}; return c;}
            );
            auto poles_v = build_poles(Q, kv_flat, v, q);
            poles.insert( poles.end(), poles_v.begin(),poles_v.end() );
        }
        // build surf
        typedef std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bs_type;
        return bs_type( invert_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
    }

    template <typename T, size_t dim,bool rational>
    auto loft(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst,const BSCurve<T, dim> &spine,size_t v_degree_max = 3)
    {
        if(bs_lst.size()<2)
        {
            throw std::length_error("loft needs at least 2 curves.");
        }
        auto bs_lst_cpy = get_BSCurves_cpy(bs_lst);

        // static auto dim_poles = dim + rational; 

        unify_degree(bs_lst_cpy);
        unify_knots(bs_lst_cpy);

        auto n_poles_v = bs_lst_cpy.size(); 
        auto n_poles_u = bs_lst_cpy.front().poles().size();
        auto ku_flat = bs_lst_cpy.front().knotsFlats();
        auto p = bs_lst_cpy.front().degree();
        auto [poles_v_params, poles_v_lst, v] = get_v_aligned_poles<T, dim, rational>(n_poles_u, n_poles_v, bs_lst_cpy);
        // // get v aligned poles
        // std::vector< points_vector<T,dim + rational> > poles_v_lst{ n_poles_u };
        // std::vector< std::vector<T> > poles_v_params{ n_poles_u };
        // for(auto i = 0 ; i < n_poles_u; i++)
        // {
        //     points_vector<T,dim + rational> poles_v{ n_poles_v };
        //     std::transform(
        //         std::execution::par,
        //         bs_lst_cpy.begin(),
        //         bs_lst_cpy.end(),
        //         poles_v.begin(),
        //         [&i](const auto &crv)
        //         {
        //             return crv.poles()[i];
        //         }
        //     );
        //     poles_v_params[i] = curve_parametrization(poles_v, KnotsCalcMode::CHORD_LENGTH, true);
        //     std::swap( poles_v_lst[i] , poles_v );
        // }
        // // build parametrization
        // std::vector<T> v(n_poles_v,0.);
        // for(auto j = 0 ; j < n_poles_v; j++)
        // {
        //     for(auto i = 0 ; i < n_poles_u; i++)
        //     {
        //         v[j] += poles_v_params[i][j];
        //     }
        //     v[j] /=  n_poles_u;
        // }
        // compute positions on spine
        std::vector<T> v_spine(bs_lst_cpy.size());
        std::transform(
                std::execution::par,
                bs_lst_cpy.begin(),
                bs_lst_cpy.end(),
                v_spine.begin(),
                [&spine](const auto &profile_)
                {
                    return extrema_CC(spine,profile_,1e-6).u1;
                }
        );
        // interpolate poles
        size_t q = fmax(fmin(v_degree_max,n_poles_v-1),1); // degree V // diff
        auto kv_flat = build_simple_mult_flat_knots<T>(v,n_poles_v * 2,q); //dif
        points_vector<T,dim + rational> poles;
        std::vector<constrType<T,dim + rational,2> > Q(n_poles_v); //diff

        for (auto pts : poles_v_lst)
        {
            auto pt_and_bs = boost::combine(pts,v_spine);
            std::transform(
                // std::execution::par, // not working with boost::combine
                pt_and_bs.begin(),
                pt_and_bs.end(),
                Q.begin(),
                [&](const auto &pt_and_bs_){
                    auto [p,v_] = pt_and_bs_;
                    return constrType<T,dim + rational,2> {p,spine.value(v_,1)};
                    }
            );
            auto poles_v = build_poles(Q, kv_flat, v, q);
            poles.insert( poles.end(), poles_v.begin(),poles_v.end() );
        }
        // build surf
        typedef std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bs_type;
        return bs_type( invert_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
    }
} // namespace gbs