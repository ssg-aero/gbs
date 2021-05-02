#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
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
        return std::make_tuple(poles_v_lst,v);
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
        
        auto [poles_v_lst,v] = get_v_aligned_poles<T,dim,rational>(n_poles_u,n_poles_v,bs_lst_cpy);

        // interpolate poles
        size_t q = fmax(fmin(v_degree_max,n_poles_v-1),1); // degree V
        auto kv_flat = build_simple_mult_flat_knots<T>(v,q);
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
        return bs_type( inverted_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
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
    auto loft(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst,const BSCurve<T, dim> &spine,size_t v_degree_max = 3)
    {
        if(bs_lst.size()<2)
        {
            throw std::length_error("loft needs at least 2 curves.");
        }
        typedef std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type bsc_type;
        std::list<bsc_type> bs_lst_cpy = get_BSCurves_cpy(bs_lst);

        // static auto dim_poles = dim + rational; 

        unify_degree(bs_lst_cpy);
        unify_knots(bs_lst_cpy);

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
        typedef std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bss_type;
        return bss_type( inverted_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
    }
} // namespace gbs