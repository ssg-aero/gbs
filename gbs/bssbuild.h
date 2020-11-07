#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
namespace gbs
{
    // TODO try to make this work
    // template <typename T, size_t dim,bool rational, template<typename> typename Container >
    // auto loft(const Container<BSCurveGeneral<T, dim,rational>*> &bs_lst,size_t v_degree_max = 3)
    template <typename T, size_t dim,bool rational>
    auto loft(const std::list<BSCurveGeneral<T, dim,rational>*> &bs_lst,size_t v_degree_max = 3)
    {
        if(bs_lst.size()<2)
        {
            throw std::length_error("loft needs at least 2 curves.");
        }
        // auto bs_lst_cpy(bs_lst);
        typedef std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type bs_type ;
        typedef std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bss_type ;
        std::list<bs_type> bs_lst_cpy(bs_lst.size());
        std::transform(
            bs_lst.begin(),
            bs_lst.end(),
            bs_lst_cpy.begin(),
            [](const auto *p_bs){return bs_type(p_bs->poles(),p_bs->knotsFlats(),p_bs->degree()); }
        );
        unify_degree(bs_lst_cpy);
        unify_knots(bs_lst_cpy);
        auto n_poles_v = bs_lst_cpy.size();
        auto n_poles_u = bs_lst_cpy.front().poles().size();
        auto ku_flat = bs_lst_cpy.front().knotsFlats();
        auto p = bs_lst_cpy.front().degree();
        // get v aligned poles
        std::vector< points_vector<T,dim> > poles_v_lst{ n_poles_u };
        std::vector< std::vector<T> > poles_v_params{ n_poles_u };
        for(auto i = 0 ; i < n_poles_u; i++)
        {
            points_vector<T,dim> poles_v{ n_poles_v };
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
        // interpolate poles
        size_t q = fmax(fmin(v_degree_max,n_poles_v-1),1); // degree V
        auto kv_flat = build_simple_mult_flat_knots<T>(v,n_poles_v,q);
        points_vector<T,dim> poles;
        std::vector<constrType<T,dim,1> > Q(n_poles_v);
        for(auto pts : poles_v_lst)
        {
            std::transform(
                std::execution::par,
                pts.begin(),
                pts.end(),
                Q.begin(),
                [](const auto &p){constrType<T,dim,1> c{p}; return c;}
            );
            auto poles_v = build_poles(Q, kv_flat, v, q);
            poles.insert( poles.end(), poles_v.begin(),poles_v.end() );
        }
        // build surf
        return BSSurface3d_d( invert_uv_poles(poles,n_poles_u),  ku_flat, kv_flat, p , q );
    }
} // namespace gbs