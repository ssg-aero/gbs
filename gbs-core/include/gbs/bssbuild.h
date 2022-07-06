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

    template <typename T, size_t dim,bool rational, typename ForwardIt>
    // auto loft(const std::vector<std::shared_ptr<BSCurveGeneral<T, dim,rational>>> &bs_lst, const std::vector<T> &v,size_t q)
    auto loft(ForwardIt first, ForwardIt last, const std::vector<T> &v,size_t q)
    {
        // auto nv = bs_lst.size();
        auto nv = std::distance(first, last);
        if(nv != v.size())
        {
            throw std::invalid_argument("incompatible sizes");
        }
        // Get curves caracteristics
        // auto [ u1, u2 ] = bs_lst.front()->bounds();
        auto [ u1, u2 ] = (*first)->bounds();
        auto max_deg_curve = *std::max_element(
            // bs_lst.begin(), bs_lst.end(),
            first, last,
            [](const auto &bsc1, const auto &bsc2){return bsc1->degree() < bsc2->degree();}
        );
        auto p = max_deg_curve->degree();
        // Get poles array with unified degree and bounds
        std::vector < std::pair<points_vector<T, dim + rational>, std::vector<T> > > poles_knots(nv);
        std::transform(
            std::execution::par,
            // bs_lst.begin(), bs_lst.end(),
            first, last,
            poles_knots.begin(),
            [u1, u2, p](const auto &bsc)
            {
                auto poles{bsc->poles()};
                auto knots{bsc->knotsFlats()};
                auto deg{bsc->degree()};
                change_bounds(u1, u2, knots);
                for(size_t i{}; i < p - deg; i++)
                {
                    increase_degree(knots, poles, deg);
                    deg++;
                }
                return std::make_pair(poles, knots);
            }
        );
        // Unify knots
        auto &pk1 = poles_knots.front();
        std::for_each(
            std::next(poles_knots.begin()),
            poles_knots.end(),
            [&pk1,p](auto &pk) { unify_knots(pk1,pk,p); }
        );
        std::for_each(
            std::next(poles_knots.begin()),
            poles_knots.end(),
            [&pk1,p](auto &pk) { unify_knots(pk,pk1,p); }
        );
        // interpolate in V direction
        auto nu =  poles_knots.front().first.size();
        std::vector< points_vector<T, dim + rational> > pole_interp(nu) ;
        auto indexes = make_range<size_t>(0, nu-1);
        auto kv_flat = build_simple_mult_flat_knots<T>(v,q);
        std::transform(
            std::execution::par,
            indexes.begin(), indexes.end(),
            pole_interp.begin(),
            [nv,q,&poles_knots,&v,&kv_flat](auto i)
            {
                points_vector<T,dim+rational> poles(nv);
                for(size_t j{}; j < nv; j++)
                {
                    poles[j] = poles_knots[j].first[i];
                }
                std::vector<gbs::constrType<T, dim+rational, 1>> Q(poles.size());
                std::transform(poles.begin(),poles.end(),Q.begin(),[](const auto &pt_){return constrType<T, dim+rational, 1>{pt_};});
                return build_poles<T,dim,1>(Q, kv_flat, v, q);
            }
        );
        // invert pole alignment
        points_vector<T, dim + rational> poles(nu*nv);
        for (size_t j{}; j < nv; j++)
        {
            for (size_t i{}; i < nu; i++)
            {
                poles[ i + nu * j] = pole_interp[i][j];
            }
        }
        // build surf
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;
        return bs_type( poles,  poles_knots.front().second, kv_flat, p , q );
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
        auto bs_lst_cpy = get_BSCurves_ptr_cpy(bs_lst);
        // static auto dim_poles = dim + rational; 

        unify_curves_degree(bs_lst_cpy);
        unify_curves_knots(bs_lst_cpy);

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
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;
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

    template <typename T, size_t dim>
    auto loft(const std::list<BSCurve<T, dim>> &bs_lst, size_t v_degree_max = 3)
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

        return gbs::loft<T, dim, false>(p_curve_lst, v_degree_max);
    }
    template <typename T, size_t dim>
    auto loft(const std::list<Curve<T, dim>*> &crv_lst, size_t v_degree_max = 3,T dev = 0.01, size_t np=100, size_t deg_approx=5)
    {
        std::list<BSCurve<T, dim>> bs_lst(crv_lst.size());
        std::transform(
            crv_lst.begin(),
            crv_lst.end(),
            bs_lst.begin(),
            [dev,np,deg_approx](const auto &crv)
            {
                auto p_bs = dynamic_cast<const gbs::BSCurve<T, dim>*>(crv);
                if(p_bs)
                {
                    return *p_bs;
                }
                else
                {
                    return approx(*crv,dev,deg_approx,gbs::KnotsCalcMode::CHORD_LENGTH,np);
                }
            }
        );
        return loft(bs_lst,v_degree_max);
    }

    template <typename T, size_t dim>
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, size_t v_degree_max = 3,T dev = 0.01, size_t np=100, size_t deg_approx=5)
    {
        std::list<BSCurve<T, dim>> bs_lst(crv_lst.size());
        std::transform(
            crv_lst.begin(),
            crv_lst.end(),
            bs_lst.begin(),
            [dev,np,deg_approx](const auto &crv)
            {
                auto p_bs = std::dynamic_pointer_cast<const gbs::BSCurve<T, dim>>(crv);
                if(p_bs)
                {
                    return *p_bs;
                }
                else
                {
                    return approx<T,dim>(*crv,dev,deg_approx,gbs::KnotsCalcMode::CHORD_LENGTH,np);
                }
            }
        );
        return loft(bs_lst,v_degree_max);
    }

    template <typename T, size_t dim>
    auto loft(const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst, const std::vector<T> &v, size_t q,T dev = 0.01, size_t np=100, size_t deg_approx=5)
    {
        std::vector<std::shared_ptr<BSCurveGeneral<T, dim, false>>> bs_lst(crv_lst.size());
        std::transform(
            crv_lst.begin(),
            crv_lst.end(),
            bs_lst.begin(),
            [dev,np,deg_approx](const auto &crv)
            {
                auto p_bs = std::dynamic_pointer_cast<BSCurve<T, dim>>(crv);
                if(p_bs)
                {
                    return p_bs;
                }
                else
                {
                    return std::make_shared< BSCurve<T, dim> >( approx<T,dim>(*crv,dev,deg_approx,gbs::KnotsCalcMode::CHORD_LENGTH,np) );
                }
            }
        );
        return loft<T,dim,false>(bs_lst.begin(), bs_lst.end() ,v, q);
    }

    
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

} // namespace gbs