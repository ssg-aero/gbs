#pragma once
#include <gbs/bssurf.h>
#include <gbs/bscapprox.h>
#include <gbs/transformpoints.h>
#include <gbs/bsctools.h>
namespace gbs 
{
    /**
     * @brief Create a surface's copy with and additional dimension. The value of the dimension can be specified, the default is 0.
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational 
     * @param srf 
     * @param val 
     * @return auto 
    **/
    template <typename T, size_t dim, bool rational>
    auto add_dimension(const BSSurfaceGeneral<T, dim, rational> &srf, T val = 0.)
    {
        points_vector<T, dim + rational + 1> poles(srf.poles().size());
        std::transform(
            std::execution::par,
            srf.poles().begin(),
            srf.poles().end(),
            poles.begin(),
            [&val](const auto &p_)
            {
                return add_dimension<T,dim+rational>(p_, val);
            });
        using bs_type = typename std::conditional<rational, BSSurfaceRational<T, dim + 1>, BSSurface<T, dim + 1>>::type;
        return bs_type(poles, srf.knotsFlatsU(), srf.knotsFlatsV(), srf.degreeU(), srf.degreeV());
    }

    /// @brief Extend surface to match curve
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param crv 
    /// @param natural_end 
    /// @param max_cont 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto extention_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        points_vector<T,dim+rational> poles_new;

        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        using bss_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;

        auto nu = srf.nPolesU();
        auto nv = srf.nPolesV();
        auto ku = srf.knotsFlatsU();
        auto kv = srf.knotsFlatsV();
        auto p = srf.degreeU();
        auto q = srf.degreeV();
        auto du = 0.;
        auto np = 100;
        auto v1 = kv.front();
        auto v2 = kv.back();
        auto u2 = ku.back();

        for(size_t i {} ; i < np ; i++)
        {
            auto v = v1 + (v2-v1)*i/(np-1.);
            auto t = srf(u2,v,1,0);
            auto dl = norm(srf(u2,v)-crv(v));
            du += dl / norm(t) / np;
        }
        auto u_new = u2 + du;

        size_t p_ext;
        std::vector<T> k_ext;
        for(size_t i {} ; i < nv ; i++)
        {
            auto poles = points_vector<T,dim+rational>{srf.poles_begin()+i*nu,srf.poles_begin()+i*nu+nu};
            auto iso  = bsc_type{poles,ku,p};
            auto crv_ext = bsc_type{extention_to_point<T, dim, rational>(
                iso,
                weight_projected_pole<T, dim, rational>(crv.poles()[i]),
                u2,
                u_new,
                natural_end,
                max_cont
            )};
            auto poles_ext = crv_ext.poles();
            poles_new.insert(poles_new.end(),poles_ext.begin(),poles_ext.end());
            if(i==0)
            {
                p_ext = crv_ext.degree();
                k_ext = crv_ext.knotsFlats();
            }
        }

        return bss_type(poles_new,k_ext,kv,p_ext,q) ;
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const BSSurfaceGeneral<T,dim,rational1> &srf1,const BSSurfaceGeneral<T,dim,rational2> &srf2)
    {
        const auto rational = rational1 || rational2;
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        bs_type srf1_ { srf1.poles(), srf1.knotsFlatsU(), srf1.knotsFlatsV(), srf1.degreeU(), srf1.degreeV() };
        bs_type srf2_ { srf2.poles(), srf2.knotsFlatsU(), srf2.knotsFlatsV(), srf2.degreeU(), srf2.degreeV() };

        while(srf1_.degreeU()<srf2_.degreeU()) srf1_.increaseDegreeU();
        while(srf2_.degreeU()<srf1_.degreeU()) srf2_.increaseDegreeU();

        auto nkv = srf2_.knotsFlatsV().size(); // TODO insert v knots
        if(nkv!=srf1_.knotsFlatsV().size())
            throw std::runtime_error("Not implemented");
        for(size_t i{}; i < nkv ; i ++)
        {
            if(std::fabs( srf2_.knotsFlatsV()[i] - srf1_.knotsFlatsV()[i]) > knot_eps<T>)
                throw std::runtime_error("Not implemented");
        }
        // join
        auto nu = srf1_.nPolesU();
        auto nv = srf1_.nPolesV();
        auto poles_ext = srf2_.poles();
        auto nu_ext = srf2_.nPolesU();
        auto ku_ext = srf2_.knotsFlatsU();
        auto ku = srf1_.knotsFlatsU();
        auto kv = srf1_.knotsFlatsV();
        auto p = srf1_.degreeU();
        auto q = srf1_.degreeV();
        auto p_ext = srf2_.degreeU();

        points_vector<T,dim+rational> poles_new;
        for(size_t i{}; i < nv; i++)
        {
            poles_new.insert(
                poles_new.end(),
                std::next( srf1_.poles().begin(), i * nu ),
                std::next( srf1_.poles().begin(), i * nu + nu)
            );
            poles_new.insert(
                poles_new.end(),
                std::next( poles_ext.begin(), i * nu_ext +1 ),
                std::next( poles_ext.begin(), i * nu_ext + nu_ext)
            );
        }
        auto dk = ku.back() - ku_ext.front();
        std::transform(
            ku_ext.begin(),
            ku_ext.end(),
            ku_ext.begin(),
            [dk](auto k){return k+dk;}
        );
        std::vector<double> k_new {ku.begin(),std::next(ku.end(),-1)};
        k_new.insert(
            k_new.end(),
            std::next(ku_ext.begin(),p_ext+1),
            ku_ext.end()
        );

        return bs_type{
            poles_new,
            k_new,
            kv,
            p,
            q
        };
    }

    /// @brief Split surface at u parameter, aka in v direction
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param u 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto split_u(const BSSurfaceGeneral<T,dim,rational> &srf, T u)
    {
        auto [u1, u2, v1, v2] = srf.bounds();
        assert ( u1<=u && u <=u2);

        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 

        bs_type srf1{srf};
        bs_type srf2{srf};

        srf1.trimU(u1, u);
        srf2.trimU(u, u2);

        return std::make_pair(
            srf1, srf2
        );
    }
    /// @brief Split surface at v parameter, aka in u direction
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param v 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto split_v(const BSSurfaceGeneral<T,dim,rational> &srf, T v)
    {
        auto [u1, u2, v1, v2] = srf.bounds();
        assert ( v1<=v && v <=v2);

        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 

        bs_type srf1{srf};
        bs_type srf2{srf};

        srf1.trimV(v1, v);
        srf2.trimV(v, v2);

        return std::make_pair(
            srf1, srf2
        );
    }

    template <typename T, size_t dim, bool rational>
    auto extended_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto srf_ext = extention_to_curve(srf,crv,natural_end,max_cont);
        return join(srf,srf_ext);
    }

    template <typename T, size_t dim, bool rational>
    auto extention(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        auto np = srf.nPolesV()*10;
        auto [u1,u2,v1,v2] = srf.bounds();
        points_vector<T,dim> pts(np);
        std::vector<T> v_crv(np);
        for(size_t i {} ; i < np ; i++)
        {
            auto v_ = v1 + (v2-v1)*i/(np-1.);
            v_crv[i] = v_;
            auto pt = srf(u2,v_);
            auto tg = srf(u2,v_,1,0);
            tg = tg / norm(tg);
            pt = pt + tg * l_ext;
            pts[i] = pt;
        }
        auto nv = srf.nPolesV();
        auto kv = srf.knotsFlatsV();
        auto q  = srf.degreeV();
        auto crv = bsc_type{approx_bound_fixed(pts,q,nv,v_crv,kv)};
        return extention_to_curve(srf,crv,natural_end,max_cont);
    }

    template <typename T, size_t dim, bool rational>
    auto extended(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto srf_ext = extention(srf,l_ext,natural_end,max_cont);
        return join(srf,srf_ext);
    }

    enum class SurfaceBound { U_START, U_END, V_START, V_END};

    template <typename T, size_t dim, bool rational>
    auto extention(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, SurfaceBound location , bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        switch (location)
        {
        case SurfaceBound::U_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.reverseU();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break;
            case SurfaceBound::V_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                srf_.reverseU();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break;        
            case SurfaceBound::V_END:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break; 
        default: // SurfaceBound::U_END
            return extention(srf,l_ext,natural_end, max_cont);
            break;
        }
    }
    template <typename T, size_t dim, bool rational>
    auto extended(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, SurfaceBound location, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        // auto srf_ext = extention(srf,l_ext,location,natural_end,max_cont);
        // return join(srf,srf_ext);
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        switch (location)
        {
        case SurfaceBound::U_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.reverseU();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.reverseU();
                return srf_ext;
            }
            break;
            case SurfaceBound::V_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                srf_.reverseU();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.reverseU();
                srf_ext.invertUV();
                return srf_ext;
            }
            break;        
            case SurfaceBound::V_END:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.invertUV();
                return srf_ext;
            }
            break; 
        default: // SurfaceBound::U_END
            return extended(srf,l_ext,natural_end, max_cont);
            break;
        }
    }
}