#pragma once
#include <gbs/bssurf.h>
#include <gbs/bscapprox.h>
namespace gbs 
{
    template <typename T, size_t dim, bool rational>
    auto extention_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        points_vector<T,dim+rational> poles_new;

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
            BSCurve iso {poles,ku,p};
            auto crv_ext = extention_to_point(iso,crv.poles()[i],u2,u_new,natural_end,max_cont);
            auto poles_ext = crv_ext.poles();
            poles_new.insert(poles_new.end(),poles_ext.begin(),poles_ext.end());
            if(i==0)
            {
                p_ext = crv_ext.degree();
                k_ext = crv_ext.knotsFlats();
            }
        }

        return BSSurface(poles_new,k_ext,kv,p_ext,q) ;
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const BSSurfaceGeneral<T,dim,rational1> &srf1,const BSSurfaceGeneral<T,dim,rational2> &srf2)
    {
        typedef std::conditional<rational1 || rational2,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bs_type; 
        bs_type srf1_ { srf1.poles(), srf1.knotsFlatsU(), srf1.knotsFlatsV(), srf1.degreeU(), srf1.degreeV() };
        bs_type srf2_ { srf2.poles(), srf2.knotsFlatsU(), srf2.knotsFlatsV(), srf2.degreeU(), srf2.degreeV() };

        while(srf1_.degreeU()<srf2_.degreeU()) srf1_.increaseDegreeU();
        while(srf2_.degreeU()<srf1_.degreeU()) srf2_.increaseDegreeU();

        auto nkv = srf2_.knotsFlatsV().size(); // TODO insert v knots
        if(nkv!=srf1_.knotsFlatsV().size())
            throw std::exception("Not implemented");
        for(size_t i{}; i < nkv ; i ++)
        {
            if(std::fabs( srf2_.knotsFlatsV()[i] - srf1_.knotsFlatsV()[i]) > knot_eps)
                throw std::exception("Not implemented");
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

        points_vector<T,dim> poles_new;
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

    template <typename T, size_t dim, bool rational>
    auto extended_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto srf_ext = extention_to_curve(srf,crv,natural_end,max_cont);
        return join(srf,srf_ext);
    }

    template <typename T, size_t dim, bool rational>
    auto extention(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto np = srf.nPolesV()*10;
        auto [u1,u2,v1,v2] = srf.bounds();
        points_vector<double,3> pts(np);
        std::vector<double> v_crv(np);
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
        auto crv = approx_bound_fixed(pts,q,nv,v_crv,kv);
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
        typedef std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type bs_type; 
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
    
}