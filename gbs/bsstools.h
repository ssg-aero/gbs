#pragma once
#include <gbs/bssurf.h>

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
}