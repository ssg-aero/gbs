#pragma once
#include <gbs/bssurf.h>

namespace gbs
{
    template <typename T, size_t dim,bool rational>
    auto Discretize(const BSSurfaceGeneral<T,dim,rational> &srf, size_t nu, size_t nv) -> gbs::points_vector<T,dim>
    {
        points_vector<T,dim> points(nu*nv);
        auto u1 = srf.knotsFlatsU().front();
        auto u2 = srf.knotsFlatsU().back();
        auto du = (u2-u1) / (nu-1);
        auto v1 = srf.knotsFlatsV().front();
        auto v2 = srf.knotsFlatsV().back();
        auto dv = (v2-v1) / (nv-1);

        auto v = v1;
        for(auto it = points.begin() ; it!= points.end() ; it=std::next(it,nu) )
        {
            std::generate(it,std::next(it,nu),[&,u=u1-du]() mutable {return srf.value(u+=du,v);});
            v+=dv;
        }

        return points;
    }
}