#include <gbs/occt/intersections.h>
#include <Extrema_ExtCS.hxx>
#include <Extrema_ExtPC.hxx>
#include <Adaptor3d_Curve.hxx>
#include <Adaptor3d_Surface.hxx>
#include <algorithm> 
namespace occt_utils
{
    Standard_EXPORT auto etrema_CS(const Adaptor3d_Curve & C, const Adaptor3d_Surface &S,double tol) -> std::list<res_CS >
    {
        Extrema_ExtCS ExtCS{C,S, tol, tol};
        std::list<res_CS > results_lst;
        size_t count = ExtCS.NbExt();
        Extrema_POnCurv P1;
        Extrema_POnSurf P2;
        for(size_t i = 1 ; i <= count ;i++)
        {
            ExtCS.Points(i, P1, P2);
            results_lst.push_back(res_CS{P1,P2});
        }
        return results_lst;
    }

    auto nearest_CS(const Adaptor3d_Curve & C, const Adaptor3d_Surface &S,double tol) -> res_CS
    {
        
        auto lst = etrema_CS(C,S,tol);
        auto compare = [&](const res_CS &res1,const res_CS &res2)
        {
            auto d1 = res1.first.Value().Distance(res1.second.Value());
            auto d2 = res2.first.Value().Distance(res2.second.Value());
            return (d1 < d2);
        };
        auto nearest = std::min_element(lst.begin(),lst.end(),compare);
        return (*nearest);
    }

    auto etrema_PC(const gp_Pnt &pt, const Adaptor3d_Curve & C,double tol)-> std::list< res_PC >
    {
        Extrema_ExtPC ExtPC{pt,C, tol};
        std::list<res_PC > results_lst;
        size_t count = ExtPC.NbExt();
        for(size_t i = 1 ; i <= count ;i++)
        {
            auto P = ExtPC.Point(i);
            results_lst.push_back(res_PC{pt,P});
        }
        return results_lst;
    }

    auto nearest_PC(const gp_Pnt &pt, const Adaptor3d_Curve & C,double tol)  -> res_PC
    {
        auto lst = etrema_PC(pt,C,tol);
        auto compare = [&](const res_PC &res1,const res_PC &res2)
        {
            auto d1 = res1.first.Distance(res1.second.Value());
            auto d2 = res2.first.Distance(res2.second.Value());
            return (d1 < d2);
        };
        auto nearest = std::min_element(lst.begin(),lst.end(),compare);
        return (*nearest);
    }
}