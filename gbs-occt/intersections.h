#pragma once

#include <utility>
#include <list>
#include <Extrema_POnCurv.hxx>
#include <Extrema_POnSurf.hxx>
class Adaptor3d_Curve;
class Adaptor3d_Surface;

typedef std::pair<Extrema_POnCurv,Extrema_POnSurf> res_CS;
typedef std::pair<gp_Pnt,Extrema_POnCurv> res_PC;

namespace occt_utils
{
    Standard_EXPORT auto etrema_CS(const Adaptor3d_Curve & C, const Adaptor3d_Surface &S,double tol=1e-6) -> std::list< res_CS >;
    Standard_EXPORT auto nearest_CS(const Adaptor3d_Curve & C, const Adaptor3d_Surface &S,double tol=1e-6) -> res_CS;
    Standard_EXPORT auto etrema_PC(const gp_Pnt &pt, const Adaptor3d_Curve & C,double tol=1e-6) -> std::list< res_PC >;
    Standard_EXPORT auto nearest_PC(const gp_Pnt &pt, const Adaptor3d_Curve & C,double tol=1e-6) -> res_PC;
}