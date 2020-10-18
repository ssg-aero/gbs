#pragma once
#include <Geom_Surface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <vector>
#include <list>
#include <algorithm>
#include <execution>
#include <gbs-occt/containers.h>
#include <gbs/bssurf.h>
#include <gbs/knotsfunctions.h>

class gp_Ax3;
namespace occt_utils
{

    auto inline BSplineSurface(  const std::vector<std::array<double,3> > &p,
                        const std::vector<double> &ku,
                        const std::vector<double> &kv,
                        const std::vector<int> &mu,
                        const std::vector<int> &mv,
                        int degU,int degV) -> Handle(Geom_BSplineSurface)
    { 
        auto nPolesV = std::reduce(std::execution::par, mv.cbegin(), mv.cend()) - degV - 1;
        auto poles_arr = occt_utils::col2_from_vec<gp_Pnt>(p,nPolesV);
	return new Geom_BSplineSurface(
		poles_arr,
		occt_utils::col_from_vec(ku),
        occt_utils::col_from_vec(kv),
		occt_utils::col_from_vec(mu),
        occt_utils::col_from_vec(mv),
        degU,degV);
	
    }

    auto inline BSplineSurface(const gbs::BSSurface<double,3> &srf){
        std::vector<double> ku;
        std::vector<double> kv;
        std::vector<int>    mu;
        std::vector<int>    mv;
        gbs::unflat_knots(srf.knotsFlatsU(),mu,ku);
        gbs::unflat_knots(srf.knotsFlatsV(),mv,kv);
        return BSplineSurface(srf.poles(),ku,kv,mu,mv,srf.degreeU(),srf.degreeV());
    }

    Standard_EXPORT auto half_planeX(const gp_Ax3 &ax, bool x_positive) -> Handle(Geom_Surface);
    Standard_EXPORT auto half_planeY(const gp_Ax3 &ax, bool y_positive) -> Handle(Geom_Surface);
    Standard_EXPORT auto loft_occ(const std::list<Handle(Geom_Curve)> &lst) -> Handle(Geom_Surface);
    Standard_EXPORT auto loft_occ(const std::vector<Handle(Geom_Curve)> &lst) -> Handle(Geom_Surface);
    Standard_EXPORT auto loft_occ(const std::vector<std::list<Handle(Geom_Curve)> > &lst) -> std::list<Handle(Geom_Surface)>;

} // namespace occt_utils