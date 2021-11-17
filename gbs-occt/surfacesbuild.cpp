#include <gbs-occt/surfacesbuild.h>
#include <gbs-occt/curvesbuild.h>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Plane.hxx>
#include <GeomFill_NSections.hxx>
#include <gbs-occt/containers.h>
#include <Geom_BSplineSurface.hxx>
#include <exception>

namespace occt_utils
{
    auto half_planeX(const gp_Ax3 &ax, bool x_positive) -> Handle(Geom_Surface)
    {
        gp_Ax3 zx_plane{ax.Location(), -ax.YDirection(), ax.XDirection()};
        Handle(Geom_Plane) pln = new Geom_Plane(zx_plane);
        double U1, U2, V1, V2;
        pln->Bounds(U1, U2, V1, V2);
        if (x_positive)
            return new Geom_RectangularTrimmedSurface(pln, 0., U2, V1, V2);
        else
            return new Geom_RectangularTrimmedSurface(pln, U1, 0., V1, V2);
    }

    auto half_planeY(const gp_Ax3 &ax, bool y_positive) -> Handle(Geom_Surface)
    {
        gp_Ax3 zx_plane{ax.Location(), ax.XDirection(), ax.YDirection()};
        Handle(Geom_Plane) pln = new Geom_Plane(zx_plane);
        double U1, U2, V1, V2;
        pln->Bounds(U1, U2, V1, V2);
        if (y_positive)
            return new Geom_RectangularTrimmedSurface(pln, U1, U2, 0., V2);
        else
            return new Geom_RectangularTrimmedSurface(pln, U1, U2, V1, 0.);
    }

    auto loft_occ(const std::list<Handle(Geom_Curve)> &lst) -> Handle(Geom_Surface)
    {
        auto S = col_from_list(to_bs_list(lst));//occt is inconsitent in behaviour, non bs curves are skiped
        return GeomFill_NSections(S).BSplineSurface();
    }

    auto loft_occ(const std::vector<Handle(Geom_Curve)> &lst) -> Handle(Geom_Surface)
    {
        return GeomFill_NSections(col_from_list(lst)).BSplineSurface();
    }
    auto loft_occ(const std::vector<std::list<Handle(Geom_Curve)>> &lst) -> std::list<Handle(Geom_Surface)>
    {
        std::list<Handle(Geom_Surface)> s_lst;
        auto nc_section = lst.front().size();
        for (auto i = 0; i < nc_section; i++)
        {
            std::list<Handle(Geom_Curve)> lst_;
            std::for_each(lst.begin(), lst.end(), [&](const auto cs) {
                if (nc_section != cs.size())
                {
                    throw std::length_error("All list must hace same length");
                }
                auto it(cs.begin());
                std::advance(it, i);
                lst_.push_back(*it);
            });
            s_lst.push_back(occt_utils::loft_occ(lst_));
        }
        return s_lst;
    }
} // namespace occt_utils