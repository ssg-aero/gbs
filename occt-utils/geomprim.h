#pragma once
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_XYZ.hxx>
#include <gp_XY.hxx>
#include <gp_Vec.hxx>
#include <gp_Vec2d.hxx>
#include <vector>
#include <array>

namespace occt_utils
{
    inline auto coord(const std::array<double, 3> &c) -> gp_XYZ
    {
        return gp_XYZ(c[0], c[1], c[2]);
    }
    inline auto coord(const std::array<double, 2> &c) -> gp_XY
    {
        return gp_XY(c[0], c[1]);
    }
    inline auto point(const std::array<double, 3> &c) -> gp_Pnt
    {
        return gp_Pnt(c[0], c[1], c[2]);
    }
    inline auto point(const std::array<double, 2> &c) -> gp_Pnt2d
    {
        return gp_Pnt2d(c[0], c[1]);
    }
    inline auto vector(const std::array<double, 3> &c) -> gp_Vec
    {
        return gp_Vec(c[0], c[1], c[2]);
    }
    inline auto vector(const std::array<double, 2> &c) -> gp_Vec2d
    {
        return gp_Vec2d(c[0], c[1]);
    }
    template <class T> inline int dimension()
    {
        T p;
        return dimension(decltype(p.Coord()));
    }

    template <>
    inline int dimension<gp_XY>()
    {
        return 2;
    }

    template <>
    inline int dimension<gp_XYZ>()
    {
        return 3;
    }
} // namespace occt_utils