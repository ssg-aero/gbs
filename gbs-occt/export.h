#pragma once

#include <vector>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <Standard_Handle.hxx>
#include <GeomAPI.hxx>
#include <gp_Pln.hxx>
#include <TopoDS_Shape.hxx>
#include <TopLoc_Location.hxx>
#include <algorithm>

namespace occt_utils
{

    inline auto add(const Handle(Geom2d_Curve) & c, IGESControl_Writer &writer, double scale = 1000) -> bool
    {
        return writer.AddGeom(GeomAPI::To3d(c, gp_Pln())->Scaled(gp_Pnt(), scale));
    }

    inline auto add(const Handle(Geom_Geometry) & c, IGESControl_Writer &writer, double scale = 1000) -> bool
    {
        return writer.AddGeom(c->Scaled(gp_Pnt(), scale));
    }

    inline auto add(const TopoDS_Shape &sh, IGESControl_Writer &writer, double scale = 1000) -> bool
    {
        gp_Trsf Trsf;
        Trsf.SetScaleFactor(scale);
        TopLoc_Location loc(Trsf);
        return writer.AddShape(sh.Moved(loc));
    }

    template <typename _InIt>
    inline auto to_iges(const _InIt &_First, const _InIt &_Last, const char *f_name, double scale = 1000) -> void
    {
        IGESControl_Controller::Init();
        IGESControl_Writer writer("MM", 0);

        auto _add_ = [&](const auto &t) { add(t, writer,scale); };
        std::for_each(_First, _Last, _add_);

        writer.Write(f_name);
    }

    template <typename Container>
    inline auto to_iges(const Container &lst, const char *f_name, double scale = 1000) -> void
    {
        to_iges(lst.cbegin(), lst.cend(), f_name, scale);
    }

    // template <typename T>
    // inline auto to_iges(const std::vector<T> &lst, const char *f_name, double scale = 1000) -> void
    // {
    //     to_iges(lst.cbegin(), lst.cend(), f_name, scale);
    // }

    // template <typename T>
    // inline auto to_iges(const std::list<T> &lst, const char *f_name, double scale = 1000) -> void
    // {
    //     to_iges(lst.cbegin(), lst.cend(), f_name, scale);
    // }

    template <typename T, size_t _Size>
    inline auto to_iges(const std::array<T, _Size> &lst, const char *f_name, double scale = 1000) -> void
    {
        to_iges(lst.cbegin(), lst.cend(), f_name, scale);
    }

} // namespace occt_utils