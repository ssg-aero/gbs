#pragma once
#include <gbs-occt/containers.h>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
namespace occt_utils
{

    inline auto to_edge(const Handle(Geom2d_Curve) & c, const gp_Pln &pln = gp_Pln()) -> TopoDS_Edge
    {
        return BRepBuilderAPI_MakeEdge(GeomAPI::To3d(c, pln)).Edge();
    }

    inline auto to_edge(const Handle(Geom_Curve) & c) -> TopoDS_Edge
    {
        return BRepBuilderAPI_MakeEdge(c).Edge();
    }

    template <class _InIt>
    auto to_sh_list(const _InIt &_First, const _InIt &_Last, TopTools_ListOfShape &L) -> void
    {
        std::for_each(_First, _Last, [&](const auto & c) {
            L.Append(to_edge(c));
        });
    }
    
    template <class _InIt>
    auto to_sh_list(const _InIt &_First, const _InIt &_Last) -> TopTools_ListOfShape
    {
        TopTools_ListOfShape L;
        to_sh_list(_First,_Last,L);
        return L;
    }

    template <class _InIt>
    auto to_wire(const _InIt &_First, const _InIt &_Last) -> TopoDS_Wire
    {
        TopTools_ListOfShape L;
        to_sh_list(_First,_Last, L);
        BRepBuilderAPI_MakeWire builder;
        builder.Add(L);
        return builder.Wire();
    }

    template <typename container>
    inline auto to_wire(const container &curves)
    {
        return to_wire(curves.begin(),curves.end());
    }

    template <typename container>
    inline auto to_face(const container &curves)
    {
        return BRepBuilderAPI_MakeFace(to_wire(curves.begin(),curves.end()), true).Face();
    }

    template <class _InIt>
    inline auto to_shell(const _InIt &_First, const _InIt &_Last,double tol) -> TopoDS_Shell
    {
        BRepBuilderAPI_Sewing Sewing(tol);
        std::for_each(_First,_Last,[&](const auto &sh){Sewing.Add(sh);});

    }

    template <typename T>
    inline TopAbs_ShapeEnum get_type()
    {
        return TopAbs_SHAPE;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Edge>()
    {
        return TopAbs_EDGE;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Vertex>()
    {
        return TopAbs_VERTEX;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Wire>()
    {
        return TopAbs_WIRE;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Face>()
    {
        return TopAbs_FACE;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Shell>()
    {
        return TopAbs_SHELL;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Solid>()
    {
        return TopAbs_SOLID;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_Compound>()
    {
        return TopAbs_COMPOUND;
    }

    template <>
    inline TopAbs_ShapeEnum get_type<TopoDS_CompSolid>()
    {
        return TopAbs_COMPSOLID;
    }

    inline auto explode(const TopoDS_Shape &sh, std::list<TopoDS_Shape> &sh_lst, TopAbs_ShapeEnum sh_type) -> void
    {
        TopTools_IndexedMapOfShape IndexedMapOfShape;
        TopExp::MapShapes(sh,sh_type,IndexedMapOfShape);
        std::for_each(IndexedMapOfShape.cbegin(),IndexedMapOfShape.cend(),
        [&](const TopoDS_Shape &sh)
        {
            sh_lst.push_back(sh);
        });
    }
    inline auto explode(const TopoDS_Shape &sh, TopAbs_ShapeEnum sh_type) -> std::list<TopoDS_Shape>
    {
        std::list<TopoDS_Shape> sh_lst;
        explode(sh, sh_lst, sh_type);
        return sh_lst;
    }

    template <typename T>
    auto explode(const TopoDS_Shape &sh, std::list<T> &sh_lst) -> void
    {
        explode(sh,sh_lst,get_type<T>());
    }

    template <typename T>
    auto explode(const TopoDS_Shape &sh) -> std::list<T>
    {
        std::list<T> sh_lst;
        explode<T>(sh, sh_lst);
        return sh_lst;
    }

} // namespace occt_utils