#include <gbs/occt/measure.h>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopoDS_Shape.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <CPnts_AbscissaPoint.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
namespace occt_utils
{
    auto area(const TopoDS_Shape &sh) -> double
    {
        GProp_GProps SProps;
        BRepGProp::SurfaceProperties(sh, SProps,false,true);
        return SProps.Mass();
    }
    auto length(const TopoDS_Shape &sh) -> double
    {
        switch (sh.ShapeType())
        {
        case TopAbs_EDGE:
            return CPnts_AbscissaPoint::Length(BRepAdaptor_Curve(TopoDS::Edge(sh)));
            break;
        case TopAbs_WIRE:
            return CPnts_AbscissaPoint::Length(BRepAdaptor_CompCurve(TopoDS::Wire(sh)));
            break;
        default:
            Standard_Failure::Raise("Not implemented.");
            break;
        }
        return -1.;
    }

    auto length(const Adaptor3d_Curve &c) -> double
    {
        return CPnts_AbscissaPoint::Length(c);
    }
    auto length(const Adaptor2d_Curve2d &c) -> double
    {
        return CPnts_AbscissaPoint::Length(c);
    }
    auto cg_lineical(const TopoDS_Shape &sh) -> gp_Pnt
    {
        GProp_GProps LProps;
        BRepGProp::LinearProperties(sh, LProps,false,true);
        return LProps.CentreOfMass();
    }
} // namespace occt_utils