#pragma once
#include <Standard.hxx>
class TopoDS_Shape;
class Adaptor3d_Curve;
class Adaptor2d_Curve2d;
class gp_Pnt;
namespace occt_utils
{
    Standard_EXPORT auto area(const TopoDS_Shape &sh) -> double;
    Standard_EXPORT auto length(const TopoDS_Shape &sh) -> double;
    Standard_EXPORT auto length(const Adaptor3d_Curve &c) -> double;
    Standard_EXPORT auto length(const Adaptor2d_Curve2d &c) -> double;
    Standard_EXPORT auto cg_lineical(const TopoDS_Shape &sh) -> gp_Pnt;
} // namespace occt_utils