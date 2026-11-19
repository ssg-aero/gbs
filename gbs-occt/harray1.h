#pragma once
/// Handled arrays of OCCT under one name for OCCT 7 and 8.
///
/// OCCT 8 deprecates the TColgp_, TColStd_, TopTools_ and AppParCurves_ typedefs for the
/// NCollection templates they name. The arrays, lists and maps were already these templates in
/// OCCT 7 and gbs-occt writes them as such (NCollection_Array1<gp_Pnt>...). The handled arrays
/// were classes of their own in OCCT 7 (DEFINE_HARRAY1): HArray1<T> names NCollection_HArray1<T>
/// with OCCT 8 and the OCCT 7 class otherwise.
#include <Standard_Version.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <AppParCurves_ConstraintCouple.hxx>
#if OCC_VERSION_MAJOR >= 8
#include <NCollection_HArray1.hxx>
#else
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_HArray1OfBoolean.hxx>
#include <AppParCurves_HArray1OfConstraintCouple.hxx>
#endif

namespace occt_utils
{
#if OCC_VERSION_MAJOR >= 8
    template <class T>
    using HArray1 = NCollection_HArray1<T>;
#else
    namespace detail
    {
        template <class T>
        struct harray1; // the element types gbs-occt uses
        template <>
        struct harray1<gp_Pnt> { using type = TColgp_HArray1OfPnt; };
        template <>
        struct harray1<gp_Pnt2d> { using type = TColgp_HArray1OfPnt2d; };
        template <>
        struct harray1<double> { using type = TColStd_HArray1OfReal; };
        template <>
        struct harray1<bool> { using type = TColStd_HArray1OfBoolean; };
        template <>
        struct harray1<AppParCurves_ConstraintCouple> { using type = AppParCurves_HArray1OfConstraintCouple; };
    } // namespace detail
    template <class T>
    using HArray1 = typename detail::harray1<T>::type;
#endif
} // namespace occt_utils
