#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/containers.h>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TColgp_HArray1OfVec.hxx>
#include <TColgp_HArray1OfVec2d.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColStd_HArray1OfBoolean.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepAdaptor_HCompCurve.hxx>
#include <Approx_Curve3d.hxx>
#include <GeomConvert.hxx>

#include <AppDef_MultiLine.hxx>
#include <AppParCurves_HArray1OfConstraintCouple.hxx>
#include <AppDef_Variational.hxx>
#include <AppParCurves_MultiBSpCurve.hxx>

namespace occt_utils
{
    auto to_bs_list(const std::list<Handle(Geom_Curve)> &lst) -> std::list<Handle(Geom_Curve)>
    {
        std::list<Handle(Geom_Curve)> bs_lst;
        std::transform(lst.begin(), lst.end(), std::back_inserter(bs_lst), [&](const Handle(Geom_Curve) & c) {
            auto bs = Handle(Geom_BSplineCurve)::DownCast(c);
            if (bs.IsNull())
            {
                return GeomConvert::CurveToBSplineCurve(c);
            }
            else
            {
                return bs;
            }
        });
        return bs_lst;
    }

    auto bscurve_c1( const TColgp_Array1OfPnt2d &pt_lst, const gp_Vec2d &t1,const gp_Vec2d &t2,bool scale_tg,double tol) -> Handle(Geom2d_BSplineCurve)
    {
        Handle( TColgp_HArray1OfPnt2d ) h_points =new TColgp_HArray1OfPnt2d(pt_lst);
        Geom2dAPI_Interpolate Interpolate(h_points,false,tol);
        Interpolate.Load(t1,t2,scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();
    }

    auto bscurve_c1( const TColgp_Array1OfPnt &pt_lst, const gp_Vec &t1,const gp_Vec &t2,bool scale_tg,double tol) -> Handle(Geom_BSplineCurve)
    {
        Handle( TColgp_HArray1OfPnt ) h_points =new TColgp_HArray1OfPnt(pt_lst);
        GeomAPI_Interpolate Interpolate(h_points,false,tol);
        Interpolate.Load(t1,t2,scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();
    }

    auto bscurve_c1(const TColgp_Array1OfPnt2d &pt_lst, const TColgp_Array1OfVec2d &tg_lst, bool scale_tg, double tol) -> Handle(Geom2d_BSplineCurve)
    {
        Handle( TColStd_HArray1OfBoolean ) flags = new TColStd_HArray1OfBoolean(1,tg_lst.Length());
        std::transform(tg_lst.begin(),tg_lst.end(),flags->begin(),[&](const gp_Vec2d &v){return v.Magnitude()>tol;});

        Geom2dAPI_Interpolate Interpolate(new TColgp_HArray1OfPnt2d(pt_lst),false,tol);
        Interpolate.Load(tg_lst,flags,scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();  
    }

    auto bscurve_c1(const TColgp_Array1OfPnt &pt_lst, const TColgp_Array1OfVec &tg_lst, bool scale_tg, double tol) -> Handle(Geom_BSplineCurve)
    {
        Handle( TColStd_HArray1OfBoolean ) flags = new TColStd_HArray1OfBoolean(1,tg_lst.Length());
        std::transform(tg_lst.begin(),tg_lst.end(),flags->begin(),[&](const gp_Vec &v){return v.Magnitude()>tol;});

        GeomAPI_Interpolate Interpolate(new TColgp_HArray1OfPnt(pt_lst),false,tol);
        Interpolate.Load(tg_lst,flags,scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();  
    }

    auto bscurve_c1(const std::vector<std::array<double, 2>> &pt_lst, const std::vector<std::array<double, 2> > &tg_lst, bool scale_tg, double tol) -> Handle(Geom2d_BSplineCurve)
    {
        return bscurve_c1(
            col_from_vec<gp_Pnt2d>(pt_lst),
            col_from_vec<gp_Vec2d>(tg_lst),
            scale_tg,
            tol
        );
    }

    auto bscurve_c1( const std::vector< std::array<double,2> > &pt_lst, const std::vector< double > &pr_lst, const std::array<double,2> &t1,const std::array<double,2> &t2,bool scale_tg,double tol) -> Handle(Geom2d_BSplineCurve)
    {
        Handle( TColgp_HArray1OfPnt2d ) points = new TColgp_HArray1OfPnt2d( col_from_vec<gp_Pnt2d>(pt_lst) );
        Handle(TColStd_HArray1OfReal) params = new TColStd_HArray1OfReal(col_from_vec(pr_lst) );
        Geom2dAPI_Interpolate Interpolate(points,params,false,tol);
        Interpolate.Load(vector(t1),vector(t2),scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();    
    }

    auto bscurve_c1(const std::vector< std::array<double, 2> > &pt_lst,double tol) -> Handle(Geom2d_BSplineCurve)
    {
        Handle( TColgp_HArray1OfPnt2d ) points = new TColgp_HArray1OfPnt2d(col_from_vec<gp_Pnt2d>(pt_lst));
        Geom2dAPI_Interpolate Interpolate(points,false,tol);
        Interpolate.Perform();
        return Interpolate.Curve();
    }

    auto bscurve_c1(const Handle(Geom2d_Curve) &back,const Handle(Geom2d_Curve) &front,bool scale_tg,double tol) -> Handle(Geom2d_BSplineCurve)
    {
        auto u1 = back->LastParameter();
        auto u2 = front->FirstParameter();
        auto t1 = back->DN(u1,1);
        auto t2 = front->DN(u2, 1);
        Handle(TColgp_HArray1OfPnt2d) points = new TColgp_HArray1OfPnt2d(col_from_vec<gp_Pnt2d>({back->Value(u1), front->Value(u2)}) );
        Geom2dAPI_Interpolate Interpolate(points, false, tol);
        Interpolate.Load(t1, t2, scale_tg);
        Interpolate.Perform();
        return Interpolate.Curve();
    }

template <typename _InIt>
    auto bs_from_list(const _InIt &_First,const _InIt &_Last,double tol) ->  Handle(Geom_BSplineCurve)
    {
        
        auto S = to_sh_list(_First,_Last);
        BRepBuilderAPI_MakeWire MakeWire;
        MakeWire.Add(S);
        Handle(BRepAdaptor_HCompCurve) h_ad = new BRepAdaptor_HCompCurve(MakeWire.Wire());
        int nSeg = std::distance(_First,_Last) *20;
        int deg = 5;
        return Approx_Curve3d(h_ad,tol,GeomAbs_C2,nSeg,deg).Curve();

    }
    auto bs_from_list(const std::list<Handle(Geom_Curve)> &lst,double tol) ->  Handle(Geom_BSplineCurve)
    {
        return bs_from_list(lst.begin(),lst.end(),tol);
    }

    auto bscurve_c2_approx( const TColgp_Array1OfPnt2d &pt_lst, const TColgp_Array1OfVec2d &tg_lst,double tol) -> Handle(Geom2d_BSplineCurve)
    {

    auto nb_csrt = Standard_Integer(pt_lst.Length());
    AppDef_MultiLine SSP(nb_csrt);

    for (int i = 1; i <= nb_csrt; i++)
    {
        if (tg_lst(i).Magnitude() > tol)
        {
            TColgp_Array1OfPnt2d tabPnt(pt_lst(i), 1, 1);
            TColgp_Array1OfVec2d tabTg(tg_lst(i), 1, 1);
            AppDef_MultiPointConstraint MultiPointConstraint(tabPnt, tabTg);
            SSP.SetValue(i, MultiPointConstraint);
        }
        else
        {
            TColgp_Array1OfPnt2d tabPnt(pt_lst(i), 1, 1);
            AppDef_MultiPointConstraint MultiPointConstraint(tabPnt);
            SSP.SetValue(i, MultiPointConstraint);
        }
    }
//Construction du paramètreage
    bool PeriodicFlag = false;
    Handle(TColStd_HArray1OfReal) ParametersPtr;
    double distance;
    int num_parameters = nb_csrt;

    if (PeriodicFlag)
        num_parameters++;
    ParametersPtr = new TColStd_HArray1OfReal(1, num_parameters);
    ParametersPtr->SetValue(1, 0.0e0);

    for (int ii = pt_lst.Lower(); ii < pt_lst.Upper(); ii++)
    {
        distance = pt_lst.Value(ii).Distance(pt_lst.Value(ii + 1));
        ParametersPtr->SetValue(ii + 1, ParametersPtr->Value(ii) + distance);
    }
    if (PeriodicFlag)
    {
        distance = pt_lst.Value(pt_lst.Upper()).Distance(pt_lst.Value(pt_lst.Lower()));
        ParametersPtr->SetValue(pt_lst.Upper() + 1, ParametersPtr->Value(pt_lst.Upper()) + distance);
    }


    Handle(AppParCurves_HArray1OfConstraintCouple) TheConstraints =
        new AppParCurves_HArray1OfConstraintCouple(1, Standard_Integer( pt_lst.Length() ) );

    int index = 1;
    std::for_each(pt_lst.begin(), pt_lst.end(), [&](const auto &pt) {

        AppParCurves_Constraint Constraint = tg_lst(index).Magnitude()>tol ? AppParCurves_TangencyPoint : AppParCurves_PassPoint;
        
        AppParCurves_ConstraintCouple ACC(index, Constraint);
        TheConstraints->SetValue(index, ACC);
        index++;
    });

    AppDef_Variational Variational(SSP,
                                   1,                               //FirstPoint
                                   Standard_Integer(nb_csrt), //LastPoint
                                   TheConstraints,
                                   5,                               // MaxDegree=14
                                   Standard_Integer(nb_csrt), //MaxSegment=100,
                                // 100,
                                   GeomAbs_C2,                      //Continuity=GeomAbs_C2,
                                   Standard_False,                  ////WithMinMax=Standard_False,
                                   Standard_True,                   //WithCutting=Standard_True,
                                   1e-3,                            //Tolerance=1.0,
                                   2                               // NbIterations=2
    );

    Variational.SetParameters(ParametersPtr);
    Variational.SetCriteriumWeight(1., 1., 1.);

    Variational.Approximate();
    AppParCurves_MultiBSpCurve TheCurve = Variational.Value();

    TColgp_Array1OfPnt2d Poles(1, TheCurve.NbPoles());

    TheCurve.Curve(1, Poles);

    return new Geom2d_BSplineCurve(Poles,
                                   TheCurve.Knots(),
                                   TheCurve.Multiplicities(),
                                   TheCurve.Degree());
    }

    //TODO: remove duplicat code
    auto bscurve_c2_approx( const TColgp_Array1OfPnt &pt_lst, const TColgp_Array1OfVec &tg_lst,double tol) -> Handle(Geom_BSplineCurve)
    {

    auto nb_csrt = Standard_Integer(pt_lst.Length());
    AppDef_MultiLine SSP(nb_csrt);

    for (int i = 1; i <= nb_csrt; i++)
    {
        if (tg_lst(i).Magnitude() > tol)
        {
            TColgp_Array1OfPnt tabPnt(pt_lst(i), 1, 1);
            TColgp_Array1OfVec tabTg(tg_lst(i), 1, 1);
            AppDef_MultiPointConstraint MultiPointConstraint(tabPnt, tabTg);
            SSP.SetValue(i, MultiPointConstraint);
        }
        else
        {
            TColgp_Array1OfPnt tabPnt(pt_lst(i), 1, 1);
            AppDef_MultiPointConstraint MultiPointConstraint(tabPnt);
            SSP.SetValue(i, MultiPointConstraint);
        }
    }
//Construction du paramètreage
    bool PeriodicFlag = false;
    Handle(TColStd_HArray1OfReal) ParametersPtr;
    double distance;
    int num_parameters = nb_csrt;

    if (PeriodicFlag)
        num_parameters++;
    ParametersPtr = new TColStd_HArray1OfReal(1, num_parameters);
    ParametersPtr->SetValue(1, 0.0e0);

    for (int ii = pt_lst.Lower(); ii < pt_lst.Upper(); ii++)
    {
        distance = pt_lst.Value(ii).Distance(pt_lst.Value(ii + 1));
        ParametersPtr->SetValue(ii + 1, ParametersPtr->Value(ii) + distance);
    }
    if (PeriodicFlag)
    {
        distance = pt_lst.Value(pt_lst.Upper()).Distance(pt_lst.Value(pt_lst.Lower()));
        ParametersPtr->SetValue(pt_lst.Upper() + 1, ParametersPtr->Value(pt_lst.Upper()) + distance);
    }


    Handle(AppParCurves_HArray1OfConstraintCouple) TheConstraints =
        new AppParCurves_HArray1OfConstraintCouple(1, Standard_Integer( pt_lst.Length() ) );

    int index = 1;
    std::for_each(pt_lst.begin(), pt_lst.end(), [&](const auto &pt) {

        AppParCurves_Constraint Constraint = tg_lst(index).Magnitude()>tol ? AppParCurves_TangencyPoint : AppParCurves_PassPoint;
        
        AppParCurves_ConstraintCouple ACC(index, Constraint);
        TheConstraints->SetValue(index, ACC);
        index++;
    });

    AppDef_Variational Variational(SSP,
                                   1,                               //FirstPoint
                                   Standard_Integer(nb_csrt), //LastPoint
                                   TheConstraints,
                                   5,                               // MaxDegree=14
                                   Standard_Integer(nb_csrt), //MaxSegment=100,
                                // 100,
                                   GeomAbs_C2,                      //Continuity=GeomAbs_C2,
                                   Standard_False,                  ////WithMinMax=Standard_False,
                                   Standard_True,                   //WithCutting=Standard_True,
                                   1e-3,                            //Tolerance=1.0,
                                   2                               // NbIterations=2
    );

    Variational.SetParameters(ParametersPtr);
    Variational.SetCriteriumWeight(1., 1., 1.);

    Variational.Approximate();
    AppParCurves_MultiBSpCurve TheCurve = Variational.Value();

    TColgp_Array1OfPnt Poles(1, TheCurve.NbPoles());

    TheCurve.Curve(1, Poles);

    return new Geom_BSplineCurve(Poles,
                                   TheCurve.Knots(),
                                   TheCurve.Multiplicities(),
                                   TheCurve.Degree());
    }

} // namespace occt_utils