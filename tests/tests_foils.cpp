#include <gtest/gtest.h>
#include <gbslib/bscbuild.h>
#include <gbslib/vecop.h>
#include <gbslib/bscinterp.h>

#include <occt-utils/curvesbuild.h>
#include <occt-utils/export.h>

#include <Eigen/Dense>

#include <GeomTools.hxx>

using gbs::operator-;
using gbs::operator+;
using gbs::operator*;
TEST(tests_foils, type1)
{
    auto r_ba = 5.;
    auto e_ba = 5.;

    auto r_bf = 0.25;
    auto e_bf = 1.;

    auto c = 25.;

    auto arc_ba = gbs::build_ellipse<double,3>(r_ba ,r_ba/ e_ba, {r_ba, 0., 0.});
    arc_ba.trim(2.5/8.,1./2.);
    arc_ba.reverse();

    auto arc_bf = gbs::build_ellipse<double,3>(r_bf ,r_bf/ e_bf, {c - r_bf, 0., 0.});
    arc_bf.trim(0.,0.5/8.);
    arc_bf.reverse();

    size_t p = 3;
    std::vector<double> u = {0.,1.};
    // std::vector<double> k = {0., 0., 0.,1./4.,1/2.,3./4, 1. ,1. ,1.};
    // std::vector<double> k = {0., 0., 0.,0.1,1/2.,0.9, 1. ,1. ,1.};
    std::vector<double> k = {0., 0., 0.,0., 0.33, 0.66 ,1., 1. ,1. ,1.};
    // std::vector<double> k = gbs::build_simple_mult_flat_knots<double>(u, 8, 5);
    auto pt0 = arc_ba.endRational();
    auto pt3 = arc_bf.beginRational();
    auto tg0 = arc_ba.endRational(1);
    auto tg3 = arc_bf.beginRational(1);
    auto cu0 = arc_ba.endRational(2);
    auto cu3 = arc_bf.beginRational(2);
    auto to0 = arc_ba.endRational(3);
    auto to3 = arc_bf.beginRational(3);

    std::vector<gbs::constrType<double, 3, 3>> Q =
        {
            {{pt0, tg0, cu0/*, to0*/}},
            {{pt3, tg3, cu3/*, to3*/}},
        };

    auto poles = gbs::build_poles(Q,k,u,p);


    auto arc = gbs::BSCurve(poles,k,p);

    GeomTools::Dump(occt_utils::BSplineCurve(arc),std::cout);

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(occt_utils::NURBSplineCurve(arc_ba)); 
    // crv_lst.push_back(occt_utils::NURBSplineCurve(arc)); 
    crv_lst.push_back(occt_utils::BSplineCurve(arc)); 
    crv_lst.push_back(occt_utils::NURBSplineCurve(arc_bf));   

    occt_utils::to_iges(crv_lst, "C:/Users/sebastien/workspace2/gbslib/tests/out/foils_type1.igs");

}