#include <gtest/gtest.h>
#include <gbs/occt/intersections.h>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Line.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_Plane.hxx>
using namespace occt_utils;
TEST(tests_intersections, nearest_CS)
{
    GeomAdaptor_Curve   lin(new Geom_Line({0,0,0},{0,0,1}));
    GeomAdaptor_Surface pln(new Geom_Plane({0,0,0},{0,0,1}));
    auto res = nearest_CS(lin,pln);
    ASSERT_DOUBLE_EQ(res.first.Parameter(),0.);
    ASSERT_DOUBLE_EQ(res.first.Value().Distance({0,0,0}),0.);
    ASSERT_DOUBLE_EQ(res.second.Value().Distance({0,0,0}),0.);
}

TEST(tests_intersections, nearest_PC)
{
    GeomAdaptor_Curve   lin(new Geom_Line({0,0,0},{0,0,1}));
    auto res = nearest_PC(gp_Pnt(),lin);
    ASSERT_DOUBLE_EQ(res.second.Parameter(),0.);
    ASSERT_DOUBLE_EQ(res.second.Value().Distance({0,0,0}),0.);
    res = nearest_PC(gp_Pnt(1,0,0),lin);
    ASSERT_DOUBLE_EQ(res.second.Parameter(),0.);
    ASSERT_DOUBLE_EQ(res.second.Value().Distance({0,0,0}),0.);

}