#include <gtest/gtest.h>
#include <gbs-io/fromjson.h>

using namespace gbs;

TEST(tests_io, json_bsc)
{

   rapidjson::Document document;
   parse_file("C:/Users/sebastien/workspace/gbslib/data/crv_bs.json",document);


   auto curve1d = make_bscurve<double,1>(document["1d_entities"][0]);
   ASSERT_EQ(curve1d.degree(), 3);
   auto curve2d = make_bscurve<double,2>(document["2d_entities"][0]);
   ASSERT_EQ(curve2d.poles().size(), 4);
   ASSERT_EQ(curve2d.knotsFlats().size(), 8);
   auto curvei2d= bscurve_interp_cn<double,2>(document["2d_entities"][1]);
   ASSERT_NEAR(curvei2d.value(0.5)[0],0.0,1e-6);
   ASSERT_NEAR(curvei2d.value(0.5)[1],0.5,1e-6);
   auto curve3d = make_bscurve<double,3>(document["3d_entities"][0]);
   ASSERT_NEAR(curve3d.poles()[2][2],1.0,1e-6);


}
