#include <gtest/gtest.h>
#include <gbs/gbslib.h>

#include <chrono>

import math;

TEST(tests_math, angle_conversion)
{
    ASSERT_DOUBLE_EQ(gbs::radians(180.),std::acos(-1.));
    ASSERT_DOUBLE_EQ(gbs::degrees(std::acos(-1.)),180.);
}

TEST(tests_math, factorial)
{
    ASSERT_EQ(gbs::factorial(4),4*3*2);
}

TEST(tests_math, binomial_law)
{
    auto C4_2 = gbs::binomial_law<int,size_t>(4,2);
    int  C4_2x= 4*3*2/(2*2); 
    ASSERT_EQ(C4_2,C4_2x);
}

TEST(tests_math, make_range)
{
    size_t n = 5;
    auto r1 =gbs::make_range(0.,1.,size_t{2*n+1});
    ASSERT_DOUBLE_EQ(r1[0],0.);
    ASSERT_DOUBLE_EQ(r1[n],0.5);
    ASSERT_DOUBLE_EQ(r1[2*n],1.);

    auto r2 =gbs::make_range(0.f,1.f,size_t{2*n+1});
    ASSERT_FLOAT_EQ(r2[0],0.f);
    ASSERT_FLOAT_EQ(r2[n],0.5f);
    ASSERT_FLOAT_EQ(r2[2*n],1.f);

    gbs::point<double,1> pt1_1{0.};
    gbs::point<double,1> pt2_1{1.};
    auto r3 =gbs::make_range(pt1_1,pt2_1,size_t{2*n+1});
    ASSERT_DOUBLE_EQ(r3[0][0],0.);
    ASSERT_DOUBLE_EQ(r3[n][0],0.5);
    ASSERT_DOUBLE_EQ(r3[2*n][0],1.);

    gbs::point<float,1> pt1_2{0.};
    gbs::point<float,1> pt2_2{1.};
    auto r4 =gbs::make_range(pt1_2,pt2_2,size_t{2*n+1});
    ASSERT_FLOAT_EQ(r4[0][0],0.);
    ASSERT_FLOAT_EQ(r4[n][0],0.5);
    ASSERT_FLOAT_EQ(r4[2*n][0],1.);


    auto r5 = gbs::make_range(0., 1., 2., 4.,3,2);
    ASSERT_DOUBLE_EQ(r5[0].first,0.0);
    ASSERT_DOUBLE_EQ(r5[1].first,0.0);
    ASSERT_DOUBLE_EQ(r5[0].second,2.0);
    ASSERT_DOUBLE_EQ(r5[1].second,4.0);
    ASSERT_DOUBLE_EQ(r5[2].first,0.5);
    ASSERT_DOUBLE_EQ(r5[3].first,0.5);
    ASSERT_DOUBLE_EQ(r5[2].second,2.0);
    ASSERT_DOUBLE_EQ(r5[3].second,4.0);
    ASSERT_DOUBLE_EQ(r5[4].first,1.0);
    ASSERT_DOUBLE_EQ(r5[5].first,1.0);
    ASSERT_DOUBLE_EQ(r5[4].second,2.0);
    ASSERT_DOUBLE_EQ(r5[5].second,4.0);
}