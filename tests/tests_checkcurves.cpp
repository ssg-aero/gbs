#include <gtest/gtest.h>
#include <gbs/curves>
#include <gbs/curvescheck.h>


TEST(tests_checkcurves, curves_bounds)
{
    gbs::BSCurve<float,2> c1{
        {
            {0.,0.},
            {1.,0.},
        },
        {0., 0., 1., 1.},
        1
    };
    gbs::BSCurve<float,2> c2{
        {
            {0.,0.},
            {1.,0.},
        },
        {0., 0., 1.01, 1.01},
        1
    };
    gbs::BSCurve<float,2> c3{
        {
            {0.,0.},
            {1.2,0.},
        },
        {0., 0., 1., 1.},
        1
    };
    auto p_c1 = std::make_shared<gbs::BSCurve<float,2>>( c1 );
    auto p_c3 = std::make_shared<gbs::BSCurve<float,2>>( c3 );

    std::vector<gbs::BSCurve<float,2>*> v1{&c1,&c2};
    std::vector<gbs::BSCurve<float,2>*> v2{&c1,&c3};
    std::vector<gbs::BSCurve<float,2>> v3{c1,c2};
    std::vector<gbs::BSCurve<float,2>> v4{c1,c3};
    std::vector<std::shared_ptr<gbs::BSCurve<float,2>>> v5{p_c1, p_c3};

    bool result;
    result = gbs::check_p_curves_bounds<float, 2>(v1.begin(), v1.end());
    ASSERT_FALSE(result);
    result = gbs::check_p_curves_bounds<float,2>(v2.begin(),v2.end());
    ASSERT_TRUE(result);
    result = gbs::check_curves_bounds<float,2>(v3.begin(),v3.end());
    ASSERT_FALSE(result);
    result = gbs::check_curves_bounds<float,2>(v4.begin(),v4.end());
    ASSERT_TRUE(result);
    result = gbs::check_p_curves_bounds<float,2>(v5.begin(),v5.end());
    ASSERT_TRUE(result);
}