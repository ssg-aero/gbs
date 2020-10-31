#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bssbuild.h>
#include <gbs-render/vtkcurvesrender.h>

const double tol = 1e-6;

using gbs::operator-;

TEST(tests_bssbuild, loft)
{
    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_3d_d poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {1.,1.,0.},
        {3.,1.,0.},
        {0.,4.,0.},
    };
    gbs::points_vector_3d_d poles2 =
    {
        {0.,0.,1.},
        {0.,1.,1.},
        {1.,1.,1.},
        {1.,1.,1.},
        {1.,1.,1.},
        {3.,1.,1.},
        {2.,4.,1.},
    };
    gbs::points_vector_3d_d poles3 =
    {
        {0.,0.,2.},
        {0.,1.,2.},
        {1.,1.,2.},
        {1.,1.,2.},
        {1.,1.,2.},
        {3.,1.,2.},
        {1.,4.,2.},
    };
    gbs::BSCurve3d_d c1(poles1,k,p);
    gbs::BSCurve3d_d c2(poles2,k,p);
    gbs::BSCurve3d_d c3(poles3,k,p);

    auto s = gbs::loft<double,3,false>(std::list<gbs::BSCurve3d_d>{{c1,c2,c3}});
    gbs::plot(s,c1,c2,c3);
}
