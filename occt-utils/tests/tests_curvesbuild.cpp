#include <gtest/gtest.h>
#include <occt-utils/curvesbuild.h>
#include <occt-utils/export.h>
#include <array>
#include <vector>

#include <GeomAPI.hxx>
// using occt_utils;

TEST(tests_curvebuild, bsc_c1)
{
    auto c1 = occt_utils::bscurve_c1<2>(
        {
            {0.,0.},
            {-.1,-.1}
        },
        {0.,5.},
        {1.,0.},
        false,
        1e-6
    );
    auto c2 = occt_utils::bscurve_c1<3>(
        {
            {0.,0.,0.},
            {-.1,-.1,-1.}
        },
        {0.,5.,0.},
        {1.,0.,0.3},
        false,
        1e-6
    );
    // auto c1 = occt_utils::bscurve_c1<4>(//won't compile
    //     {
    //         {0.,0.},
    //         {-.1,-.1}
    //     },
    //     {0.,5.},
    //     {1.,0.},
    //     false,
    //     1e-6
    // );
    // occt_utils::to_iges({c1},"tests/out/c1.igs");
    occt_utils::to_iges(std::vector<Handle(Geom_Geometry)>{GeomAPI::To3d(c1, gp_Pln()),c2},"tests/out/bsc_c1.igs");
}

TEST(tests_curvebuild, bscurve_c2_approx)
{
    // auto c1 = occt_utils::bscurve_c2_approx<2>(
    //     {
    //         {0., 0.},
    //         {-.1, -.1}
    //     },
    //     {
    //         {0., 5.},
    //         {1., 0.}
    //     },
    //     1e-6);
    auto pt = occt_utils::col_from_vec<gp_Pnt2d>({std::array<double, 2>{0., 0.},
                                                  std::array<double, 2>{0.3, 0.1},
                                                  std::array<double, 2>{1., 1.}});
    auto tg = occt_utils::col_from_vec<gp_Vec2d>({std::array<double, 2>{1., 0.},
                                                  std::array<double, 2>{0., 0.},
                                                  std::array<double, 2>{0., 1.}});
    auto c1 = occt_utils::bscurve_c2_approx(pt, tg, 1e-6);
    occt_utils::to_iges(std::vector<Handle(Geom_Geometry)>{GeomAPI::To3d(c1, gp_Pln())}, "tests/out/bscurve_c2_approx.igs");
}