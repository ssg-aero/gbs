#include <gtest/gtest.h>
#include <gbs-occt/geomprim.h>
#include <gbs-occt/containers.h>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>

// using occt_utils;

TEST(tests_geomprim, point)
{
    auto pt1 = occt_utils::point({1., 2., 3});
    ASSERT_DOUBLE_EQ(pt1.X(), 1.);
    ASSERT_DOUBLE_EQ(pt1.Y(), 2.);
    ASSERT_DOUBLE_EQ(pt1.Z(), 3.);
}

template <size_t d>
inline auto pt(const std::array<double, d> &coord)
{
    std::conditional<d == 2, gp_Pnt2d, gp_Pnt>::type p;
    for (int i = 1; i <= d; i++)
        p.SetCoord(i, coord[i - 1]);
    return p;
}

template <size_t d>
inline auto pt_vec(const std::vector< std::array<double, d> > &coord)
{
    std::conditional<d == 2, gp_Pnt2d, gp_Pnt>::type p;

    // for (int i = 1; i <= d; i++)
    //     p.SetCoord(i, coord[i - 1]);
    // return p;
}

TEST(tests_geomprim, template_dim)
{
    auto pt2d = pt<2>({1, 2});
    ASSERT_DOUBLE_EQ(pt2d.X(), 1.);
    ASSERT_DOUBLE_EQ(pt2d.Y(), 2.);

    auto pt3d = pt<3>({1, 2, 3});
    ASSERT_DOUBLE_EQ(pt3d.X(), 1.);
    ASSERT_DOUBLE_EQ(pt3d.Y(), 2.);
    ASSERT_DOUBLE_EQ(pt3d.Z(), 3.);

    //     auto pt3d_err = pt({1,2,3}); // won't compile
    auto pt3d_ok = pt(std::array<double, 3>({1, 2., 3.}));
    ASSERT_DOUBLE_EQ(pt3d_ok.X(), 1.);
    ASSERT_DOUBLE_EQ(pt3d_ok.Y(), 2.);
    ASSERT_DOUBLE_EQ(pt3d_ok.Z(), 3.);

    // auto pt2d_vec = occt_utils::col_from_vec< std::conditional<true, gp_Pnt2d, gp_Pnt>::type >(
    // {
    //     {1.,2.},
    //     {3.,4.},
    // }
    // );

    // std::conditional<d == 2,gp_Pnt2d,gp_Pnt> pt_type;
    //     return bscurve_c1(
    //         col_from_vec<pt_type>(pt_lst),
}