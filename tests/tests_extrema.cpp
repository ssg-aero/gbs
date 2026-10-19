#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include <gbs/vecop.ixx>
#endif
#include <gbs/bssinterp.h>
#include <gbs/bscbuild.h>
#include <gbs/extrema.h>
#include <gbs/bselementary.h>
#include <gbs/bssbuild.h>
#include <gbs-brep/offset.h>

#include <numbers>

#include <gbs-render/vtkGbsRender.h>

using gbs::operator-;



using namespace gbs;



    template <typename... T, std::size_t... I>
    auto subtuple_(const std::tuple<T...>& t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<I>(t)...);
    }

    template <int Trim, typename... T>
    auto subtuple(const std::tuple<T...>& t) {
    return subtuple_(t, std::make_index_sequence<sizeof...(T) - Trim>());
    }

    // template <typename... Targs>
    //     auto bracket(Targs... Fargs)
    //     {
    //         auto tuple = std::tie(Fargs...);
    //         if
    //     }
    


TEST(tests_extrema, PC)
{
    std::vector<gbs::constrType<double, 3, 1>> Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
            {{1., 1., 2.}},
            {{0., 1., 1.}},
            {{0., -1., 4.}},
        };

    gbs::points_vector<double,3> pts;

    auto crv = gbs::interpolate(Q,2,gbs::KnotsCalcMode::CHORD_LENGTH);

    {
        auto u = 0.3;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
    }
    {
        auto u = 0.7;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
    }

}

TEST(tests_extrema, PC_f)
{
    std::vector<gbs::constrType<float, 3, 1>> Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
            {{1., 1., 2.}},
            {{0., 1., 1.}},
            {{0., -1., 4.}},
        };

    gbs::points_vector<float,3> pts;

    auto crv = gbs::interpolate(Q,2,gbs::KnotsCalcMode::CHORD_LENGTH);

    {
        auto u = 0.3;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point<float,3>(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 1e-5);
    }
    {
        auto u = 0.7;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point<float,3>(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 1e-5);
    }

}

TEST(tests_extrema, PS)
{
    std::vector<double> ku_flat = {0.,0.,0.,5.,5.,5.};
    std::vector<double> kv_flat = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;

    //Pij avec j inner loop
    const std::vector<std::array<double,3> > poles =
    {
        {0,0,0},{1,0,1},{2,0,2},
        {0,1,0},{1.2,1,1},{2,1,1},
        {0,2,1},{1.2,2,1},{2,2,2},
        {0,3,1},{1.1,3,1},{2,3,3},
        {0,4,0},{1,4,1},{2,4,0},
    };
    

    gbs::BSSurface<double,3> srf(poles,ku_flat,kv_flat,p,q);
    auto u = 2.5;
    auto v = 1.;

    auto pt = srf.value(u,v);

    // gbs::plot(
    //     srf,
    //     gbs::points_vector<double,3>{srf(res.u,res.v),pt}
    // );


    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt(srf, pt, 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 5e-6);
        ASSERT_NEAR(res_v, v, 5e-6);
    }


    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};

    std::vector<std::array<double,4> > poles_t = {
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,6,8,2}, {12,24,12,6}, {4,6,0,2},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1}
                                                    };

    //Pij avec i inner loop
    std::vector<std::array<double,4> > poles_r(poles_t.size());
    int ni = 5 , nj =8;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            poles_r[j + nj * i] = poles_t[i + ni * j];
        }
    }

    gbs::BSSurfaceRational<double,3> srfNURBS(poles_r,ku,kv,p,q);

    pt = srfNURBS.value(u,v);
    // gbs::plot(
    //     srfNURBS,
    //     gbs::points_vector<double,3>{pt}
    // );

    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt(srfNURBS, pt, 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
        ASSERT_NEAR(res_v, v, 1e-6);
    }

}

TEST(tests_extrema, PS_f)
{
    std::vector<float> ku_flat = {0.,0.,0.,5.,5.,5.};
    std::vector<float> kv_flat = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;

    //Pij avec j inner loop
    const std::vector<std::array<float,3> > poles =
    {
        {0,0,0},{1,0,1},{2,0,2},
        {0,1,0},{1.2,1,1},{2,1,1},
        {0,2,1},{1.2,2,1},{2,2,2},
        {0,3,1},{1.1,3,1},{2,3,3},
        {0,4,0},{1,4,1},{2,4,0},
    };
    

    gbs::BSSurface<float,3> srf(poles,ku_flat,kv_flat,p,q);
    auto u = 2.5;
    auto v = 1.;

    auto pt = srf.value(u,v);

    // gbs::plot(
    //     srf,
    //     gbs::points_vector<double,3>{srf(res.u,res.v),pt}
    // );


    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt<float,3>(srf, pt, 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 5e-5);
        ASSERT_NEAR(res_v, v, 5e-5);
    }
}

TEST(tests_extrema, CS)
{
    std::vector<double> ku_flat = {0.,0.,2.,2.};
    std::vector<double> kv_flat = {0.,0.,1.,1.};
    size_t p = 1;
    size_t q = 1;
    const std::vector<std::array<double,3> > polesS =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,0},
    };
    const std::vector<std::array<double,3> > polesC =
    {
        {0.5,0.5,-1. },{0.5,0.5,1.},
    };
    gbs::BSSurface<double,3> srf(polesS,ku_flat,kv_flat,p,q);
    gbs::BSCurve<double,3> crv(polesC,ku_flat,p);

    // gbs::plot(srf,crv);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve(srf, crv, 1.e-6);
        ASSERT_LT(res_d, 1e-6);
        ASSERT_NEAR(res_u, 1., 1e-6);
        ASSERT_NEAR(res_v, 0.5, 1e-6);
        ASSERT_NEAR(res_uc, 1., 1e-6);
    }

    std::vector<double> ku_flat2 = {0.,0.,0.,1.,1.,1.};
    std::vector<double> kv_flat2 = {0.,0.,0.,1.,1.,1.};
    p = 2;
    q = 2;
    const std::vector<std::array<double,4> > polesS2 =
    {
        {0,0,0  ,1},{0.5,0,0   ,1},{1,0,0  ,1},
        {0,0.5,0,1},{5*0.75,5*0.5,5*1.,5},{1,0.5,0,1},
        {0,1,0  ,1},{0.5,1.,0   ,1},{1,1,0 ,1},
    };
        const std::vector<std::array<double,3> > polesC2 =
    {
        {0.5,0.5,-.5 },{0.,0.5,0.},{0.5,0.5,0.75},
    };
    gbs::BSSurfaceRational<double,3> srf2(polesS2,ku_flat2,kv_flat2,p,q);
    gbs::BSCurve<double,3> crv2(polesC2,ku_flat2,p);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve(srf2, crv2, 1.e-6);
        auto I = crv2.value(res_uc);

        ASSERT_LT(res_d, 2e-5);
    }

    // auto colors = vtkSmartPointer<vtkNamedColors>::New();
    // gbs::plot(
    //     srf2,
    //     crv2,
    //     gbs::make_actor(gbs::points_vector<double,3>{I},30.,true,colors->GetColor4d("Blue").GetData())
    //     );

}
TEST(tests_extrema, CS_f)
{
    std::vector<float> ku_flat = {0.,0.,1.,1.};
    std::vector<float> kv_flat = {0.,0.,1.,1.};
    size_t p = 1;
    size_t q = 1;
    const std::vector<std::array<float,3> > polesS =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,0},
    };
    const std::vector<std::array<float,3> > polesC =
    {
        {0.5,0.5,-1. },{0.5,0.5,1.},
    };
    gbs::BSSurface<float,3> srf(polesS,ku_flat,kv_flat,p,q);
    gbs::BSCurve<float,3> crv(polesC,ku_flat,p);

    // gbs::plot(srf,crv);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve<float,3>(srf, crv, 1.e-6);
        ASSERT_LT(res_d, 1e-6);
        ASSERT_NEAR(res_u, 0.5, 1e-6);
        ASSERT_NEAR(res_v, 0.5, 1e-6);
        ASSERT_NEAR(res_uc, 0.5, 1e-6);
    }

}
TEST(tests_extrema, CC)
{
    auto c1 = gbs::build_circle<double,3>(1.);
    auto c2 = gbs::build_segment<double,3>({0.,0.,0.},{1.,1.,0.});
    auto [res_u1,res_u2,res_d] = gbs::extrema_curve_curve(c1,c2,1.e-6);

    ASSERT_LT(res_d,5e-6);
    ASSERT_NEAR(res_u2,1.,5e-6);
}

TEST(tests_extrema, CC_f)
{
    auto c1 = gbs::build_circle<float,3>(1.);
    auto c2 = gbs::build_segment<float,3>({0.,0.,0.},{1.,1.,0.});
    auto [res_u1,res_u2,res_d] = gbs::extrema_curve_curve<float,3>(c1,c2,1.e-6);

    ASSERT_LT(res_d,5e-6);
    ASSERT_NEAR(res_u2,1.,5e-6);
}

TEST(tests_extrema, max_coordinate)
{
    gbs::BSCurve2d_d crv{
        {
            {0.0,0.},
            {0.1,0.},
            {0.4,1.},
            {0.6,1.},
            {0.9,0.},
            {1.0,0.},
        },
        {0., 0.4, 0.6, 1.},
        {4,1,1,4},
        3
    };
    auto [u,r] = gbs::max_coordinate(crv, 1, 0.5);
    ASSERT_NEAR(u,0.5,1e-6);
}
// The extrema return the true distance (they returned its square), stay in the bounds and
// find the closest point on either side of the seam of closed curves and surfaces.
namespace
{
    template <typename T>
    const gbs::ax2<T, 3> z_axis{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};

    template <typename T>
    auto polar(T r, T a, T z = T(0)) -> gbs::point<T, 3>
    {
        return {r * std::cos(a), r * std::sin(a), z};
    }

    template <typename T>
    T tol_of() { return std::is_same_v<T, float> ? T(2e-5) : T(1e-10); }

    // angles on both sides of the seam at 0 and in the middle
    template <typename T>
    const std::array<T, 5> angles{T(0), T(1e-3), T(2) * std::numbers::pi_v<T> - T(1e-3), T(1), std::numbers::pi_v<T>};

    template <typename T>
    void point_curve_true_distance()
    {
        const T tol = tol_of<T>();
        const auto circle = gbs::build_circle<T, 3>(T(1));
        for (const T a : angles<T>)
        {
            // in the plane, outside and inside, and above the circle
            for (const auto &[p, d] : std::array{std::pair{polar(T(1.25), a), T(0.25)}, std::pair{polar(T(0.6), a), T(0.4)},
                                                 std::pair{polar(T(1), a, T(0.3)), T(0.3)}})
            {
                const auto [u, dist] = gbs::extrema_curve_point(circle, p, T(1e-6));
                ASSERT_NEAR(dist, d, tol);
                ASSERT_LT(gbs::norm(circle(u) - polar(T(1), a)), 100 * tol);
                // local refinement from a start near the solution
                const auto [u_l, dist_l] = gbs::extrema_curve_point(circle, p, std::clamp(u + T(0.05), circle.bounds()[0], circle.bounds()[1]), T(1e-6));
                ASSERT_NEAR(dist_l, d, tol);
            }
        }
        // projection beyond the range of an arc: its start
        const auto arc = gbs::build_circle_arc<T>(T(1), T(0), std::numbers::pi_v<T> / T(2), z_axis<T>);
        const auto p = polar(T(1), -std::numbers::pi_v<T> / T(4));
        const auto [u, dist] = gbs::extrema_curve_point(arc, p, T(1e-6));
        ASSERT_NEAR(u, arc.bounds()[0], tol);
        ASSERT_NEAR(dist, gbs::norm(p - polar(T(1), T(0))), tol);
    }

    template <typename T>
    void point_surface_true_distance()
    {
        const T tol = tol_of<T>();
        const auto cylinder = gbs::build_cylinder<T>(T(1), z_axis<T>, T(0), T(2));
        const auto sphere = gbs::build_sphere<T>(T(2), z_axis<T>);
        const auto torus = gbs::build_torus<T>(T(3), T(1), z_axis<T>);
        for (const T a : angles<T>)
        {
            {
                const auto [u, v, d] = gbs::extrema_surf_pnt(cylinder, polar(T(1.5), a, T(1)), T(1e-6));
                ASSERT_NEAR(d, T(0.5), tol);
                ASSERT_LT(gbs::norm(cylinder(u, v) - polar(T(1), a, T(1))), 100 * tol);
            }
            {
                // beyond the top of the cylinder: v on its bound, u exact (the bound is active)
                const auto [u, v, d] = gbs::extrema_surf_pnt(cylinder, polar(T(1.5), a, T(3)), T(1e-6));
                ASSERT_NEAR(v, T(2), tol);
                ASSERT_NEAR(d, std::sqrt(T(1.25)), tol);
                ASSERT_LT(gbs::norm(cylinder(u, v) - polar(T(1), a, T(2))), 100 * tol);
            }
            {
                const auto dir = polar(std::cos(T(0.4)), a, std::sin(T(0.4)));
                const auto [u, v, d] = gbs::extrema_surf_pnt(sphere, T(3) * dir, T(1e-6));
                ASSERT_NEAR(d, T(1), tol);
                ASSERT_LT(gbs::norm(sphere(u, v) - T(2) * dir), 100 * tol);
            }
            {
                // near the seams of both parameters of the torus
                const T b = -a;
                const auto n = polar(std::cos(b), a, std::sin(b));
                const auto on = polar(T(3), a) + n;
                const auto [u, v, d] = gbs::extrema_surf_pnt(torus, on + T(0.3) * n, T(1e-6));
                ASSERT_NEAR(d, T(0.3), tol);
                ASSERT_LT(gbs::norm(torus(u, v) - on), 100 * tol);
            }
        }
    }

    template <typename T>
    void curve_curve_true_distance()
    {
        const T tol = tol_of<T>();
        const T h = T(0.3);
        const auto s1 = gbs::build_segment<T, 3>({-1., 0., 0.}, {1., 0., 0.});
        const auto s2 = gbs::build_segment<T, 3>({0., -1., h}, {0., 1., h});
        {
            const auto [u1, u2, d] = gbs::extrema_curve_curve(s1, s2, T(1e-6));
            ASSERT_NEAR(d, h, tol);
            ASSERT_LT(gbs::norm(s1(u1) - gbs::point<T, 3>{0., 0., 0.}), 100 * tol);
            ASSERT_LT(gbs::norm(s2(u2) - gbs::point<T, 3>{0., 0., h}), 100 * tol);
        }
        {
            // a circle and a segment above its seam
            const auto circle = gbs::build_circle<T, 3>(T(1));
            const auto s = gbs::build_segment<T, 3>({0.5, 0., h}, {2., 0., h});
            const auto [u1, u2, d] = gbs::extrema_curve_curve(circle, s, T(1e-6));
            ASSERT_NEAR(d, h, tol);
            ASSERT_LT(gbs::norm(circle(u1) - gbs::point<T, 3>{1., 0., 0.}), 100 * tol);
            const auto [u1_l, u2_l, d_l] = gbs::extrema_curve_curve(circle, s, circle.bounds()[1] - T(0.1), s.bounds()[0], T(1e-6));
            ASSERT_NEAR(d_l, h, tol);
            // the list version
            std::vector<const gbs::Curve<T, 3> *> lst{&s1, &s2, &s};
            const auto res = gbs::extrema_curve_curve_lst(circle, lst.begin(), lst.end(), T(1e-6));
            ASSERT_NEAR(res[0][2], T(0), tol);
            ASSERT_NEAR(res[1][2], h, tol);
            ASSERT_NEAR(res[2][2], h, tol);
        }
    }

    template <typename T>
    void surface_curve_true_distance()
    {
        const T tol = tol_of<T>();
        const auto sphere = gbs::build_sphere<T>(T(1), z_axis<T>);
        for (const T a : angles<T>)
        {
            const auto dir = polar(T(1), a);
            const auto ortho = polar(T(1), a + std::numbers::pi_v<T> / T(2));
            const auto s = gbs::build_segment<T, 3>(T(1.5) * dir - T(2) * ortho, T(1.5) * dir + T(2) * ortho);
            const auto [u, v, uc, d] = gbs::extrema_surf_curve(sphere, s, T(1e-6));
            ASSERT_NEAR(d, T(0.5), tol);
            ASSERT_LT(gbs::norm(sphere(u, v) - dir), 100 * tol);
            ASSERT_LT(gbs::norm(s(uc) - T(1.5) * dir), 100 * tol);
            const auto [u_l, v_l, uc_l, d_l] = gbs::extrema_surf_curve(sphere, s, uc + T(0.1), u, v + T(0.1), T(1e-6));
            ASSERT_NEAR(d_l, T(0.5), tol);
        }
    }
}

TEST(tests_extrema, true_distance_PC)
{
    point_curve_true_distance<double>();
    point_curve_true_distance<float>();
}

TEST(tests_extrema, true_distance_PS)
{
    point_surface_true_distance<double>();
    point_surface_true_distance<float>();
}

TEST(tests_extrema, true_distance_CC)
{
    curve_curve_true_distance<double>();
    curve_curve_true_distance<float>();
}

TEST(tests_extrema, true_distance_CS)
{
    surface_curve_true_distance<double>();
    surface_curve_true_distance<float>();
}

TEST(tests_extrema, build_intersections_threshold)
{
    // the curves are 1e-4 apart: farther than the tolerance (the squared distance was 1e-8)
    const auto crv = gbs::build_segment<double, 3>({0., 0., 0.}, {1., 0., 0.});
    std::vector<gbs::BSCurve<double, 3>> lst{gbs::build_segment<double, 3>({0.5, -1., 1e-4}, {0.5, 1., 1e-4})};
    ASSERT_THROW(gbs::build_intersections(lst.begin(), lst.end(), crv, 1e-6), std::invalid_argument);
    lst[0] = gbs::build_segment<double, 3>({0.5, -1., 1e-8}, {0.5, 1., 1e-8});
    const auto u = gbs::build_intersections(lst.begin(), lst.end(), crv, 1e-6);
    ASSERT_NEAR(u[0], 0.5, 1e-9);
}

// Cases of the review: all samples at the same distance (pole, axis), no second derivatives,
// Gauss-Newton and diagonal fallbacks, a Newton step pointing out of the box.
TEST(tests_extrema, pole_and_axis_starts)
{
    const ax2<double, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    const auto sphere = gbs::build_sphere<double>(1., ax);
    {
        // on the axis above the pole: every sample of the edge v = pi/2 is the pole
        const auto [u, v, d] = gbs::extrema_surf_pnt(sphere, gbs::point<double, 3>{0., 0., 1.7}, 1e-6);
        ASSERT_NEAR(d, 0.7, 1e-12);
    }
    {
        // 0.039 rad from the pole: the two starts were two copies of the pole
        const double a = 3.11284, b = 1.53228, r = 1.35336;
        const gbs::point<double, 3> p{r * std::cos(b) * std::cos(a), r * std::cos(b) * std::sin(a), r * std::sin(b)};
        const auto [u, v, d] = gbs::extrema_surf_pnt(sphere, p, 1e-6);
        ASSERT_NEAR(d, r - 1., 1e-12);
    }
    {
        // on the axis of a cylinder: every sample of a row is at the same distance
        const auto cylinder = gbs::build_cylinder<double>(1., ax, 0., 2.);
        const auto [u, v, d] = gbs::extrema_surf_pnt(cylinder, gbs::point<double, 3>{0., 0., 1.3}, 1e-6);
        ASSERT_NEAR(d, 1., 1e-12);
        ASSERT_NEAR(v, 1.3, 1e-9);
    }
}

TEST(tests_extrema, no_second_derivatives)
{
    // an offset surface has first derivatives only: Gauss-Newton, without exceptions
    const ax2<double, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    auto cylinder = std::make_shared<gbs::BSSurfaceRational<double, 3>>(gbs::build_cylinder<double>(1., ax, 0., 2.));
    const gbs::brep::OffsetSurface<double> offset{cylinder, 0.5};
    ASSERT_FALSE(offset.has_second_derivatives());
    ASSERT_THROW(offset.value(0.3, 1., 2, 0), std::runtime_error);
    const auto q = offset.value(0.3, 1.);
    const double R = std::hypot(q[0], q[1]); // 1.5 or 0.5 according to the orientation
    for (const double a : {0., 1., 2.5, 6.2})
    {
        const gbs::point<double, 3> p{3. * std::cos(a), 3. * std::sin(a), 0.7};
        const auto [u, v, d] = gbs::extrema_surf_pnt(offset, p, 1e-6);
        ASSERT_NEAR(d, 3. - R, 1e-8);
    }
}

TEST(tests_extrema, gauss_newton_and_diagonal_fallbacks)
{
    // a quarter of circle with weights 1, 7.07, 100 (same curve, strongly varying speed):
    // far from the curve the Hessian is not positive definite
    const double w = std::sqrt(0.5), l = 10.;
    const std::vector<std::array<double, 3>> poles{{1., 0., 0.}, {1., 1., 0.}, {0., 1., 0.}};
    const std::vector<double> knots{0., 0., 0., 1., 1., 1.}, weights{1., l * w, l * l};
    const gbs::BSCurveRational<double, 3> arc(poles, knots, weights, 2);
    for (const auto &p : {gbs::point<double, 3>{3., 4., 0.}, gbs::point<double, 3>{30., 2., 5.}, gbs::point<double, 3>{0.1, 0.05, 0.}})
    {
        const double rho = std::hypot(p[0], p[1]);
        const auto [u, d] = gbs::extrema_curve_point(arc, p, 1e-6);
        ASSERT_NEAR(d, std::hypot(rho - 1., p[2]), 1e-12);
        const auto [u_l, d_l] = gbs::extrema_curve_point(arc, p, 0.999, 1e-6); // from the far end
        ASSERT_NEAR(d_l, std::hypot(rho - 1., p[2]), 1e-12);
    }
    // parallel segments: JᵀJ is singular, every pair of the overlap is a minimum
    const auto s1 = gbs::build_segment<double, 3>({0., 0., 0.}, {1., 0., 0.});
    const auto s2 = gbs::build_segment<double, 3>({0.2, 0.3, 0.}, {1.5, 0.3, 0.});
    const auto [u1, u2, d] = gbs::extrema_curve_curve(s1, s2, 1e-6);
    ASSERT_NEAR(d, 0.3, 1e-12);
    const auto [u1_l, u2_l, d_l] = gbs::extrema_curve_curve(s1, s2, 0., 1., 1e-6);
    ASSERT_NEAR(d_l, 0.3, 1e-12);
}

TEST(tests_extrema, newton_step_out_of_the_box)
{
    // a sheared planar patch S(u, v) = u X + v (0.9, 1, 0): the minimum lies on the edge u = 0,
    // where the coupled Newton step points out of the box
    const std::vector<std::array<double, 3>> poles{{0., 0., 0.}, {1., 0., 0.}, {0.9, 1., 0.}, {1.9, 1., 0.}};
    const gbs::BSSurface<double, 3> patch(poles, {0., 0., 1., 1.}, {0., 0., 1., 1.}, 1, 1);
    const gbs::point<double, 3> e{0.9, 1., 0.}, p{-0.5, 1.2, 0.3};
    const double t = (p * e) / (e * e);
    const double d_exact = gbs::norm(p - t * e);
    {
        const auto [u, v, d] = gbs::extrema_surf_pnt(patch, p, 1e-6);
        ASSERT_NEAR(d, d_exact, 1e-12);
        ASSERT_NEAR(u, 0., 1e-12);
        ASSERT_NEAR(v, t, 1e-9);
    }
    for (const auto &[u0, v0] : {std::pair{0., 0.95}, std::pair{0., 0.05}, std::pair{0.5, 0.5}, std::pair{1., 1.}})
    {
        const auto [u, v, d] = gbs::extrema_surf_pnt(patch, p, u0, v0, 1e-6);
        ASSERT_NEAR(d, d_exact, 1e-12);
    }
}

TEST(tests_extrema, global_by_subdivision)
{
    // points where the two nearest of 30 samples lie in the basin of a local minimum: the
    // subdivision finds the global one (checked against 200 000 samples)
    const std::vector<gbs::constrType<double, 3, 1>> Q{{{0., 0., 0.}}, {{1., 0., 0.}}, {{1., 1., 0.}}, {{1., 1., 2.}}, {{0., 1., 1.}}, {{0., -1., 4.}}};
    const auto crv = gbs::interpolate(Q, 3, gbs::KnotsCalcMode::CHORD_LENGTH);
    const std::array<std::pair<gbs::point<double, 3>, double>, 2> cases{
        std::pair{gbs::point<double, 3>{0.35054339991977651, -0.092326199527380637, 1.8786730147085242}, 1.05915626261},
        std::pair{gbs::point<double, 3>{0.9900054766472135, -0.20054087438863033, 2.9136017688174327}, 1.40501174382}};
    for (const auto &[p, d_min] : cases)
    {
        const auto [u, d] = gbs::extrema_curve_point(crv, p, 1e-6);
        ASSERT_NEAR(d, d_min, 1e-10);
        const auto sampled = gbs::project_point_curve(crv, p, crv.bounds(), 30, 1e-6);
        ASSERT_GT(sampled[1], d_min + 1e-2); // the sampled version misses it
    }
}
