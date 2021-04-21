#include <gtest/gtest.h>
#include <gbs/transform.h>
#include <gbs/gbslib.h>
#include <gbs/bscbuild.h>
#include <gbs-io/print.h>

using namespace gbs;
const double tol = 1e-6;
const auto PI =acos(-1);

TEST(tests_transform, transform_point)
{
    {
        gbs::point<double, 3> x{0, 0, 0};
        gbs::translate(x, {1., 0., 0.});
        ASSERT_DOUBLE_EQ(x[0], 1.);
    }

    {
        gbs::point<double, 2> x{1, 0};
        gbs::rotate(x, PI / 2);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
    }

    {
        gbs::point<double, 3> x{1, 0, 0};
        gbs::rotate(x, PI / 2, {0.,0.,1.});
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
    }

    {
        gbs::point<double, 3> x{1, 0, 0};
        gbs::rotate(x, PI / 2, {1.,0.,0.});
        ASSERT_NEAR(x[0], 1.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
    }

    {
        gbs::point<double, 3> x{1, 0, 0};
        gbs::rotate(x, PI / 2, {0.,1.,0.});
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        ASSERT_NEAR(x[2],-1.,tol);
    }

    {
        gbs::point<double, 3> x{1, 0, 0};
        gbs::scale(x,2.);
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
    }

    {
        gbs::point<double, 3> x{1, 1, 0};
        gbs::scale(x,2.,0);
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
    }

}


TEST(tests_transform, transform_BSCurve)
{
    {
        auto crv = gbs::build_segment<double,3>({0.,0.,0.},{1.,0.,0.});
        gbs::translate(crv,{1.,1.,0.});
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 1.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
        x = crv(1.);
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
        ASSERT_NEAR(x[2], 0.,tol);
    }

    {
        auto crv = gbs::build_segment<double,2>({0.,0.},{1.,0.});
        gbs::rotate(crv,PI / 2);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        x = crv(1.);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
    }

    {
        auto crv = gbs::build_segment<double,2>({0.,0.},{1.,0.});
        gbs::scale(crv,2.);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        x = crv(1.);
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
    }

    {
        auto crv = gbs::build_segment<double,2>({0.,0.},{1.,1.});
        gbs::scale(crv,2.,0);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        x = crv.end();
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
    }
}

TEST(tests_transform, transform_BSCurveRational)
{
    {
        auto crv = gbs::build_ellipse<double,3>(1.,2.);
        auto x0= crv(0.0);
        auto x1= crv(0.5);
        std::array<double,3> t = {1.,1.,0.};
        gbs::translate(crv,t);
        auto x = crv(0.);
        ASSERT_NEAR(norm(x-x0-t), 0.,tol);

        x = crv(0.5);
        ASSERT_NEAR(norm(x-x1-t), 0.,tol);
    }

    {
        auto crv = gbs::build_ellipse<double,2>(1.,2.);
        auto x0= crv(0.0);
        auto x1= crv(0.25);
        gbs::rotate(crv,PI / 2);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 1.,tol);
        x = crv(0.25);
        ASSERT_NEAR(x[0],-2.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
    }
    
    {
        auto crv = gbs::build_ellipse<double,2>(1.,2.);
        auto x0= crv(0.0);
        auto x1= crv(0.25);
        gbs::scale(crv,2.);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 2.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        x = crv(0.25);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 4.,tol);
    }
    {
        auto crv = gbs::build_ellipse<double,2>(1.,2.);
        auto x0= crv(0.0);
        auto x1= crv(0.25);
        gbs::scale(crv,2.,1);
        auto x = crv(0.);
        ASSERT_NEAR(x[0], 1.,tol);
        ASSERT_NEAR(x[1], 0.,tol);
        x = crv(0.25);
        ASSERT_NEAR(x[0], 0.,tol);
        ASSERT_NEAR(x[1], 4.,tol);
    }
}
/*
using namespace gbs;

template<typename T>
auto build_trf_loc_and_matrix(const ax2<T,3> &R)
{
    Matrix3<double> P;
    Vector3<double> Loc;
    auto O  = R[0];
    auto P1 = R[2] / norm(R[2]);
    auto P3 = R[1] / norm(R[1]);
    auto P2 = P3 ^ P1;
    std::copy(O.begin(),O.end(),Loc.data());
    std::copy(P1.begin(),P1.end(),P.col(0).data());
    std::copy(P2.begin(),P2.end(),P.col(1).data());
    std::copy(P3.begin(),P3.end(),P.col(2).data());
    return std::make_pair(Loc,P);
}
template<typename T>
auto build_trf_matrix(const Vector3<T> &O, const Matrix3<T> &P) -> Matrix4<T>
{
    Matrix3<double> M;
    M.setZero();
    for(size_t i {} ; i < 3; i++)
    {
        for (size_t j {}; j < 3; j++)
        {
            M(i, j) = P(i, j);
        }
        M(i,3) = O(i);
    }
    M(3,3) = 1.;
    return M;
}

template<typename T>
auto build_trf_matrix(const ax2<T,3> &R) -> Matrix4<T>
{
    auto [O, P] = build_trf_loc_and_matrix(R);
    return build_trf_matrix(O,P);
}


template <typename T>
auto build_trf_loc_and_matrix(const ax2<T,3> &R_from,const ax2<T,3> &R_to)
{
    auto [O_to, P_to] = build_trf_loc_and_matrix(R_to);
    auto [O_from, P_from] = build_trf_loc_and_matrix(R_from);

    auto P = P_to;
    Vector3<double> Loc   = -1. * P * O_to;
    // Loc.reverse();
    O_from =  P * O_from;
    Loc = Loc + O_from;
    P = P * P_from;
    return std::make_pair(Loc,P);
}

template <typename T>
auto build_trf_matrix(const ax2<T,3> &R_from,const ax2<T,3> &R_to) ->Matrix4<T>
{
    auto [O, P] = build_trf_loc_and_matrix(R_from, R_to);
    return build_trf_matrix(O,P);
}

template<typename T>
auto transform(point<T,3> &pt, const Vector3<T> &Loc, const Matrix3<T> &M) -> void
{
    Vector3<double> P;
    std::copy(pt.begin(),pt.end(),P.data());
    P = M * P;
    P = P + Loc;

    for(size_t i {}; i < 3; i++)
    {
        pt[i] = P(i);
    }
}

template<typename T>
auto transformed(const point<T,3> &pt, const Vector3<T> &Loc, const Matrix3<T> &M) -> point<T,3>
{
    point<T,3> pt_moved {pt};
    transform(pt_moved,Loc,M);
    return pt_moved;
}

template<typename T>
auto transform(point<T,3> &pt, const Matrix4<T> &M) -> void
{
    Vector4<double> P;
    std::copy(pt.begin(),pt.end(),P.data());
    P(3) = 1.;
    P = M * P;

    for(size_t i {}; i < 3; i++)
    {
        pt[i] = P(i);
    }
}

template<typename T>
auto transformed(const point<T,3> &pt, const Matrix4<T> &M) -> point<T,3>
{
    point<T,3> pt_moved {pt};
    transform(pt_moved,M);
    return pt_moved;
}
*/

TEST(tests_transform,base_change)
{

    {
        ax2<double,3> R_from {
            {
                {0.,0.,0.},
                {0.,0.,1.},
                {1.,0.,0.}
            }
        };

        ax2<double, 3> R_to{
            {
                {0., 0., 0.},
                {0., 1., 0.}, // z
                {0., 0., 1.}  // x
            }
        };
        auto [Loc, P] = build_trf_loc_and_matrix<double>(R_from, R_to);
        ASSERT_DOUBLE_EQ(0., gbs::norm(transformed(R_from[2], Loc, P) - R_to[2]));
        ASSERT_DOUBLE_EQ(0., gbs::norm(transformed(R_from[1], Loc, P) - R_to[1]));
    }


    {
        ax2<double,3> R_from {
            {
                {0.,0.,0.},
                {0.,0.,1.},
                {1.,0.,0.}
            }
        };

        ax2<double, 3> R_to{
            {
                {0., 1., 0.},
                {0., 0., 1.}, // z
                {0., 1., 0.}  // x
            }
        };
        auto [Loc, P] = build_trf_loc_and_matrix<double>(R_from, R_to);
        auto [Loc_to, P_to] = build_trf_loc_and_matrix<double>(R_to);
        auto M = build_trf_matrix<double>(R_from, R_to);
        std::cerr<<Loc<<std::endl<<std::endl;
        {
            {
                auto pt = transformed({0., 0., 0.}, Loc, P);
                // std::cerr << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
                ASSERT_DOUBLE_EQ(0., gbs::norm(pt - gbs::point<double, 3>{1., 0., 0.}));
            }
            {
                auto pt = transformed({1., 0., 0.}, Loc, P);
                // std::cerr << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
                ASSERT_DOUBLE_EQ(0., gbs::norm(pt - gbs::point<double, 3>{1., 1., 0.}));
            }
            {
                auto pt = transformed({1., 0., 0.}, M);
                // std::cerr << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
                ASSERT_DOUBLE_EQ(0., gbs::norm(pt - gbs::point<double, 3>{1., 1., 0.}));
            }
        }
    }
}

// TODO test with surfaces