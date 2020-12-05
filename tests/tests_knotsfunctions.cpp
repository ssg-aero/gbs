#include <gtest/gtest.h>

#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bscanalysis.h>

#include <gbs-occt/export.h>
#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/discetize.h>

#include <Eigen/Dense>

#include <GeomTools.hxx>

#include <algorithm>
const double tol = 1e-7;

using namespace gbs;
TEST(tests_knotsfunctions, insert_knot)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k,p);
    auto c1_3d_dp_cp(c1_3d_dp);

    c1_3d_dp.insertKnot(0.5,2);

    c1_3d_dp.insertKnot(1.5,10);
    std::vector<int> mult;
    std::vector<double> knots;
    gbs::unflat_knots(c1_3d_dp.knotsFlats(), mult, knots);
    ASSERT_EQ(mult[3],2);

    c1_3d_dp.insertKnot(4.5,2);
    for( int i = 0 ; i < 100; i++)
    {
        auto u = i / 99. * 5.;
        auto d = distance(c1_3d_dp.value(u),c1_3d_dp_cp.value(u));
        ASSERT_LT(d,tol);
    }


    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp    ));
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp_cp ));

    // for( auto c : crv_lst) GeomTools::Dump(c,std::cout);

    std::vector<std::array<double,4> > polesW =
    {
        {0.,0.,0.,1.5},
        {0.,1.,0.,1.},
        {1.,1.,0.,1.},
        {1.,1.,1.,1.},
        {1.,1.,2.,1.},
        {3.,1.,1.,1.},
        {0.,4.*0.5,1.*0.5,0.5},
    };

    auto c1_3d_dp_w = gbs::BSCurveRational<double,3>(polesW,k,p);
    auto c1_3d_dp_w_cp(c1_3d_dp_w);
    c1_3d_dp_w.insertKnot(2.5,1);
    for( int i = 0 ; i < 100; i++)
    {
        auto u = i / 99. * 5.;
        auto d = distance(c1_3d_dp_w.value(u),c1_3d_dp_w_cp.value(u));
        ASSERT_LT(d,tol);
    }

    crv_lst.push_back( occt_utils::NURBSplineCurve( c1_3d_dp_w    ));
    crv_lst.push_back( occt_utils::NURBSplineCurve( c1_3d_dp_w_cp ));

    occt_utils::to_iges(crv_lst,"insert_knot.igs");
}

TEST(tests_knotsfunctions, refine)
{
    std::vector<double> k = { 0. , 0.1 , 0.4 , 0.8 , 1.};
    std::vector<size_t> m = {6,1,3,2,6};
    auto k_flat = flat_knots(k,m);

    std::vector<double> k_;
    unflat_knots(k_flat,m,k_);
    auto p = m.front()-1;
    auto np= k_flat.size()-p-1;
    std::vector<std::array<double,3> > poles(np);
    std::generate(poles.begin(),poles.end(),[&,a=0.]() mutable
    {
        auto p = std::array<double,3>{1.*cos(a),1.*sin(a),0.3*a};
        a+=2*M_PI/(np-1.);
        return p; }
    );
    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k_flat,p);


    auto mult_max = *std::max_element(std::next(m.begin()),std::next(m.end(),-1));

    for(auto i = 1 ; i < k.size()-1; i++)
    {
        if(m[i]!=mult_max)
        {
            for(auto j = 0 ; j < mult_max-m[i]; j++) 
            {
                insert_knot(k[i],p,k_flat,poles);
            }
        }
    }

    ASSERT_EQ(poles.size(), k_flat.size()-p-1);

    auto c2_3d_dp = gbs::BSCurve<double,3>(poles,k_flat,p);

    // auto n = 30;
    // auto u = gbs::make_range(k.front(),k.back(),n);
    // std::vector<double> k_new(k_);
    // k_new.front()=k_front();
    // std::transform(
    //     u.begin(),u.end(),
    //     k_new.begin(),
    //     std::next(k_)
    // )
    std::vector<Handle_Geom_Geometry> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp    ));
    crv_lst.push_back( occt_utils::BSplineCurve( c2_3d_dp    ));


    occt_utils::to_iges(crv_lst,"tests/out/refine.igs");
}

TEST(tests_knotsfunctions, changeBounds)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k,p);
    gbs::change_bounds(1.,2.,k);
    auto c2_3d_dp = gbs::BSCurve<double,3>(poles,k,p);

    ASSERT_NEAR(k.front(),1.,knot_eps);
    ASSERT_NEAR(k.back(),2.,knot_eps);

    auto pts = gbs::discretize(c2_3d_dp,20);
    auto deviation = gbs::dev_from_points(pts,c1_3d_dp);
    ASSERT_LT(deviation.d_max,1.5e-5);
}

TEST(tests_knotsfunctions, remove_knot)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k,p);
    auto c1_3d_dp_cp(c1_3d_dp);

    c1_3d_dp.insertKnot(0.5,2); // a priori occt teste si mult > deg, límplémentation gbs semble passer dans ce cas et donne la valeur correcte.

    c1_3d_dp.removeKnot(1.,tol);

    for( int i = 0 ; i < 100; i++)
    {
        auto u = i / 99. * 5.;
        auto d = distance(c1_3d_dp.value(u),c1_3d_dp_cp.value(u));
        std::cout << d <<std::endl;
        // ASSERT_LT(d,tol);
    }


    std::vector<Handle_Geom_Geometry> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp    ));
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp_cp ));


    occt_utils::to_iges(crv_lst,"tests/out/remove_knot.igs");
}

TEST(tests_knotsfunctions, reparam1)
{
    std::vector<double> k1 = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    auto c1_3d_dp = gbs::BSCurve<double,3>(poles1,k1,p);

    auto n = poles1.size();

    std::vector<double> k2 = {0., 0., 0., 0.5, 1., 2.5, 4, 5., 5., 5.};

    auto u = gbs::make_range(0., 5., n);

    Eigen::MatrixXd N(n, n);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            {
                N(i, j) = gbs::basis_function(u[i], j, p,0, k2);
            }
        }
    }
    auto N_inv = N.partialPivLu(); //TODO solve banded system

    VectorX<double> b(n);
    std::vector<std::array<double, 3>> poles2(n);
    for (int d = 0; d < 3; d++)
    {
        for (int i = 0; i < n; i++)
        {
            b(i) = c1_3d_dp.value(u[i])[d];//Pas top 
        }

        auto x = N_inv.solve(b);
        
        for (int i = 0; i < n; i++)
        {
            poles2[i][d] = x(i);
        }
    }

    auto c2_3d_dp = gbs::BSCurve<double,3>(poles2,k2,p);
    std::vector<Handle_Geom_Geometry> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp ));
    crv_lst.push_back( occt_utils::BSplineCurve( c2_3d_dp ));

    occt_utils::to_iges(crv_lst,"tests/out/reparam1.igs");
}

TEST(tests_knotsfunctions, reparam2)
{
    std::vector<double> k1 = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    auto c1_3d_dp = gbs::BSCurve<double,3>(poles1,k1,p);


    std::vector<double> k2 = {0., 0., 0., 0.5, 1., 2.5, 4, 5., 5., 5.};

    auto n = 1000;
    auto np = poles1.size();
    auto u = gbs::make_range(0., 5., n);

    Eigen::MatrixXd N(n, np);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < np; j++)
        {
            {
                N(i, j) = gbs::basis_function(u[i], j, p,0, k2);
            }
        }
    }
    // auto N_inv = (N.transpose() * N ).partialPivLu(); //TODO solve banded system
    auto N_inv = N.colPivHouseholderQr();

    VectorX<double> b(n);
    std::vector<std::array<double, 3>> poles2(np);

    for (int d = 0; d < 3; d++)
    {
        for (int i = 0; i < n; i++)
        {
            b(i) = c1_3d_dp.value(u[i])[d];
        }

        auto x = N_inv.solve(b);
        
        for (int i = 0; i < np; i++)
        {
            poles2[i][d] = x(i);
        }
    }

    auto c2_3d_dp = gbs::BSCurve<double,3>(poles2,k2,p);
    std::vector<Handle_Geom_Geometry> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp ));
    crv_lst.push_back( occt_utils::BSplineCurve( c2_3d_dp ));

    occt_utils::to_iges(crv_lst,"tests/out/reparam2.igs");
}

TEST(tests_knotsfunctions, reparam3)
{
    size_t p = 3;
    std::vector<double> k1 = {0., 0., 0., 0., 1, 2, 3, 4, 5.,5., 5., 5.};

    std::vector<std::array<double,3> > poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {0.3,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    auto np = poles1.size();
    ASSERT_EQ(k1.size()-p-1,np);

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles1,k1,p);

    // std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1., 2.5, 4, 5., 5., 5., 5.};
    std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1.,2., 2.5,3., 4, 5., 5., 5., 5.};
    // std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1.,1.5, 2.5,3.1, 4, 5., 5., 5., 5.};
    np = k2.size()-p-1;

    auto n = 300;
    auto u = gbs::make_range(0., 5., n);

    Eigen::MatrixXd N(n-2, np-2);
    for (size_t i = 1; i < n-1; i++)
    {
        for (size_t j = 1; j < np-1; j++)
        {
            {
                N(i-1, j-1) = gbs::basis_function(u[i], j, p,0, k2);
            }
        }
    }
    // auto N_inv = (N.transpose() * N ).partialPivLu(); //TODO solve banded system
    auto N_inv = N.colPivHouseholderQr();

    VectorX<double> b(n-2);
    std::vector<std::array<double, 3>> poles2(np);
    poles2.front() = poles1.front();
    poles2.back() = poles1.back();
    for (int d = 0; d < 3; d++)
    {
        for (int i = 1; i < n-1; i++)
        {
            b(i-1) = c1_3d_dp.value(u[i])[d] - gbs::basis_function(u[i], 0, p,0, k2) * (poles2.front()[d]) - gbs::basis_function(u[i], np-1, p,0, k2) * (poles2.back()[d]);
        }

        auto x = N_inv.solve(b);
        
        for (int i = 1; i < np-1; i++)
        {
            poles2[i][d] = x(i-1);
        }
    }
    auto c2_3d_dp = gbs::BSCurve<double,3>(poles2,k2,p);

    auto d = distance(c1_3d_dp.value(1.),c2_3d_dp.value(0.5));

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp ));
    crv_lst.push_back( occt_utils::BSplineCurve( c2_3d_dp ));

    for(auto c : crv_lst) GeomTools::Dump(c,std::cout);

    occt_utils::to_iges(crv_lst,"tests/out/reparam3.igs");
}

// TEST(tests_knotsfunctions, reparam4)
// {
//     size_t p = 5;
//     std::vector<double> k1 = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1. };

//     std::vector<std::array<double,3> > poles1 =
//     {
//         {0.,0.,0.},
//         {0.,1.,0.},
//         {0.3,1.,0.},
//         {1.,1.,0.},
//         {1.,1.,1.},
//         {1.,1.,2.},
//         {3.,1.,1.},
//         {1.,2.,1.2},
//         {0.,4.,1.},
//     };
//     auto np = poles1.size();
//     ASSERT_EQ(k1.size()-p-1,np);

//     auto c1_3d_dp = gbs::BSCurve(poles1,k1,p);
// /*
//     // std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1., 2.5, 4, 5., 5., 5., 5.};
//     std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1.,2., 2.5,3., 4, 5., 5., 5., 5.};
//     // std::vector<double> k2 = {0., 0., 0., 0., 0.5, 1.,1.5, 2.5,3.1, 4, 5., 5., 5., 5.};
//     np = k2.size()-p-1;

//     auto n = 300;
//     auto u = gbs::make_range(0., 5., n);

//     Eigen::MatrixXd N(n-2, np-2);
//     for (size_t i = 1; i < n-1; i++)
//     {
//         for (size_t j = 1; j < np-1; j++)
//         {
//             {
//                 N(i-1, j-1) = gbs::basis_function(u[i], j, p,0, k2);
//             }
//         }
//     }
//     // auto N_inv = (N.transpose() * N ).partialPivLu(); //TODO solve banded system
//     auto N_inv = N.colPivHouseholderQr();

//     Eigen::VectorX<double> b(n-2);
//     std::vector<std::array<double, 3>> poles2(np);
//     poles2.front() = poles1.front();
//     poles2.back() = poles1.back();
//     for (int d = 0; d < 3; d++)
//     {
//         for (int i = 1; i < n-1; i++)
//         {
//             b(i-1) = c1_3d_dp.value(u[i])[d] - gbs::basis_function(u[i], 0, p,0, k2) * (poles2.front()[d]) - gbs::basis_function(u[i], np-1, p,0, k2) * (poles2.back()[d]);
//         }

//         auto x = N_inv.solve(b);
        
//         for (int i = 1; i < np-1; i++)
//         {
//             poles2[i][d] = x(i-1);
//         }
//     }
//     auto c2_3d_dp = gbs::BSCurve(poles2,k2,p);
// */
//     std::vector<Handle_Geom_Curve> crv_lst;
//     crv_lst.push_back( occt_utils::BSplineCurve( c1_3d_dp ));
//     // crv_lst.push_back( occt_utils::BSplineCurve( c2_3d_dp ));

//     // for(auto c : crv_lst) GeomTools::Dump(c,std::cout);

//     occt_utils::to_iges(crv_lst,"tests/out/reparam4.igs");
// }