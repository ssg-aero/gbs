#include <gtest/gtest.h>

#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bscanalysis.h>
#include <gbs/bssbuild.h>

#include <gbs-render/vtkGbsRender.h>

#include <Eigen/Dense>

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

    auto c2_3d_dp = gbs::BSCurve<double,3>(poles,k,p);
    auto pt = c2_3d_dp(2.33);
    auto ik = c2_3d_dp.insertKnot(2.33,p);

    gbs::plot(
        gbs::crv_dsp<double,3,false>{
            .c =&c2_3d_dp,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {0.,1.,0.},
            .col_ctrl = {0.,0.,0.},
            .show_curvature=true,
            } // c++20
            ,
        gbs::crv_dsp<double,3,false>{
            .c =&c1_3d_dp,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {0.,1.,0.},
            .col_ctrl = {0.,0.,0.},
            .show_curvature=true,
            } // c++20
            
            
        );

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

    ASSERT_NEAR(k.front(),1.,knot_eps<double>);
    ASSERT_NEAR(k.back(),2.,knot_eps<double>);

    auto pts = gbs::discretize(c2_3d_dp,20);
    auto [u_max, d_max, d_avg] = gbs::dev_from_points(pts,c1_3d_dp);
    ASSERT_LT(d_max,1.5e-5);
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

    c1_3d_dp.removeKnot(1.,knot_eps<double>);

    for( int i = 0 ; i < 100; i++)
    {
        auto u = i / 99. * 5.;
        auto d = distance(c1_3d_dp.value(u),c1_3d_dp_cp.value(u));
        std::cout << d <<std::endl;
        // ASSERT_LT(d,tol);
    }

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

    gbs::VectorX<double> b(n);
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

    gbs::VectorX<double> b(n);
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

    gbs::VectorX<double> b(n-2);
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

TEST(tests_knotsfunctions, surf_knot_insert)
{

    size_t p = 2;
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    gbs::points_vector_3d_d poles1 =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {0.7,1.,0.},
        {1.,1.3,0.},
        {1.8,1.,0.3},
        {3.,1.,0.},
        {0.,4.,0.},
    };
    gbs::points_vector_3d_d poles2 =
    {
        {0.,0.,1.},
        {0.,1.,1.},
        {1.,1.,1.},
        {1.3,0.4,1.},
        {1.5,0.5,1.},
        {3.,1.,1.5},
        {2.,4.,1.},
    };
    gbs::points_vector_3d_d poles3 =
    {
        {0.,0.,2.},
        {0.,1.,2.},
        {0.5,1.,2.},
        {1.,1.,2.5},
        {1.5,1.,2.5},
        {3.,1.,2.},
        {1.,4.,2.},
    };
    gbs::BSCurve3d_d c1(poles1,k,p);
    gbs::BSCurve3d_d c2(poles2,k,p);
    gbs::BSCurve3d_d c3(poles3,k,p);
    gbs::BSCurve3d_d sp({{1.5,1.5,0.},{1.5,1.5,3.}},{0.,0.,3.,3.},1);

    std::list<gbs::BSCurveGeneral<double,3,false>*> bs_lst = {&c1,&c2,&c3};
    auto s = gbs::loft( bs_lst , sp);
    auto s1{s};

    auto u = 2.3;
    auto v = 0.33;

    auto pt_check = s(u,v);

    // s.insertKnotU(u,2); 
    // s.insertKnotV(v);

    auto pt = s(u,v);
    // gbs::plot(s,s1);
    BSCurve<double, 3> c{s.polesU(2), s.knotsFlatsU(), s.degreeU()};
    // BSCurve<double, 3> c{s.polesV(2), s.knotsFlatsV(), s.degreeV()};
    // c.changeBounds(0, 1);
    // std::cerr << c(0.5)[0] << " " << c(0.5)[1] << " " << c(0.5)[2];
    gbs::plot(c,s,s1,gbs::make_actor(points_vector<double,3>{pt_check},25.));
    // gbs::plot(
    //     gbs::crv_dsp<double, 3, false>{
    //         .c = &c});

    // ASSERT_LT(distance(pt,pt_check ),tol);

}

TEST(tests_knotsfunctions, increase_degree)
{
    using T = double;
    using namespace gbs;
    using namespace std;
    const size_t dim{3};

    vector<T> k = {0., 0., 0., 1, 2, 2, 4, 5., 5., 5.};
    vector<array<T,dim> > poles =
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
    auto crv1 = gbs::BSCurve<double,3>(poles,k,p);

    // Split bezier
    auto bezier_seg = bezier_segments(k, poles, p);

    vector<BSCurve<T,dim>> c_lst(bezier_seg.size());
    transform(bezier_seg.begin(), bezier_seg.end(),c_lst.begin(),
        [p](const auto &pa){
            return BSCurve<T,dim>{pa.first, {pa.second.first, pa.second.second},{p+1,p+1},p};
        }
    );
    
    // Increase Bezier degree
    size_t t {3}; // increase 3 times
    transform(
        bezier_seg.begin(), 
        bezier_seg.end(), 
        bezier_seg.begin(), 
        [p,t](const pair< vector<array<T,dim>>, pair<T,T> > &pa)
        {
            return make_pair(
                increase_bezier_degree(pa.first, p, t),
                pa.second
            );
        }
    );

    // Join poles and knots
    vector<array<T,dim>> poles_new{bezier_seg.front().first.begin(), bezier_seg.front().first.end()};
    vector<size_t> mults_new{p+t+1};
    vector<T> knots{bezier_seg.front().second.first,bezier_seg.front().second.second};
    for(size_t s{1}; s<bezier_seg.size(); s++)
    {
        const auto & seg = bezier_seg[s];
        for(size_t i{1}; i <= p+t; i++)
        {
            poles_new.push_back(seg.first[i]);
        }
        mults_new.push_back(p+t);
        knots.push_back(seg.second.second);
    }
    mults_new.push_back(p+t+1);
    BSCurve<T,dim> crv_new{poles_new, knots, mults_new,p+t};
    // Remove knots

    crv_dsp<T,dim,false> d_c1{
        .c =&crv1,
        .col_crv = {1.,0.,0.},
        .poles_on = true,
        .col_poles = {1.,0.,0.},
    };

    crv_dsp<T,dim,false> d_c_new{
        .c =&crv_new,
        .col_crv = {1.,0.,0.},
        .poles_on = true,
        .col_poles = {1.,0.,0.},
    };
   


    vector<BSCurve<T,dim>> c_lst_new(bezier_seg.size());
    transform(bezier_seg.begin(), bezier_seg.end(),c_lst_new.begin(),
        [p,t](const auto &pa){
            return BSCurve<T,dim>{pa.first, {pa.second.first, pa.second.second},{p+t+1,p+t+1},p+t};
        }
    );
      
    vector<crv_dsp<T,dim,false>> d_lst(c_lst.size());
    srand( (unsigned)time( NULL ) );
    transform(c_lst.begin(), c_lst.end(),d_lst.begin(),
    [](const auto &c){
        return crv_dsp<T,dim,false>{
            .c =&c,
            .col_crv = {(T) rand()/RAND_MAX,(T) rand()/RAND_MAX,(T) rand()/RAND_MAX},
            .poles_on = true,
            .col_poles = {0.,1.,0.},};
    });

    vector<crv_dsp<T,dim,false>> d_lst_new(c_lst_new.size());
    srand( (unsigned)time( NULL ) );
    transform(c_lst_new.begin(), c_lst_new.end(),d_lst_new.begin(),
    [](const auto &c){
        return crv_dsp<T,dim,false>{
            .c =&c,
            .col_crv = {(T) rand()/RAND_MAX,(T) rand()/RAND_MAX,(T) rand()/RAND_MAX},
            .poles_on = true,
            .col_poles = {0.,1.,0.},};
    });

    plot(  
        // d_lst_new
        d_c_new
        );
}

// TEST(tests_knotsfunctions, elevate_degree)
// {
//     using T = double;
//     using namespace gbs;
//     using namespace Eigen;
//     const size_t dim{2};

//     std::vector<T> knots= {0.0,0.10653709580151404,0.13170091743633459,0.20458556663833828};
//     std::vector<size_t> mults = {2,1,1,2};
//     auto knots_flats = flat_knots(knots,mults);
//     std::vector<std::array<T, dim>>poles =  {{0.3611771432346219,0.10511112605663975},{0.5,0.135},{0.53,0.15},{0.6179529455155457,0.19125644692155786}};
//     size_t p{1};

//     elevate_degree(knots_flats, poles, 3);
// }