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

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

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
    std::vector<int> mult;
    std::vector<double> knots;
    gbs::unflat_knots(c1_3d_dp.knotsFlats(), mult, knots);
    ASSERT_EQ(mult[1],2);

    c1_3d_dp.insertKnot(1.5,10);

    ASSERT_EQ(mult[1],2);

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


    if(PLOT_ON) gbs::plot(
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

TEST(tests_knotsfunctions, increase_deg_Pieg94)
{
    using namespace gbs;
    using T = double;
    const size_t dim{2};

    T ua =0.3;
    T ub=0.7;
    size_t ma = 1;
    size_t mb = 1;
    std::vector<T> U = {0., 0., 0., 0., ua, ub, 1., 1., 1., 1.};
    std::vector<std::array<T,dim> > Pw =
    {
        { 0.0, 0.0},
        {-0.2, 0.5},
        { 0.2, 1.0},
        { 0.8, 1.0},
        { 1.2, 0.5},
        { 1.0, 0.0}
    };
    size_t p = 3;  
    
    auto c1 = gbs::BSCurve<T,dim>(Pw,U,p);

    auto bezier_seg = bezier_segments(c1.knotsFlats(), c1.poles(), c1.degree());

    size_t t{2};
    bezier_seg[0].first = increase_bezier_degree(bezier_seg[0].first,p,t);
    bezier_seg[1].first = increase_bezier_degree(bezier_seg[1].first,p,t);
    bezier_seg[2].first = increase_bezier_degree(bezier_seg[2].first,p,t);

    p+=t;

    auto bez0 = gbs::BSCurve<T,dim>(bezier_seg[0].first,{bezier_seg[0].second.first,bezier_seg[0].second.second},{p+1,p+1},p);
    auto bez1 = gbs::BSCurve<T,dim>(bezier_seg[1].first,{bezier_seg[1].second.first,bezier_seg[1].second.second},{p+1,p+1},p);
    auto bez2 = gbs::BSCurve<T,dim>(bezier_seg[2].first,{bezier_seg[2].second.first,bezier_seg[2].second.second},{p+1,p+1},p);

    std::vector<std::array<T,dim> > Pw_new;
    Pw_new.insert(Pw_new.end(),bezier_seg[0].first.begin(),bezier_seg[0].first.end());
    Pw_new.insert(Pw_new.end(),std::next(bezier_seg[1].first.begin()),bezier_seg[1].first.end());
    Pw_new.insert(Pw_new.end(),std::next(bezier_seg[2].first.begin()),bezier_seg[2].first.end());

    U = flat_knots<T,size_t>({0.,ua,ub,1.},{p+1,p,p,p+1});
    remove_knot(ua, p-(ma+t), U, Pw_new, p);
    remove_knot(ub, p-(mb+t), U, Pw_new, p);
    auto c2 = gbs::BSCurve<T,dim>(Pw_new,U,p);

    auto c3{c1};
    c3.increaseDegree(t);

    if(PLOT_ON)
    {
        plot(
            crv_dsp<T,dim,false>{.c=&c1  ,.poles_on=true, .col_poles={1,0,0}},
            // crv_dsp<T,dim,false>{.c=&bez0,.poles_on=true, .col_poles={0,1,0}},
            // crv_dsp<T,dim,false>{.c=&bez1,.poles_on=true, .col_poles={0,1,0}},
            // crv_dsp<T,dim,false>{.c=&bez2,.poles_on=true, .col_poles={0,1,0}},
            crv_dsp<T,dim,false>{.c=&c3  ,.poles_on=true, .col_poles={0,1,}},
            crv_dsp<T,dim,false>{.c=&c2  ,.poles_on=true, .col_poles={0,0,1}}
        );
    }

}

TEST(tests_knotsfunctions, remove_knot_algo)
{
    GTEST_SKIP() << "Skipping this test for now.";
    using namespace gbs;
    using T = double;
    const size_t dim{3};

    std::vector<T> U = {0., 0., 0., 0.,1, 2, 2, 3, 4, 5., 5., 5., 5.};
    std::vector<std::array<T,dim> > Pw =
    {
        {0.,0.,0.},
        {0.,0.,1.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {2.,1.,2.},
        {0.,4.,1.},
    };
    size_t p = 3;

    size_t num{2};
    T u{1.5};
    for(size_t i{}; i<num; i++)
        insert_knot(u, p, U, Pw);

    auto r = std::distance( U.begin(), std::ranges::upper_bound(U, u) ) -1;
    auto s = multiplicity(u,U);
    auto m = U.size();
    auto n = Pw.size();
    auto first = r-p;
    auto last  = r-s;


    size_t t{0};
    std::list<std::array<T,dim>> Pi;
    std::list<std::array<T,dim>> Pj;
    for( ; t < num; t++)
    {
        for(size_t i{}; i < first; i++)
            ASSERT_DOUBLE_EQ(basis_function(u, i, p, 0, U), 0.);
        for(size_t i= last+1; i < n; i++)
            ASSERT_DOUBLE_EQ(basis_function(u, i, p, 0, U), 0.);
        auto i{first};
        auto j{last};
        Pi = {Pw[first-1]};
        Pj = {Pw[last+1]};
        while(j>i+t)
        {
            auto ai = (u-U[i  ])/(U[i+p+1+t]-U[i  ]);
            auto aj = (u-U[j-t])/(U[j+p+1  ]-U[j-t]);

            Pi.push_back( (Pw[i]-(1-ai)*Pi.back())/ ai);
            Pj.push_front((Pw[j]-  aj*Pj.front())/ (1-aj));
            i++;
            j--;
        }
        bool remove_flag=false;
        if(j<i+t)
        {
            if(distance(Pi.back(), Pj.front())<tol)
            {
                remove_flag=true;
                Pj.pop_front();
            }
        }
        else
        {
            auto ai = (u-U[i  ])/(U[i+p+1+t]-U[i  ]);
            if(distance(Pw[i], ai*Pi.back()+(1-ai)*Pj.front())<tol)
            {
                remove_flag=true;
                Pj.pop_front();
                Pi.pop_back();
            }
        }

        if(!remove_flag)
        {
            break;
        }
        else{
            auto i{first};
            auto j{last};

        }
        first--; last++;

    }
    if(t>0){
        auto Pw_beg=Pw.begin();
        std::vector<std::array<T,dim>> head{Pw_beg,std::next(Pw_beg,first)};
        std::vector<std::array<T,dim>> tail{std::next(Pw_beg,last+1), Pw.end()};
        Pw = std::move(head);
        Pw.insert(Pw.end(), Pi.begin(), Pi.end());
        Pw.insert(Pw.end(), Pj.begin(), Pj.end());
        Pw.insert(Pw.end(), tail.begin(), tail.end());

        U.erase(std::next(U.begin(),r-t+1), std::next(U.begin(),r+1));
    }

    ASSERT_EQ(U.size()-p-1, Pw.size());
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

    double u = 0.5;

    auto insert_func = [&](size_t t)
    {
        for (size_t i{}; i < t; i++)
            insert_knot(u, p, k, poles);
    };

    insert_func(p+2);

    size_t t{2};

    ASSERT_EQ( remove_knot(u,p,t,k,poles,tol), p );

    auto c2_3d_dp = gbs::BSCurve<double,3>(poles,k,p);

    int n = 1000;
    for( int i = 0 ; i < n; i++)
    {
        auto u = i / (n-1.) * 5.;
        auto d =gbs::distance(c1_3d_dp.value(u),c2_3d_dp.value(u));
        ASSERT_LT(d,tol);
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
    if(PLOT_ON) gbs::plot(c,s,s1,gbs::make_actor(points_vector<double,3>{pt_check},25.));
    // gbs::plot(
    //     gbs::crv_dsp<double, 3, false>{
    //         .c = &c});

    // ASSERT_LT(distance(pt,pt_check ),tol);

}

TEST(tests_knotsfunctions, bezier_segments)
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
    auto crv = gbs::BSCurve<double,3>(poles,k,p);

    // Split bezier
    auto bezier_seg = bezier_segments(k, poles, p);

    vector<BSCurve<T,dim>> c_lst(bezier_seg.size());
    transform(bezier_seg.begin(), bezier_seg.end(),c_lst.begin(),
        [p](const auto &pa){
            return BSCurve<T,dim>{pa.first, {pa.second.first, pa.second.second},{p+1,p+1},p};
        }
    );

    ranges::for_each( c_lst, [&crv](const auto &seg){
        auto [u1, u2] = seg.bounds();
        ASSERT_LT(norm(crv(u1)-seg(u1)), tol);
        ASSERT_LT(norm(crv(u2)-seg(u2)), tol);
    });

    if(PLOT_ON)
    {
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
        plot(d_lst);
    }

}

TEST(tests_knotsfunctions, increase_bezier_degree)
{
    using T = double;
    using namespace gbs;
    using namespace std;
    const size_t dim{3};

    vector<array<T,dim> > poles =
    {
        {0.,0.,0.2},
        {0.,1.,0.2},
        {1.,1.,0.2},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
    };
    size_t p = 5;
    size_t t = 3;

    auto poles_t = increase_bezier_degree(poles, p, t);

    BSCurve<T,dim> c1{poles,{0.,1.},{p+1,p+1},p};
    BSCurve<T,dim> c2{poles_t,{0.,1.},{p+t+1,p+t+1},p+t};
    BSCurveRational<T,dim-1> c1r{poles,{0.,1.},{p+1,p+1},p};
    BSCurveRational<T,dim-1> c2r{poles_t,{0.,1.},{p+t+1,p+t+1},p+t};

    for(int i{}; i < 100; i++)
    {
        T u = i / 99.;
        ASSERT_LT(norm(c1(u)-c2(u)), tol);
        ASSERT_LT(norm(c1r(u)-c2r(u)), tol);
    }

    if(PLOT_ON)
    {
        crv_dsp<T,dim,false> d_c1{
            .c =&c1,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {1.,0.,0.},
        };
        crv_dsp<T,dim,false> d_c2{
            .c =&c2,
            .col_crv = {0.,0.,1.},
            .poles_on = true,
            .col_poles = {0.,0.,1.},
        };
        plot(d_c1, d_c2);

        crv_dsp<T,dim-1,true> d_c1r{
            .c =&c1r,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {1.,0.,0.},
        };
        crv_dsp<T,dim-1,true> d_c2r{
            .c =&c2r,
            .col_crv = {0.,0.,1.},
            .poles_on = true,
            .col_poles = {0.,0.,1.},
        };
        plot(d_c1r, d_c2r);
    }
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
    auto crv2 = gbs::BSCurve<double,3>(poles,k,p);
    size_t t = 3;
    crv1.increaseDegree(t);

    ASSERT_EQ(crv1.degree(),p+t);

    size_t np = 1000;
    for(size_t i{}; i < np; i++)
    {
        T u = 5. * i / (np-1.);
        ASSERT_LT( distance(crv1(u),crv2(u)), std::numeric_limits<T>::epsilon()*100);
    }
  
    if(PLOT_ON){
        crv_dsp<T,dim,false> d_c1{
            .c =&crv1,
            .col_crv = {1.,0.,0.},
            .poles_on = true,
            .col_poles = {1.,0.,0.},
        };
        crv_dsp<T,dim,false> d_c2{
            .c =&crv2,
            .col_crv = {0.,0.,1.},
            .poles_on = true,
            .col_poles = {0.,0.,1.},
        };
        plot( d_c1, d_c2 );
    }
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