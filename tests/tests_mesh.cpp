#include <gtest/gtest.h>
#include <gbs/bscbuild.h>
#include <gbs/mesh/mshedge.h>
#include <gbs/mesh/tfi.h>
#include <gbs/mesh/smoothing.h>
#include <gbs/render/vtkcurvesrender.h>
#include <gbs/render/vtkgridrender.h>

TEST(tests_mesh, msh_ed)
{
    auto r_cir = 1.2f;
    auto cir = gbs::build_circle<float,2>(r_cir);
    auto p_cir = std::make_shared<gbs::BSCurveRational<float,2>>(cir);
    auto np = 11;
    gbs::msh_edge<float,2> ed1{p_cir};
    ed1.set_n_points(np);
    ed1.compute_pnts();
    ASSERT_FLOAT_EQ( r_cir,ed1.points()[0][0]);
    ASSERT_FLOAT_EQ(-r_cir,ed1.points()[5][0]);
    ASSERT_FLOAT_EQ( r_cir,ed1.points()[10][0]);
}

TEST(tests_mesh, tfi_blend_functions)
{
    std::vector<float> ksi_i_f = {0., 0.5, 1.};

    const auto n = 2;
    auto a_n = gbs::build_tfi_blend_function_with_derivatives<float, n>(ksi_i_f);
    for (auto d = 0; d < n; d++)
    {
        for (auto d_ = 0; d_ < n; d_++)
        {
            for (auto i = 0; i < ksi_i_f.size(); i++)
            {
                for (auto i_ = 0; i_ < ksi_i_f.size(); i_++)
                {
                    ASSERT_NEAR(
                        a_n[i][d].value(ksi_i_f[i_], d_),
                        gbs::kronecker<float>(i, i_) * gbs::kronecker<float>(d, d_),
                        1e-8
                        );
                }
            }
        }
    }
}

inline auto make_msh_crv()
{
    const size_t dim = 2;
    using T = double;
    using bsc2d = gbs::BSCurve<T,dim>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using points= std::vector<std::array<T,dim>>;
    using reals = std::vector<T>;
    using uints = std::vector<size_t>;

    reals knots{0., 1.};
    uints mults{3, 3};
    size_t p{2};

    auto u_crv1 = std::make_shared<bsc2d>(
        points{
            {0.00,0.50},
            {1.00,0.75},
            {1.00,1.50},
        },
        knots,
        mults,
        p
    );
    auto u_crv2 = std::make_shared<bsc2d>(
        points{
            {0.00, 1.00},
            {0.50, 1.00},
            {0.50, 1.50},
        },
        knots,
        mults,
        p
    );

    auto v_crv1 = std::make_shared<bsc2d>(
        points{
            {0.00,0.50},
            {0.10,0.75},
            // {0.00,0.75},
            {0.00,1.00},
        },
        knots,
        mults,
        p
    );
    auto v_crv2 = std::make_shared<bsc2d>(
        points{
            {1.00, 1.50},
            {0.75, 1.40},
            // {0.75, 1.50},
            {0.50, 1.50},
        },
        knots,
        mults,
        p
    );

    bsc2d_lst u_lst { u_crv1, u_crv2};
    bsc2d_lst v_lst { v_crv1, v_crv2};
    return std::make_tuple(u_lst, v_lst);
}

TEST(tests_mesh, msh_curves_set_sizes)
{
    const size_t dim = 2;
    using T = double;
    using bsc2d = gbs::BSCurve<T,dim>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using points= std::vector<std::array<T,dim>>;
    using reals = std::vector<T>;
    using uints = std::vector<size_t>;

    auto [u_lst, v_lst] = make_msh_crv();

    reals u{0.,  1.};
    reals v{0.,  1.};

    const size_t P = 1;
    const size_t Q = 1;

    size_t nu = 30;
    size_t nv = 15;
    auto nui = gbs::msh_curves_set_sizes(u_lst,u,nu);
    auto nvj = gbs::msh_curves_set_sizes(v_lst,v,nv);

    ASSERT_EQ(nui[0], nu);
    ASSERT_EQ(nvj[0], nv);

}
TEST(tests_mesh, msh_curves_set)
{
    const size_t dim = 2;
    using T = double;
    using bsc2d = gbs::BSCurve<T,dim>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using points= std::vector<std::array<T,dim>>;
    using reals = std::vector<T>;
    using uints = std::vector<size_t>;

    auto [u_lst, v_lst] = make_msh_crv();


    reals u{0.,  1.};
    reals v{0.,  1.};

    const size_t P = 1;
    const size_t Q = 1;

    size_t nu = 30;
    size_t nv = 15;
    auto nui = gbs::msh_curves_set_sizes(u_lst,u,nu);
    auto nvj = gbs::msh_curves_set_sizes(v_lst,v,nv);

    auto X_u = gbs::msh_curves_set<T, dim, P>(u_lst, nui, u);
    auto X_v = gbs::msh_curves_set<T, dim, Q>(v_lst, nvj, v);
    T tol = 1.e-6;

    for( size_t i {} ; i < u_lst.size(); i ++)
    {
        const auto &crv = *u_lst[i];
        for( const auto &pt : X_u )
        {
            auto d = gbs::extrema_curve_point<T,dim>(crv,pt[i][0],tol)[1];
            ASSERT_TRUE( d < tol );
        }
    }
    for( size_t i {} ; i < v_lst.size(); i ++)
    {
        const auto &crv = *v_lst[i];
        for( const auto &pt : X_v )
        {
            auto d = gbs::extrema_curve_point<T,dim>(crv,pt[i][0],tol)[1];
            ASSERT_TRUE( d < tol );
        }
    }   

}

TEST(tests_mesh, msh_curves_set_non_compatible_knots)
{
    size_t p = 3;
    using T = double;
    using namespace gbs;
    // Set 2 curves with diferent parametrization
    points_vector<T,2> crv1_pts{
        {0.00045, 0.01000}, // i0
        {0.00601, 0.01000},
        {0.02664, 0.01632},
        {0.03090, 0.02746},
        {0.03591, 0.03414},
        {0.03900, 0.03600}, //i1
        {0.05082, 0.04138},
        {0.06000, 0.04550}, //i2
        {0.06557, 0.04677},
        {0.07500, 0.04700}, //i3
        {0.07999, 0.04700},
        {0.10500, 0.04700}, //i4
        {0.11014, 0.04699},
        {0.13825, 0.04493},
        {0.16003, 0.03984}  // i5
    };
    auto crv1 = interpolate(
        crv1_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );

    points_vector<T,2> crv2_pts{
        {0.00010, 0.08550}, // i0
        {0.03900, 0.08550}, // i1
        {0.06000, 0.08550}, // i2
        {0.08000, 0.08550}, // i3
        {0.11000, 0.08550}, // i4
        {0.16000, 0.08550}  // i5
    };
    auto crv2 = interpolate(
        crv2_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );
    // set hard points
    std::vector<size_t> i_crv_1 {0, 5, 7, 9, 11, crv1_pts.size() - 1 };
    std::vector<size_t> i_crv_2 {0, 1, 2, 3, 4 , crv2_pts.size() - 1};
    // get hard points coordinates
    std::vector<std::vector<T>> u_lst(2);
    for( auto i : i_crv_1)
    {
        u_lst[0].push_back( extrema_curve_point(crv1, crv1_pts[i],1e-6)[0] );
    }
    for( auto i : i_crv_2)
    {
        u_lst[1].push_back( extrema_curve_point(crv2, crv2_pts[i],1e-6)[0] );
    }
    // build sizes
    std::vector<std::shared_ptr<Curve<T, 2>>> crv_lst{
        std::make_shared<BSCurve<T,2>>(crv1),
        std::make_shared<BSCurve<T,2>>(crv2)
    };
    size_t n{60};
    auto nui = msh_curves_set_sizes<T,2>( crv_lst, u_lst, n);

    auto X_ksi = msh_curves_set<T,2,1>(crv_lst, nui, u_lst);

    ASSERT_TRUE( X_ksi[0].size() == X_ksi[1].size() );

    auto ni = X_ksi.size();
    points_vector<T,2> crv1_msh( ni );
    points_vector<T,2> crv2_msh( ni );

    for(size_t i{}; i < ni ; i++)
    {
        crv1_msh[i] = X_ksi[i][0][0];
        crv2_msh[i] = X_ksi[i][1][0];
    }

    ASSERT_LT( norm( crv1_msh[ 0 ] -  crv1_pts[ 0 ]), 1e-7);
    ASSERT_LT( norm( crv2_msh[ 0 ] -  crv2_pts[ 0 ]), 1e-7);

    size_t i_msh{};
    for(size_t i{}; i < nui.size(); i++)
    {
        i_msh += nui[i] - 1;
        ASSERT_LT( norm( crv1_msh[ i_msh ] -  crv1_pts[ i_crv_1[i+1] ]), 1e-7);
        ASSERT_LT( norm( crv2_msh[ i_msh ] -  crv2_pts[ i_crv_2[i+1] ]), 1e-7);
    }

    // plot(
    //     crv1
    //     , crv2
    //     , crv1_msh
    //     , crv2_msh
    // );
}

TEST(tests_mesh, msh_curves_lattice_non_compatible_knots)
{
    size_t p = 3;
    using T = double;
    using namespace gbs;
    // Set 2 curves with diferent parametrization
    points_vector<T,2> crv1_pts{
        {0.00045, 0.01000}, // i0
        {0.00601, 0.01000},
        {0.02664, 0.01632},
        {0.03090, 0.02746},
        {0.03591, 0.03414},
        {0.03900, 0.03600}, //i1
        {0.05082, 0.04138},
        {0.06000, 0.04550}, //i2
        {0.06557, 0.04677},
        {0.07500, 0.04700}, //i3
        {0.07999, 0.04700},
        {0.10500, 0.04700}, //i4
        {0.11014, 0.04699},
        {0.13825, 0.04493},
        {0.16003, 0.03984}  // i5
    };
    auto crv1 = interpolate(
        crv1_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );

    points_vector<T,2> crv2_pts{
        {0.00010, 0.08550}, // i0
        {0.03900, 0.08550}, // i1
        {0.06000, 0.08550}, // i2
        {0.08000, 0.08550}, // i3
        {0.11000, 0.08550}, // i4
        {0.16000, 0.08550}  // i5
    };
    auto crv2 = interpolate(
        crv2_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );
    // set hard points
    std::vector<size_t> i_crv_1 {0, 5, 7, 9, 11, crv1_pts.size() - 1 };
    std::vector<size_t> i_crv_2 {0, 1, 2, 3, 4 , crv2_pts.size() - 1};
    // get hard points coordinates
    std::vector<std::vector<T>> u_lst(2);
    for( auto i : i_crv_1)
    {
        u_lst[0].push_back( extrema_curve_point(crv1, crv1_pts[i],1e-6)[0] );
    }
    for( auto i : i_crv_2)
    {
        u_lst[1].push_back( extrema_curve_point(crv2, crv2_pts[i],1e-6)[0] );
    }
    // Build opposite curves set
    std::vector<std::shared_ptr<Curve<T, 2>>> iso_u;
    auto n_hard_pts = i_crv_1.size();
    std::vector<std::vector<T>> v_lst;
    for(size_t i{}; i <n_hard_pts; i++)
    {
        iso_u.push_back(
            std::make_shared<BSCurve<T,2>>(build_segment(
                crv1_pts[i_crv_1[i]],
                crv2_pts[i_crv_2[i]]
            ))
        );
        v_lst.push_back(
            {
                iso_u.back()->bounds()[0],
                iso_u.back()->bounds()[1]
            }
        );
    }
    // build sizes
    std::vector<std::shared_ptr<Curve<T, 2>>> iso_v{
        std::make_shared<BSCurve<T,2>>(crv1),
        std::make_shared<BSCurve<T,2>>(crv2)
    };
    size_t ni{120};
    size_t nj{21};
    // size_t ni{40};
    // size_t nj{7};
    auto nvi = msh_curves_set_sizes<T,2>( iso_v, u_lst, ni);
    auto nui = msh_curves_set_sizes<T,2>( iso_u, v_lst, nj);

    const size_t P = 1;
    const size_t Q = 1;
    const bool slope_ctrl = true;

    auto [X_ksi, X_eth, X_ksi_eth, ksi, eth] = msh_curves_lattice<T,2,P,Q>(iso_u,iso_v,u_lst,v_lst,nui,nvi);
    auto pts = tfi_mesh_2d<T,2,P,Q,slope_ctrl>(
        X_ksi, 
        X_eth, 
        X_ksi_eth,
        make_range<T>(0,1,v_lst.size()),
        make_range<T>(0,1,u_lst.size()),
        ksi,
        eth
    );

    
    points_vector<T,2> eth_msh, ksi_msh, ksi_eth_msh;
    for(size_t i{}; i < X_eth.size() ; i++)
    {
        for(size_t j{}; j < X_eth[i].size() ; j++)
        {
            eth_msh.push_back( X_eth[i][j][0] );
        }
    }
    for(size_t i{}; i < X_ksi.size() ; i++)
    {
        for(size_t j{}; j < X_ksi[i].size() ; j++)
        {
            ksi_msh.push_back( X_ksi[i][j][0] );
        }
    }
    for(size_t i{}; i < X_ksi_eth.size() ; i++)
    {
        for(size_t j{}; j < X_ksi_eth[i].size() ; j++)
        {
            ksi_eth_msh.push_back( X_ksi_eth[i][j][0][0] );
        }
    }

    std::vector<size_t> vi(nvi.size()+1);
    vi[0] = 0;
    std:transform(
        nvi.begin(), nvi.end(),
        vi.begin(),
        std::next(vi.begin()),
        []( auto nvi_, auto i_prev ){return i_prev + nvi_-1;}
    );

    for( size_t i{1}; i < vi.size(); i++)
    {
        auto [it, err_max] = elliptic_structured_smoothing(pts,X_ksi.size(),vi[i-1],vi[i],0, nui[0]-1,100, 1e-5);
        printf("iterations : %i, error max: %.3e\n", int(it), err_max);
    }

    plot(
        pts,
        iso_u,
        iso_v
        // eth_msh,
        // ksi_msh,
        // ksi_eth_msh
    );
}

TEST(tests_mesh, msh_curves_lattice_non_compatible_knots_final)
{
    size_t p = 3;
    using T = double;
    using namespace gbs;
    // Set 2 curves with diferent parametrization
    points_vector<T,2> crv1_pts{
        {0.00045, 0.01000}, // i0
        {0.00600, 0.01000},
        {0.02500, 0.01500},
        {0.03000, 0.02000},
        {0.03500, 0.03200},
        {0.03900, 0.03600}, //i1
        {0.05000, 0.04200},
        {0.06000, 0.04500}, //i2
        {0.06600, 0.04800},
        {0.07500, 0.05000}, //i3
        {0.08000, 0.05000},
        {0.10500, 0.04800}, //i4
        {0.11000, 0.04700},
        {0.14000, 0.04500},
        {0.16000, 0.04000}  // i5
    };
    auto crv1 = interpolate(
        crv1_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );

    points_vector<T,2> crv2_pts{
        {0.00010, 0.08000}, // i0
        {0.03900, 0.09000}, // i1
        {0.06000, 0.08700}, // i2
        {0.08000, 0.08500}, // i3
        {0.11000, 0.08000}, // i4
        {0.16000, 0.08000}  // i5
    };
    auto crv2 = interpolate(
        crv2_pts ,   
        p,
        KnotsCalcMode::CHORD_LENGTH
    );
    // set hard points
    std::vector<size_t> i_crv_1 {0, 5, 7, 9, 11, crv1_pts.size() - 1 };
    std::vector<size_t> i_crv_2 {0, 1, 2, 3, 4 , crv2_pts.size() - 1};
    // Build opposite curves set
    std::vector<std::shared_ptr<Curve<T, 2>>> iso_u;
    auto n_hard_pts = i_crv_1.size();
    for(size_t i{}; i <n_hard_pts; i++)
    {
        iso_u.push_back(
            std::make_shared<BSCurve<T,2>>(build_segment(
                crv1_pts[i_crv_1[i]],
                crv2_pts[i_crv_2[i]]
            ))
        );
    }
    std::vector<std::shared_ptr<Curve<T, 2>>> iso_v{
        std::make_shared<BSCurve<T,2>>(crv1),
        std::make_shared<BSCurve<T,2>>(crv2)
    };

    size_t nu{120};
    size_t nv{21};
    // size_t nu{40};
    // size_t nv{11};
    const size_t P = 1;
    const size_t Q = 1;
    const bool slope_ctrl = true;

    auto [pts, nv_, nu_, n_iso_u, n_iso_v]  = tfi_mesh_2d<T,2,P,Q,slope_ctrl>(
        iso_u,
        iso_v, 
        nv,
        nu,
        1e-6
    );

    std::vector<size_t> vi(n_iso_v.size()+1);
    vi[0] = 0;
    std:transform(
        n_iso_v.begin(), n_iso_v.end(),
        vi.begin(),
        std::next(vi.begin()),
        []( auto nvi_, auto i_prev ){return i_prev + nvi_-1;}
    );

    for( size_t i{1}; i < vi.size(); i++)
    {
        auto [it, err_max] = elliptic_structured_smoothing(pts,nv_,vi[i-1],vi[i],0, n_iso_u[0]-1,100, 1e-5);
        printf("iterations : %i, error max: %.3e\n", int(it), err_max);
    }
    auto grid_actor = make_structuredgrid_actor(pts,nv_, nu_);
    plot(
        pts
        ,grid_actor
        ,iso_u
        ,iso_v
    );


}

inline auto make_msh_srf()
{
    const size_t dim = 2;
    using T = double;
    using bss2d = gbs::BSSurface<T,dim>;
    using reals = std::vector<T>;
    using points= std::vector<std::array<T,dim>>;

    reals knots{0., 0., 0., 1., 1., 1.};
    size_t p{2};

    return std::make_shared<bss2d>(
        points{
            {0.0,0.0}, {2.5,0.0}, {5.5,2.0},
            {1.0,2.0}, {3.0,2.0}, {5.5,3.5},
            {0.5,3.5}, {3.0,3.5}, {5.5,5.0},
        },
        knots,
        knots,
        p,
        p
    );

}

TEST(tests_mesh, tfi_mesh_2d_srf_opt)
{
    const size_t dim = 2;
    using T = double;
    using bss = gbs::BSSurface<T,dim>;
    using bsc = gbs::BSCurve<T,dim>;
    using bsc2d = gbs::BSCurve<T,2>;
    using reals = std::vector<T>;
    using point = std::array<T,dim>;
    using points= std::vector<point>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using uints = std::vector<size_t>;

    using gbs::operator+=;
    using gbs::operator+;
    using gbs::operator-;
    using gbs::operator*;

    reals knots{0., 0., 0., 1., 1., 1.};
    size_t p{2};

    auto srf = bss{
        points{
            {0.0,0.0}, {2.5,2.0}, {5.5,2.0},
            {1.0,2.0}, {3.0,2.0}, {5.5,3.5},
            {0.5,3.5}, {3.0,3.5}, {5.5,5.0},
        },
        knots,
        knots,
        p,
        p
    };

    // reals ksi_i{0.0,0.3, 0.5,0.8,1.};
    // reals eth_j{0.0,0.4,1.};
    // const size_t P = 2;
    // const size_t Q = 2;
    reals ksi_i{0.0, 1.};
    reals eth_j{0.0, 1.};
    const size_t P = 1;
    const size_t Q = 1;

    const bool slope_ctrl = true;

    size_t nu = 30;
    size_t nv = 20;
    auto p_srf = std::make_shared<bss>(srf);
    auto [ pts, ni, nj, n_iso_eth, n_iso_ksi ] = gbs::tfi_mesh_2d<T, dim, P, Q, slope_ctrl>(p_srf,ksi_i, eth_j, nu, nv);
    gbs::project_points_on_surface(*p_srf,pts);
    bsc2d_lst iso_eth, iso_ksi;
    for( auto eth : eth_j)
    {
        auto p_uv_crv = std::make_shared<bsc2d>( gbs::build_segment(point{ksi_i.front(),eth},point{ksi_i.back(),eth}) );
        iso_eth.push_back( std::make_shared<gbs::CurveOnSurface<T,dim>>(p_uv_crv, p_srf) );
    }
    for( auto ksi : ksi_i)
    {
        auto p_uv_crv = std::make_shared<bsc2d>( gbs::build_segment(point{ksi,eth_j.front()},point{ksi,eth_j.back()}) );
        iso_ksi.push_back( std::make_shared<gbs::CurveOnSurface<T,dim>>(p_uv_crv, p_srf) );
    }

    ASSERT_TRUE( check_curve_lattice(iso_ksi, iso_eth, ksi_i, eth_j, 1e-6) );

    gbs::plot( srf, pts);
    
    
}

TEST(tests_mesh, tfi_mesh_2d_no_hard_vtx_opt)
{
    const size_t dim = 2;
    using T = double;
    using bsc2d = gbs::BSCurve<T,dim>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using points= std::vector<std::array<T,dim>>;
    using reals = std::vector<T>;
    using uints = std::vector<size_t>;

    auto [iso_eth ,iso_ksi] = make_msh_crv();

    size_t nu = 30;
    size_t nv = 15;

    reals ksi_i{0.,  1.};
    reals eth_j{0.,  1.};

    auto ni = ksi_i.size();
    auto nj = eth_j.size();

    for(size_t i{}; i< ni; i++)
    {
        for(size_t j{}; j< nj; j++)
        {
            if(gbs::distance(iso_eth[j]->value(ksi_i[i]), iso_ksi[i]->value(eth_j[j])) > 1e-6)
            {
                auto pt1 = iso_eth[j]->value(ksi_i[i]);
                auto pt2 = iso_ksi[i]->value(eth_j[j]);
                auto k   = ksi_i[i];
                auto e   = eth_j[i];
            }
        }
    }

    const size_t P = 1;
    const size_t Q = 1;

    // uints n_iso_ksi{nu};
    // uints n_iso_eth{nv};
    auto n_iso_ksi = gbs::msh_curves_set_sizes(iso_eth,ksi_i,nv);
    auto n_iso_eth = gbs::msh_curves_set_sizes(iso_ksi,eth_j,nu);

    T tol = 1.e-6;

    auto [X_ksi, X_eth, X_ksi_eth, ksi, eth] = gbs::msh_curves_lattice<T,dim,P,Q>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth);
    auto pts = gbs::tfi_mesh_2d(X_ksi, X_eth, X_ksi_eth, ksi_i, eth_j, ksi, eth);
    size_t i{}, j{};

    for(i = 0 ; i < nv ; i ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nv * j],X_ksi[i][0][0]) < tol);
    }
    j = nu -1;
    for(i = 0 ; i < nv ; i ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nv * j],X_ksi[i][1][0]) < tol);
    }
    i = 0;
    for(j = 0 ; j < nu ; j ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nv * j],X_eth[j][0][0]) < tol);
    }
    i = nv-1;
    for(j = 0 ; j < nu ; j ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nv * j],X_eth[j][1][0]) < tol);
    }


    gbs::plot( iso_eth, iso_ksi, pts);
}
 