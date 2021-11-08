#include <gtest/gtest.h>
#include <gbs/bscbuild.h>
#include <gbs-mesh/mshedge.h>
#include <gbs-mesh/tfi.h>
#include <gbs-render/vtkcurvesrender.h>

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
                    ASSERT_FLOAT_EQ(
                        a_n[i][d].value(ksi_i_f[i_], d_),
                        gbs::kronecker<float>(i, i_) * gbs::kronecker<float>(d, d_));
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

TEST(tests_mesh, sweep_mesh)
{
    const size_t dim = 2;
    using T = double;
    using bsc2d = gbs::BSCurve<T,dim>;
    using bsc2d_lst = std::vector<std::shared_ptr<gbs::Curve<T,dim>>>;
    using points= std::vector<std::array<T,dim>>;
    using reals = std::vector<T>;
    using uints = std::vector<size_t>;

    auto [u_lst, v_lst] = make_msh_crv();

    size_t nu = 30;
    size_t nv = 15;

    reals u{0.,  1.};
    reals v{0.,  1.};

    const size_t P = 1;
    const size_t Q = 1;

    uints nui{nu};
    uints nvj{nv};

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

    reals ksi_i{0.0,0.3, 0.5,0.8,1.};
    reals eth_j{0.0,0.4,1.};
    uints n_iso_eth{10, 7, 7, 5};
    uints n_iso_ksi{6, 12 };
    const size_t P = 2;
    const size_t Q = 2;
    const bool slope_ctrl = true;
    // reals ksi_i{0.0, 1.};
    // reals eth_j{0.0, 1.};
    // uints n_iso_eth{30};
    // uints n_iso_ksi{15 };
    // const size_t P = 1;
    // const size_t Q = 1;

    auto L   = ksi_i.size();
    auto M   = eth_j.size();

    auto p_srf = std::make_shared<bss>(srf);

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

    auto [X_ksi, X_eth, X_ksi_eth, ksi, eth] = gbs::msh_curves_lattice<T,dim,P,Q>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth, p_srf);

    auto pts = gbs::tfi_mesh_2d<T,dim,P,Q, slope_ctrl>(X_ksi, X_eth, X_ksi_eth, ksi_i, eth_j, ksi, eth);


    // gbs::plot(srf, iso_eth, iso_ksi, X_eth[0], X_ksi[0], X_ksi_eth[0][0]);
    // gbs::plot(srf, iso_eth, iso_ksi, X_ksi_eth[0][0]);
    gbs::plot( iso_eth, iso_ksi, pts);
    
    
}

TEST(test_mesh, tfi_mesh_2d_no_hard_vtx_opt)
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

    uints n_iso_ksi{nu};
    uints n_iso_eth{nv};

    T tol = 1.e-6;

    auto [X_ksi, X_eth, X_ksi_eth, ksi, eth] = gbs::msh_curves_lattice<T,dim,P,Q>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth);
    auto pts = gbs::tfi_mesh_2d(X_ksi, X_eth, X_ksi_eth, ksi_i, eth_j, ksi, eth);
    size_t i{}, j{};

    for(i = 0 ; i < nu ; i ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nu * j],X_ksi[i][0][0]) < tol);
    }
    j = nv -1;
    for(i = 0 ; i < nu ; i ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nu * j],X_ksi[i][1][0]) < tol);
    }
    i = 0;
    for(j = 0 ; j < nv ; j ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nu * j],X_eth[j][0][0]) < tol);
    }
    i = nu-1;
    for(j = 0 ; j < nv ; j ++)
    {
        ASSERT_TRUE(gbs::distance(pts[i + nu * j],X_eth[j][1][0]) < tol);
    }


    gbs::plot( iso_eth, iso_ksi, pts);
}
 