#pragma once
#include <vector>
#include <execution>
#include <gbs/gbslib.h>
#include <gbs/curvescheck.h>
#include <gbs/extrema.h>
#include <gbs/bscbuild.h>
namespace gbs
{

    template <typename T, size_t P>
    auto build_tfi_blend_function_with_derivatives(const std::vector<T> &ksi_i) -> std::vector<std::array<BSCfunction<T>, P>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<std::array<BSCfunction<T>, P>> alpha_i(n_ksi_i);

        for (int i = 0; i < n_ksi_i; i++)
        {
            for (auto n = 0; n < P; n++)
            {
                std::vector<constrType<T, 1, P + 1>> dji{n_ksi_i, {0.}};
                dji[i][n] = {1.};
                alpha_i[i][n] = BSCfunction(interpolate(dji, ksi_i));
            }
        }
        return alpha_i;
    }

    template <typename T, size_t P, bool slope_control>
    auto build_tfi_blend_function_with_derivatives(const std::vector<T> &ksi_i) -> std::vector<std::array<BSCfunction<T>, P>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<std::array<BSCfunction<T>, P>> alpha_i(n_ksi_i);

        for (int i = 0; i < n_ksi_i; i++)
        {
            for (auto n = 0; n < P; n++)
            {
                std::vector<constrType<T, 1, P + slope_control>> dji{n_ksi_i, {0.}};
                dji[i][n] = {1.};
                alpha_i[i][n] = BSCfunction(interpolate(dji, ksi_i));
            }
        }
        return alpha_i;
    }

    template <typename T>
    auto build_tfi_blend_function(const std::vector<T> &ksi_i, bool slope_control) -> std::vector<BSCfunction<T>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<BSCfunction<T>> alpha_i;

        for (int i = 0; i < n_ksi_i; i++)
        {
            if (slope_control)
            {
                std::vector<constrType<T, 1, 2>> dji{n_ksi_i, {0.}};
                dji[i][0] = {1.};
                alpha_i.push_back(interpolate(dji, ksi_i));
            }
            else
            {
                if (n_ksi_i == 2)
                {
                    std::vector<constrType<T, 1, 1>> dji{n_ksi_i, {0.}};
                    dji[i][0] = {1.};
                    alpha_i.push_back(interpolate(dji, ksi_i));
                }
                else
                {
                    points_vector<T, 1> dji{n_ksi_i, {0.}};
                    dji[i] = {1.};
                    alpha_i.push_back(interpolate(dji, ksi_i, 2));
                }
            }
        }
        return alpha_i;
    }

    template <typename T, size_t dim>
    inline auto plus(T l_, const std::shared_ptr<Curve<T, dim>> &crv)
    {
        return l_ + length(*crv);
    }

    template <typename T, size_t dim>
    inline auto plus(const std::shared_ptr<Curve<T, dim>> &crv, T l_)
    {
        return plus<T,dim>(l_, crv);
    }

    template <typename T, size_t dim>
    inline auto plus(const std::shared_ptr<Curve<T, dim>> &crv1, const std::shared_ptr<Curve<T, dim>> &crv2)
    {
        return length(*crv1) + length(*crv2);
    }

    template <typename T, size_t dim>
    inline auto plus(T l_, const std::shared_ptr<Curve<T, dim>> &crv, T u1, T u2)
    {
        return l_ + length(*crv, u1, u2);
    }

    template <typename T, size_t dim>
    inline auto plus(const std::shared_ptr<Curve<T, dim>> &crv, T l_, T u1, T u2)
    {
        return l_ + length(*crv, u1, u2);
    }

    template <typename T>
    inline auto plus(T l1, T l2)
    {
        return l1 + l2;
    }

    template <typename T>
    inline auto plus(T l1, T l2, T u1, T u2)
    {
        return l1 + l2;
    }

    template <typename T, size_t dim>
    inline auto plus(const std::shared_ptr<Curve<T, dim>> &crv1, const std::shared_ptr<Curve<T, dim>> &crv2,  T u1, T u2)
    {
        return length(*crv1, u1, u2) + length(*crv2, u1, u2);
    }

    template <typename T, size_t dim>
    auto msh_curves_set_sizes(
        const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst,
        const std::vector<T> &u,
        T dm
        ) -> std::vector<size_t>
    {
        std::vector<size_t> nui(u.size()-1);
        std::transform(
            u.begin(),
            std::next(u.end(), -1),
            std::next(u.begin(), 1),
            nui.begin(),
            [&crv_lst,dm](auto u1, auto u2)
            {
                auto l = std::reduce(
                    crv_lst.begin(),crv_lst.end(),T{},
                    // https://en.cppreference.com/w/cpp/algorithm/reduce
                    // The behavior is non-deterministic if binary_op is not associative or not commutative. 
                    // [u1, u2](T l_, const std::shared_ptr<Curve<T, dim>> &crv)
                    // {
                    //     return l_ + length(*crv, u1, u2);
                    // }
                    [u1, u2](auto l_, auto crv)
                    {
                        return plus(l_,crv,u1,u2);
                    }
                );
                l /= crv_lst.size();
                auto n = static_cast<size_t>( round(l / dm)) + 1; 
                return n;
            }
        );
        return nui;
    }

    template <typename T, size_t dim>
    auto msh_curves_set_sizes(
        const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst,
        const std::vector<T> &u,
        size_t n
        ) -> std::vector<size_t>
    {
        auto l = std::reduce(
            crv_lst.begin(),crv_lst.end(),T{},
            // https://en.cppreference.com/w/cpp/algorithm/reduce
            // The behavior is non-deterministic if binary_op is not associative or not commutative. 
            // [](T l_, const std::shared_ptr<Curve<T, dim>> &crv)
            // {
            //     return l_ + length(*crv);
            // }
            [](auto l_, auto crv)
            {
                return plus(l_, crv);
            }
        );
        l /= crv_lst.size();
        auto dm = l / ( n - 1 );
        return msh_curves_set_sizes(crv_lst, u, dm);
    }

    template <typename T, size_t dim, size_t P>
    auto msh_curves_set(
        const std::vector<std::shared_ptr<Curve<T, dim>>> &crv_lst,
        const std::vector<size_t> &nui,
        const std::vector<T> &u) -> std::vector<std::vector<std::array<point<T,dim> , P>>>
    {
        if(nui.size() != (u.size()-1))
        {
            throw std::out_of_range("Incorrect sizes." + std::to_string(nui.size()) + " " + std::to_string(u.size()-1));
        }

        if(!check_p_curves_bounds(crv_lst.begin(), crv_lst.end()))
        {
            throw std::invalid_argument("Curves must have the same bounds.");
        }

        size_t ni{crv_lst.size()};
        std::vector<std::vector<T>> params(ni);
        size_t nu{u.size()};
        for (size_t i{}; i < ni ; i++)
        {
            auto crv = crv_lst[i];
            for (size_t ui{}; ui < nu - 1 ; ui++) // mesh between hard points
            {
                // mesh segment
                auto u1 = u[ui];
                auto u2 = u[ui + 1];
                auto params_crv_seg = uniform_distrib_params(*crv, u1, u2, nui[ui]);
                if (ui != u.size() - 2) // pop end for all but last segment
                    params_crv_seg.pop_back();
                params[i].insert(params[i].end(), params_crv_seg.begin(), params_crv_seg.end());
            }
        }
        size_t ni_{params[0].size()};

        std::vector<std::vector<std::array<point<T,dim> , P>>> X_ksi(ni_);
        
        for( size_t i_{} ; i_ < ni_; i_++)
        {
            X_ksi[i_] = std::vector<std::array<point<T,dim> , P>>(ni);
            for (size_t i{}; i < ni ; i++)
            {
                auto crv = crv_lst[i];
                for(size_t n{} ; n < P; n++)
                {
                    X_ksi[i_][i][n] = crv->value(params[i][i_], n);
                }
            }
        }

        return X_ksi;
    }

    /**
     * @brief Check if lattice is compatible for TFI and Gordon surface
     * 
     * @tparam T 
     * @tparam dim 
     * @param iso_ksi 
     * @param iso_eth 
     * @param ksi_i 
     * @param eth_j 
     * @param tol 
     * @return true 
     * @return false 
     */
    template <typename T, size_t dim>
    auto check_curve_lattice(
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_ksi,
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_eth,
        const std::vector<T> &ksi_i,
        const std::vector<T> &eth_j,
        T tol
    ) -> bool
    {
        auto ni = ksi_i.size();
        auto nj = eth_j.size();
        for(size_t i{}; i< ni; i++)
        {
            for(size_t j{}; j< nj; j++)
            {
                if(distance(iso_eth[j]->value(ksi_i[i]),iso_ksi[i]->value(eth_j[j])) > tol)
                    return false;
            }
        }
        return true;
    }

    template <typename T, size_t dim, size_t P, size_t Q>
    auto msh_curves_lattice(
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_ksi,
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_eth,
        const std::vector<T> &ksi_i,
        const std::vector<T> &eth_j,
        const std::vector<size_t> &n_iso_ksi,
        const std::vector<size_t> &n_iso_eth,
        const std::shared_ptr<Surface<T,dim>> p_srf = nullptr
    )
    {
        auto X_ksi = msh_curves_set<T,dim,P>(iso_ksi, n_iso_ksi, eth_j);
        auto X_eth = msh_curves_set<T,dim,Q>(iso_eth, n_iso_eth, ksi_i);

        auto L   = ksi_i.size();
        auto M   = eth_j.size();

        std::vector<std::vector<std::array<std::array<point<T, dim>,Q>,P>>> X_ksi_eth(L);
        for(size_t i{};  i < L; i++)
        {
            X_ksi_eth[i] = std::vector<std::array<std::array<point<T, dim>,Q>,P>>(M);
            for(size_t j{};  j < M; j++)
            {
                for(size_t n{} ; n < P; n++)
                {
                    for(size_t m{}; m < Q ; m++)
                    {
                        if(p_srf)
                        {
                            X_ksi_eth[i][j][n][m] = p_srf->value(ksi_i[i],eth_j[j],m,n); // warning m n, why inverted
                        }
                        else
                        {
                            if(m==0)
                                X_ksi_eth[i][j][n][m] = iso_ksi[i]->value(eth_j[j],n);
                            else if(n==0)
                                X_ksi_eth[i][j][n][m] = iso_eth[j]->value(ksi_i[i],m);
                        }
                    }
                }
            }
        }
        std::vector<T> arr_i{0}, arr_j{0};
        for (auto index : n_iso_eth)
            arr_i.push_back(arr_i.back() + index - 1);
        for (auto index : n_iso_ksi)
            arr_j.push_back(arr_j.back() + index - 1);
        
        auto ksi = interpolate(ksi_i,arr_i,1);
        auto eth = interpolate(eth_j,arr_j,1);
        return std::make_tuple(X_ksi, X_eth, X_ksi_eth, ksi, eth);
    }

    template <typename T, size_t dim, size_t P, size_t Q, bool slope_ctrl = false>
    auto tfi_mesh_2d(
        const std::vector<std::vector<std::array<point<T,dim> , P>>> &X_ksi, 
        const std::vector<std::vector<std::array<point<T,dim> , Q>>> &X_eth, 
        const std::vector<std::vector<std::array<std::array<point<T,dim>,Q>,P>>> &X_ksi_eth,
        const std::vector<T> &ksi_i, 
        const std::vector<T> &eth_j,
        const BSCfunction<T> &ksi,
        const BSCfunction<T> &eth
    )
    {
        auto alpha_i = build_tfi_blend_function_with_derivatives<T, P, slope_ctrl>(ksi_i);
        auto beta_j = build_tfi_blend_function_with_derivatives<T, Q, slope_ctrl>(eth_j);
        auto ni_ = X_eth.size();
        auto nj_ = X_ksi.size();
        auto L   = ksi_i.size();
        auto M   = eth_j.size();
        
        points_vector<T, dim> pts(ni_ * nj_);
        auto i_range = make_range<size_t>(0, ni_ - 1);

        std::for_each(
            std::execution::par,
            i_range.begin(), i_range.end(),
            [&](size_t i_)
            {
                for (size_t j_{}; j_ < nj_; j_++)
                {
                    point<T, dim> U{}, V{}, UV{};
                    auto ksi_i_{ksi(i_)};
                    auto eth_j_{eth(j_)};
                    for(size_t i{} ; i < L; i++)
                    {
                        for(size_t n{} ; n < P; n++)
                        {
                            U += alpha_i[i][n](ksi_i_) * X_ksi[j_][i][n];
                        }
                    }
                    for(size_t j{} ; j < M; j++)
                    {
                        for(size_t m{} ; m < Q; m++)
                        {
                            V += beta_j[j][m](eth_j_) * X_eth[i_][j][m];
                        }
                    }
                    for(size_t i{} ; i < L; i++)
                    {
                        for(size_t j{} ; j < M; j++)
                        {
                            for(size_t n{} ; n < P; n++)
                            {
                                for(size_t m{} ; m < Q; m++)
                                {
                                    UV += alpha_i[i][n](ksi_i_) * beta_j[j][m](eth_j_) * X_ksi_eth[i][j][n][m];
                                }
                            }
                        }
                    }
                    pts[j_ + nj_ * i_] = U + V - UV;
                }
            }
        );
        return pts;
    }

    template <typename T, size_t dim>
    void project_points_on_surface(const Surface<T,dim> &srf, points_vector<T, dim> &pts, T tol_projection =1e-6 )
    {
        std::transform(
            std::execution::par,
            pts.begin(),
            pts.end(),
            pts.begin(),
            [&srf,tol_projection](const auto &pt)
            {
                auto [u,v, d] =  extrema_surf_pnt(srf,pt,tol_projection);
                return srf(u,v);
            }
        );
    }

    template <typename T, size_t dim, size_t P, size_t Q, bool slope_ctrl = false>
    auto tfi_mesh_2d( 
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_ksi,
        const std::vector<std::shared_ptr<Curve<T, dim>>> &iso_eth,
        const std::vector<T> &ksi_i ,
        const std::vector<T> &eth_j, 
        size_t n_ksi, 
        size_t n_eth,
        const std::shared_ptr<Surface<T,dim>> p_srf = nullptr
    )
    {
        auto n_iso_eth = msh_curves_set_sizes(iso_eth,ksi_i,n_eth);
        auto n_iso_ksi = msh_curves_set_sizes(iso_ksi,eth_j,n_ksi);
        auto [X_ksi, X_eth, X_ksi_eth, ksi, eth] = msh_curves_lattice<T,dim,P,Q>(iso_ksi, iso_eth, ksi_i, eth_j, n_iso_ksi, n_iso_eth, p_srf);
        return std::make_tuple(
            tfi_mesh_2d<T,dim,P,Q, slope_ctrl>(X_ksi, X_eth, X_ksi_eth, ksi_i, eth_j, ksi, eth),
            X_ksi.size(),
            X_eth.size(),
            n_iso_ksi,
            n_iso_eth
        );
    }


    template <typename T, size_t dim, size_t P, size_t Q, bool slope_ctrl = false>
    auto tfi_mesh_2d( const std::shared_ptr<Surface<T,dim>> &p_srf, const std::vector<T> &ksi_i ,const std::vector<T> &eth_j, size_t n_ksi, size_t n_eth)
    {
        std::vector<std::shared_ptr<Curve<T,dim>>> iso_eth, iso_ksi;
        for( auto eth : eth_j)
        {
            auto p_uv_crv = std::make_shared<BSCurve<T,2>>( build_segment(point{ksi_i.front(),eth},point{ksi_i.back(),eth}) );
            iso_eth.push_back( std::make_shared<CurveOnSurface<T,dim>>(p_uv_crv, p_srf) );
        }
        for( auto ksi : ksi_i)
        {
            auto p_uv_crv = std::make_shared<BSCurve<T,2>>( build_segment(point{ksi,eth_j.front()},point{ksi,eth_j.back()}) );
            iso_ksi.push_back( std::make_shared<CurveOnSurface<T,dim>>(p_uv_crv, p_srf) );
        }
        return tfi_mesh_2d<T, dim, P, Q, slope_ctrl>(iso_ksi, iso_eth, ksi_i, eth_j, n_ksi, n_eth, p_srf);
    }


}