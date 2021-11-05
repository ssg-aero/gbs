#pragma once
#include <vector>
#include <execution>
#include <gbs/gbslib.h>
#include <gbs/curvescheck.h>
namespace gbs
{

    template <typename T, size_t P>
    auto build_tfi_blend_function_with_derivatives(const std::vector<T> &ksi_i) -> std::vector<std::array<gbs::BSCfunction<T>, P>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<std::array<gbs::BSCfunction<T>, P>> alpha_i(n_ksi_i);

        for (int i = 0; i < n_ksi_i; i++)
        {
            for (auto n = 0; n < P; n++)
            {
                std::vector<gbs::constrType<T, 1, P + 1>> dji{n_ksi_i, {0.}};
                dji[i][n] = {1.};
                alpha_i[i][n] = gbs::BSCfunction(gbs::interpolate(dji, ksi_i));
            }
        }
        return alpha_i;
    }

    template <typename T, size_t P, bool slope_control>
    auto build_tfi_blend_function_with_derivatives(const std::vector<T> &ksi_i) -> std::vector<std::array<gbs::BSCfunction<T>, P>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<std::array<gbs::BSCfunction<T>, P>> alpha_i(n_ksi_i);

        for (int i = 0; i < n_ksi_i; i++)
        {
            for (auto n = 0; n < P; n++)
            {
                std::vector<gbs::constrType<T, 1, P + slope_control>> dji{n_ksi_i, {0.}};
                dji[i][n] = {1.};
                alpha_i[i][n] = gbs::BSCfunction(gbs::interpolate(dji, ksi_i));
            }
        }
        return alpha_i;
    }

    template <typename T>
    auto build_tfi_blend_function(const std::vector<T> &ksi_i, bool slope_control) -> std::vector<gbs::BSCfunction<T>>
    {
        auto n_ksi_i = ksi_i.size();
        std::vector<gbs::BSCfunction<T>> alpha_i;

        for (int i = 0; i < n_ksi_i; i++)
        {
            if (slope_control)
            {
                std::vector<gbs::constrType<T, 1, 2>> dji{n_ksi_i, {0.}};
                dji[i][0] = {1.};
                alpha_i.push_back(gbs::interpolate(dji, ksi_i));
            }
            else
            {
                if (n_ksi_i == 2)
                {
                    std::vector<gbs::constrType<T, 1, 1>> dji{n_ksi_i, {0.}};
                    dji[i][0] = {1.};
                    alpha_i.push_back(gbs::interpolate(dji, ksi_i));
                }
                else
                {
                    gbs::points_vector<T, 1> dji{n_ksi_i, {0.}};
                    dji[i] = {1.};
                    alpha_i.push_back(gbs::interpolate(dji, ksi_i, 2));
                }
            }
        }
        return alpha_i;
    }

    template <typename T, size_t dim, size_t P>
    auto msh_curves_set(
        const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &crv_lst,
        const std::vector<size_t> &nui,
        const std::vector<T> &u) -> std::vector<std::vector<std::array<gbs::point<T,dim> , P>>>
    {
        if(nui.size() != (u.size()-1))
        {
            throw std::exception("Incorrect sizes.");
        }

        if(!check_p_curves_bounds(crv_lst.begin(), crv_lst.end()))
        {
            throw std::exception("Curves must have the same bounds.");
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
                auto params_crv_seg = gbs::uniform_distrib_params(*crv, u1, u2, nui[ui]);
                if (ui != u.size() - 2) // pop end for all but last segment
                    params_crv_seg.pop_back();
                params[i].insert(params[i].end(), params_crv_seg.begin(), params_crv_seg.end());
            }
        }
        size_t ni_{params[0].size()};

        std::vector<std::vector<std::array<gbs::point<T,dim> , P>>> X_ksi(ni_);
        
        for( size_t i_{} ; i_ < ni_; i_++)
        {
            X_ksi[i_] = std::vector<std::array<gbs::point<T,dim> , P>>(ni);
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
        const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &iso_ksi,
        const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &iso_eth,
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
                if(gbs::distance(iso_eth[j]->value(ksi_i[i]),iso_ksi[i]->value(eth_j[j])) > tol)
                    return false;
            }
        }
        return true;
    }

    template <typename T, size_t dim, size_t P, size_t Q>
    auto msh_curves_lattice(
        const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &iso_ksi,
        const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &iso_eth,
        const std::vector<T> &ksi_i,
        const std::vector<T> &eth_j,
        const std::vector<size_t> &n_iso_ksi,
        const std::vector<size_t> &n_iso_eth,
        const std::shared_ptr<Surface<T,dim>> p_srf = nullptr
    )
    {
        auto X_ksi = gbs::msh_curves_set<T,dim,P>(iso_ksi, n_iso_ksi, eth_j);
        auto X_eth = gbs::msh_curves_set<T,dim,Q>(iso_eth, n_iso_eth, ksi_i);

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
        
        auto ksi = gbs::interpolate(ksi_i,arr_i,1);
        auto eth = gbs::interpolate(eth_j,arr_j,1);
        return std::make_tuple(X_ksi, X_eth, X_ksi_eth, ksi, eth);
    }

    template <typename T, size_t dim, size_t P, size_t Q, bool slope_ctrl = false>
    auto tfi_mesh_2d(
        const std::vector<std::vector<std::array<gbs::point<T,dim> , P>>> &X_ksi, 
        const std::vector<std::vector<std::array<gbs::point<T,dim> , Q>>> &X_eth, 
        const std::vector<std::vector<std::array<std::array<point<T,dim>,Q>,P>>> &X_ksi_eth,
        const std::vector<T> &ksi_i, 
        const std::vector<T> &eth_j,
        const auto &ksi,
        const auto &eth
    )
    {
        auto alpha_i = gbs::build_tfi_blend_function_with_derivatives<T, P, slope_ctrl>(ksi_i);
        auto beta_j = gbs::build_tfi_blend_function_with_derivatives<T, Q, slope_ctrl>(eth_j);
        auto ni_ = X_eth.size();
        auto nj_ = X_ksi.size();
        auto L   = ksi_i.size();
        auto M   = eth_j.size();
        
        gbs::points_vector<T, dim> pts(ni_ * nj_);
        auto i_range = gbs::make_range<size_t>(0, ni_ - 1);

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

}