#pragma once
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <exception>
namespace gbs
{

    // ( u, v, x, du, dv)
    template < typename T, size_t dim>
    using bss_constraint = std::tuple<T,T,point<T,dim>,size_t,size_t>; 

    template < typename T, size_t dim>
    auto get_constraints_bounds(const std::vector<bss_constraint<T, dim>> &Q)
    {
        auto u_min = std::get<0>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<0>(Q1) < std::get<0>(Q2);}));
        auto u_max = std::get<0>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<0>(Q1) > std::get<0>(Q2);}));
        auto v_min = std::get<1>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<1>(Q1) < std::get<1>(Q2);}));
        auto v_max = std::get<1>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<1>(Q1) > std::get<1>(Q2);}));
        return std::array<T,4>{u_min, u_max, v_min, v_max};
    }

    template <typename Container>
    auto extract_U(size_t i_V,const Container &Q, size_t n_col) -> Container
    {
        if (Q.size() % n_col)
        {
            throw std::length_error("size error");
        }
        Container row( n_col );
        std::copy(std::next(Q.begin(),i_V*n_col),std::next(Q.begin(),(i_V+1)*n_col),row.begin());
        return row;
    }

    template <typename Container>
    auto extract_V(size_t i_U,const Container &Q, size_t n_col) -> Container
    {
        if (Q.size() % n_col)
        {
            throw std::length_error("size error");
        }
        auto n_row = Q.size() / n_col;
        Container col( n_row );
        for(auto j = 0 ; j < n_row; j++)
        {
            col[j] = Q[i_U+n_col*j];
        }
        return col;
    }

    template <typename T>
    auto fill_poles_matrix(MatrixX<T> &N,const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q)
    {
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = k_flat_u.size() - (p+1); 
        auto n_poles_v = k_flat_v.size() - (q+1); 
        size_t iRow, iCol;
        for (size_t iv = 0; iv < n_params_v; iv++)
        {
            for (size_t iu = 0; iu < n_params_u; iu++)
            {
                iRow = iu + n_params_u * iv;
                for (size_t j = 0; j < n_poles_v; j++)
                {
                    auto Nv = basis_function(v[iv], std::next(k_flat_v.begin(), j), q, k_flat_v.end());
                    for (size_t i = 0; i < n_poles_u; i++)
                    {
                        iCol = i + n_poles_u * j;
                        N(iRow, iCol) = basis_function(u[iu], std::next(k_flat_u.begin(), i), p, k_flat_u.end()) * Nv;
                    }
                }
            }
        }
    }

    template <typename T, size_t dim>
    auto build_poles(const auto &N_inv, const std::vector<std::array<T, dim>> &Q)
    {
        auto n_poles = N_inv.cols();
        auto n_pt = Q.size();
        VectorX<T> b(n_pt);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt; i++)
            {
                b(i) = Q[i][d]; //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return poles;
    }

    template <typename T, size_t dim>
    auto build_poles(const std::vector<std::array<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        auto n_poles = Q.size();
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = n_params_u;
        auto n_poles_v = n_params_v;
        if (n_poles != n_params_u * n_params_v)
        {
            throw std::length_error("size error");
        }
        MatrixX<T> N(n_poles, n_poles);
        fill_poles_matrix(N,k_flat_u,k_flat_v,u,v,p,q);

        return build_poles(N.partialPivLu(),Q);//TODO solve banded system

    }

    template <typename T, size_t dim>
    auto interpolate(const std::vector<std::array<T, dim>> &Q, size_t n_poles_v, size_t p, size_t q, KnotsCalcMode mode) -> BSSurface<T, dim>
    {

        if (Q.size() % n_poles_v)
        {
            throw std::length_error("size error");
        }

        size_t n_poles_u = Q.size() / n_poles_v;

        auto avg_p = [&](auto extract_f, size_t ni, size_t nj) mutable {
            
            std::vector<T> params(nj , T(0));
            for (auto i = 0; i < ni; i++)
            {
                auto pts = extract_f(i);
                auto params_i = curve_parametrization(pts, mode, true);
                std::transform(
                    std::execution::par,
                    params_i.begin(),
                    params_i.end(),
                    params.begin(),
                    params.begin(),
                    [&](const auto p_new, const auto p_) {
                        return p_ + p_new / ni;
                    });
            }
            return params;
        };

        auto u = avg_p(
            [&](size_t i) {
                return extract_U(i, Q, n_poles_u);
            },
            n_poles_v, n_poles_u);
        adimension(u);
        auto ku = build_simple_mult_flat_knots<T>(u, p);
        auto v = avg_p(
            [&](size_t i) {
                return extract_V(i, Q, n_poles_u);
            },
            n_poles_u, n_poles_v);
        adimension(v);
        auto kv = build_simple_mult_flat_knots<T>(v, q);
        
        return BSSurface<T, dim>(build_poles(Q,ku,kv,u,v,p,q),ku,kv,p,q);

    }
    /**
     * @brief Internal use for template <typename T, size_t dim>
     * auto interpolate(const std::vector<bss_constraint<T,dim>> &Q, size_t n_polesv, size_t p, size_t q, const std::array<T,4> &bounds = {0.,1.,0.,1.}, const std::array<size_t,2> &mults = {1,1} ) -> BSSurface<T, dim>
     * 
     * Use it at your own risks!
     * 
     * @tparam T 
     * @tparam dim 
     * @param Q 
     * @param k_flat_u 
     * @param k_flat_v 
     * @param p 
     * @param q 
     * @return std::vector<std::array<T, dim>> 
     */
    template <typename T, size_t dim>
    auto build_poles(const std::vector<bss_constraint<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        // TODO sort Q by contrain order to get a block systen
        auto n_poles = Q.size();
        auto n = k_flat_u.size() - p - 1;
        auto m = k_flat_v.size() - q - 1;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);
        for (int k = 0; k < n_poles; k++) // Eigen::ColMajor is default
        {
            auto [u, v, x, du, dv] = Q[k];
            for (int j{}; j < m; j++)
            {
                for (int i{}; i < n; i++)
                {
                    auto l = i + n * j;
                    N(k, l) = basis_function(u, i, p, du, k_flat_u) *
                              basis_function(v, j, q, dv, k_flat_v);
                }
            }
        }

        auto N_inv = N.partialPivLu();

        VectorX<T> b(n_poles);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_poles; i++)
            {
                b(i) = std::get<2>(Q[i])[d];
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return poles;
    }
    /**
     * @brief Build BSSurface<T,dim> satisfying the constraints vector provided: 
     *          * constraints parameters u,v are within bounds
     *          * Q.size() % n_polesv == 0 ( constraints have to be properly shaped aka Q.size() = n_polesv* n_polesv)
     *          * (n_polesu-p-1) % mu == 0 and (n_polesv-q-1) % mv == 0
     * 
     * @tparam T 
     * @tparam dim 
     * @param Q 
     * @param n_polesv 
     * @param p : degree U
     * @param q : degree V
     * @param bounds 
     * @param mults 
     * @return BSSurface<T, dim> 
     */
    template <typename T, size_t dim>
    auto interpolate(const std::vector<bss_constraint<T,dim>> &Q, size_t n_polesv, size_t p, size_t q, const std::array<T,4> &bounds = {0.,1.,0.,1.}, const std::array<size_t,2> &mults = {1,1} ) -> BSSurface<T, dim>
    {
        auto [u1, u2, v1, v2] = bounds;
        auto [mu, mv        ] = mults;

        auto [u_min, u_max, v_min, v_max] = get_constraints_bounds(Q);
        if(u_min < u1 || u_max > u2 || v_min < v1 || v_max > v2)
        {
            throw std::invalid_argument("constains must be specified within bounds");
        }

        if( Q.size() % n_polesv !=0)
        {
            throw std::length_error("interpolate, Q.size() = n_polesv * n_polesu required");
        }
        auto n_polesu = Q.size() / n_polesv;
        auto ku = build_mult_flat_knots<T>(u1, u2, 2 + (n_polesu-p-1)/mu, p, mu);
        auto kv = build_mult_flat_knots<T>(v1, v2, 2 + (n_polesv-q-1)/mv, q, mv);
        if(ku.size() - p - 1 != n_polesu)
        {
            throw std::invalid_argument("(n_polesu-p-1) % mu == 0 required");
        }
        if(kv.size() - q - 1 != n_polesv)
        {
            throw std::invalid_argument("(n_polesv-q-1) % mv == 0 required");
        }
        return BSSurface<T, dim>(build_poles(Q,ku,kv,p,q),ku,kv,p,q);

    }
} // namespace gbs