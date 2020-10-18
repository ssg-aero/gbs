#pragma once
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <exception>
namespace gbs
{

    template <typename Container>
    auto extract_U(size_t i_V,const Container &Q, size_t n_col) -> Container
    {
        if (Q.size() % n_col)
        {
            std::exception("size error");
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
            std::exception("size error");
        }
        auto n_row = Q.size() / n_col;
        Container col( n_row );
        for(auto j = 0 ; j < n_row; j++)
        {
            col[j] = Q[i_U+n_col*j];
        }
        return col;
    }

    template <typename T, size_t dim>
    auto build_poles(const std::vector<std::array<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        auto n_pt = Q.size();
        auto n_poles = Q.size();
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = n_params_u; //ajouter nc
        auto n_poles_v = n_params_v; //ajouter nc
        if (n_poles != n_params_u * n_params_v)
        {
            std::exception("size error");
        }
        MatrixX<T> N(n_poles, n_poles);

        size_t iRow, iCol;
        for (size_t iv = 0; iv < n_params_v; iv++)
        {
            for (size_t iu = 0; iu < n_params_u; iu++)
            {
                iRow = iu + n_params_u * iv;
                for (size_t j = 0; j < n_poles_v; j++)
                {
                    auto Nv = gbs::basis_function(v[iv], std::next(k_flat_v.begin(), j), q, k_flat_v.end());
                    for (size_t i = 0; i < n_poles_u; i++)
                    {
                        iCol = i + n_poles_u * j;
                        N(iRow, iCol) = gbs::basis_function(u[iu], std::next(k_flat_u.begin(), i), p, k_flat_u.end()) * Nv;
                    }
                }
            }
        }

        auto N_inv = N.partialPivLu(); //TODO solve banded system

        VectorX<T> b(n_poles);
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
    auto interpolate(const std::vector<std::array<T, dim>> &Q, size_t n_poles_v, size_t p, size_t q, gbs::KnotsCalcMode mode) -> gbs::BSSurface<T, dim>
    {

        if (Q.size() % n_poles_v)
        {
            std::exception("size error");
        }

        size_t n_poles_u = Q.size() / n_poles_v;

        auto avg_p = [&](auto extract_f, size_t ni, size_t nj) mutable {
            
            std::vector<T> params(nj , T(0));
            for (auto i = 0; i < ni; i++)
            {
                auto pts = extract_f(i);
                auto params_i = gbs::curve_parametrization(pts, mode, true);
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
        auto ku = build_simple_mult_flat_knots<T>(u, n_poles_u, p);
        auto v = avg_p(
            [&](size_t i) {
                return extract_V(i, Q, n_poles_u);
            },
            n_poles_u, n_poles_v);
        adimension(v);
        auto kv = build_simple_mult_flat_knots<T>(v, n_poles_v, q);
        
        auto poles = build_poles(Q,ku,kv,u,v,p,q);

        return gbs::BSSurface<T, dim>(poles,ku,kv,p,q);

    }
} // namespace gbs