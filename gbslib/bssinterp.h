#pragma once
#include <gbslib/bscinterp.h>
#include <gbslib/bssurf.h>
#include <exception>
namespace gbs
{
    template <typename T, size_t dim>
    auto build_poles(const std::vector<std::array<T, dim>> Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
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
        Eigen::MatrixX<T> N(n_poles, n_poles);

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

        Eigen::VectorX<T> b(n_poles);
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
} // namespace gbs