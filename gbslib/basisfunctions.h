#pragma once
#include <vector>
#include <utility>
#include <gbslib/gbslib.h>
#include <execution>
#include <algorithm>

namespace gbs
{
    template <class _InIt, typename T>
    T basis_function(T u, const _InIt &it, size_t p, const _InIt &_Last)
    {
        auto u_last = *std::next(_Last, -1);

        T ui = *it;
        T ui1 = *std::next(it, 1);
        if (p == 0)
        {
            return ((ui <= u) && (u < ui1)) || (fabs(ui1 - u_last) < knot_eps && fabs(u - u_last) < knot_eps)
                       ? T(1.)
                       : T(0.);
        }
        else
        {
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps)
            {
                C1 = (u - ui) / C1;
                C1 *= basis_function(u, it, p - 1, _Last);
            }
            if (C2 > knot_eps)
            {
                C2 = (ui1p - u) / C2;
                C2 *= basis_function(u, std::next(it, 1), p - 1, _Last);
            }
            return C1 + C2;
        }
    }

    template <template <class> class Container, class T>
    T basis_function(T u, size_t i, size_t p, const Container<T> &k)
    {
        return basis_function(u, std::next(k.begin(), i), p, k.end());
    }

    template <class _InIt, typename T>
    T basis_function(T u, const _InIt &it, size_t p, size_t d, const _InIt &_Last)
    {
        if (d == 0)
        {
            return basis_function(u, it, p, _Last);
        }
        else if (d >= p)
        {
            return 0.;
        }
        else
        {
            T ui = *it;
            T ui1 = *std::next(it, 1);
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps)
            {
                C1 = basis_function(u, it, p - 1, d - 1, _Last) / C1;
            }
            if (C2 > knot_eps)
            {
                C2 = basis_function(u, std::next(it, 1), p - 1, d - 1, _Last) / C2;
            }
            return p * (C1 - C2);
        }
    }

    template <typename T>
    auto find_span(size_t n, size_t p, T u, const std::vector<T> &U)
    {
        return --std::lower_bound(std::next(U.begin(), p), std::next(U.begin(), n + 1), u);
    }
    
    template <typename T>
    auto find_span2(size_t n, size_t p, T u, const std::vector<T> &U)
    {
        if (u >= U[n + 1])
            return std::next(U.begin(), n);
        auto low = p, high = n + 1, mid = (low + high) / 2;
        while (u < U[mid] || u >= U[mid + 1])
        {
            if (u < U[mid])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }
        return std::next(U.begin(), mid);
    }

    // template <template <class> class Container, class T>
    // T basis_function(T u, size_t i, size_t p, size_t d, const Container<T> &k)
    template <typename T>
    T basis_function(T u, size_t i, size_t p, size_t d, const std::vector<T> &k)
    {
        return basis_function(u, std::next(k.begin(), i), p, d, k.end());
    }
    template <class T, size_t dim>
    std::array<T, dim> eval_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p, size_t d = 0)
    //template <template <class> class Container, class T, size_t dim>
    //std::array<T, dim> eval_value_simple(const T &u, const Container<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p,size_t d=0)
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        auto i_max = std::min(find_span(n_poles,p,u,k) - k.begin() + d , k.size()-p-2) ;
        // printf("d: %zi\n",d);
        // /*
        for (size_t i = std::max(int(0),int(i_max-p)); i <= i_max; i++)
        // for (size_t i = 0; i < n_poles; i++) //Reducing span for few pole makes things worst
        {
            auto N = basis_function(u, i, p, d, k);
            for (size_t d_ = 0; d_ < dim; d_++)
            {
                pt[d_] += N * poles[i][d_];
            }
        }
        return pt;
    }
    template <class T, size_t dim>
    std::array<T, dim> eval_value_simple(T u,T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t du = 0, size_t dv = 0)
        {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_polesU = ku.size()-p-1;
        size_t n_polesV = kv.size()-q-1;

        auto i_max = std::min(find_span(n_polesU,p,u,ku) - ku.begin() + du , ku.size()-p-2) ;
        auto j_max = std::min(find_span(n_polesV,q,v,kv) - kv.begin() + dv , kv.size()-q-2) ;
        T Nu,Nv;

            // for (size_t j = 0; j < n_polesV; j++)
        for (size_t j = std::max(int(0), int(j_max - q)); j <= j_max; j++)
        {
            Nv = basis_function(v, j, q, dv, kv);
            for (size_t i = std::max(int(0), int(i_max - p)); i <= i_max; i++)
            {
                Nu = basis_function(u, i, p, du, ku);
                for (size_t d_ = 0; d_ < dim; d_++)
                {
                    pt[d_] += Nu * Nv * poles[i + n_polesU * j][d_];
                }
            }
        }
        return pt;
    }
} // namespace gbs