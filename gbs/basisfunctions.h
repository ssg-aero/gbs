#pragma once
#include <vector>
#include <utility>
#include <gbs/gbslib.h>
#include <execution>
#include <algorithm>

namespace gbs
{
    template <class _InIt, typename T>
    T BasisFunction(T u, const _InIt &it, size_t p, const _InIt &_Last)
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
                C1 *= BasisFunction(u, it, p - 1, _Last);
            }
            if (C2 > knot_eps)
            {
                C2 = (ui1p - u) / C2;
                C2 *= BasisFunction(u, std::next(it, 1), p - 1, _Last);
            }
            return C1 + C2;
        }
    }

    template <template <class> class Container, class T>
    T BasisFunction(T u, size_t i, size_t p, const Container<T> &k)
    {
        return BasisFunction(u, std::next(k.begin(), i), p, k.end());
    }

    template <class _InIt, typename T>
    T BasisFunction(T u, const _InIt &it, size_t p, size_t d, const _InIt &_Last)
    {
        if (d == 0)
        {
            return BasisFunction(u, it, p, _Last);
        }
        else if (d > p)
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
                C1 = BasisFunction(u, it, p - 1, d - 1, _Last) / C1;
            }
            if (C2 > knot_eps)
            {
                C2 = BasisFunction(u, std::next(it, 1), p - 1, d - 1, _Last) / C2;
            }
            return p * (C1 - C2);
        }
    }

    template <typename T>
    auto FindSpan(size_t n, size_t p, T u, const std::vector<T> &U)
    {
        return --std::lower_bound(std::next(U.begin(), p), std::next(U.begin(), n + 1), u);
    }
    
    template <typename T>
    auto FindSpan2(size_t n, size_t p, T u, const std::vector<T> &U)
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
    // T BasisFunction(T u, size_t i, size_t p, size_t d, const Container<T> &k)
    template <typename T>
    T BasisFunction(T u, size_t i, size_t p, size_t d, const std::vector<T> &k)
    {
        return BasisFunction(u, std::next(k.begin(), i), p, d, k.end());
    }
    template <class T, size_t dim>
    std::array<T, dim> EvalValueRecursive(T u, const std::vector<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p, size_t d = 0,bool use_span_reduction =true)
    //template <template <class> class Container, class T, size_t dim>
    //std::array<T, dim> EvalValueRecursive(const T &u, const Container<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p,size_t d=0)
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        size_t i_max{n_poles-1}, i_min{0};
        if (use_span_reduction)//Reducing span for few pole makes things worst
        {
            i_max = FindSpan(n_poles, p, u, k) - k.begin();
            i_max = std::min(i_max, k.size() - p - 2);
            i_min = std::max(int(0),int(i_max-p));
        }
        // printf("d: %zi\n",d);
        // /*
        for (auto i = i_min; i <= i_max; i++)
        {
            auto N = BasisFunction(u, i, p, d, k);
            std::transform(poles[i].begin(), poles[i].end(), pt.begin(), pt.begin(), [&](const auto val_, const auto tot_) { return tot_ + N * val_; });
        }
        return pt;
    }

    template <class T, size_t dim>
    std::array<T, dim> EvalValueRecursive(T u, T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t du = 0, size_t dv = 0)
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_polesU = ku.size() - p - 1;
        size_t n_polesV = kv.size() - q - 1;

        size_t i_max = FindSpan(n_polesU, p, u, ku) - ku.begin();
        size_t j_max = FindSpan(n_polesV, q, v, kv) - kv.begin();
        i_max = std::min( i_max, ku.size() - p - 2);
        j_max = std::min( j_max, kv.size() - q - 2);

        T Nu, Nv;

        // for (size_t j = 0; j < n_polesV; j++)
        for (size_t j = std::max(int(0), int(j_max - q)); j <= j_max; j++)
        {
            Nv = BasisFunction(v, j, q, dv, kv);
            // for (size_t i = 0; i < n_polesU; i++)
            for (size_t i = std::max(int(0), int(i_max - p)); i <= i_max; i++)
            {
                Nu = BasisFunction(u, i, p, du, ku);
                
                std::transform(poles[i + n_polesU * j].begin(),poles[i + n_polesU * j].end(),pt.begin(),pt.begin(),[&](const auto val_,const auto tot_){return tot_+Nu*Nv*val_;});
            }
        }
        return pt;
    }
} // namespace gbs