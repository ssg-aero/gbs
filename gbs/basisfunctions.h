#pragma once
#include <vector>
#include <utility>
#include <gbs/gbslib.h>
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
        auto pos = std::lower_bound(std::next(U.begin(), p), std::next(U.begin(), n + 1), u);
        if( pos != U.begin())
            return --pos;
        else
            return pos;
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
    std::array<T, dim> eval_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p, size_t d = 0,bool use_span_reduction =true)
    //template <template <class> class Container, class T, size_t dim>
    //std::array<T, dim> eval_value_simple(const T &u, const Container<T> &k, const std::vector<std::array<T, dim>> &poles, size_t p,size_t d=0)
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        size_t i_max{n_poles-1}, i_min{0};
        if (use_span_reduction)//Reducing span for few pole makes things worst
        {
            i_max = find_span(n_poles, p, u, k) - k.begin();
            i_max = std::min(i_max, k.size() - p - 2);
            i_min = std::max(int(0),int(i_max-p));
        }
        // printf("d: %zi\n",d);
        // /*
        for (auto i = i_min; i <= i_max; i++)
        {
            auto N = basis_function(u, i, p, d, k);
            std::transform(poles[i].begin(), poles[i].end(), pt.begin(), pt.begin(), [&](const auto val_, const auto tot_) { return tot_ + N * val_; });
        }
        return pt;
    }

    template <class T, size_t dim>
    std::array<T, dim> eval_value_simple(T u, T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t du = 0, size_t dv = 0)
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_polesU = ku.size() - p - 1;
        size_t n_polesV = kv.size() - q - 1;

        size_t i_max = find_span(n_polesU, p, u, ku) - ku.begin();
        size_t j_max = find_span(n_polesV, q, v, kv) - kv.begin();
        i_max = std::min( i_max, ku.size() - p - 2);
        j_max = std::min( j_max, kv.size() - q - 2);

        T Nu, Nv;

        // for (size_t j = 0; j < n_polesV; j++)
        for (size_t j = std::max(int(0), int(j_max - q)); j <= j_max; j++)
        {
            Nv = basis_function(v, j, q, dv, kv);
            // for (size_t i = 0; i < n_polesU; i++)
            for (size_t i = std::max(int(0), int(i_max - p)); i <= i_max; i++)
            {
                Nu = basis_function(u, i, p, du, ku);
                
                std::transform(poles[i + n_polesU * j].begin(),poles[i + n_polesU * j].end(),pt.begin(),pt.begin(),[&](const auto val_,const auto tot_){return tot_+Nu*Nv*val_;});
            }
        }
        return pt;
    }

    /**
     * @brief Utility function, add weights to poles array
     * 
     * @tparam T      : precision of curves
     * @tparam dim    : space dimension of curve (aka 1D, 2D, 3D,...)
     * @param poles   : non rational poles
     * @param weights : weights
     * @return std::vector<std::array<T, dim + 1>> : the unified poles
     */
    template <typename T, size_t dim>
    auto merge_weights(const std::vector<std::array<T, dim>> &poles, const std::vector<T> &weights) -> std::vector<std::array<T, dim + 1>>
    {
        std::vector<std::array<T, dim + 1>> p(poles.size());
        std::transform(
            std::execution::par,
            poles.begin(), poles.end(), weights.begin(), p.begin(),
            [](const auto &v, const auto &c) {
                std::array<T, dim + 1> n;
                // std::copy(v.begin(), v.end(), n.begin());
                std::transform(v.begin(), v.end(), n.begin(),[&](const auto &v_){return v_*c;});
                n.back() = c;
            });
        return p;
    }
    /**
     * @brief Project point with rational definition
     * 
     * @tparam T : precision of curves
     * @tparam dim : space dimension of curve (aka 1D, 2D, 3D,...)
     * @param pt : point coordinates to project
     * @return std::array<T, dim - 1> 
     */
    template <typename T, size_t dim>
    auto weight_projection(const std::array<T, dim> &pt) -> std::array<T, dim - 1>
    {
        std::array<T, dim - 1> r;
        std::transform(pt.begin(), std::next(pt.end(), -1), r.begin(), [&pt](const auto &pt_) { return pt_ / pt.back(); });
        return r;
    }

    template <typename T, size_t dim> 
    auto separate_weights(const std::vector<std::array<T, dim+1>> &poles_and_weights,std::vector<std::array<T, dim>> &poles, std::vector<T> &weights) ->void
    {
        poles.resize(poles_and_weights.size());
        weights.resize(poles_and_weights.size());

        std::transform(
            std::execution::par,
            poles_and_weights.begin(),
            poles_and_weights.end(),
            poles.begin(),
            [](const auto &pw_) { return weight_projection<T, dim+1>(pw_); });

        std::transform(
            std::execution::par,
            poles_and_weights.begin(),
            poles_and_weights.end(),
            weights.begin(),
            [](const auto &po_) {return po_.back(); } );
    }

    template <typename T, size_t dim>
    auto add_weight(std::array<T, dim> pole, T w) -> std::array<T, dim+1>
    {
        std::array<T, dim+1> new_pole;
        std::transform(
            pole.begin(),
            pole.end(),
            new_pole.begin(),
            [&w](const auto &x_){return x_*w;}
        );
        new_pole.back()=w;
        return new_pole;
    }

    template <typename T, size_t dim> 
    auto add_weights_coord(const std::vector<std::array<T, dim>> &poles) -> std::vector<std::array<T, dim+1>>
    {
        std::vector<std::array<T, dim+1>> poles_with_weights(poles.size());
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles_with_weights.begin(),
            [](const auto p_)
            {

                return add_weight(p_,T(1.));
            }
        );
        return poles_with_weights;
    }

    template <typename T, size_t dim>
    auto scale_poles(std::vector<std::array<T, dim>> &poles, T scale) -> void
    {
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            [&scale](const auto &p_) { return p_ * scale; });
    }

    template <typename T, size_t dim>
    auto scale_weights(std::vector<std::array<T, dim>> &poles, T mean_value = T(1.0)) -> void
    {
        // auto mean_value_ = std::reduce(
        //     std::execution::par,
        //     poles.begin(),
        //     poles.end(),
        //     T(0.),
        //     [](const auto &s_, const auto &v_) { s_.back() + v_.back(); });
        T mean_value_ = 0.;
        std::for_each(
            // std::execution::par,
            poles.begin(),
            poles.end(),
            [&mean_value_] (const auto &v_) mutable { mean_value_ += v_.back(); });
        scale_poles(poles, mean_value / mean_value_);
    }

    template <typename T, size_t dim>
    auto scaled_poles(const std::vector<std::array<T, dim>> &poles, T scale) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_poles(new_poles);
        return new_poles;
    }

    template <typename T, size_t dim>
    auto scaled_weights(const std::vector<std::array<T, dim>> &poles, T mean_value = T(1.0)) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_weights(new_poles);
        return new_poles;
    }

    template <class T, size_t dim>
    auto eval_rational_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim+1>> &poles, size_t p, size_t d = 0) -> std::array<T, dim>
    {
        if (d == 0)
        {
            return weight_projection(gbs::eval_value_simple(u, k, poles , p, d, true));
        }
        else
        {
            auto wu = eval_value_simple<T,dim+1>(u, k, poles, p, 0, false).back();
            auto Ckw = eval_value_simple<T,dim+1>(u, k, poles, p, d, false);
            Ckw.back() = 1.;                  //
            auto Ak = weight_projection(Ckw); // not real projection just drop last coord
            std::array<T, dim> sum{Ak};
            for (int i = 1; i <= d; i++)
            {
                auto wi = eval_value_simple<T,dim+1>(u, k, poles, p, i, false).back();
                auto C = eval_rational_value_simple<T,dim>(u, k, poles, p, d - i);
                sum = sum - binomial_law<T>(d, i) * wi * C;
            }
            sum = sum / wu;
            return sum;
        }
    }

} // namespace gbs