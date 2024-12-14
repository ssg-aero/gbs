#pragma once
#include <vector>
#include <utility>
#include <gbs/gbslib.h>
import vecop;
#include <gbs/maths.h>
#include <execution>
#include <algorithm>
#include <cmath>

namespace gbs
{
    /**
     * @brief Basis function used to compute BSpline
     *
     * @tparam InIt
     * @tparam T
     * @param u  Parameter on BSpline object
     * @param it Flat knots iterator
     * @param p  BSpline object degree
     * @param _Last Container's end
     * @return T
     */
    template <std::input_iterator InIt, std::floating_point T>
    T basis_function(T u, const InIt &it, size_t p, const InIt &last)
    {
        auto u_last = *std::next(last, -1);

        T ui = *it;
        T ui1 = *std::next(it, 1);
        if (p == 0)
        {
            return ((ui <= u) && (u < ui1)) || (std::abs(ui1 - u_last) < knot_eps<T> && std::abs(u - u_last) < knot_eps<T>)
                       ? T(1.)
                       : T(0.);
        }
        else
        {
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps<T>)
            {
                C1 = (u - ui) / C1;
                C1 *= basis_function(u, it, p - 1, last);
            }
            if (C2 > knot_eps<T>)
            {
                C2 = (ui1p - u) / C2;
                C2 *= basis_function(u, std::next(it, 1), p - 1, last);
            }
            return C1 + C2;
        }
    }

    /**
     * @brief Basis function used to compute BSpline's derivatives
     *
     * @tparam InIt
     * @tparam T
     * @param u   Parameter on BSpline object
     * @param it  Flat knots iterator
     * @param p   BSpline object degree
     * @param d   Derivative order
     * @param _Last
     * @return T
     */
    template <std::input_iterator InIt, std::floating_point T>
    T basis_function(T u, const InIt &it, size_t p, size_t d, const InIt &last)
    {
        if (d == 0)
        {
            return basis_function(u, it, p, last);
        }
        else if (d > p)
        {
            return T(0.);
        }
        else
        {
            T ui = *it;
            T ui1 = *std::next(it, 1);
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps<T>)
            {
                C1 = basis_function(u, it, p - 1, d - 1, last) / C1;
            }
            if (C2 > knot_eps<T>)
            {
                C2 = basis_function(u, std::next(it, 1), p - 1, d - 1, last) / C2;
            }
            return p * (C1 - C2);
        }
    }

    template <typename T>
    auto basis_funcs(size_t i, size_t p, T u, const std::vector<T> &k, std::vector<T>& N)
    {
        
        std::vector<T> left(p+1), right(p+1);
        T saved{}, temp{};
        N[0] = static_cast<T>(1.);
        for(size_t j{1}; j <= p; j++)
        {
            left[j]  = u - k[i+1-j];
            right[j] = k[i+j] - u;
            saved = 0.0;
            for( size_t r{0}; r <j; r++)
            {
                temp = N[r] / (right[r+1]+left[j-r]);
                N[r] = saved + right[r+1]*temp;
                saved = left[j-r]*temp;
            }
            N[j] = saved;
        }
    }

    template <typename T>
    auto find_span(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        auto first = std::next(k.begin(), p);
        auto last  = std::next(k.begin(), n + 1);
        auto pos = std::lower_bound(first, last, u);
        if( pos != first)
            return std::next(pos,-1);
        else
            return pos;
    }
    
    template <typename T>
    auto find_span2(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        if (u >= k[n + 1])
            return std::next(k.begin(), n);
        auto low = p, high = n + 1, mid = (low + high) / 2;
        while (u < k[mid] || u >= k[mid + 1])
        {
            if (u < k[mid])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }
        return std::next(k.begin(), mid);
    }

    template <typename T>
    auto find_span_range(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        size_t i2 = find_span(n, p, u, k) - k.begin();
        i2 = std::min(static_cast<long long> (i2), static_cast<long long> (k.size() - p - 2));
        size_t i1 = std::max(static_cast<long long>(0), static_cast<long long>(i2-p));

        return std::make_pair(i1,i2);
    }

    // template <typename T>
    // auto eval_index_span(size_t n, size_t p, T u, const std::vector<T> &k, size_t d, bool use_span_reduction, size_t &i_min, size_t &i_max)
    // {
    //     if (use_span_reduction)//Reducing span for few pole makes things worst
    //     {
    //         i_max = find_span(n, p, u, k) - k.begin();
    //         i_max = std::min(i_max, k.size() - p - 2);
    //         i_min = std::max(int(0),int(i_max-p));
    //     }
    // }

    template <typename T>
    T basis_function(T u, size_t i, size_t p, size_t d, const std::vector<T> &k)
    {
        return basis_function(u, std::next(k.begin(), i), p, d, k.end());
    }

    /**
     * @brief BSpline curve evaluation using simple recursive basis functions
     * 
     * @tparam T 
     * @tparam dim 
     * @param u parameter on curve
     * @param k flat knots
     * @param poles poles
     * @param p degree
     * @param d derivative order
     * @param use_span_reduction 
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    auto eval_value_decasteljau(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> std::array<T, dim>
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        size_t i_max{n_poles-1}, i_min{0};
        if (use_span_reduction && d == 0 )//Reducing span for few pole makes things worst // TODO fix for d != 0
        {
            i_max = find_span(n_poles, p, u, k) - k.begin();
            i_max = std::min(i_max, k.size() - p - 2);
            i_min = std::max(int(0),int(i_max-p));
        }

        for (auto i = i_min; i <= i_max; i++)
        {
            auto N = basis_function(u, i, p, d, k);

            std::transform(
                // std::execution::par,
                poles[i].begin(),
                poles[i].end(),
                pt.begin(),
                pt.begin(),
                [&] (const auto val_, const auto tot_)
                {
                    return tot_ + N * val_;
                }
            );
        }
        
        return pt;
    }


    // TODO Seems dead code
    /**
     * @brief BSpline curve evaluation using simple recursive basis functions
     * 
     * @tparam T 
     * @tparam dim 
     * @param u parameter on curve
     * @param k flat knots
     * @param poles poles
     * @param p degree
     * @param d derivative order
     * @param use_span_reduction 
     * @return T
     */
    template <typename T, size_t dim>
    auto eval_value_decasteljau(T u, const std::vector<T> &k, const std::vector<T> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> T
    {
        T pt{};
        size_t n = poles.size();
        size_t i_max, i_min;
        eval_index_span(n, p, u, k, d, use_span_reduction, i_min, i_max);

        for (auto i = i_min; i <= i_max; i++)
        {
            auto N = basis_function(u, i, p, d, k);
            pt += N * poles[i];
        }
        
        return pt;
    }


    template <typename T, size_t dim>
    auto eval_value_par(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> std::array<T, dim>
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        size_t i_max{n_poles-1}, i_min{0};
        if (use_span_reduction && d == 0 )//Reducing span for few pole makes things worst // TODO fix for d != 0
        {
            i_max = find_span(n_poles, p, u, k) - k.begin();
            i_max = std::min(i_max, k.size() - p - 2);
            i_min = std::max(int(0),int(i_max-p));
        }
        // The benifit of paralelization is not obvious
        std::vector<T> Ni(i_max+1-i_min);
        auto indexes = make_range(i_min,i_max);
        std::transform(
            std::execution::par,
            indexes.begin(),
            indexes.end(),
            Ni.begin(),
            [u,p,d,&k](size_t i) {
                return basis_function(u, i, p, d, k);
            });
        
        points_vector<T, dim> poles_in_span(poles.begin()+i_min,poles.begin()+i_max+1);
        std::transform(
            std::execution::par,
            Ni.begin(),
            Ni.end(),
            poles_in_span.begin(),
            poles_in_span.begin(),
            [](auto N_, const auto &pole_)
            {
                return pole_ * N_;
            }
        );

        // // return point<T,dim>{};        
        
        return std::reduce(
            std::execution::par,
            poles_in_span.cbegin(),
            poles_in_span.cend()
        );

        // // auto indexes = make_range(i_min,i_max);

        // for (auto i = i_min; i <= i_max; i++)
        // {
        //     auto N = basis_function(u, i, p, d, k);
        //     // auto pole = poles[i];
        //     // for( size_t d {} ; d < dim ; d++ )
        //     // {
        //     //     pt[d] += N * pole[d];
        //     // }
        //     std::transform(
        //         // std::execution::par,
        //         poles[i].begin(),
        //         poles[i].end(),
        //         pt.begin(),
        //         pt.begin(),
        //         [&] (const auto val_, const auto tot_)
        //         {
        //             return tot_ + N * val_;
        //         }
        //     );
        // }
        
        // return pt;
    }

    template <typename T, size_t dim>
    auto eval_value_deboor_cox(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p)
    {
        auto n_poles = poles.size();

        auto i =  std::min<size_t>( find_span2(n_poles, p, u, k) - k.begin(), k.size() - p - 2);

        std::vector<T> N(p+1);
        basis_funcs( i, p, u, k, N );

        point<T,dim> pt{0., 0., 0.};
        
        std::for_each(N.begin(), N.end(), [&] (const auto &N_)
        {
            std::transform(
                poles[i-p].begin(),
                poles[i-p].end(),
                pt.begin(),
                pt.begin(),
                [N_] (const auto pole_coordinate, const auto pt_coordinate)
                {
                    return pt_coordinate + N_ * pole_coordinate;
                }
            );
            i++;
        });
        return pt;
    }
    /**
     * @brief BSpline surface evaluation using simple recursive basis functions
     * 
     * @tparam T 
     * @tparam dim 
     * @param u first parameter on curve
     * @param v second parameter on curve
     * @param ku first flats knots
     * @param kv second flats knots
     * @param poles poles
     * @param p u degree
     * @param q v degree
     * @param du derivative order on u
     * @param dv derivative order on v
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    std::array<T, dim> eval_value_decasteljau(T u, T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t du = 0, size_t dv = 0)
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
                return n;
            });
        return p;
    }
    /**
     * @brief Project point with rational definition if last coord is null -> return origin point
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
        if(pt.back()!=T{0}) 
        {
            std::transform(pt.begin(), std::next(pt.end(), -1), r.begin(), [&pt](const auto &pt_)
                           { return pt_ / pt.back(); });
        }
        else
        {
            r.fill(T{0.});
        }
        
        return r;
    }
    template <typename T, size_t dim, bool rational>
    inline auto weight_projected_pole(const point<T,dim+rational> &p) -> point<T,dim>
    {
        std::array<T, dim> r;
        if( rational )
        {
            if(p.back()!=T{0}) 
            {
                std::transform(p.begin(), std::next(p.end(), -1), r.begin(), [&p](const auto &p_)
                            { return p_ / p.back(); });
            }
            else
            {
                r.fill(T{0.});
            }
            
        }
        else
        {
            std::copy(p.begin(), p.end(), r.begin());
        }
        return r;
    }

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles_and_weights compact poles and weights storage
     * @param poles vector receiving projected poles 
     * @param weights vector receiving weights
     */
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
    /**
     * @brief Add weight coordinates to a pole, the other coordinates are scaled
     * 
     * @tparam T 
     * @tparam dim 
     * @param pole 
     * @param w 
     * @return std::array<T, dim+1> 
     */
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
    /**
     * @brief Add weight corrdinate colum, the other coordinates are scaled
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @return std::vector<std::array<T, dim+1>> 
     */
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

    /**
     * @brief Builds a pole vector containing weights
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param weights 
     * @return std::vector<std::array<T, dim+1>> 
     */
    template <typename T, size_t dim> 
    auto add_weights_coord(const std::vector<std::array<T, dim>> &poles, const std::vector<T> &weights) -> std::vector<std::array<T, dim+1>>
    {
        if (poles.size() != weights.size())
        {
            throw std::length_error("BSCurveGeneral: wrong pole vector length.");
        }
        std::vector<std::array<T, dim+1>> poles_with_weights(poles.size());
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            weights.begin(),
            poles_with_weights.begin(),
            [](const auto p_, const auto w_)
            {

                return add_weight(p_,w_);
            }
        );
        return poles_with_weights;
    }

    /**
     * @brief Scale poles with origin {0. ... 0.}
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param scale 
     */
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
        mean_value_ /= poles.size();
        scale_poles(poles, mean_value / mean_value_);
    }

    template <typename T, size_t dim>
    auto scaled_poles(const std::vector<std::array<T, dim>> &poles, T scale) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_poles(new_poles,scale);
        return new_poles;
    }

    template <typename T, size_t dim>
    auto scaled_weights(const std::vector<std::array<T, dim>> &poles, T mean_value = T(1.0)) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_weights(new_poles);
        return new_poles;
    }
    /**
     * @brief Rational evaluation os a pol set and knots the last coordinate stands for weight and the others ares scales according to this weight
     * 
     * @tparam T 
     * @tparam dim 
     * @param u 
     * @param k 
     * @param poles 
     * @param p 
     * @param d 
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    auto eval_rational_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim+1>> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> std::array<T, dim>
    {
        if (d == 0)
        {
            return weight_projection(eval_value_decasteljau(u, k, poles , p, d, use_span_reduction));
        }
        else
        {
            auto wu = eval_value_decasteljau<T,dim+1>(u, k, poles, p, 0, false).back();
            auto Ckw = eval_value_decasteljau<T,dim+1>(u, k, poles, p, d, false);
            Ckw.back() = 1.;                  //
            auto Ak = weight_projection(Ckw); // not real projection just drop last coord
            std::array<T, dim> sum{Ak};
            for (size_t i = 1; i <= d; i++)
            {
                auto wi = eval_value_decasteljau<T,dim+1>(u, k, poles, p, i, false).back();
                auto C = eval_rational_value_simple<T,dim>(u, k, poles, p, d - i);
                sum = sum - binomial_law<T>(d, i) * wi * C;
            }
            sum = sum / wu;
            return sum;
        }
    }

} // namespace gbs