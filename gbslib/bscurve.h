#pragma once
#include <gbslib/basisfunctions.h>
#include <gbslib/knotsfunctions.h>
#include <gbslib/bsctools.h>
#include <gbslib/math.h>
#include <gbslib/vecop.h>

#include <vector>
#include <array>
#include <any>
namespace gbs
{
    template <typename T, size_t dim>
    auto merge_weights(const std::vector<std::array<T, dim>> &poles, const std::vector<T> &weigths) -> std::vector<std::array<T, dim + 1>>
    {
        std::vector<std::array<T, dim + 1>> p(poles.size());
        std::transform(
            std::execution::par,
            poles.begin(), poles.end(), weigths.begin(), p.begin(),
            [](const auto &v, const auto &c) {
                std::array<T, dim + 1> n;
                std::copy(v.begin(), v.end(), n.begin());
                n.back() = c;
            });
        return p;
    }

    template <typename T, size_t dim>
    auto weight_projection(const std::array<T, dim> &p) -> std::array<T, dim - 1>
    {
        std::array<T, dim - 1> r;
        std::transform(p.begin(), std::next(p.end(), -1), r.begin(), [&p](const auto &p_) { return p_ / p.back(); });
        return r;
    }

    template <typename T, size_t dim>
    class BSCurve
    {
        bool m_rational;
        size_t m_deg;
        std::vector<std::array<T, dim>> m_poles;
        std::vector<T> m_knotsFlats;

    public:
        BSCurve(const std::vector<std::array<T, dim>> &poles,
                const std::vector<T> &knots,
                const std::vector<size_t> &mult,
                size_t deg) : m_poles(poles),
                              m_deg(deg),
                              m_rational(false),
                              m_knotsFlats(flat_knots(knots, mult))
        {
        }

        BSCurve(const std::vector<std::array<T, dim>> &poles,
                const std::vector<T> &weights,
                const std::vector<T> &knots,
                const std::vector<size_t> &mult,
                size_t deg) : m_poles(merge_weights(poles, weights)),
                              m_deg(deg),
                              m_rational(true),
                              m_knotsFlats(flat_knots(knots, mult))
        {
        }

        BSCurve(const std::vector<std::array<T, dim>> &poles,
                const std::vector<T> &knots_flats,
                size_t deg) : m_poles(poles),
                              m_knotsFlats(knots_flats),
                              m_rational(false),
                              m_deg(deg) 
        {
        }
        auto value(T u, size_t d = 0) const -> std::array<T, dim> 
        {
            return gbs::eval_value_simple(u, m_knotsFlats, m_poles, m_deg, d);
            
        }

        auto begin(size_t d = 0) const -> std::array<T, dim> 
        {
            return value(m_knotsFlats.front(),d);
        }

        auto end(size_t d = 0) const -> std::array<T, dim> 
        {
            return value(m_knotsFlats.back(),d);
        }

        auto valueRationnal(T u, size_t d = 0) const -> std::array<T, dim-1> 
        {
            if (d == 0)
            {
                return weight_projection(value(u));
            }
            else
            {
                /*
                auto wu = value(u).back();
                auto Ckw= value(u,d);
                Ckw.back() = 1.;
                auto Ak = weight_projection(Ckw); // not real projection just drop last coord
                std::vector<std::array<T, dim-1>> v{d};
                std::transform(v.begin(), v.end(), v.begin(),
                               [&, i = 1] (const auto v_) mutable {
                                   auto wi = value(u, i).back();
                                   auto C = valueRationnal(u, d - i);
                                   auto res = binomial_law<T>(d,i) * wi * C;
                                    i++;
                                    return res;
                               });
                std::array<T, dim-1> sum{};
                // sum = std::reduce(v.begin(),v.end(), sum,gbs::operator+<T, dim-1>);
                for(auto v_ : v)
                {
                    sum = sum +v_;
                }
                return (Ak-sum) / wu;
                */
                auto wu = value(u).back();
                auto Ckw= value(u,d);
                Ckw.back() = 1.;
                auto Ak = weight_projection(Ckw); // not real projection just drop last coord
                std::array<T, dim-1> sum{Ak};
                for (int i = 1; i <= d; i++)
                {
                    auto wi = value(u, i).back();
                    auto C = valueRationnal(u, d - i);
                    sum = sum - binomial_law<T>(d, i) * wi * C;
                }
                sum = sum / wu;
                return sum;
            }
        }

        auto beginRationnal(size_t d = 0) const -> std::array<T, dim-1> 
        {
            return valueRationnal(m_knotsFlats.front(),d);
        }

        auto endRationnal(size_t d = 0) const -> std::array<T, dim-1> 
        {
            return valueRationnal(m_knotsFlats.back(),d);
        }

        auto degree() const noexcept -> size_t
        {
            return m_deg;
        }

        auto knotsFlats() const noexcept -> const std::vector<T> &
        {
            return m_knotsFlats;
        }

        auto insertKnot(T u, size_t m = 1) -> void //Fail saife, i.e. if fails, curve stays in precedent state
        {
            for (auto i = 0; i < m; i++)
                insert_knot(u, m_deg, m_knotsFlats, m_poles);
        }

        auto removeKnot(T u, T tol, size_t m = 1) -> void //Fail saife, i.e. if fails, curve stays in precedent state
        {
            for (auto i = 0; i < m; i++)
                remove_knot(u, m_deg, m_knotsFlats, m_poles, tol);
        }

        auto poles() const noexcept -> const std::vector<std::array<T, dim>> &
        {
            return m_poles;
        }

        auto isRational() const -> bool { return m_rational; }
        auto setRational(bool rational) -> void { m_rational = rational; }
        auto reverse() -> void
        {
            std::reverse(m_poles.begin(), m_poles.end());
            auto k1 = m_knotsFlats.front();
            auto k2 = m_knotsFlats.back();
            std::reverse(m_knotsFlats.begin(), m_knotsFlats.end());
            std::transform(m_knotsFlats.begin(), m_knotsFlats.end(),
                           m_knotsFlats.begin(),
                           [&](const auto k_) {
                               return k1 + k2 - k_;
                           });
        }

        auto trim(T u1, T u2) -> void
        {
            gbs::trim(m_deg, m_knotsFlats,m_poles, u1, u2);
        }
    };


} // namespace gbs