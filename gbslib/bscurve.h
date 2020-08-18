#pragma once
#include <gbslib/basisfunctions.h>
#include <gbslib/knotsfunctions.h>
#include <vector>
#include <array>
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

        BSCurve(const std::vector<std::array<T, dim - 1>> &poles,
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
                              m_deg(deg) ,
                              m_rational(false)
        {
        }

        auto value(T u, size_t d = 0) const -> std::array<T, dim> 
        {
            return gbs::eval_value_simple(u, m_knotsFlats, m_poles, m_deg, d);
        }

        auto valueRationnal(T u, size_t d = 0) const -> std::array<T, dim-1> 
        {
            return weight_projection(value(u, d));
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

        auto isRational() -> bool {return m_rational;}
    };
} // namespace gbs