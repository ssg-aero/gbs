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
                std::copy(v.begin(), v.end(), n.begin());
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
    /**
 * @brief Géneral BSpline curve class, any kind of precision, space dimension with rational definition capability
 * 
 * @tparam T    : curve precision
 * @tparam dim  : space dimension of curve (aka 1D, 2D, 3D,...)
 */
    template <typename T, size_t dim, bool rational>
    class BSCurveGeneral
    {
        bool m_rational;
        size_t m_deg;
        std::vector<std::array<T, dim + rational>> m_poles;
        std::vector<T> m_knotsFlats;

    public:
        /**
     * @brief Construct a new BSCurve object, non rational definition
     * 
     * @param poles : array of poles
     * @param knots : array of knots
     * @param mult  : array of knots multiplicity
     * @param deg   : curve's degree
     */
        BSCurveGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                       const std::vector<T> &knots,
                       const std::vector<size_t> &mult,
                       size_t deg) : m_poles(poles),
                                     m_deg(deg),
                                     m_rational(false),
                                     m_knotsFlats(flat_knots(knots, mult))
        {
        }
        /**
         * @brief Construct a new BSCurve object, rational definition
         * 
         * @param poles   : array of poles
         * @param weights : array of poles' weights
         * @param knots   : array of knots
         * @param mult    : array of knots multiplicity
         * @param deg     : curve's degree
         */
        BSCurveGeneral(const std::vector<std::array<T, dim>> &poles,
                       const std::vector<T> &weights,
                       const std::vector<T> &knots,
                       const std::vector<size_t> &mult,
                       size_t deg) : m_poles(merge_weights(poles, weights)),
                                     m_deg(deg),
                                     m_rational(true),
                                     m_knotsFlats(flat_knots(knots, mult))
        {
        }
        /**
         * @brief Construct a new BSCurve object
         * 
         * @param poles       : array of poles
         * @param knots_flats : flat knots 
         * @param deg         : curve's degree
         */
        BSCurveGeneral(const std::vector<std::array<T, dim + rational>> &poles,
                       const std::vector<T> &knots_flats,
                       size_t deg) : m_poles(poles),
                                     m_knotsFlats(knots_flats),
                                     m_rational(false),
                                     m_deg(deg)
        {
        }
        /**
         * @brief Non rational curve evaluation
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> = 0;

        // auto value(T u, size_t d = 0) const -> std::array<T, dim>  {}
        /**
         * @brief Non rational curve's begin
         * 
         * @param d : derivative order
         * @return std::array<T, dim> 
         */
        auto begin(size_t d = 0) const -> std::array<T, dim>
        {
            return value(m_knotsFlats.front(), d);
        }
        /**
         * @brief Non rational curve's end
         * 
         * @param d : derivative order
         * @return std::array<T, dim> 
         */
        auto end(size_t d = 0) const -> std::array<T, dim>
        {
            return value(m_knotsFlats.back(), d);
        }

        /**
         * @brief curve's degree
         * 
         * @return size_t 
         */
        auto degree() const noexcept -> size_t
        {
            return m_deg;
        }
        /**
         * @brief curves's flat knots
         * 
         * @return const std::vector<T>& 
         */
        auto knotsFlats() const noexcept -> const std::vector<T> &
        {
            return m_knotsFlats;
        }
        /**
         * @brief Insert knot with the given multiplicity
         * 
         * @param u : knot value
         * @param m : knot's multiplicity
         */
        auto insertKnot(T u, size_t m = 1) -> void //Fail safe, i.e. if fails, curve stays in precedent state
        {
            for (auto i = 0; i < m; i++)
                insert_knot(u, m_deg, m_knotsFlats, m_poles);
        }
        /**
         * @brief Try to remove m times the given knot
         * 
         * @param u   : knot value
         * @param tol : tolerance on curve
         * @param m   : knot removal occurrences
         */
        auto removeKnot(T u, T tol, size_t m = 1) -> void //Fail safe, i.e. if fails, curve stays in precedent state
        {
            for (auto i = 0; i < m; i++)
                remove_knot(u, m_deg, m_knotsFlats, m_poles, tol);
        }
        /**
         * @brief Curve's poles
         * 
         * @return const std::vector<std::array<T, dim>>& 
         */
        auto poles() const noexcept -> const std::vector<std::array<T, dim + rational>> &
        {
            return m_poles;
        }
        /**
         * @brief 
         * 
         * @return bool
         */
        auto isRational() const -> bool { return m_rational; }
        /**
         * @brief Set the Rational object
         * 
         * @param rational 
         */
        auto setRational(bool rational) -> void { m_rational = rational; }
        /**
         * @brief reverse curve orientation
         * 
         */
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
        /**
         * @brief Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones
         * 
         * @param u1 
         * @param u2 
         */
        auto trim(T u1, T u2) -> void
        {
            gbs::trim(m_deg, m_knotsFlats, m_poles, u1, u2);
        }
    };

    template <typename T, size_t dim>
    class BSCurveRational : public BSCurveGeneral<T, dim, true>
    {
    public:
        BSCurveRational(const std::vector<std::array<T, dim + 1>> &poles,
                        const std::vector<T> &knots_flats,
                        size_t deg) : BSCurveGeneral<T, dim, true>(poles, knots_flats, deg) {}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            if (d == 0)
            {
                return weight_projection(gbs::eval_value_simple(u, knotsFlats(), poles(), degree(), d, true));
            }
            else
            {
                // auto wu = value(u).back();
                auto wu = gbs::eval_value_simple(u, knotsFlats(), poles(), degree(), 0, false).back();
                // auto Ckw= value(u,d);
                auto Ckw = gbs::eval_value_simple(u, knotsFlats(), poles(), degree(), d, false);
                Ckw.back() = 1.;
                auto Ak = weight_projection(Ckw); // not real projection just drop last coord
                std::array<T, dim> sum{Ak};
                for (int i = 1; i <= d; i++)
                {
                    // auto wi = value(u, i).back();
                    auto wi = gbs::eval_value_simple(u, knotsFlats(), poles(), degree(), i, false).back();
                    auto C = value(u, d - i);
                    sum = sum - binomial_law<T>(d, i) * wi * C;
                }
                sum = sum / wu;
                return sum;
            }
        }
    };

    template <typename T, size_t dim>
    class BSCurve : public BSCurveGeneral<T, dim, false>
    {
    public:
        BSCurve(const std::vector<std::array<T, dim>> &poles,
                const std::vector<T> &knots_flats,
                size_t deg) : BSCurveGeneral<T, dim, false>(poles, knots_flats, deg) {}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            return gbs::eval_value_simple(u, knotsFlats(), poles(), degree(), d);
        }
    };

    using BSCurve2d_f = BSCurve<float,2>;
    using BSCurve3d_f = BSCurve<float,3>;
    using BSCurve2d_d = BSCurve<double,2>;
    using BSCurve3d_d = BSCurve<double,3>;
    using BSCurveRational2d_f = BSCurveRational<float,2>;
    using BSCurveRational3d_f = BSCurveRational<float,3>;
    using BSCurveRational2d_d = BSCurveRational<double,2>;
    using BSCurveRational3d_d = BSCurveRational<double,3>;
} // namespace gbs