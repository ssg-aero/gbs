#pragma once
#include <memory>
#include <gbs/bscurve.h>
namespace gbs {
/**
 * @brief A curve that reverses the parameterization of an underlying curve.
 *
 * This class is useful for creating curves that need to be evaluated in reverse,
 * such as when dealing with control points that are defined in reverse order.
 */
template <typename T, size_t dim> class CurveReversed : public Curve<T, dim> 
{
private:
    /**
     * @brief The underlying curve to be extended.
     */
    std::shared_ptr<Curve<T, dim>> m_curve;

public:
    /**
     * @brief Construct a new CurveReversed object.
     *
     * @param crv Shared pointer to the underlying curve.
     */
    CurveReversed(std::shared_ptr<Curve<T, dim>> crv)
        : m_curve(std::move(crv))
    {
    }

    /**
     * @brief Evaluate the curve (and optionally its derivatives) at the parameter \p u.
     *
     * This method reverses the parameterization of the underlying curve.
     *
     * @param u The parameter value at which to evaluate the curve.
     * @param d The order of derivative to compute (0 for value, 1 for first derivative).
     * @return std::array<T, dim> The evaluated value or derivative at \p u.
     */
    virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
    {
        auto u1 = m_curve->bounds()[0];
        auto u2 = m_curve->bounds()[1];

        u = u2 - (u - u1); // Reverse the parameterization
        T sign = (d % 2 == 0) ? 1 : -1;

        return sign * m_curve->value(u, d);
    }

    virtual auto bounds() const -> std::array<T, 2> override
    {
        return m_curve->bounds();
    }

    const auto & basisCurve() const
    {
        return m_curve;
    }
};
} // namespace gbs
