#pragma once

#include <limits>
#include <memory>
#include <array>
#include <gbs/bscurve.h>

namespace gbs
{

/**
 * @brief Extends a given curve to an infinite domain.
 *
 * This class wraps an existing curve whose domain is \([u_1, u_2]\). For values of
 * \p u less than \p u_1 or greater than \p u_2, it either clamps to the boundary
 * value or performs a linear extrapolation based on the curve’s first derivative at
 * that boundary.
 *
 * @tparam T   Floating-point type used for the curve parameter and coordinates.
 * @tparam dim Dimension of the curve.
 */
template <typename T, size_t dim>
class CurveExtended : public Curve<T, dim>
{
private:
    /**
     * @brief The underlying curve to be extended.
     */
    std::shared_ptr<Curve<T, dim>> m_curve;

    /**
     * @brief Whether to clamp (true) or linearly extrapolate (false) outside [u1, u2].
     */
    bool m_clamped;

public:
    /**
     * @brief Construct a new CurveExtended object.
     *
     * @param crv      Shared pointer to the underlying curve.
     * @param clamped  If true, values of \p u outside the domain [u1, u2] are clamped
     *                 to the respective boundary value. If false, linear extrapolation
     *                 is used outside [u1, u2].
     */
    CurveExtended(std::shared_ptr<Curve<T, dim>> crv, bool clamped = false)
        : m_curve(std::move(crv))
        , m_clamped(clamped)
    {
    }

    /**
     * @brief Evaluate the curve (and optionally its derivatives) at the parameter \p u.
     *
     * - If \p u is within [u1, u2], this method calls the underlying curve’s value().
     * - If \p u is less than u1 or greater than u2, the result depends on \p m_clamped:
     *   - If clamped, it returns the boundary value/derivative at u1 or u2.
     *   - If not clamped, it performs a linear extrapolation using the boundary value
     *     plus the boundary derivative times \c (u - boundary).
     *   - For derivative orders (\p d) higher than 1, the linear extrapolation returns
     *     zero beyond the first derivative.
     *
     * @param u The parameter at which to evaluate the curve or its derivative.
     * @param d The order of the derivative. Defaults to 0 (position).
     * @return std::array<T, dim> The curve value (or derivative) at \p u.
     */
    virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
    {
        // Extract the original bounds of the underlying curve
        auto [u1, u2] = m_curve->bounds();

        // If u is before the underlying domain
        if (u < u1)
        {
            if (m_clamped)
            {
                // Simply clamp to the lower boundary
                return m_curve->value(u1, d);
            }
            else
            {
                // Perform linear extrapolation (only up to first derivative)
                if (d == 0)
                {
                    auto val   = m_curve->value(u1, 0);  // boundary value
                    auto deriv = m_curve->value(u1, 1);  // boundary derivative
                    std::array<T, dim> result{};
                    for (size_t i = 0; i < dim; ++i)
                    {
                        result[i] = val[i] + deriv[i] * (u - u1);
                    }
                    return result;
                }
                else if (d == 1)
                {
                    // First derivative is constant outside
                    return m_curve->value(u1, 1);
                }
                else
                {
                    // Higher-order derivatives are zero for a linear extension
                    return {};
                }
            }
        }
        // If u is after the underlying domain
        else if (u > u2)
        {
            if (m_clamped)
            {
                // Clamp to the upper boundary
                return m_curve->value(u2, d);
            }
            else
            {
                // Linear extrapolation (only up to first derivative)
                if (d == 0)
                {
                    auto val   = m_curve->value(u2, 0);
                    auto deriv = m_curve->value(u2, 1);
                    std::array<T, dim> result{};
                    for (size_t i = 0; i < dim; ++i)
                    {
                        result[i] = val[i] + deriv[i] * (u - u2);
                    }
                    return result;
                }
                else if (d == 1)
                {
                    // First derivative is constant outside
                    return m_curve->value(u2, 1);
                }
                else
                {
                    // Higher-order derivatives are zero in linear extension
                    return {};
                }
            }
        }

        // Otherwise, u is within [u1, u2], so just delegate to the underlying curve
        return m_curve->value(u, d);
    }

    /**
     * @brief Returns the effective domain of the extended curve.
     *
     * Since this curve is extended to be valid everywhere, it returns the full range
     * of representable values for type \p T.
     *
     * @return std::array<T, 2> The extreme numeric limits for \p T.
     */
    virtual auto bounds() const -> std::array<T, 2> override
    {
        return {
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::max()
        };
    }
};

} // namespace gbs