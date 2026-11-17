#pragma once
/**
 * @file curveoffset.h
 * @brief Offset curves: C(u) + f(u) n(u), the offset f a function of the parameter.
 *
 * Notation. C the basis, n = M / |M| with M = J C' in the plane (J the rotation by +90 degrees:
 * n on the left of C) and M = D x C' in space (D a fixed direction: n normal to D and to C').
 * M^(j) = J C^(j+1) or D x C^(j+1); with w = |M|, Leibniz's rule on w² = M . M and M = w n gives
 *
 *   w^(k) = ( (M . M)^(k) - sum_{0<i<k} C(k,i) w^(i) w^(k-i) ) / 2w
 *   n^(k) = ( M^(k) - sum_{0<i<=k} C(k,i) w^(i) n^(k-i) ) / w
 *   P^(k) = C^(k) + sum_{i<=k} C(k,i) f^(i) n^(k-i)
 *
 * exact at any order k when C has the derivatives of order k + 1 and f those of order k (R. T.
 * Farouki, C. A. Neff, *Analytic properties of plane offset curves*, CAGD 7, 1990; T. Maekawa,
 * *An overview of offset curves and surfaces*, CAD 31(3), 1999). The same recurrences give the
 * offset surfaces of `offsets.h`.
 */
#include <gbs/bscurve.h>
#include <gbs/bscanalysis.h>

#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

namespace gbs
{
    namespace offset_detail
    {
        template <std::floating_point T>
        constexpr auto binomial(std::size_t n, std::size_t k) -> T
        {
            T r{1};
            for (std::size_t i = 1; i <= k; ++i)
                r = r * T(n - k + i) / T(i);
            return r;
        }

        /// Derivatives n^(0..k) of n = M / |M| from M^(0..k) (file header); nullopt where |M| <= small.
        template <std::floating_point T, std::size_t dim>
        auto unit_derivatives(const std::vector<point<T, dim>> &M, T small) -> std::optional<std::vector<point<T, dim>>>
        {
            const std::size_t q = M.size();
            const T w0 = norm(M[0]);
            if (!(w0 > small) || !std::isfinite(w0))
                return std::nullopt;
            std::vector<T> w(q, T(0));
            std::vector<point<T, dim>> n(q, point<T, dim>{});
            w[0] = w0;
            n[0] = M[0] / w0;
            for (std::size_t k = 1; k < q; ++k)
            {
                T g{}, ww{};
                for (std::size_t i = 0; i <= k; ++i)
                {
                    g += binomial<T>(k, i) * (M[i] * M[k - i]);
                    if (i > 0 && i < k)
                        ww += binomial<T>(k, i) * w[i] * w[k - i];
                }
                w[k] = (g - ww) / (T(2) * w0);
                auto r = M[k];
                for (std::size_t i = 1; i <= k; ++i)
                    r = r - (binomial<T>(k, i) * w[i]) * n[k - i];
                n[k] = r / w0;
            }
            return n;
        }

        /// P^(d) = C^(d) + sum C(d,i) f^(i) n^(d-i), M^(j) = normal(C^(j+1)).
        template <typename T, std::size_t dim, typename Func, typename Normal>
        auto offset_value(const Curve<T, dim> &c, const Func &f, T u, std::size_t d, Normal &&normal) -> point<T, dim>
        {
            std::vector<point<T, dim>> M(d + 1);
            for (std::size_t j = 0; j <= d; ++j)
                M[j] = normal(c.value(u, j + 1));
            const auto n = unit_derivatives<T, dim>(M, T(0));
            if (!n)
                throw std::domain_error("offset curve: no normal (null tangent, or tangent along the direction)");
            auto p = c.value(u, d);
            for (std::size_t i = 0; i <= d; ++i)
                p = p + (binomial<T>(d, i) * f(u, i)) * (*n)[d - i];
            return p;
        }
    }

    template <typename T, size_t dim, typename Func>
    class CurveOffset : public Curve<T, dim>
    {
        const std::shared_ptr<Curve<T, dim>> p_crv_;
        const std::shared_ptr<Func> f_offset_;

    public:
        CurveOffset(const CurveOffset<T, dim, Func> &crv) = default;
        CurveOffset(const std::shared_ptr<Curve<T, dim>> &crv, const Func &f_offset) : p_crv_{crv}, f_offset_{std::make_shared<Func>(f_offset)}
        {
        }
        CurveOffset(const std::shared_ptr<Curve<T, dim>> &crv, const std::shared_ptr<Func> &f_offset) : p_crv_{crv}, f_offset_{f_offset}
        {
        }
        CurveOffset(const BSCurve<T, dim> &crv, const Func &f_offset) : p_crv_{std::make_shared<BSCurve<T, dim>>(crv)},
                                                                        f_offset_{std::make_shared<Func>(f_offset)}
        {
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_crv_->bounds();
        }
        auto changeBounds(const std::array<T, 2> &b) -> void
        {
            p_crv_->changeBounds(b);
        }
        auto basisCurve() const -> const Curve<T, dim> &
        {
            return *p_crv_;
        }
        auto offset() const -> const Func &
        {
            return *f_offset_;
        }
    };

    /**
     * @brief An offset curve in the plane: C + f n, n the unit normal on the left of C
     * (derivatives of any order: file header).
     *
     * @tparam T Floating-point type
     * @tparam Func A callable object representing the offset function
     */
    template <typename T, typename Func>
    class CurveOffset2D : public CurveOffset<T, 2, Func>
    {
    public:
        /**
         * @brief Constructor that takes a shared_ptr of a Curve, a shared_ptr of an offset function, and a specified direction
         *
         * @param crv Shared pointer to the base curve
         * @param f_offset Shared pointer to the offset function
         */
        CurveOffset2D(const std::shared_ptr<Curve<T, 2>> &crv, const std::shared_ptr<Func> &f_offset)
            : CurveOffset<T, 2, Func>{crv, f_offset}
        {
        }
        CurveOffset2D(const std::shared_ptr<Curve<T, 2>> &crv, const Func &f_offset)
            : CurveOffset<T, 2, Func>{crv, f_offset}
        {
        }
        CurveOffset2D(const BSCurve<T, 2> &crv, const Func &f_offset)
            : CurveOffset<T, 2, Func>{crv, f_offset}
        {
        }
        CurveOffset2D(const CurveOffset2D<T, Func> &crv) = default;
        /** @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, 2>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 2> override;
    };

    template <typename T, typename Func>
    auto CurveOffset2D<T, Func>::value(T u, size_t d) const -> std::array<T, 2>
    {
        return offset_detail::offset_value(this->basisCurve(), this->offset(), u, d,
                                           [](const point<T, 2> &v) { return point<T, 2>{-v[1], v[0]}; });
    }

    /**
    * @brief An offset curve in space: C + f n, n = D x C' / |D x C'| for a fixed direction D
    * (derivatives of any order: file header).
    *
    * @tparam T Floating-point type
    * @tparam Func A callable object representing the offset function
    */
    template <typename T, typename Func>
    class CurveOffset3D : public CurveOffset<T, 3, Func>
    {
    public:
        /**
         * @brief Constructor that takes a shared_ptr of a Curve, a shared_ptr of an offset function, and a specified direction
         *
         * @param crv Shared pointer to the base curve
         * @param f_offset Shared pointer to the offset function
         * @param direction Direction vector of the offset in 3D space
         */
        CurveOffset3D(const std::shared_ptr<Curve<T, 3>> &crv, const std::shared_ptr<Func> &f_offset, const std::array<T, 3> &direction)
            : CurveOffset<T, 3, Func>{crv, f_offset}, direction_{normalized( direction )}
        {
        }
        CurveOffset3D(const std::shared_ptr<Curve<T, 3>> &crv, const Func &f_offset, const std::array<T, 3> &direction)
            : CurveOffset<T, 3, Func>{crv, f_offset}, direction_{normalized( direction )}
        {
        }
        CurveOffset3D(const CurveOffset3D<T, Func> &crv) = default;
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 3> override;
        const std::array<T,3>& direction() const {return direction_;}
    private:
        std::array<T, 3> direction_;               ///< Specified direction vector of the offset in 3D space
    };


    template <typename T, typename Func>
    auto CurveOffset3D<T, Func>::value(T u, size_t d) const -> std::array<T, 3>
    {
        return offset_detail::offset_value(this->basisCurve(), this->offset(), u, d,
                                           [this](const point<T, 3> &v) { return cross(direction_, v); });
    }
}
