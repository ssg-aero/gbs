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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
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

        /**
         * Derivatives of n = M / |M| on the rectangle a <= A, b <= B of orders in two variables
         * (row-major, B + 1 per row; one variable: B = 0) from those of M, by the recurrences of
         * the file header, |M| > 0 checked by the caller; w receives the derivatives of |M|. Each
         * (a, b) needs only the orders below it in both variables: the rectangle is closed.
         */
        template <std::floating_point T, std::size_t dim>
        void unit_derivatives(std::span<const point<T, dim>> M, std::size_t A, std::size_t B, std::span<T> w, std::span<point<T, dim>> n)
        {
            const std::size_t q = B + 1;
            const T w0 = norm(M[0]);
            w[0] = w0;
            n[0] = M[0] / w0;
            for (std::size_t a = 0; a <= A; ++a)
                for (std::size_t b = 0; b <= B; ++b)
                {
                    if (a == 0 && b == 0)
                        continue;
                    T g{}, ww{};
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                        {
                            const T c = binomial<T>(a, i) * binomial<T>(b, j);
                            g += c * (M[i * q + j] * M[(a - i) * q + (b - j)]);
                            if (!(i == 0 && j == 0) && !(i == a && j == b))
                                ww += c * w[i * q + j] * w[(a - i) * q + (b - j)];
                        }
                    w[a * q + b] = (g - ww) / (T(2) * w0);
                    auto r = M[a * q + b];
                    for (std::size_t i = 0; i <= a; ++i)
                        for (std::size_t j = 0; j <= b; ++j)
                            if (!(i == 0 && j == 0))
                                r = r - (binomial<T>(a, i) * binomial<T>(b, j) * w[i * q + j]) * n[(a - i) * q + (b - j)];
                    n[a * q + b] = r / w0;
                }
        }

        /**
         * P^(d) = C^(d) + sum C(d,i) f^(i) n^(d-i), M^(j) = normal(C^(j+1)). Where M vanishes
         * (|M| below 1e-10 of |C'| and of |C''| times the length of the bounds: a null tangent, or a
         * tangent along the direction in space), the point at an end of the curve whose tangent
         * vanishes (equal end poles) takes the normal of its limit from inside, along C'' at the
         * start and -C'' at the end; elsewhere, and for the derivatives, std::domain_error.
         */
        template <typename T, std::size_t dim, typename Func, typename Normal>
        auto offset_value(const Curve<T, dim> &c, const Func &f, T u, std::size_t d, Normal &&normal) -> point<T, dim>
        {
            std::vector<point<T, dim>> M(d + 1), n(d + 1);
            std::vector<T> w(d + 1);
            for (std::size_t j = 0; j <= d; ++j)
                M[j] = normal(c.value(u, j + 1));
            const auto [u1, u2] = c.bounds();
            const T l = std::isfinite(u2 - u1) && u2 > u1 ? u2 - u1 : T(1);
            const auto c1 = c.value(u, 1), c2 = c.value(u, 2);
            if (!(norm(M[0]) * l > T(1e-10) * std::max(norm(c1) * l, norm(c2) * l * l)))
            {
                if (d == 0 && norm(c1) * l <= T(1e-10) * norm(c2) * l * l && (u <= u1 || u >= u2))
                {
                    const auto m = normal(u <= u1 ? c2 : T(-1) * c2);
                    if (norm(m) > T(0))
                        return c.value(u) + f(u, 0) * (m / norm(m));
                }
                throw std::domain_error("offset curve: no normal (null tangent, or tangent along the direction)");
            }
            unit_derivatives<T, dim>(M, d, 0, w, n);
            auto p = c.value(u, d);
            for (std::size_t i = 0; i <= d; ++i)
                p = p + (binomial<T>(d, i) * f(u, i)) * n[d - i];
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
