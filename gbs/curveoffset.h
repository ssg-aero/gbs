#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscanalysis.h>

namespace gbs
{

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
     * @brief A class representing an offset curve in 2D space
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
        const auto &crv = this->basisCurve();
        const auto &off = this->offset();
        switch (d)
        {
        case 0:
            return crv.value(u) + normal_direction(crv, u) * off(u, 0);
            break;
        case 1:
        {

            auto d1 = crv.value(u, 1);
            auto d2 = crv.value(u, 2);
            auto r = sq_norm(d1);
            auto sqrt_r = std::sqrt(r);
            auto drdu = 2. * d1 * d2;
            auto o0 = off(u, 0);
            auto o1 = off(u, 1);

            auto x1 = d1[0];
            x1 -= o0 * d2[1] / sqrt_r;
            x1 -= d1[1] * (o1 * r - 0.5 * o0 * drdu) / (r * sqrt_r);

            auto y1 = d1[1];
            y1 += o0 * d2[0] / sqrt_r;
            y1 += d1[0] * (o1 * r - 0.5 * o0 * drdu) / (r * sqrt_r);

            return {x1, y1};
        }
        break;
        case 2:
        {
            auto d1 = crv.value(u, 1);
            auto d2 = crv.value(u, 2);
            auto d3 = crv.value(u, 3);
            auto o0 = off(u, 0);
            auto o1 = off(u, 1);
            auto o2 = off(u, 2);
            auto r = sq_norm(d1);
            auto sqrt_r = std::sqrt(r);
            auto drdu = 2. * d1 * d2;
            auto d2rdu2 = 2. * ((d2 * d2) + d1 * d3);

            auto x2 = d2[0];
            x2 -= (o2 * d1[1] + 2. * o1 * d2[1] + o0 * d3[1]) / sqrt_r;
            x2 += (o1 * d1[1] * drdu + o0 * d2[1] * drdu + o0 * d1[1] * d2rdu2 / 2.) / sqrt_r / r;
            x2 -= 3. / 4. * o0 * d1[1] * drdu * drdu / sqrt_r / r / r;

            auto y2 = d2[1];
            y2 += (o2 * d1[0] + 2. * o1 * d2[0] + o0 * d3[0]) / sqrt_r;
            y2 -= (o1 * d1[0] * drdu + o0 * d2[0] * drdu + o0 * d1[0] * d2rdu2 / 2.) / sqrt_r / r;
            y2 += 3. / 4. * o0 * d1[0] * drdu * drdu / sqrt_r / r / r;

            return {x2, y2};
        }
        break;
        default:
            throw std::runtime_error("Not implemented yet.");
            break;
        }
    }

    /**
    * @brief A class representing an offset curve in 3D space
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
        const auto &crv = this->basisCurve();
        const auto &off = this->offset();

        switch (d)
        {
        case 0:
            {
                auto crv_val = crv.value(u);
                auto crv_d1 = crv.value(u, 1);
                auto tangent = normalized(crv_d1);
                auto normal = normalized(cross(direction_, tangent));
                return crv_val + normal * off(u, 0);
            }
            break;
        case 1:
            {
                auto crv_d1 = crv.value(u, 1);
                auto crv_d2 = crv.value(u, 2);
                auto tangent = normalized(crv_d1);
                auto normal = normalized(cross(direction_, tangent));
                auto curvature = cross(crv_d1, crv_d2) / (std::pow(norm(crv_d1), 3));
                auto normal_derivative = normalized(curvature * off(u, 0));
                auto off_d1 = off(u, 1);
                return crv_d1 + normal_derivative * off_d1;
            }
            break;
        case 2:
            {
                auto crv_d1 = crv.value(u, 1);
                auto crv_d2 = crv.value(u, 2);
                auto crv_d3 = crv.value(u, 3);
                auto tangent = normalized(crv_d1);
                auto normal = normalized(cross(direction_, tangent));
                auto curvature = cross(crv_d1, crv_d2) / (std::pow(norm(crv_d1), 3));
                auto normal_derivative = normalized(curvature * off(u, 0));

                auto d2t = (crv_d2 - ( crv_d2 * tangent) * tangent) / norm(crv_d1);
                auto d2n = cross(d2t, direction_) + cross(tangent, cross(crv_d2, direction_));
                auto off_d1 = off(u, 1);
                auto off_d2 = off(u, 2);

                return crv_d2 + d2n * off_d1 + normal_derivative * off_d2;
            }
            break;
        default:
            throw std::runtime_error("Not implemented yet.");
            break;
        }
    }

}