#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscanalysis.h>

namespace gbs
{

    template <typename T, size_t dim>
    class CurveOffset : public Curve<T, dim>
    {
        const std::shared_ptr<Curve<T, dim>> p_crv_;
        BSCfunction<T> f_offset_;

    public:
        CurveOffset(const std::shared_ptr<Curve<T, dim>> &crv, const BSCfunction<T> &f_offset) : p_crv_{crv}, f_offset_{f_offset}
        {
            f_offset_.changeBounds(p_crv_->bounds());
        }
        /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            switch (d)
            {
            case 0:
                return p_crv_->value(u) + normal_direction(*(p_crv_), u) * f_offset_(u);
                break;
            default:
                throw std::exception("Not implemented yet.");
                break;
            }
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_crv_->bounds();
        }
        auto basisCurve() const -> const Curve<T, dim> &
        {
            return *p_crv_;
        }
    };

    template <typename T>
    class CurveOffset<T, 2>: public Curve<T, 2>
    {
        const std::shared_ptr<Curve<T, 2>> p_crv_;
        BSCfunction<T> f_offset_;
        public:
        CurveOffset(const std::shared_ptr<Curve<T, 2>> &crv, const BSCfunction<T> &f_offset) : p_crv_{crv}, f_offset_{f_offset}
        {
            f_offset_.changeBounds(p_crv_->bounds());
        }
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 2> override;
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_crv_->bounds();
        }
        auto basisCurve() const -> const Curve<T, 2> &
        {
            return *p_crv_;
        }
    };
    /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
    template <typename T>
    auto CurveOffset<T, 2>::value(T u, size_t d) const -> std::array<T, 2>
    {
        const auto &crv = this->basisCurve();
        switch (d)
        {
        case 0:
            return crv.value(u) + normal_direction(crv, u) * f_offset_(u);
            break;
        case 1:
        {

            // auto d0 = crv.value(u, 0);
            auto d1 = crv.value(u, 1);
            auto d2 = crv.value(u, 2);
            auto r = sq_norm(d1);
            auto sqrt_r= std::sqrt(r);
            auto n1 = point<T, 2>{ d1[1],-d1[0]};
            auto n2 = point<T, 2>{ d2[1],-d2[0]};
            return d1 +
                   n1 / sqrt_r * f_offset_(u, 1) +
                   (n2 / sqrt_r - n1 * (d1 * d2) / sqrt_r / r) * f_offset_(u);
        }
        break;
        default:
            throw std::exception("Not implemented yet.");
            break;
        }
    }
}