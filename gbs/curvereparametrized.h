#pragma once
#include <gbs/bscurve.h>

namespace gbs{

    template <typename T, size_t dim>
    class CurveReparametrized : public Curve<T,dim>
    {
        const std::shared_ptr<Curve<T,dim>>   p_crv_;
        const BSCfunction<T> f_u_;
    public:
        CurveReparametrized(const std::shared_ptr<Curve<T,dim>> &crv, const BSCfunction<T> &f_u) : p_crv_{crv}, f_u_{f_u} {}
        /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim>  override
        {
            switch (d)
            {
            case 0:
                return p_crv_->value(f_u_(u));
                break;
            case 1:
                return p_crv_->value(f_u_(u),1)*f_u_(u,1);
                break;
            case 2:
            {
                auto fd0 = f_u_(u);
                auto fd1 = f_u_(u, 1);
                return p_crv_->value(fd0, 1) * f_u_(u, 2) + p_crv_->value(fd0, 2) * fd1 * fd1;
            }
                break;
            default:
                throw std::exception("Not implemented yet.");
                break;
            }

        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return f_u_.bounds();
        }
    };
    }