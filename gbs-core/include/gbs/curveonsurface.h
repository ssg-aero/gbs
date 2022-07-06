#pragma once
#include <gbs/curves>
#include <gbs/bssurf.h>

namespace gbs
{

    template <typename T, size_t dim>
    class CurveOnSurface : public Curve<T, dim>
    {
        const std::shared_ptr<Curve<T, 2>> p_crv_;
        const std::shared_ptr<Surface<T, dim>> p_srf_;

    public:
        CurveOnSurface(const CurveOnSurface<T,dim> &crv) = default;
        CurveOnSurface(const std::shared_ptr<Curve<T, 2>> &crv, const std::shared_ptr<Surface<T, dim>>& srf) : p_crv_{crv}, p_srf_{srf}
        {
            
        }
        CurveOnSurface(const BSCurve<T, 2> &crv, const BSSurface<T, dim>& srf) : p_crv_{std::make_shared<BSCurve<T, 2>>( crv)}, p_srf_{std::make_shared<BSSurface<T, dim>>( srf )}
        {
            
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
            auto [u_, v_] = p_crv_->value(u);
            switch (d)
            {
            case 0:
                return p_srf_->value(u_,v_);
                break;
            case 1:
                {
                    auto [du_, dv_] = p_crv_->value(u,1);
                    auto D1U = p_srf_->value(u_,v_,1,0);
                    auto D1V = p_srf_->value(u_,v_,0,1);
                    return du_*D1U + dv_*D1V;
                }
                break;
            case 2:
                {
                    auto [du_, dv_] = p_crv_->value(u,1);
                    auto [d2u_, d2v_] = p_crv_->value(u,2);
                    auto D1U = p_srf_->value(u_,v_,1,0);
                    auto D1V = p_srf_->value(u_,v_,0,1);
                    auto D2U = p_srf_->value(u_,v_,2,0);
                    auto D2V = p_srf_->value(u_,v_,0,2);
                    auto D2UV = p_srf_->value(u_,v_,1,1);
                    auto V2 = d2u_ * D1U + d2u_ * D1V + 2 * du_ * dv_ * D2UV;
                    return du_ * du_ * D2U + dv_ * dv_ * D2V + V2;
                }
                break;
            default:
                throw std::runtime_error("Not implemented yet.");
                break;
            }
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_crv_->bounds();
        }
        auto basisCurve() const -> const Curve<T, 2> &
        {
            return *p_crv_;
        }
        auto basisSurface() const -> const Surface<T, dim> &
        {
            return *p_srf_;
        }
    };
}