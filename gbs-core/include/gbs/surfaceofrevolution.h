#pragma once
#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/transform.h>
#include <numbers>

namespace gbs
{
    template <typename , size_t>
    class CurveOnSurface;

    template <typename T>
    class SurfaceOfRevolution : public Surface<T,3>
    {
        const std::shared_ptr<Curve<T, 2>> p_crv_;
        const ax1<T,3> ax_;
        const ax2<T,3> ax2_;
        std::array<T,2> angle_span_;
        Matrix4<T> M_;

    public:
        SurfaceOfRevolution(const SurfaceOfRevolution<T> &srf) = default;
        SurfaceOfRevolution(const std::shared_ptr<Curve<T, 2>> &crv, const ax2<T, 3> &ax, T a1 = 0., T a2 = 2. * std::numbers::pi) : p_crv_{crv}, ax2_{ax}, ax_{{ax[0],ax[1]/norm(ax[1])}}, angle_span_{{a1, a2}}
        {
            T i = 1.;
            T o = 0.;
            M_ = build_trf_matrix<T>(
                {{
                    {o,o,o},
                    {o,i,o},
                    {o,o,i}
                }}
                    ,
                ax
            );
        }
        SurfaceOfRevolution(const BSCurve<T, 2> &crv, const ax2<T, 3> &ax, T a1 = 0., T a2 = 2. * std::numbers::pi) : 
            SurfaceOfRevolution{std::make_shared<BSCurve<T, 2>>(crv), ax, a1, a2}
        {        }
        SurfaceOfRevolution(const CurveOnSurface<T, 2> &crv, const ax2<T, 3> &ax, T a1 = 0., T a2 = 2. * std::numbers::pi) : 
            SurfaceOfRevolution{std::make_shared<CurveOnSurface<T, 2>>(crv), ax, a1, a2}
        {        }
        /**
         * @brief  Surface evaluation at parameters {u,v}
         * 
         * @param u  : u parameter on surface
         * @param v  : v tangential position
         * @param du : u derivative order
         * @param dv : v derivative order
         * @return point<T, dim> const 
         */
        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            point<T,3> pt;
            v = fmod(v,2*pi);
            if(dv==0)
            {
                auto pt2d = p_crv_->value(u,du) ;
                pt = add_dimension<T,2>(pt2d,0.);
                transform(pt,M_);
            }
            else
            {
                auto pt2d = p_crv_->value(u, du);
                pt = add_dimension<T, 2>(pt2d, 0.);
                transform(pt,M_);
                if (du == 0)
                {
                    pt = pt - ax_[0];
                }

                pt = ax_[1] ^ pt;
                rotate<T>(pt,(dv-1)*pi/2.,ax_[1]);

            }
            return rotated<T>( pt , v, ax_[1]);
            
        }
        /**
         * @brief return surface's bounds {U1,U2,V1,V2}
         * 
         * @return point<T, dim>
         */
        virtual auto bounds() const -> std::array<T, 4> override
        {
            auto [u1,u2] = p_crv_->bounds();
            auto [v1,v2] = angle_span_;
            return { u1, u2, v1, v2};
        }

        auto operator()(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> { return value(u, v, du, dv); };

        const auto basisCurve() const noexcept {return p_crv_;}
        const auto & transformation() const noexcept {return M_;}
        const auto & axis() const noexcept {return ax_;}
        const auto & axis2() const noexcept {return ax2_;}

    };
}