#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscanalysis.h>

namespace gbs
{

    template <typename T, size_t dim, typename Func>
    class CurveOffset : public Curve<T, dim>
    {
        const std::shared_ptr<Curve<T, dim>> p_crv_;
        const std::shared_ptr< Func > f_offset_;

    public:
        CurveOffset(const CurveOffset<T,dim,Func>& crv) = default;
        CurveOffset(const std::shared_ptr<Curve<T, dim>> &crv, const Func&f_offset) : p_crv_{crv}, f_offset_{std::make_shared<Func>( f_offset )}
        {
            // f_offset_.changeBounds(p_crv_->bounds());
        }
        CurveOffset(const BSCurve<T, dim> &crv, const Func &f_offset) : 
            p_crv_{std::make_shared<BSCurve<T,dim>>(crv)}, 
            f_offset_{std::make_shared<Func>( f_offset )}
        {
            // f_offset_.changeBounds(p_crv_->bounds());
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
                // return p_crv_->value(u) + normal_direction(*(p_crv_), u) * f_offset_(u);
                return p_crv_->value(u) + normal_direction(*(p_crv_), u) * (*f_offset_)(u);
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
        auto changeBounds(const std::array<T, 2> &b) -> void
        {
            p_crv_->changeBounds(b);
        }
        auto basisCurve() const -> const Curve<T, dim> &
        {
            return *p_crv_;
        }
        
    };

    template <typename T, typename Func>
    class CurveOffset<T, 2,Func>: public Curve<T, 2>
    {
        const std::shared_ptr<Curve<T, 2>> p_crv_;
        const std::shared_ptr< Func > f_offset_;
        public:
        CurveOffset(const std::shared_ptr<Curve<T, 2>> &crv, const std::shared_ptr< Func > &f_offset) : p_crv_{crv}, f_offset_{ f_offset }
        { }
        CurveOffset(const std::shared_ptr<Curve<T, 2>> &crv, const Func &f_offset) : p_crv_{crv}, f_offset_{ std::make_shared<Func>( f_offset )}
        { }
        CurveOffset(const BSCurve<T, 2> &crv, const Func &f_offset) : 
            p_crv_{std::make_shared<BSCurve<T,2>>(crv)}, 
            f_offset_{std::make_shared<Func>( f_offset )}
        {
            // f_offset_.changeBounds(p_crv_->bounds());
        }
        virtual auto value(T u, size_t d = 0) const -> std::array<T, 2> override;
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_crv_->bounds();
        }
        auto changeBounds(const std::array<T, 2> &b) -> void
        {
            p_crv_->changeBounds(b);
        }
        auto basisCurve() const -> const Curve<T, 2> &
        {
            return *p_crv_;
        }
        auto offset() const -> const Func &
        {
            return *f_offset_;
        }
    };
    /**
         * @brief Curve evaluation at parameter u
         *
         * @param u : parameter on curve
         * @param d : derivative order
         * @return std::array<T, dim>
         */
    template <typename T, typename Func>
    auto CurveOffset<T, 2,Func>::value(T u, size_t d ) const -> std::array<T, 2>
    {
        const auto &crv = this->basisCurve();
        const auto &off = this->offset();
        switch (d)
        {
        case 0:
            return crv.value(u) + normal_direction(crv, u) * off(u,0);
            break;
        case 1:
            {

                auto d1 = crv.value(u, 1);
                auto d2 = crv.value(u, 2);
                auto r = sq_norm(d1);
                auto sqrt_r= std::sqrt(r);
                auto drdu = 2. * d1 * d2;
                auto o0 = off(u,0);
                auto o1 = off(u,1);

                auto x1 = d1[0];
                x1 -=  o0 * d2[1]/sqrt_r;
                x1 -= d1[1]*( o1*r - 0.5 * o0 * drdu ) / ( r * sqrt_r);

                auto y1 = d1[1];
                y1 +=  o0 * d2[0]/sqrt_r;
                y1 += d1[0]*( o1*r - 0.5 * o0 * drdu ) / ( r * sqrt_r);

                return {x1,y1};

            }
            break;
        case 2:
            {
                auto d1 = crv.value(u, 1);
                auto d2 = crv.value(u, 2);
                auto d3 = crv.value(u, 3);
                auto o0 = off(u,0);
                auto o1 = off(u,1);
                auto o2 = off(u,2);
                auto r = sq_norm(d1);
                auto sqrt_r= std::sqrt(r);
                auto drdu = 2.* d1 * d2;
                auto d2rdu2 = 2.* ((d2 * d2) + d1 * d3);

                auto x2 = d2[0];
                x2 -= (o2 * d1[1] + 2. * o1 * d2[1] + o0*d3[1]) / sqrt_r;
                x2 += (o1 * d1[1] * drdu + o0 * d2[1] * drdu + o0 * d1[1] * d2rdu2 /2.) / sqrt_r / r;
                x2 -= 3. / 4. * o0 * d1[1] * drdu * drdu / sqrt_r / r / r;

                auto y2 = d2[1];
                y2 += (o2 * d1[0] + 2. * o1 * d2[0] + o0* d3[0]) / sqrt_r;
                y2 -= (o1 * d1[0] * drdu + o0 * d2[0] * drdu + o0 * d1[0] * d2rdu2 /2.) / sqrt_r / r;
                y2 += 3. / 4. * o0 * d1[0] * drdu * drdu / sqrt_r / r / r;

                return {x2,y2};
            }
            break;
        default:
            throw std::exception("Not implemented yet.");
            break;
        }
    }
}