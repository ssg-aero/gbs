#pragma once
#include <gbs/bscurve.h>


namespace gbs
{
    template <typename T, size_t dim>
    class CurveTrimmed : public Curve<T, dim>
    {
        const std::shared_ptr<Curve<T, dim>> p_crv_;
        std::array<T,2> bounds_;
        public:
        CurveTrimmed(const std::shared_ptr<Curve<T, dim>> &p_crv, T u1, T u2) : p_crv_{p_crv}, bounds_{u1,u2} {}
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            if (u < this->bounds()[0] - knot_eps<T>|| u > this->bounds()[1] + knot_eps<T>)
            {
                throw OutOfBoundsCurveEval(u,this->bounds());
            }
            return p_crv_->value(u,d);
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return bounds_;
        }
        const auto basisCurve() const noexcept {return p_crv_;}
    };
}