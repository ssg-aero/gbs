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
            assert(u>=bounds_[0] && u<=bounds_[1]);
            return p_crv_->value(u,d);
        }
        virtual auto bounds() const -> std::array<T, 2>
        {
            return bounds_;
        }
    };
}