#pragma once

#include <gbs/bscurve.h>
#include <algorithm>
namespace gbs
{
    template <typename T, size_t dim>
    class CurveComposite : public Curve<T, dim>
    {
        std::vector<std::shared_ptr<Curve<T,dim>>> crv_lst_;
        std::vector<T> u_;
        public:
        CurveComposite(const std::vector<std::shared_ptr<Curve<T,dim>>> &crv_lst) : crv_lst_{crv_lst} 
        {
            assert(crv_lst.size());
            auto n = crv_lst.size();
            u_ = std::vector<T>(n+1,crv_lst.front()->bounds()[0]);
            for( size_t i {}; i < n ; i++)
            {
                auto [u1,u2] = crv_lst_[i]->bounds();
                u_[i+1] = u_[i] + (u2-u1);
            }
        }
        virtual auto bounds() const -> std::array<T, 2> override
        {
            return {u_.front(),u_.back()};
        }
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim>
        {
            auto pos = std::upper_bound(u_.begin(), u_.end(), u);
            if (pos == u_.end())
            {
                pos = std::next(pos, -2);
            }
            else
            {
                pos = std::next(pos, -1);
            }
            auto i = std::distance(u_.begin(), pos);
            const auto &crv_loc = crv_lst_[i];
            auto u_loc = u - *pos + crv_loc->bounds()[0];
            return crv_loc->value(u_loc, d);
        }
        const auto & curves() const {return crv_lst_;}
    };
}
