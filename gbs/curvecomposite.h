#pragma once

#include <cassert>
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
        virtual auto value(T u, size_t d = 0) const -> std::array<T, dim> override
        {
            // the curve whose piece holds u: the first one below the start (a parameter rounded
            // under it), the last one at the end and past it; never outside the list
            const auto n = static_cast<std::ptrdiff_t>(crv_lst_.size());
            const auto i = std::clamp<std::ptrdiff_t>(std::distance(u_.begin(), std::upper_bound(u_.begin(), u_.end(), u)) - 1, 0, n - 1);
            const auto &crv_loc = crv_lst_[static_cast<std::size_t>(i)];
            return crv_loc->value(u - u_[static_cast<std::size_t>(i)] + crv_loc->bounds()[0], d);
        }
        const auto & curves() const {return crv_lst_;}
    };
}
