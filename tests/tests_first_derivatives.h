#pragma once
/// A surface seen through its first derivatives only (`has_second_derivatives()` false, value()
/// throws beyond them), for the tests of the paths taken without second derivatives.
#include <gbs/bssurf.h>

#include <memory>
#include <stdexcept>

namespace gbs_tests
{
    template <std::floating_point T>
    class FirstDerivativesOnly : public gbs::Surface<T, 3>
    {
        std::shared_ptr<gbs::Surface<T, 3>> s_;

    public:
        explicit FirstDerivativesOnly(std::shared_ptr<gbs::Surface<T, 3>> s) : s_{std::move(s)} {}
        auto bounds() const -> std::array<T, 4> override { return s_->bounds(); }
        auto has_second_derivatives() const -> bool override { return false; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> gbs::point<T, 3> override
        {
            if (du + dv > 1)
                throw std::runtime_error("FirstDerivativesOnly: only first derivatives");
            return s_->value(u, v, du, dv);
        }
    };
}
