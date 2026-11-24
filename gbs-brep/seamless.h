#pragma once

/**
 * @file seamless.h
 * @brief Native BREP core: a surface closing across a seam seen past it (`SeamlessSurface`).
 *
 * Made by the fillets round such a face (gbs-brep/fillet.h, architecture note
 * docs/sources/design/brep4_pr07_fillet_corners.md), and stored in the model as the surface of
 * their contact curves (a `CurveOnSurface`): a class of its own, outside the fillets, so that
 * what reads and writes a model (gbs-io/brep_json.h) knows it without the whole kernel.
 */

#include <gbs/bssurf.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>

namespace gbs::brep
{
    /**
     * A surface closing on itself across a seam without being periodic (a closed B-spline:
     * the skinned sections of a blade) evaluated past its seam: a parameter out of its bounds
     * brought into them by the period, the width of the bounds. The balls of a fillet along
     * an edge round such a face cross the seam: the centre line, the stations and the
     * contacts are computed on this view, their parameters running on continuously.
     */
    template <std::floating_point T>
    class SeamlessSurface : public Surface<T, 3>
    {
        std::shared_ptr<Surface<T, 3>> s_;
        std::array<T, 2> per_;

        auto in(T u, T v) const -> point<T, 2>
        {
            const auto b = s_->bounds();
            point<T, 2> w{u, v};
            for (std::size_t k = 0; k < 2; ++k)
                if (const T p = per_[k]; p > T(0) && (w[k] < b[2 * k] || w[k] > b[2 * k + 1]))
                {
                    w[k] = b[2 * k] + std::fmod(w[k] - b[2 * k], p);
                    if (w[k] < b[2 * k])
                        w[k] += p;
                }
            return w;
        }

    public:
        SeamlessSurface(std::shared_ptr<Surface<T, 3>> s, const std::array<T, 2> &periods) : s_{std::move(s)}, per_{periods} {}
        [[nodiscard]] auto basis() const -> const std::shared_ptr<Surface<T, 3>> & { return s_; }
        [[nodiscard]] auto periods() const -> std::array<T, 2> { return per_; }
        auto bounds() const -> std::array<T, 4> override { return s_->bounds(); }
        auto has_second_derivatives() const -> bool override { return s_->has_second_derivatives(); }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const auto w = in(u, v);
            return s_->value(w[0], w[1], du, dv);
        }
        auto derivatives(T u, T v, size_t nu, size_t nv, std::array<T, 3> *SKL) const -> void override
        {
            const auto w = in(u, v);
            s_->derivatives(w[0], w[1], nu, nv, SKL);
        }
    };
}
