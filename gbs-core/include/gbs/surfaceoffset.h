#pragma once
#include <gbs/bssurf.h>

namespace gbs
{

    template <typename T, size_t dim, typename Func>
    class SurfaceOffset : public Surface<T, dim>
    {
        const std::shared_ptr<Surface<T, dim>> p_srf_;
        const std::shared_ptr<Func> f_offset_;

    public:
        SurfaceOffset(const SurfaceOffset<T, dim, Func> &srf) = default;

        SurfaceOffset(const std::shared_ptr<Surface<T, dim>> &srf, const Func &f_offset) : p_srf_{srf},
                                                                                           f_offset_{std::make_shared<Func>(f_offset)}
        {
        }

        // CurveOffset(const BSSurface<T, dim> &srf, const Func &f_offset) : p_srf_{std::make_shared<BSSurface<T, dim>>(srf)},
        //                                                                 f_offset_{std::make_shared<Func>(f_offset)}
        // {
        // }

        virtual auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> std::array<T, dim> override
        {
            // switch (d)
            // {
            // case 0:
                return p_srf_->value(u) + normal_direction(*(p_srf_), u, v) * (*f_offset_)(u, v);
            //     break;
            // default:
            //     throw std::runtime_error("Not implemented yet.");
            //     break;
            // }
        }

        virtual auto bounds() const -> std::array<T, 2> override
        {
            return p_srf_->bounds();
        }

        auto changeBounds(const std::array<T, 2> &b) -> void
        {
            p_srf_->changeBounds(b);
        }

        auto basisSurface() const -> const Surface<T, dim> &
        {
            return *p_srf_;
        }
    };
}