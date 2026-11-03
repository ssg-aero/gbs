#pragma once
#include <gbs/bssanalysis.h>
#include <gbs/bssurf.h>
#include <memory>
#include <stdexcept>
namespace gbs
{
    /**
     * @brief The surface at a distance f_offset(u, v) from a surface of 3D, along its normal
     * (normal_direction: at a pole, the limit of the normal). Points only: the derivatives throw.
     */
    template <typename T, size_t dim, typename Func>
    class SurfaceOffset : public Surface<T, dim>
    {
        static_assert(dim == 3, "SurfaceOffset: the normal of a surface is defined in 3D");
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
            if (du != 0 || dv != 0)
                throw std::runtime_error("SurfaceOffset: derivatives not implemented yet");
            return p_srf_->value(u, v) + normal_direction(*p_srf_, u, v) * (*f_offset_)(u, v);
        }

        virtual auto has_second_derivatives() const -> bool override { return false; }

        virtual auto bounds() const -> std::array<T, 4> override
        {
            return p_srf_->bounds();
        }

        auto basisSurface() const -> const Surface<T, dim> &
        {
            return *p_srf_;
        }
    };
}
