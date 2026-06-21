#pragma once
#include <gbs/bssapprox.h>
#include <gbs/bssurf.h>
#include <gbs/execution.h>
namespace gbs
{
    template <typename T>
    auto offset_points(const BSSfunction<T> &f_ep, const Surface<T, 3> &cl_srf, const std::vector<T> &u, const std::vector<T> &v)
    {
        const size_t nu = u.size();
        const size_t nv = v.size();
        // Flatten the (u,v) grid in the existing U-first, V-outer order, then run
        // the offset kernel (value + du + dv + normal + offset) as one size-gated
        // parallel transform writing each node to its own slot (#89). Each node is
        // an independent surface eval, so the result is bit-identical to the old
        // nested push_back loop — just no longer serialized by a growing container.
        std::vector<std::array<T, 2>> uv(nu * nv);
        for (size_t j = 0; j < nv; ++j)
            for (size_t i = 0; i < nu; ++i)
                uv[j * nu + i] = {u[i], v[j]};

        points_vector<T, 3> pts(nu * nv);
        transform_threshold(
            uv.begin(), uv.end(), pts.begin(),
            [&f_ep, &cl_srf](const std::array<T, 2> &uv_) -> point<T, 3>
            {
                const T u_ = uv_[0], v_ = uv_[1];
                auto pt = cl_srf(u_, v_);
                auto tu = cl_srf(u_, v_, 1, 0);
                auto tv = cl_srf(u_, v_, 0, 1);
                auto n = tu ^ tv;
                n = n / norm(n);
                return pt + n * f_ep(u_, v_);
            });
        return pts;
    }

    template <typename T, bool rational>
    auto offset_approx(const BSSfunction<T> &f_ep, const BSSurfaceGeneral<T, 3, rational> &cl_srf,size_t nu , size_t nv )
    {
        auto [u1, u2, v1, v2] = cl_srf.bounds();
        auto u = make_range<T>(u1, u2, nu);
        auto v = make_range<T>(v1, v2, nv);

        auto pts = offset_points(f_ep,cl_srf,u,v);

        auto ku_flat = cl_srf.knotsFlatsU();
        auto kv_flat = cl_srf.knotsFlatsV();
        auto p = cl_srf.degreeU();
        auto q = cl_srf.degreeV();
        auto poles = approx(pts, ku_flat, kv_flat, u, v, p, q);

        using bs_type = typename std::conditional<rational, BSSurfaceRational<T, 3>, BSSurface<T, 3>>::type;

        return bs_type{
            poles,
            ku_flat,
            kv_flat,
            p,
            q};
    }
}
//occt 2d
  // P(u) = p(u) + Offset * Ndir / R
  // with R = || p' ^ Z|| and Ndir = P' ^ Z

  // P'(u)  = p'(u) + (Offset / R**2) * (DNdir/DU * R -  Ndir * (DR/R))

  // P"(u)  = p"(u) + (Offset / R) * (D2Ndir/DU - DNdir * (2.0 * Dr/ R**2) +
  //          Ndir * ( (3.0 * Dr**2 / R**4) - (D2r / R**2)))

  // P"'(u) = p"'(u) + (Offset / R) * (D3Ndir - (3.0 * Dr/R**2 ) * D2Ndir -
  //          (3.0 * D2r / R2) * DNdir) + (3.0 * Dr * Dr / R4) * DNdir -
  //          (D3r/R2) * Ndir + (6.0 * Dr * Dr / R4) * Ndir +
  //          (6.0 * Dr * D2r / R4) * Ndir - (15.0 * Dr* Dr* Dr /R6) * Ndir
//occt 3d
  // P(u) = p(u) + Offset * Ndir / R
  // with R = || p' ^ V|| and Ndir = P' ^ direction (local normal direction)

  // P'(u) = p'(u) + (Offset / R**2) * (DNdir/DU * R -  Ndir * (DR/R))

  // P"(u) = p"(u) + (Offset / R) * (D2Ndir/DU - DNdir * (2.0 * Dr/ R**2) +
  //         Ndir * ( (3.0 * Dr**2 / R**4) - (D2r / R**2)))

  //P"'(u) = p"'(u) + (Offset / R) * (D3Ndir - (3.0 * Dr/R**2) * D2Ndir -
  //         (3.0 * D2r / R2) * DNdir + (3.0 * Dr * Dr / R4) * DNdir -
  //         (D3r/R2) * Ndir + (6.0 * Dr * Dr / R4) * Ndir +
  //         (6.0 * Dr * D2r / R4) * Ndir - (15.0 * Dr* Dr* Dr /R6) * Ndir