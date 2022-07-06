#pragma once
#include <gbs/bssapprox.h>
#include <gbs/bssurf.h>
namespace gbs
{
    template <typename T>
    auto offset_points(const gbs::BSSfunction<T> &f_ep, const gbs::Surface<T, 3> &cl_srf, const std::vector<T> &u, const std::vector<T> &v)
    {
        size_t nu = u.size();
        size_t nv = v.size();
        gbs::points_vector<T, 3> pts;
        T u_, v_;
        for (auto j = 0; j < nv; j++)
        {
            v_ = v[j];
            for (auto i = 0; i < nu; i++)
            {
                u_ = u[i];
                auto pt = cl_srf(u_, v_);
                auto tu = cl_srf(u_, v_, 1, 0);
                auto tv = cl_srf(u_, v_, 0, 1);
                auto n = tu ^ tv;
                n = n / gbs::norm(n);
                pt = pt + n * f_ep(u_, v_);
                pts.push_back(pt);
            }
        }
        return pts;
    }

    template <typename T, bool rational>
    auto offset_approx(const gbs::BSSfunction<T> &f_ep, const gbs::BSSurfaceGeneral<T, 3, rational> &cl_srf,size_t nu , size_t nv )
    {
        auto [u1, u2, v1, v2] = cl_srf.bounds();
        auto u = gbs::make_range<T>(u1, u2, nu);
        auto v = gbs::make_range<T>(v1, v2, nv);

        auto pts = offset_points(f_ep,cl_srf,u,v);

        auto ku_flat = cl_srf.knotsFlatsU();
        auto kv_flat = cl_srf.knotsFlatsV();
        auto p = cl_srf.degreeU();
        auto q = cl_srf.degreeV();
        auto poles = gbs::approx(pts, ku_flat, kv_flat, u, v, p, q);

        using bs_type = typename std::conditional<rational, gbs::BSSurfaceRational<T, 3>, gbs::BSSurface<T, 3>>::type;

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