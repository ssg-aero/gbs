#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <gbs/bssurf.h>
#include <gbs/bscapprox.h>
#include <gbs/transformpoints.h>
#include <gbs/bsctools.h>
namespace gbs 
{
    template <std::floating_point T, size_t dim>
    void unify_knots(std::vector<std::pair<T, size_t>> &km_out_set, std::vector<T> &k_out_set, std::vector<std::vector<std::array<T, dim>>> &poles_out_set, size_t degree, const std::vector<std::pair<T, size_t>> &km_in)
    {
        std::vector<std::pair<T, size_t>> km_out;
        std::vector<T> k_out;
        for(auto &poles_out: poles_out_set )
        {
            // Copy to preserve
            km_out = std::vector<std::pair<T, size_t>>{km_out_set};
            k_out  = std::vector<T>{k_out_set};
            unify_knots(km_out, k_out, poles_out, degree,km_in);
        }
        // Use last updated
        // std::ranges::copy(km_out_set, km_out.begin());
        // std::ranges::copy(k_out_set,  k_out.begin() );

        km_out_set = km_out;
        k_out_set = k_out;

    }

    template <std::floating_point T, size_t dim>
    auto unify_knots(std::vector<BSSurfaceInfo<T, dim >> &surfaces_info) -> void {
        // Assume that the first surface's info is used as a reference
        auto &[polesUV0, knotsU0, knotsV0, degU0, degV0] = surfaces_info.front();

        // Determine the start and end of the knot vectors in both U and V directions
        auto kU_start = knotsU0.front();
        auto kU_end = knotsU0.back();
        auto kV_start = knotsV0.front();
        auto kV_end = knotsV0.back();

        // Convert the first surface's knot vectors to a multiplicity format (not shown, similar to unflat_knots)
        auto kmU0 = unflat_knots(knotsU0);
        auto kmV0 = unflat_knots(knotsV0);

        // Adjust bounds and unify knots for each surface with the first surface
        std::for_each(
            std::next(surfaces_info.begin()),
            surfaces_info.end(),
            [kU_start, kU_end, kV_start, kV_end, degU0, degV0, &knotsU0, &kmU0, &knotsV0, &kmV0, &polesUV0](auto &surface_info)
            {
                auto &[polesUV, knotsU, knotsV, degU, degV] = surface_info;
                change_bounds(kU_start, kU_end, knotsU);
                change_bounds(kV_start, kV_end, knotsV);

                auto kmU = unflat_knots(knotsU);
                auto kmV = unflat_knots(knotsV);

                // Unify knots in U and V directions (not shown, similar to unify_knots)
                unify_knots(kmU0, knotsU0, polesUV0, degU0, kmU);
                // TODO transpose to polesVU0
                auto polesVU0 = transpose_poles(polesUV0);
                unify_knots(kmV0, knotsV0, polesVU0, degV0, kmV);
                polesUV0 = transpose_poles(polesVU0);
            });

        // Apply the unified knot vectors to all surfaces
        std::for_each(
            std::next(surfaces_info.begin()),
            surfaces_info.end(),
            [&](auto &surface_info)
            {
                auto &[polesUV, knotsU, knotsV, degU, degV] = surface_info;
                auto kmU = unflat_knots(knotsU);
                auto kmV = unflat_knots(knotsV);

                unify_knots(kmU, knotsU, polesUV, degU, kmU0);
                auto polesVU = transpose_poles(polesUV);
                unify_knots(kmV, knotsV, polesVU, degV, kmV0);
                polesUV = transpose_poles(polesVU);
            });
    }

    /**
     * @brief Create a surface's copy with and additional dimension. The value of the dimension can be specified, the default is 0.
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational 
     * @param srf 
     * @param val 
     * @return auto 
    **/
    template <typename T, size_t dim, bool rational>
    auto add_dimension(const BSSurfaceGeneral<T, dim, rational> &srf, T val = 0.)
    {
        points_vector<T, dim + rational + 1> poles(srf.poles().size());
        std::transform(
            GBS_PAR_EXEC
            srf.poles().begin(),
            srf.poles().end(),
            poles.begin(),
            [&val](const auto &p_)
            {
                if constexpr (rational)
                    return add_dimension_homogeneous<T, dim + rational>(p_, val); // the weight stays the last coordinate
                else
                    return add_dimension<T, dim>(p_, val);
            });
        using bs_type = typename std::conditional<rational, BSSurfaceRational<T, dim + 1>, BSSurface<T, dim + 1>>::type;
        return bs_type(poles, srf.knotsFlatsU(), srf.knotsFlatsV(), srf.degreeU(), srf.degreeV());
    }

    /// @brief Extend surface to match curve
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param crv 
    /// @param natural_end 
    /// @param max_cont 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto extention_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        points_vector<T,dim+rational> poles_new;

        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        using bss_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type;

        auto nu = srf.nPolesU();
        auto nv = srf.nPolesV();
        auto ku = srf.knotsFlatsU();
        auto kv = srf.knotsFlatsV();
        auto p = srf.degreeU();
        auto q = srf.degreeV();
        auto du = 0.;
        auto np = 100;
        auto v1 = kv.front();
        auto v2 = kv.back();
        auto u2 = ku.back();

        for(size_t i {} ; i < np ; i++)
        {
            auto v = v1 + (v2-v1)*i/(np-1.);
            auto t = srf(u2,v,1,0);
            auto dl = norm(srf(u2,v)-crv(v));
            du += dl / norm(t) / np;
        }
        auto u_new = u2 + du;

        size_t p_ext;
        std::vector<T> k_ext;
        for(size_t i {} ; i < nv ; i++)
        {
            auto poles = points_vector<T,dim+rational>{srf.poles_begin()+i*nu,srf.poles_begin()+i*nu+nu};
            auto iso  = bsc_type{poles,ku,p};
            auto crv_ext = bsc_type{extention_to_point<T, dim, rational>(
                iso,
                weight_projected_pole<T, dim, rational>(crv.poles()[i]),
                u2,
                u_new,
                natural_end,
                max_cont
            )};
            auto poles_ext = crv_ext.poles();
            poles_new.insert(poles_new.end(),poles_ext.begin(),poles_ext.end());
            if(i==0)
            {
                p_ext = crv_ext.degree();
                k_ext = crv_ext.knotsFlats();
            }
        }

        return bss_type(poles_new,k_ext,kv,p_ext,q) ;
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const BSSurfaceGeneral<T,dim,rational1> &srf1,const BSSurfaceGeneral<T,dim,rational2> &srf2)
    {
        const auto rational = rational1 || rational2;
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        bs_type srf1_ { srf1.poles(), srf1.knotsFlatsU(), srf1.knotsFlatsV(), srf1.degreeU(), srf1.degreeV() };
        bs_type srf2_ { srf2.poles(), srf2.knotsFlatsU(), srf2.knotsFlatsV(), srf2.degreeU(), srf2.degreeV() };

        while(srf1_.degreeU()<srf2_.degreeU()) srf1_.increaseDegreeU();
        while(srf2_.degreeU()<srf1_.degreeU()) srf2_.increaseDegreeU();

        auto nkv = srf2_.knotsFlatsV().size(); // TODO insert v knots
        if(nkv!=srf1_.knotsFlatsV().size())
            throw std::runtime_error("Not implemented");
        for(size_t i{}; i < nkv ; i ++)
        {
            if(std::fabs( srf2_.knotsFlatsV()[i] - srf1_.knotsFlatsV()[i]) > knot_eps<T>)
                throw std::runtime_error("Not implemented");
        }
        // join
        auto nu = srf1_.nPolesU();
        auto nv = srf1_.nPolesV();
        auto poles_ext = srf2_.poles();
        auto nu_ext = srf2_.nPolesU();
        auto ku_ext = srf2_.knotsFlatsU();
        auto ku = srf1_.knotsFlatsU();
        auto kv = srf1_.knotsFlatsV();
        auto p = srf1_.degreeU();
        auto q = srf1_.degreeV();
        auto p_ext = srf2_.degreeU();

        points_vector<T,dim+rational> poles_new;
        for(size_t i{}; i < nv; i++)
        {
            poles_new.insert(
                poles_new.end(),
                std::next( srf1_.poles().begin(), i * nu ),
                std::next( srf1_.poles().begin(), i * nu + nu)
            );
            poles_new.insert(
                poles_new.end(),
                std::next( poles_ext.begin(), i * nu_ext +1 ),
                std::next( poles_ext.begin(), i * nu_ext + nu_ext)
            );
        }
        auto dk = ku.back() - ku_ext.front();
        std::transform(
            ku_ext.begin(),
            ku_ext.end(),
            ku_ext.begin(),
            [dk](auto k){return k+dk;}
        );
        std::vector<double> k_new {ku.begin(),std::next(ku.end(),-1)};
        k_new.insert(
            k_new.end(),
            std::next(ku_ext.begin(),p_ext+1),
            ku_ext.end()
        );

        return bs_type{
            poles_new,
            k_new,
            kv,
            p,
            q
        };
    }

    /// @brief Split surface at u parameter, aka in v direction
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param u 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto split_u(const BSSurfaceGeneral<T,dim,rational> &srf, T u)
    {
        auto [u1, u2, v1, v2] = srf.bounds();
        assert ( u1<=u && u <=u2);

        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 

        bs_type srf1{srf};
        bs_type srf2{srf};

        srf1.trimU(u1, u);
        srf2.trimU(u, u2);

        return std::make_pair(
            srf1, srf2
        );
    }
    /// @brief Split surface at v parameter, aka in u direction
    /// @tparam T 
    /// @tparam dim 
    /// @tparam rational 
    /// @param srf 
    /// @param v 
    /// @return 
    template <typename T, size_t dim, bool rational>
    auto split_v(const BSSurfaceGeneral<T,dim,rational> &srf, T v)
    {
        auto [u1, u2, v1, v2] = srf.bounds();
        assert ( v1<=v && v <=v2);

        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 

        bs_type srf1{srf};
        bs_type srf2{srf};

        srf1.trimV(v1, v);
        srf2.trimV(v, v2);

        return std::make_pair(
            srf1, srf2
        );
    }

    template <typename T, size_t dim, bool rational>
    auto extended_to_curve(const BSSurfaceGeneral<T,dim,rational> &srf,const BSCurveGeneral<T,dim,rational> &crv, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto srf_ext = extention_to_curve(srf,crv,natural_end,max_cont);
        return join(srf,srf_ext);
    }

    template <typename T, size_t dim, bool rational>
    auto extention(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        using bsc_type = typename std::conditional<rational,BSCurveRational<T, dim>,BSCurve<T, dim>>::type;
        auto np = srf.nPolesV()*10;
        auto [u1,u2,v1,v2] = srf.bounds();
        points_vector<T,dim> pts(np);
        std::vector<T> v_crv(np);
        for(size_t i {} ; i < np ; i++)
        {
            auto v_ = v1 + (v2-v1)*i/(np-1.);
            v_crv[i] = v_;
            auto pt = srf(u2,v_);
            auto tg = srf(u2,v_,1,0);
            tg = tg / norm(tg);
            pt = pt + tg * l_ext;
            pts[i] = pt;
        }
        auto nv = srf.nPolesV();
        auto kv = srf.knotsFlatsV();
        auto q  = srf.degreeV();
        auto crv = bsc_type{approx_bound_fixed(pts,q,nv,v_crv,kv)};
        return extention_to_curve(srf,crv,natural_end,max_cont);
    }

    template <typename T, size_t dim, bool rational>
    auto extended(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        auto srf_ext = extention(srf,l_ext,natural_end,max_cont);
        return join(srf,srf_ext);
    }

    enum class SurfaceBound { U_START, U_END, V_START, V_END};

    template <typename T, size_t dim, bool rational>
    auto extention(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, SurfaceBound location , bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        switch (location)
        {
        case SurfaceBound::U_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.reverseU();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break;
            case SurfaceBound::V_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                srf_.reverseU();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break;        
            case SurfaceBound::V_END:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                return extention(srf_,l_ext,natural_end, max_cont);
            }
            break; 
        default: // SurfaceBound::U_END
            return extention(srf,l_ext,natural_end, max_cont);
            break;
        }
    }
    template <typename T, size_t dim, bool rational>
    auto extended(const BSSurfaceGeneral<T,dim,rational> &srf,T l_ext, SurfaceBound location, bool natural_end , std::optional<size_t> max_cont = std::nullopt)
    {
        // auto srf_ext = extention(srf,l_ext,location,natural_end,max_cont);
        // return join(srf,srf_ext);
        using bs_type = typename std::conditional<rational,BSSurfaceRational<T, dim>,BSSurface<T, dim>>::type; 
        switch (location)
        {
        case SurfaceBound::U_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.reverseU();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.reverseU();
                return srf_ext;
            }
            break;
            case SurfaceBound::V_START:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                srf_.reverseU();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.reverseU();
                srf_ext.invertUV();
                return srf_ext;
            }
            break;        
            case SurfaceBound::V_END:
            {
                bs_type srf_ {srf.poles(),srf.knotsFlatsU(),srf.knotsFlatsV(),srf.degreeU(),srf.degreeV()};
                srf_.invertUV();
                auto srf_ext = extended(srf_,l_ext,natural_end, max_cont);
                srf_ext.invertUV();
                return srf_ext;
            }
            break; 
        default: // SurfaceBound::U_END
            return extended(srf,l_ext,natural_end, max_cont);
            break;
        }
    }

    namespace detail
    {
        /// The polar form (blossom) of the polynomial of the span [k[s], k[s + 1]] of the B-spline of
        /// poles `p`, knots `k` and degree q, at the q values `t`: de Boor's algorithm with the value
        /// t[r - 1] at its level r.
        template <typename T, size_t d>
        auto blossom(const std::vector<std::array<T, d>> &p, const std::vector<T> &k, size_t q, size_t s,
                     const std::vector<T> &t) -> std::array<T, d>
        {
            std::vector<std::array<T, d>> b(p.begin() + std::ptrdiff_t(s - q), p.begin() + std::ptrdiff_t(s + 1));
            for (size_t r = 1; r <= q; ++r)
                for (size_t j = q; j >= r; --j)
                {
                    const size_t i = s - q + j;
                    const T a = (t[r - 1] - k[i]) / (k[i + q + 1 - r] - k[i]);
                    for (size_t c = 0; c < d; ++c)
                        b[j][c] = (T(1) - a) * b[j - 1][c] + a * b[j][c];
                }
            return b[q];
        }

        /// The clamped knots `k` (degree q) with those of their end (`at_end`) or of their start moved to `e`.
        template <typename T>
        auto moved_end(std::vector<T> k, size_t q, bool at_end, T e) -> std::vector<T>
        {
            if (at_end)
                std::fill(k.end() - std::ptrdiff_t(q + 1), k.end(), e);
            else
                std::fill(k.begin(), k.begin() + std::ptrdiff_t(q + 1), e);
            return k;
        }

        /// The poles `p` of the clamped B-spline (p, k, q) made those of the same curve on the knots
        /// `moved` (`moved_end`): the q poles next to that end become the blossom of the polynomial
        /// of the end span at the moved knots, the others are kept. The curve is unchanged within its
        /// bounds, and its end span goes on as the same polynomial: no knot is inserted.
        template <typename T, size_t d>
        void continued_end(std::vector<std::array<T, d>> &p, const std::vector<T> &k, const std::vector<T> &moved,
                           size_t q, bool at_end)
        {
            const size_t n = p.size();
            const size_t s = at_end ? n - 1 : q;         // the end span [k[s], k[s + 1]]
            const size_t first = at_end ? n - q : 0;     // the poles [first, first + q) that change
            std::vector<std::array<T, d>> changed(q);
            std::vector<T> t(q);
            for (size_t i = first; i < first + q; ++i)
            {
                // pole i is the blossom at the knots i + 1 .. i + q
                std::copy_n(moved.begin() + std::ptrdiff_t(i + 1), q, t.begin());
                changed[i - first] = blossom(p, k, q, s, t);
            }
            std::ranges::copy(changed, p.begin() + std::ptrdiff_t(first));
        }

        /// The smallest value of f on [a, b]: at n + 1 points, then by golden section between the
        /// neighbours of the smallest.
        template <typename T, typename F>
        auto smallest_on(const F &f, T a, T b, size_t n) -> T
        {
            if (n == 0)
                throw std::invalid_argument("smallest_on: at least one interval.");
            size_t best = 0;
            T smallest = std::numeric_limits<T>::max();
            for (size_t i = 0; i <= n; ++i)
            {
                const T fi = f(a + (b - a) * T(i) / T(n));
                if (fi < smallest)
                {
                    smallest = fi;
                    best = i;
                }
            }
            T lo = a + (b - a) * T(best > 0 ? best - 1 : 0) / T(n);
            T hi = a + (b - a) * T(std::min(best + 1, n)) / T(n);
            const T g = (std::sqrt(T(5)) - T(1)) / T(2);
            T x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
            T f1 = f(x1), f2 = f(x2);
            for (int it = 0; it < 64; ++it)
            {
                if (f1 < f2)
                {
                    hi = x2, x2 = x1, f2 = f1;
                    x1 = hi - g * (hi - lo), f1 = f(x1);
                }
                else
                {
                    lo = x1, x1 = x2, f1 = f2;
                    x2 = lo + g * (hi - lo), f2 = f(x2);
                }
            }
            return std::min({smallest, f1, f2});
        }

        template <typename S>
        struct bs_surface_traits;
        template <typename T, size_t dim>
        struct bs_surface_traits<BSSurface<T, dim>>
        {
            using value_type = T;
        };
        template <typename T, size_t dim>
        struct bs_surface_traits<BSSurfaceRational<T, dim>>
        {
            using value_type = T;
        };
    } // namespace detail

    /**
     * @brief `srf` gone on past its ends across v by the parameter `dv`: before v1 by dv[0], after
     * v2 by dv[1] (a dv <= 0 leaves that end). Each end span across v goes on as the same
     * polynomial: its end knots are moved and the q poles next to them made the blossom of that
     * polynomial at the moved knots, the other poles and knots kept, no knot inserted. The surface
     * is unchanged within its bounds; a rational one is continued in its homogeneous poles.
     *
     * Each column of poles is continued on its own: surfaces sharing a boundary (a column of
     * poles) and their v-knots, continued by the same `dv`, still share it.
     * @throw std::invalid_argument a dv not finite, v-knots not clamped
     */
    template <typename T, size_t dim, bool rational>
    auto extrapolated_v(const BSSurfaceGeneral<T, dim, rational> &srf, const std::array<T, 2> &dv)
        -> std::conditional_t<rational, BSSurfaceRational<T, dim>, BSSurface<T, dim>>
    {
        if (!std::isfinite(dv[0]) || !std::isfinite(dv[1]))
            throw std::invalid_argument("extrapolated_v: the parameters must be finite.");
        const size_t q = srf.degreeV(), nu = srf.nPolesU(), nv = srf.nPolesV();
        const auto &k0 = srf.knotsFlatsV();
        const auto clamped = [&](auto first) { return std::all_of(first, first + std::ptrdiff_t(q + 1), [&](T k) { return k == *first; }); };
        if (!clamped(k0.begin()) || !clamped(k0.end() - std::ptrdiff_t(q + 1)))
            throw std::invalid_argument("extrapolated_v: the v-knots must be clamped.");
        const bool after = dv[1] > T(0), before = dv[0] > T(0);
        const auto k1 = after ? detail::moved_end(k0, q, true, k0.back() + dv[1]) : k0;
        const auto k2 = before ? detail::moved_end(k1, q, false, k1.front() - dv[0]) : k1;
        auto poles = srf.poles();
        std::vector<std::array<T, dim + rational>> col(nv);
        for (size_t i = 0; i < nu; ++i)
        {
            for (size_t j = 0; j < nv; ++j)
                col[j] = poles[i + j * nu];
            if (after)
                detail::continued_end(col, k0, k1, q, true);
            if (before)
                detail::continued_end(col, k1, k2, q, false);
            for (size_t j = 0; j < nv; ++j)
                poles[i + j * nu] = col[j];
        }
        return {poles, srf.knotsFlatsU(), k2, srf.degreeU(), q};
    }

    /**
     * @brief The smallest speed |dS/dv| along the iso-v line `v` of `srf`: at `n` + 1 points across
     * u, then refined between the neighbours of the slowest.
     * @throw std::invalid_argument n == 0
     */
    template <typename T, size_t dim>
    auto slowest_across_v(const Surface<T, dim> &srf, T v, size_t n = 64) -> T
    {
        const auto [u1, u2, v1, v2] = srf.bounds();
        return detail::smallest_on([&](T u) { return norm(srf.value(u, v, 0, 1)); }, u1, u2, n);
    }

    namespace detail
    {
        /**
         * The parameter by which `surfaces` (sharing their v-knots) go on past their end
         * (`at_end`) or their start for every point of it to get `length` away at least: from the
         * slowest point (length / its speed), raised while the nearest end point of the continued
         * surfaces falls short. Throws when a point of the end has no speed across v, when the
         * continuation does not get that far, or when it turns back on the way.
         */
        template <typename S, typename T>
        auto v_extension_for_length(const std::vector<S> &surfaces, bool at_end, T length) -> T
        {
            constexpr size_t n = 64, m = 8; // points across u, along the continuation
            const T sign = at_end ? T(1) : T(-1);
            const auto end_of = [&](const S &s) { return s.bounds()[at_end ? 3 : 2]; };
            const auto continued = [&](const S &s, T dv) {
                std::array<T, 2> d{};
                d[at_end ? 1 : 0] = dv;
                return extrapolated_v(s, d);
            };
            T slowest = std::numeric_limits<T>::max(), fastest{};
            for (const auto &s : surfaces)
            {
                const auto [u1, u2, v1, v2] = s.bounds();
                slowest = std::min(slowest, slowest_across_v(s, end_of(s), n));
                for (size_t i = 0; i <= n; ++i)
                    fastest = std::max(fastest, norm(s.value(u1 + (u2 - u1) * T(i) / T(n), end_of(s), 0, 1)));
            }
            if (!(slowest > std::numeric_limits<T>::epsilon() * fastest))
                throw std::invalid_argument("extrapolated_v_by_lengths: an end with a point without speed across v cannot go on by a length.");
            T dv = length / slowest;
            for (int it = 0; it < 100; ++it)
            {
                T reach = std::numeric_limits<T>::max();
                for (const auto &s : surfaces)
                {
                    const auto c = continued(s, dv);
                    const auto [u1, u2, v1, v2] = s.bounds();
                    const T ve = end_of(s);
                    reach = std::min(reach, smallest_on([&](T u) { return distance(c.value(u, ve + sign * dv), s.value(u, ve)); }, u1, u2, n));
                }
                if (reach >= length)
                {
                    // getting away all along the way: (S - S(ve)) . dS/dv along the continuation
                    for (const auto &s : surfaces)
                    {
                        const auto c = continued(s, dv);
                        const auto [u1, u2, v1, v2] = s.bounds();
                        const T ve = end_of(s);
                        for (size_t i = 0; i <= n; ++i)
                        {
                            const T u = u1 + (u2 - u1) * T(i) / T(n);
                            const auto p0 = s.value(u, ve);
                            for (size_t j = 1; j <= m; ++j)
                            {
                                const T v = ve + sign * dv * T(j) / T(m);
                                if (!(sign * ((c.value(u, v) - p0) * c.value(u, v, 0, 1)) > T(0)))
                                    throw std::invalid_argument("extrapolated_v_by_lengths: the continuation of an end turns back on the way.");
                            }
                        }
                    }
                    return dv;
                }
                if (!(reach > T(0)))
                    throw std::invalid_argument("extrapolated_v_by_lengths: the continuation of an end turns back on the way.");
                dv *= length / reach * (T(1) + T(1e-9));
            }
            throw std::invalid_argument("extrapolated_v_by_lengths: the continuation of an end does not get that far.");
        }
    } // namespace detail

    /**
     * @brief `surfaces` (sharing their v-knots) gone on together past their ends across v, every
     * point of an end at least `lengths` away from where it was (before v1, after v2; 0: that end
     * left): extrapolated_v by the same parameter, from the one taking the slowest point of their
     * ends (slowest_across_v) that far, raised until the nearest point gets there (a continuation
     * slowing down). Surfaces sharing a boundary still share it.
     * @throw std::invalid_argument surfaces whose v-knots differ, a length not finite, an end with
     * a point without speed across v, a continuation that does not get that far or turns back on
     * the way
     */
    template <typename S>
    auto extrapolated_v_by_lengths(std::vector<S> surfaces,
                                   const std::array<typename detail::bs_surface_traits<S>::value_type, 2> &lengths)
        -> std::vector<S>
    {
        using T = typename detail::bs_surface_traits<S>::value_type;
        for (const auto &s : surfaces)
            if (s.knotsFlatsV() != surfaces.front().knotsFlatsV() || s.degreeV() != surfaces.front().degreeV())
                throw std::invalid_argument("extrapolated_v_by_lengths: the surfaces must share their v-knots.");
        if (!std::isfinite(lengths[0]) || !std::isfinite(lengths[1]))
            throw std::invalid_argument("extrapolated_v_by_lengths: the lengths must be finite.");
        std::array<T, 2> dv{};
        for (size_t end = 0; end < 2; ++end)
            if (lengths[end] > T(0) && !surfaces.empty())
                dv[end] = detail::v_extension_for_length(surfaces, end == 1, lengths[end]);
        for (auto &s : surfaces)
            s = extrapolated_v(s, dv);
        return surfaces;
    }
}
