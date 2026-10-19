#pragma once

/**
 * @file pcurve.h
 * @brief Native BREP core, stage 1: pcurve of an edge on a surface, either
 * extracted exactly (the edge is a `CurveOnSurface` on that very surface) or
 * projected point by point and interpolated.
 *
 * Design: docs/sources/design/brep_core.md, section 6.2; architecture note
 * docs/sources/design/brep_pr05_face_wire.md.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <expected>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include <gbs-brep/model.h>
#include <gbs-brep/closure.h>
#include <gbs/bscinterp.h>
#include <gbs/curveonsurface.h>
#include <gbs/intersect/intersection_curve.h>
#include <gbs/extrema.h>
#include <gbs/primitives.h>
#include <gbs/transform.h>

namespace gbs::brep
{
    /// Parameters of the face-from-wire builders.
    template <std::floating_point T>
    struct MakeFaceOptions
    {
        T tol = brep_default_tolerance<T>;          ///< face tolerance; edges are raised to at least this
        T pcurve_tol = brep_pcurve_approx_tol<T>;   ///< accepted 3D deviation of a projected pcurve (and of the edge from the surface)
        std::size_t pcurve_degree = 3;              ///< degree of the interpolated pcurves
        std::size_t n_samples_min = 9;              ///< first sampling of an edge
        std::size_t n_samples_max = 1025;           ///< sampling refined (2n - 1) until pcurve_tol is met or this is reached
    };

    using gbs::PCurveFit;

    /// Why a pcurve could not be built.
    enum class PCurveErrc : std::uint8_t
    {
        OffSurface,    ///< a point of the edge is farther than pcurve_tol from the surface
        CrossesSeam,   ///< the edge goes across the seam of a closed surface (split it first)
        NotConverged,  ///< pcurve_tol not reached with n_samples_max samples
    };

    namespace detail
    {
        /**
         * Gauss-Newton projection of `p` on `srf` from the seed `uv`, clamped to the
         * parametric rectangle. Returns {u, v, distance}. Converges quadratically for
         * points on the surface; at a singular point (pole) it stops where it is.
         */
        template <std::floating_point T>
        auto project_point(const Surface<T, 3> &srf, const point<T, 3> &p, point<T, 2> uv) -> std::array<T, 3>
        {
            const auto [u1, u2, v1, v2] = srf.bounds();
            const T eps = std::numeric_limits<T>::epsilon() * 16;
            for (int it = 0; it < 30; ++it)
            {
                const auto r = srf.value(uv[0], uv[1]) - p;
                const auto su = srf.value(uv[0], uv[1], 1, 0);
                const auto sv = srf.value(uv[0], uv[1], 0, 1);
                const T a = (su * su), b = (su * sv), c = (sv * sv);
                const T g1 = (su * r), g2 = (sv * r);
                const T det = a * c - b * b;
                if (!(det > eps * a * c) || !std::isfinite(det))
                    break; // singular (pole) or degenerate metric
                const T du = -(c * g1 - b * g2) / det;
                const T dv = -(a * g2 - b * g1) / det;
                const point<T, 2> next{std::clamp(uv[0] + du, u1, u2), std::clamp(uv[1] + dv, v1, v2)};
                const T step = std::abs(next[0] - uv[0]) / (u2 - u1) + std::abs(next[1] - uv[1]) / (v2 - v1);
                uv = next;
                if (step < T(1e-15))
                    break;
            }
            return {uv[0], uv[1], distance(srf.value(uv[0], uv[1]), p)};
        }

        /**
         * Global projection: bracketed derivative-free search, then Gauss-Newton refinement.
         * The search may stop in a local minimum on the border of the parametric rectangle (a
         * point near the seam of a closed surface, on the other side): if the result is not on
         * the surface, Gauss-Newton is restarted from the nearest node of a 32 x 32 grid.
         */
        template <std::floating_point T>
        auto project_point_global(const Surface<T, 3> &srf, const point<T, 3> &p, T tol) -> std::array<T, 3>
        {
            // a primitive inverts the point in closed form
            if (auto e = dynamic_cast<const ElementarySurface<T> *>(&srf))
            {
                const auto [u1, u2, v1, v2] = srf.bounds();
                auto uv = e->parameters(p);
                const auto r = project_point(srf, p, point<T, 2>{std::clamp(uv[0], u1, u2), std::clamp(uv[1], v1, v2)});
                if (r[2] <= tol)
                    return r;
            }
            auto [u, v, d] = extrema_surf_pnt(srf, p, tol);
            auto best = project_point(srf, p, point<T, 2>{u, v});
            if (best[2] > d)
                best = {u, v, d};
            if (best[2] <= tol)
                return best;
            const auto [u1, u2, v1, v2] = srf.bounds();
            constexpr int n = 32;
            point<T, 2> seed{u1, v1};
            T dmin = std::numeric_limits<T>::max();
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= n; ++j)
                {
                    const point<T, 2> uv{u1 + (u2 - u1) * T(i) / T(n), v1 + (v2 - v1) * T(j) / T(n)};
                    const T dij = distance(srf.value(uv[0], uv[1]), p);
                    if (dij < dmin)
                        dmin = dij, seed = uv;
                }
            const auto grid = project_point(srf, p, seed);
            return grid[2] < best[2] ? grid : best;
        }

        /// True if the parametric point is a singular point of the surface for the coordinate k (dS/dk ~ 0).
        template <std::floating_point T>
        bool free_coordinate(const Surface<T, 3> &srf, const point<T, 2> &uv, std::size_t k)
        {
            const auto su = srf.value(uv[0], uv[1], 1, 0);
            const auto sv = srf.value(uv[0], uv[1], 0, 1);
            const T nu = norm(su), nv = norm(sv);
            const T scale = std::max(nu, nv);
            return scale > T(0) && (k == 0 ? nu : nv) < T(1e-7) * scale;
        }

        /**
         * Makes the parametric samples of an edge continuous on a closed surface:
         * a sample lying on the seam (u = u1 or u = u2) takes the value consistent
         * with its resolved neighbours; if the whole edge lies on the seam, the
         * value nearest to `hint` (where the previous co-edge ended) is taken.
         * Samples at a singular point copy the free coordinate of a neighbour.
         * Returns false if two consecutive samples are more than half a period
         * apart (the edge crosses the seam).
         */
        template <std::floating_point T>
        bool make_continuous(const Surface<T, 3> &srf, const SurfaceClosure &cl, std::vector<point<T, 2>> &uv,
                             const std::optional<point<T, 2>> &hint)
        {
            const auto [u1, u2, v1, v2] = srf.bounds();
            const std::array<T, 2> lo{u1, v1}, hi{u2, v2};
            const std::array<bool, 2> closed{cl.closed_u, cl.closed_v};
            const std::size_t n = uv.size();

            for (std::size_t k = 0; k < 2; ++k)
            {
                const T period = hi[k] - lo[k];
                const T seam_eps = T(1e-7) * period;
                std::vector<std::uint8_t> resolved(n, 1);
                for (std::size_t i = 0; i < n; ++i)
                {
                    const bool on_seam = closed[k] && (std::abs(uv[i][k] - lo[k]) < seam_eps || std::abs(uv[i][k] - hi[k]) < seam_eps);
                    if (on_seam || free_coordinate(srf, uv[i], k))
                        resolved[i] = 0;
                }
                auto pick = [&](std::size_t i, T ref) {
                    if (free_coordinate(srf, uv[i], k))
                        uv[i][k] = ref; // singular point: any value is the same 3D point
                    else
                        uv[i][k] = std::abs(lo[k] - ref) <= std::abs(hi[k] - ref) ? lo[k] : hi[k];
                    resolved[i] = 1;
                };
                const auto first = std::ranges::find(resolved, std::uint8_t{1});
                if (first == resolved.end())
                {
                    const T ref = hint ? (*hint)[k] : lo[k];
                    for (std::size_t i = 0; i < n; ++i)
                        pick(i, ref);
                }
                else
                {
                    const auto f = static_cast<std::size_t>(first - resolved.begin());
                    for (std::size_t i = f; i-- > 0;)
                        pick(i, uv[i + 1][k]);
                    for (std::size_t i = f + 1; i < n; ++i)
                        if (!resolved[i])
                            pick(i, uv[i - 1][k]);
                }
                if (closed[k])
                    for (std::size_t i = 1; i < n; ++i)
                        if (std::abs(uv[i][k] - uv[i - 1][k]) > period / 2)
                            return false;
            }
            return true;
        }
    } // namespace detail

    /**
     * @brief Exact pcurve of an edge on `srf`: the 2D curve of the edge's
     * `CurveOnSurface` when it lies on that very surface (same `shared_ptr`),
     * nullptr otherwise. Its parameter is the edge's parameter.
     */
    template <std::floating_point T>
    [[nodiscard]] auto extract_pcurve(const Edge<T> &e, const std::shared_ptr<Surface<T, 3>> &srf) -> std::shared_ptr<Curve<T, 2>>
    {
        auto cos = std::dynamic_pointer_cast<CurveOnSurface<T, 3>>(e.curve);
        if (cos && cos->p_basisSurface() == srf)
            return cos->p_basisCurve();
        return nullptr;
    }

    namespace detail
    {
        /**
         * The iso line of a B-spline surface at the parameter `fixed` of v (`along_u`) or of u:
         * the row or column of poles on a bound (within rounding of it), a knot insertion inside;
         * nullopt beyond the bounds.
         */
        template <typename S, std::floating_point T>
        auto row_or_iso(const S &srf, bool along_u, T fixed) -> std::optional<decltype(srf.isoU(T{}))>
        {
            using C = decltype(srf.isoU(T{}));
            const auto [u1, u2, v1, v2] = srf.bounds();
            const T lo = along_u ? v1 : u1, hi = along_u ? v2 : u2;
            const T eps = T(1e-12) * (T(1) + hi - lo);
            if (fixed < lo - eps || fixed > hi + eps)
                return std::nullopt;
            if (std::abs(fixed - lo) <= eps)
                return along_u ? C{srf.polesU(0), srf.knotsFlatsU(), srf.degreeU()} : C{srf.polesV(0), srf.knotsFlatsV(), srf.degreeV()};
            if (std::abs(fixed - hi) <= eps)
                return along_u ? C{srf.polesU(srf.nPolesV() - 1), srf.knotsFlatsU(), srf.degreeU()}
                               : C{srf.polesV(srf.nPolesU() - 1), srf.knotsFlatsV(), srf.degreeV()};
            try
            {
                return along_u ? srf.isoV(fixed) : srf.isoU(fixed);
            }
            catch (const std::exception &)
            {
                return std::nullopt;
            }
        }

        /**
         * The B-spline of a curve, parametrized like it: the curve itself, or the iso line of a
         * B-spline surface that a `CurveOnSurface` follows (its pcurve a segment at constant u or
         * v, linear in the parameter), as the edges of natural faces are. nullptr otherwise.
         */
        template <std::floating_point T>
        auto as_bspline(const std::shared_ptr<Curve<T, 3>> &c) -> std::shared_ptr<Curve<T, 3>>
        {
            if (dynamic_cast<const BSCurve<T, 3> *>(c.get()) || dynamic_cast<const BSCurveRational<T, 3> *>(c.get()))
                return c;
            const auto *cos = dynamic_cast<const CurveOnSurface<T, 3> *>(c.get());
            if (!cos)
                return nullptr;
            const auto *seg = dynamic_cast<const BSCurve<T, 2> *>(cos->p_basisCurve().get());
            if (!seg || seg->degree() != 1 || seg->poles().size() != 2)
                return nullptr;
            const auto p0 = seg->poles()[0], p1 = seg->poles()[1];
            const auto [t0, t1] = seg->bounds();
            auto iso = [&](const auto &srf) -> std::shared_ptr<Curve<T, 3>> {
                const auto [u1, u2, v1, v2] = srf.bounds();
                const T eu = T(1e-12) * (T(1) + u2 - u1), ev = T(1e-12) * (T(1) + v2 - v1);
                const bool along_u = std::abs(p1[1] - p0[1]) <= ev && std::abs(p1[0] - p0[0]) > eu;
                const bool along_v = std::abs(p1[0] - p0[0]) <= eu && std::abs(p1[1] - p0[1]) > ev;
                if (!along_u && !along_v)
                    return nullptr;
                const int k = along_u ? 0 : 1;
                auto curve = row_or_iso(srf, along_u, along_u ? p0[1] : p0[0]);
                if (!curve)
                    return nullptr;
                // the parameter w of the iso is a t + b along the edge: reversed when a < 0 (its
                // parameter becomes k1 + k2 - w), then the knots mapped to t
                T a = (p1[k] - p0[k]) / (t1 - t0), b = p0[k] - a * t0;
                if (a < T(0))
                {
                    const T k1 = curve->knotsFlats().front(), k2 = curve->knotsFlats().back();
                    curve->reverse();
                    a = -a;
                    b = k1 + k2 - b;
                }
                auto knots = curve->knotsFlats();
                for (auto &x : knots)
                    x = (x - b) / a;
                using C = std::remove_cvref_t<decltype(*curve)>;
                return std::make_shared<C>(curve->poles(), knots, curve->degree());
            };
            if (const auto *bs = dynamic_cast<const BSSurface<T, 3> *>(cos->p_basisSurface().get()))
                return iso(*bs);
            if (const auto *bs = dynamic_cast<const BSSurfaceRational<T, 3> *>(cos->p_basisSurface().get()))
                return iso(*bs);
            return nullptr;
        }
    }

    /**
     * @brief Exact pcurve of a B-spline edge (rational or not, or an iso line of a B-spline
     * surface, `detail::as_bspline`) on a plane: the curve trimmed to the range of the edge, its
     * poles mapped into the parameters of the plane, S(u, v) = O + u X + v Y, homogeneous
     * coordinates kept for a rational curve, same knots, so parametrized like the edge. Its
     * deviation is the largest distance of a pole to the plane (the curve lies in the convex hull
     * of its poles).
     *
     * @return nullopt when `srf` is not a plane, the edge is not a B-spline, or a pole is farther
     * than `tol` from the plane
     */
    template <std::floating_point T>
    [[nodiscard]] auto plane_pcurve(const Edge<T> &e, const Surface<T, 3> &srf, T tol) -> std::optional<PCurveFit<T>>
    {
        const auto *pl = dynamic_cast<const Plane<T> *>(&srf);
        if (!pl || !e.curve)
            return std::nullopt;
        const auto &f = pl->frame();
        const auto crv = detail::as_bspline(e.curve);
        if (!crv)
            return std::nullopt;
        T dev{};
        // (trimmed to the edge: the curve may leave the plane beyond it)
        if (const auto *b0 = dynamic_cast<const BSCurve<T, 3> *>(crv.get()))
        {
            auto bs = *b0;
            bs.trim(e.u1, e.u2);
            points_vector<T, 2> poles;
            for (const auto &q : bs.poles())
            {
                const auto d = q - f.o;
                dev = std::max(dev, std::abs(d * f.z));
                poles.push_back({d * f.x, d * f.y});
            }
            if (dev > tol)
                return std::nullopt;
            return PCurveFit<T>{std::make_shared<BSCurve<T, 2>>(poles, bs.knotsFlats(), bs.degree()), dev, true};
        }
        if (const auto *b0 = dynamic_cast<const BSCurveRational<T, 3> *>(crv.get()))
        {
            auto bs = *b0;
            bs.trim(e.u1, e.u2);
            // homogeneous poles (w P, w): (w (P - O) . X, w (P - O) . Y, w)
            std::vector<std::array<T, 3>> poles;
            for (const auto &q : bs.poles())
            {
                const T w = q[3];
                if (!(w > T(0)))
                    return std::nullopt;
                const point<T, 3> wd{q[0] - w * f.o[0], q[1] - w * f.o[1], q[2] - w * f.o[2]};
                dev = std::max(dev, std::abs(wd * f.z) / w);
                poles.push_back({wd * f.x, wd * f.y, w});
            }
            if (dev > tol)
                return std::nullopt;
            return PCurveFit<T>{std::make_shared<BSCurveRational<T, 2>>(poles, bs.knotsFlats(), bs.degree()), dev, true};
        }
        return std::nullopt;
    }

    /**
     * @brief A pcurve moved by `shift` in (u, v) (a whole period across a seam): a copy of the
     * B-spline, rational or not, translated (`gbs::translate`); the pcurve itself for a zero
     * shift; nullptr for another curve moved.
     */
    template <std::floating_point T>
    [[nodiscard]] auto shifted_pcurve(const std::shared_ptr<Curve<T, 2>> &pc, const point<T, 2> &shift) -> std::shared_ptr<Curve<T, 2>>
    {
        if (shift[0] == T(0) && shift[1] == T(0))
            return pc;
        if (const auto *bs = dynamic_cast<const BSCurve<T, 2> *>(pc.get()))
        {
            auto c = std::make_shared<BSCurve<T, 2>>(*bs);
            translate(*c, shift);
            return c;
        }
        if (const auto *bs = dynamic_cast<const BSCurveRational<T, 2> *>(pc.get()))
        {
            auto c = std::make_shared<BSCurveRational<T, 2>>(*bs);
            translate(*c, shift);
            return c;
        }
        return nullptr;
    }

    /**
     * @brief Pcurve of an edge on `srf` by projection (exact for a B-spline on a plane,
     * `plane_pcurve`): the edge is sampled at
     * n parameters t_i, each C(t_i) is projected on the surface (global search
     * for the first sample, Gauss-Newton seeded by the previous one after),
     * the parametric samples are made continuous across the seam and the
     * poles, and a B-spline of degree `pcurve_degree` interpolates (t_i, uv_i):
     * the pcurve is parametrized like the edge. The 3D deviation is measured
     * at the midpoints; the sampling is refined until it is <= `pcurve_tol`.
     *
     * @param hint parametric point where the previous co-edge of the wire ends,
     *             used only when the whole edge lies on a seam.
     */
    template <std::floating_point T>
    [[nodiscard]] auto project_pcurve(const Edge<T> &e, const Surface<T, 3> &srf, const SurfaceClosure &cl,
                                      const std::optional<point<T, 2>> &hint, const MakeFaceOptions<T> &opts)
        -> std::expected<PCurveFit<T>, PCurveErrc>
    {
        // a B-spline on a plane: exact, whatever its parametrization (a rational one is only C1
        // at its knots, and an interpolated pcurve would need hundreds of samples)
        if (auto exact = plane_pcurve(e, srf, opts.pcurve_tol))
            return *exact;
        const auto &crv = *e.curve;
        std::size_t n = std::max<std::size_t>(opts.n_samples_min, 2);
        for (;;)
        {
            std::vector<T> t(n);
            std::vector<point<T, 2>> uv(n);
            T dev = 0;
            for (std::size_t i = 0; i < n; ++i)
            {
                t[i] = e.u1 + (e.u2 - e.u1) * T(i) / T(n - 1);
                const auto p = crv.value(t[i]);
                auto r = i == 0 ? detail::project_point_global(srf, p, opts.pcurve_tol * T(1e-3))
                                : detail::project_point(srf, p, uv[i - 1]);
                if (i > 0 && r[2] > opts.pcurve_tol) // seed fell on the wrong branch: restart globally
                    r = detail::project_point_global(srf, p, opts.pcurve_tol * T(1e-3));
                if (r[2] > opts.pcurve_tol)
                    return std::unexpected(PCurveErrc::OffSurface);
                uv[i] = {r[0], r[1]};
                dev = std::max(dev, r[2]);
            }
            if (!detail::make_continuous(srf, cl, uv, hint))
                return std::unexpected(PCurveErrc::CrossesSeam);

            std::vector<constrType<T, 2, 1>> q(n);
            for (std::size_t i = 0; i < n; ++i)
                q[i] = {uv[i]};
            auto ip = interpolate<T, 2>(q, t, std::min(opts.pcurve_degree, n - 1));
            // keep the pcurve inside the parametric rectangle (an edge on its border is interpolated with
            // rounding errors beyond it)
            auto pc = std::make_shared<BSCurve<T, 2>>(clamp_poles(ip, srf.bounds()));

            for (std::size_t i = 0; i + 1 < n; ++i)
            {
                const T tm = (t[i] + t[i + 1]) / 2;
                const auto w = pc->value(tm);
                dev = std::max(dev, distance(srf.value(w[0], w[1]), crv.value(tm)));
            }
            if (dev <= opts.pcurve_tol)
                return PCurveFit<T>{std::move(pc), dev, false};
            if (n >= opts.n_samples_max)
                return std::unexpected(PCurveErrc::NotConverged);
            n = std::min(2 * n - 1, opts.n_samples_max);
        }
    }

} // namespace gbs::brep
