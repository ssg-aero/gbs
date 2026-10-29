#pragma once
/**
 * @file sweep.h
 * @brief Linear extrusion (prism) and revolution of edges, wires and faces into faces, shells and
 * solids.
 *
 * Each edge of the profile gives a natural face of its swept surface:
 * - extrusion: a plane for a segment normal to the direction, a cylinder for a circle normal to
 *   it, the exact NURBS extrusion otherwise (`build_extrusion`);
 * - revolution: a cylinder or a cone for a segment in a half plane of the axis, a sphere or a
 *   torus for a circle in a plane of the axis, the exact NURBS revolution otherwise
 *   (`build_revolution`; a segment normal to the axis sweeps a planar annulus, which a natural
 *   face of a plane cannot bound); an edge lying on the axis sweeps nothing.
 * A curve which is neither a primitive nor a B-spline is approximated within `approx_tol` (the
 * interpolation refined until then) and counted. The faces of a face profile are closed by copies
 * of the profile (moved to the end of the sweep), and everything is sewn (`sew`, within the
 * deviation of the approximations) into a shell, closed into a solid for a face. On failure the
 * model is left as it was.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/placement.h>
#include <gbs-brep/sew.h>
#include <gbs/bselementary.h>
#include <gbs/bsctools.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <functional>
#include <span>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    template <std::floating_point T>
    struct SweepOptions
    {
        T tol{brep_default_tolerance<T>};  ///< tolerance of the faces and of the sewing
        T approx_tol{T(1e-5)};             ///< approximation of the curves that are neither primitives nor B-splines
    };

    template <std::floating_point T>
    struct SweepResult
    {
        ShapeId shape;                 ///< face (edge profile), shell (wire profile) or solid (face profile)
        std::size_t approximated{};    ///< curves approximated by B-splines
        T deviation{};                 ///< largest deviation of these approximations (within approx_tol)
    };

    namespace sweep_detail
    {
        template <std::floating_point T>
        auto unit(const point<T, 3> &v) -> point<T, 3> { return v / norm(v); }

        /// Is the edge a straight segment (within tol)? Its ends.
        template <std::floating_point T>
        auto segment(const Edge<T> &e, T tol) -> std::optional<std::array<point<T, 3>, 2>>
        {
            const auto a = e.curve->value(e.u1), b = e.curve->value(e.u2);
            const auto ab = b - a;
            const T l = norm(ab);
            if (!(l > tol))
                return std::nullopt;
            for (int k = 1; k < 16; ++k)
            {
                const auto p = e.curve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(16));
                if (norm(cross(p - a, ab)) / l > tol)
                    return std::nullopt;
            }
            return std::array<point<T, 3>, 2>{a, b};
        }

        /**
         * NURBS of the edge curve on its range: exact for B-splines and primitives; otherwise a cubic
         * B-spline through points of the curve, on the parameter of the edge, their number doubled
         * until it is within approx_tol of the curve (`interpolation_gap`: the middles, the quarters
         * and the top of the parabola through them). Throws BRepError when 4096 spans do not reach
         * it; the deviation reached raises `deviation`.
         */
        template <std::floating_point T, typename Visit>
        void with_nurbs(const Edge<T> &e, T approx_tol, std::size_t &approximated, T &deviation, Visit &&visit)
        {
            if (const auto *bs = dynamic_cast<const BSCurve<T, 3> *>(e.curve.get()))
            {
                auto c = *bs;
                c.trim(e.u1, e.u2);
                return visit(c);
            }
            if (const auto *bsr = dynamic_cast<const BSCurveRational<T, 3> *>(e.curve.get()))
            {
                auto c = *bsr;
                c.trim(e.u1, e.u2);
                return visit(c);
            }
            if (const auto *el = dynamic_cast<const ElementaryCurve<T> *>(e.curve.get()))
            {
                const auto n = el->to_nurbs();
                const T a = el->nurbs_parameter(e.u1), b = el->nurbs_parameter(e.u2);
                if (const auto *r = dynamic_cast<const BSCurveRational<T, 3> *>(n.get()))
                {
                    auto c = *r;
                    c.trim(std::min(a, b), std::max(a, b));
                    return visit(c);
                }
                auto c = dynamic_cast<const BSCurve<T, 3> &>(*n);
                c.trim(std::min(a, b), std::max(a, b));
                return visit(c);
            }
            // any other curve: a cubic B-spline through its points, on the parameter of the edge
            ++approximated;
            for (std::size_t n = 64; n <= 4096; n *= 2)
            {
                std::vector<T> t(n + 1);
                std::vector<constrType<T, 3, 1>> q(n + 1);
                for (std::size_t k = 0; k <= n; ++k)
                {
                    t[k] = e.u1 + (e.u2 - e.u1) * T(k) / T(n);
                    q[k] = {e.curve->value(t[k])};
                }
                auto c = interpolate<T, 3>(q, t, std::size_t{3});
                const T dev = detail::interpolation_gap(std::span<const T>{t}, [&](T s) { return distance(c.value(s), e.curve->value(s)); }, T(0), approx_tol);
                if (dev <= approx_tol)
                {
                    deviation = std::max(deviation, dev);
                    return visit(c);
                }
            }
            throw BRepError("sweep: a curve that is neither a primitive nor a B-spline cannot be approximated within approx_tol");
        }

        /// Rotation of angle a around the axis (o, z), as a placement.
        template <std::floating_point T>
        auto rotation(const point<T, 3> &o, const point<T, 3> &z, T a) -> Placement<T>
        {
            const T c = std::cos(a), s = std::sin(a), t = T(1) - c;
            const T x = z[0], y = z[1], w = z[2];
            Placement<T> p;
            p.r = {{{t * x * x + c, t * x * y - s * w, t * x * w + s * y},
                    {t * x * y + s * w, t * y * y + c, t * y * w - s * x},
                    {t * x * w - s * y, t * y * w + s * x, t * w * w + c}}};
            point<T, 3> ro;
            for (int i = 0; i < 3; ++i)
                ro[i] = p.r[i][0] * o[0] + p.r[i][1] * o[1] + p.r[i][2] * o[2];
            p.t = o - ro;
            return p;
        }

        template <std::floating_point T>
        auto frame(const point<T, 3> &o, const point<T, 3> &z, const point<T, 3> &x) -> Frame3<T>
        {
            Frame3<T> f;
            f.o = o;
            f.z = unit<T>(z);
            f.x = unit<T>(x - (x * f.z) * f.z);
            f.y = cross(f.z, f.x);
            return f;
        }

        /// Surface swept by the edge along d.
        template <std::floating_point T>
        auto prism_surface(const Edge<T> &e, const point<T, 3> &d, const SweepOptions<T> &o, std::size_t &approximated, T &deviation)
            -> std::shared_ptr<Surface<T, 3>>
        {
            const T L = norm(d);
            const auto dn = d / L;
            if (const auto seg = segment(e, o.tol))
            {
                const auto x = (*seg)[1] - (*seg)[0];
                if (std::abs(unit<T>(x) * dn) <= T(1e-12))
                    return std::make_shared<Plane<T>>(frame<T>((*seg)[0], cross(x, d), x), std::array<T, 4>{T(0), norm(x), T(0), L});
            }
            if (const auto *el = dynamic_cast<const EllipseCurve<T> *>(e.curve.get()); el && el->is_circle() && norm(cross(el->frame().z, dn)) <= T(1e-12))
            {
                const T s = el->frame().z * dn > T(0) ? T(1) : T(-1);
                return std::make_shared<CylindricalSurface<T>>(el->radius1(), el->frame(),
                                                               std::array<T, 4>{e.u1, e.u2, s > 0 ? T(0) : -L, s > 0 ? L : T(0)});
            }
            std::shared_ptr<Surface<T, 3>> r;
            with_nurbs(e, o.approx_tol, approximated, deviation, [&](const auto &c) {
                using S = decltype(build_extrusion(c, d, T(0), T(1)));
                r = std::make_shared<S>(build_extrusion(c, d, T(0), T(1)));
            });
            return r;
        }

        /// Surface swept by the edge around the axis (O, Z) from 0 to angle; nullptr for an edge on the axis.
        template <std::floating_point T>
        auto revolution_surface(const Edge<T> &e, const point<T, 3> &O, const point<T, 3> &Z, T angle, const SweepOptions<T> &o,
                                std::size_t &approximated, T &deviation) -> std::shared_ptr<Surface<T, 3>>
        {
            auto h = [&](const point<T, 3> &p) { return (p - O) * Z; };
            auto radial = [&](const point<T, 3> &p) { return p - O - h(p) * Z; };
            // on the axis: nothing
            bool on_axis = true;
            for (int k = 0; k <= 8; ++k)
                on_axis = on_axis && norm(radial(e.curve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(8)))) <= o.tol;
            if (on_axis)
                return nullptr;
            if (const auto seg = segment(e, o.tol))
            {
                const auto [a, b] = *seg;
                const auto ra = radial(a), rb = radial(b);
                const auto ref = norm(ra) >= norm(rb) ? ra : rb;
                // in a half plane of the axis: the radial directions are the same
                if (norm(cross(ra, rb)) <= T(1e-9) * (T(1) + norm(ra) * norm(rb)) && ra * rb >= T(0))
                {
                    const T h1 = h(a), h2 = h(b), r1 = norm(ra), r2 = norm(rb);
                    if (std::abs(r1 - r2) <= o.tol)
                        return std::make_shared<CylindricalSurface<T>>(r1, frame<T>(O, Z, ref), std::array<T, 4>{T(0), angle, std::min(h1, h2), std::max(h1, h2)});
                    if (std::abs(h1 - h2) > o.tol)
                    {
                        // cone: radius r1 at the height h1, semi-angle from the slope
                        const T k = (r2 - r1) / (h2 - h1);
                        return std::make_shared<ConicalSurface<T>>(r1, std::atan(k), frame<T>(O + h1 * Z, Z, ref),
                                                                   std::array<T, 4>{T(0), angle, std::min(T(0), h2 - h1), std::max(T(0), h2 - h1)});
                    }
                }
            }
            if (const auto *el = dynamic_cast<const EllipseCurve<T> *>(e.curve.get()); el && el->is_circle())
            {
                const auto &f = el->frame();
                const auto c = f.o;
                // in a plane of the axis, on one side of it
                if (std::abs(f.z * Z) <= T(1e-12) && std::abs((c - O) * f.z) <= o.tol)
                {
                    const T r = el->radius1(), R = norm(radial(c));
                    const auto pm = e.curve->value((e.u1 + e.u2) / T(2));
                    const auto ref = R > o.tol ? radial(c) : radial(pm);
                    // latitude (sphere) or tube angle (torus) along the arc, unwrapped
                    auto angle_of = [&](const point<T, 3> &p) {
                        const T ph = h(p) - h(c), pr = radial(p) * unit<T>(ref);
                        return R > o.tol ? std::atan2(ph, pr - R) : std::atan2(ph, pr);
                    };
                    std::vector<T> vs;
                    for (int k = 0; k <= 32; ++k)
                    {
                        T v = angle_of(e.curve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(32)));
                        if (!vs.empty())
                            v = vs.back() + std::remainder(v - vs.back(), T(2) * std::numbers::pi_v<T>);
                        vs.push_back(v);
                    }
                    const T v1 = std::min(vs.front(), vs.back()), v2 = std::max(vs.front(), vs.back());
                    bool same_side = true;
                    for (int k = 0; k <= 8; ++k)
                        same_side = same_side && radial(e.curve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(8))) * ref >= -o.tol;
                    if (same_side && R <= o.tol && v1 >= -std::numbers::pi_v<T> / 2 - T(1e-9) && v2 <= std::numbers::pi_v<T> / 2 + T(1e-9))
                        return std::make_shared<SphericalSurface<T>>(r, frame<T>(O + h(c) * Z, Z, ref),
                                                                     std::array<T, 4>{T(0), angle, std::max(v1, -std::numbers::pi_v<T> / 2), std::min(v2, std::numbers::pi_v<T> / 2)});
                    if (R > o.tol)
                        return std::make_shared<ToroidalSurface<T>>(R, r, frame<T>(O + h(c) * Z, Z, ref), std::array<T, 4>{T(0), angle, v1, v2});
                }
            }
            std::shared_ptr<Surface<T, 3>> r;
            with_nurbs(e, o.approx_tol, approximated, deviation, [&](const auto &c) {
                r = std::make_shared<BSSurfaceRational<T, 3>>(build_revolution(c, ax1<T, 3>{{O, Z}}, T(0), angle));
            });
            return r;
        }

        /// Faces of the sweep of the profile edges, then the caps of a face profile, sewn (on failure, `sweep` rolls back).
        template <std::floating_point T, typename SurfaceOf>
        auto build(Model<T> &m, const ShapeId &profile, const SweepOptions<T> &o, SurfaceOf surface_of,
                   const std::optional<Placement<T>> &end_cap) -> BuildResult<SweepResult<T>>
        {
            const bool is_face = std::holds_alternative<FaceId>(profile);
            SweepResult<T> res;
            std::vector<FaceId> faces;
            for (auto eid : explore<EdgeId>(m, profile))
            {
                const auto &e = m.edge(eid);
                if (e.degenerate || !e.curve)
                    continue;
                auto srf = surface_of(e, res.approximated, res.deviation);
                if (!srf)
                    continue;
                auto f = make_face(m, srf, o.tol);
                if (!f)
                    return std::unexpected(f.error());
                faces.push_back(*f);
            }
            if (faces.empty())
                return detail::fail(BuildErrc::EmptyInput, "sweep: the profile sweeps no surface");
            if (std::holds_alternative<EdgeId>(profile))
            {
                res.shape = faces.front();
                return res;
            }
            if (is_face && end_cap) // a complete revolution has no caps
            {
                faces.push_back(std::get<FaceId>(copy_shape(m, profile, Placement<T>{})));
                faces.push_back(std::get<FaceId>(copy_shape(m, profile, *end_cap)));
            }
            // the sides of approximated curves meet the caps within their deviation
            auto sewn = sew(m, faces, SewOptions<T>{std::max(o.tol, res.deviation)});
            if (!sewn)
                return std::unexpected(sewn.error());
            if (!is_face)
            {
                if (sewn->shells.size() != 1)
                    return detail::fail(BuildErrc::Disconnected, "sweep: the swept faces do not form one shell");
                res.shape = sewn->shells.front();
                return res;
            }
            // a face with holes turned completely: the holes sweep cavities, inside the largest shell
            std::vector<std::pair<T, ShellId>> by_volume;
            for (auto sh : sewn->shells)
                by_volume.emplace_back(std::abs(signed_volume(m, sh, 64)), sh);
            std::ranges::sort(by_volume, std::greater<>{});
            std::vector<ShellId> voids;
            for (std::size_t k = 1; k < by_volume.size(); ++k)
                voids.push_back(by_volume[k].second);
            auto solid = make_solid(m, by_volume.front().second, voids);
            if (!solid)
                return std::unexpected(solid.error());
            res.shape = *solid;
            return res;
        }

        /// `build`, the model left as it was on failure (also when a step throws BRepError).
        template <std::floating_point T, typename SurfaceOf>
        auto sweep(Model<T> &m, const ShapeId &profile, const SweepOptions<T> &o, SurfaceOf surface_of,
                   std::type_identity_t<std::optional<Placement<T>>> end_cap) -> BuildResult<SweepResult<T>>
        {
            if (!m.alive(profile))
                return detail::fail(BuildErrc::InvalidId, "sweep: dead or invalid profile");
            if (!std::holds_alternative<FaceId>(profile) && !std::holds_alternative<WireId>(profile) && !std::holds_alternative<EdgeId>(profile))
                return detail::fail(BuildErrc::InvalidId, "sweep: the profile must be an edge, a wire or a face");
            const auto size = detail::arena_size(m);
            try
            {
                auto r = build(m, profile, o, surface_of, end_cap);
                if (!r)
                    detail::rollback(m, size);
                return r;
            }
            catch (const BRepError &e)
            {
                detail::rollback(m, size);
                return detail::fail(BuildErrc::PCurveApproximation, e.what());
            }
        }
    }

    /**
     * @brief Linear extrusion of a profile along the vector d: an edge gives a face, a wire a
     * shell, a face a solid (closed by the profile and its copy moved by d).
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_prism(Model<T> &m, const ShapeId &profile, const std::type_identity_t<point<T, 3>> &d,
                                  std::type_identity_t<SweepOptions<T>> o = {}) -> BuildResult<SweepResult<T>>
    {
        if (!(norm(d) > o.tol))
            return detail::fail(BuildErrc::InvalidBounds, "make_prism: null direction");
        return sweep_detail::sweep(m, profile, o, [&](const Edge<T> &e, std::size_t &n, T &dev) { return sweep_detail::prism_surface(e, d, o, n, dev); },
                                   Placement<T>{.t = d});
    }

    /**
     * @brief Revolution of a profile around the axis {origin, direction} by `angle` (radians, in
     * (0, 2 pi]): an edge gives a face, a wire a shell, a face a solid (closed by the profile and
     * its turned copy unless the turn is complete). Edges on the axis sweep nothing.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_revolution(Model<T> &m, const ShapeId &profile, const ax1<T, 3> &axis, std::type_identity_t<T> angle = T(2) * std::numbers::pi_v<T>,
                                       std::type_identity_t<SweepOptions<T>> o = {}) -> BuildResult<SweepResult<T>>
    {
        if (!(angle > T(0)) || angle > T(2) * std::numbers::pi_v<T> * (T(1) + T(1e-12)))
            return detail::fail(BuildErrc::InvalidBounds, "make_revolution: angle out of (0, 2 pi]");
        if (!(norm(axis[1]) > T(0)))
            return detail::fail(BuildErrc::InvalidBounds, "make_revolution: null axis");
        const auto Z = sweep_detail::unit<T>(axis[1]);
        const bool full = angle >= T(2) * std::numbers::pi_v<T> * (T(1) - T(1e-12));
        std::optional<Placement<T>> cap;
        if (!full)
            cap = sweep_detail::rotation<T>(axis[0], Z, angle);
        return sweep_detail::sweep(
            m, profile, o, [&](const Edge<T> &e, std::size_t &n, T &dev) { return sweep_detail::revolution_surface(e, axis[0], Z, angle, o, n, dev); },
            full ? std::optional<Placement<T>>{} : cap);
    }
}
