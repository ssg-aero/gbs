#pragma once
/**
 * @file pipe.h
 * @brief Sweep of a profile along a path (pipe): an edge gives a face per path edge, a wire a
 * shell, a face a solid.
 *
 * The profile is carried by a moving frame of the path (`gbs/movingframe.h`): rotation-minimizing
 * by default, Frenet or with a fixed binormal on request. It is taken as it is: the frame at the
 * start of the path is the reference, the profile need not lie in the normal plane nor touch the
 * path. Each edge of the profile and each edge of the path give a face:
 * - along a segment, where the frame does not turn: the prism of the placed profile edge (plane,
 *   cylinder, exact NURBS extrusion, as `make_prism`);
 * - along an arc of circle, where the frame turns with the arc (any planar path for a
 *   rotation-minimizing or Frenet frame): the revolution of the placed profile edge (cylinder,
 *   cone, sphere, torus, exact NURBS revolution, as `make_revolution`);
 * - otherwise a B-spline: the profile edge (its NURBS) placed by the frame at stations of the path
 *   edge and skinned (`loft`, cubic across), the stations doubled until the surface is within
 *   `approx_tol` of the exact sweep (the profile placed by the frame at the same parameter),
 *   measured between the stations. Counted in `SweepResult::approximated`.
 *
 * The faces are sewn into a shell, closed by two copies of a face profile (at the start, and
 * placed by the frame at the end) into a solid. A closed path gives no caps: the frame must bring
 * the profile back onto itself (the rotation-minimizing one does, its residual angle spread along
 * the path).
 *
 * Refused (`NonGeneric`): a path that is not G1 (a corner whose angle, times the size of the
 * profile, exceeds `tol`); a profile beyond the radius of curvature of the path, where the swept
 * surface folds; a frame undefined on the path (Frenet where the curvature vanishes, a fixed
 * binormal along the tangent). Design note: docs/sources/design/brep4_pr02_pipe.md.
 */
#include <gbs-brep/sweep.h>
#include <gbs/loft/loftGeneric.h>
#include <gbs/movingframe.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    using PipeFrame = gbs::FrameLaw;

    template <std::floating_point T>
    struct PipeOptions
    {
        T tol{brep_default_tolerance<T>}; ///< tolerance of the faces, of the sewing and of the corners of the path
        T approx_tol{T(1e-5)};            ///< approximated surfaces: deviation from the exact sweep
        PipeFrame frame{PipeFrame::RotationMinimizing};
        point<T, 3> binormal{};           ///< the direction kept by FixedBinormal
    };

    namespace pipe_detail
    {
        template <std::floating_point T>
        struct Piece
        {
            EdgeId edge;
            std::shared_ptr<Curve<T, 3>> curve;
            T u1, u2;
            bool reversed;
        };

        /// The path as one curve: the piece k on [k, k + 1], in the direction of the wire.
        template <std::floating_point T>
        class PathCurve : public Curve<T, 3>
        {
            std::vector<Piece<T>> p_;

        public:
            explicit PathCurve(std::vector<Piece<T>> p) : p_{std::move(p)} {}
            [[nodiscard]] auto pieces() const -> const std::vector<Piece<T>> & { return p_; }
            /// Value of the piece k at s (in [k, k + 1], or slightly beyond).
            auto at(std::size_t k, T s, size_t d = 0) const -> std::array<T, 3>
            {
                const auto &q = p_[k];
                const T w = s - T(k), l = q.u2 - q.u1;
                const T u = q.reversed ? q.u2 - w * l : q.u1 + w * l;
                auto v = q.curve->value(u, d);
                const T f = std::pow(q.reversed ? -l : l, T(d));
                return d == 0 ? v : f * v;
            }
            auto value(T s, size_t d = 0) const -> std::array<T, 3> override
            {
                const auto k = std::size_t(std::clamp(std::floor(s), T(0), T(p_.size() - 1)));
                return at(k, s, d);
            }
            auto bounds() const -> std::array<T, 2> override { return {T(0), T(p_.size())}; }
        };

        /// The rigid motion carrying the frame at s0 to the frame at s.
        template <std::floating_point T>
        auto motion(const MovingFrame<T> &f, T s0, T s) -> Placement<T>
        {
            const auto a = f(s0), b = f(s);
            return Placement<T>::frame(b.o, b.x, b.y, b.z) * Placement<T>::frame(a.o, a.x, a.y, a.z).inverse();
        }

        template <std::floating_point T>
        bool same(const Placement<T> &a, const Placement<T> &b, T eps, T size)
        {
            for (std::size_t i = 0; i < 3; ++i)
            {
                if (std::abs(a.t[i] - b.t[i]) > eps * size)
                    return false;
                for (std::size_t j = 0; j < 3; ++j)
                    if (std::abs(a.r[i][j] - b.r[i][j]) > eps)
                        return false;
            }
            return true;
        }

        /// The profile edge placed by p.
        template <std::floating_point T>
        auto placed(const Edge<T> &e, const Placement<T> &p) -> Edge<T>
        {
            auto r = e;
            r.curve = moved_curve(e.curve, p);
            return r;
        }

        /**
         * B-spline sweep of the profile edge along the piece k: its NURBS placed at stations and
         * skinned, the stations doubled until the surface is within approx_tol of the exact sweep.
         */
        template <std::floating_point T>
        auto skinned(const Edge<T> &e, const MovingFrame<T> &f, std::size_t k, const PipeOptions<T> &o, std::size_t &approximated,
                     T &deviation) -> std::shared_ptr<Surface<T, 3>>
        {
            std::shared_ptr<Surface<T, 3>> r;
            sweep_detail::with_nurbs(e, o.approx_tol, approximated, deviation, [&](const auto &c) {
                using C = std::remove_cvref_t<decltype(c)>;
                // parameters of the profile where the gap is measured: inside each span
                std::vector<T> us;
                const auto kn = c.knots();
                for (std::size_t i = 0; i + 1 < kn.size(); ++i)
                    for (int j = 0; j <= 4; ++j)
                        us.push_back(kn[i] + (kn[i + 1] - kn[i]) * T(j) / T(4));
                const T s0 = T(k);
                for (std::size_t n = 4; n <= 4096; n *= 2)
                {
                    std::vector<T> v(n + 1);
                    std::vector<C> sections;
                    sections.reserve(n + 1);
                    for (std::size_t i = 0; i <= n; ++i)
                    {
                        v[i] = s0 + T(i) / T(n);
                        auto ci = c;
                        const auto p = motion(f, T(0), v[i]);
                        gbs::transform(ci, [&p](const point<T, 3> &x) { return p.apply(x); });
                        sections.push_back(std::move(ci));
                    }
                    auto srf = loft_generic<T, 3>(sections, v, std::min<std::size_t>(3, n));
                    T dev{};
                    for (std::size_t i = 0; i < n && dev <= o.approx_tol; ++i)
                        for (int j = 1; j < 4; ++j)
                        {
                            const T s = v[i] + (v[i + 1] - v[i]) * T(j) / T(4);
                            const auto p = motion(f, T(0), s);
                            for (T u : us)
                                dev = std::max(dev, distance(srf.value(u, s), p.apply(c.value(u))));
                        }
                    if (dev <= o.approx_tol)
                    {
                        ++approximated;
                        deviation = std::max(deviation, dev);
                        r = std::make_shared<decltype(srf)>(std::move(srf));
                        return;
                    }
                }
                throw BRepError("make_pipe: the swept surface cannot be approximated within approx_tol");
            });
            return r;
        }
    }

    /**
     * @brief Sweep of a profile (an edge, a wire or a face) along a path (an edge or a G1 wire),
     * carried by a moving frame: an edge gives a face (a shell for a path of several edges), a wire
     * a shell, a face a solid (closed by the profile and its copy at the end, unless the path is
     * closed). On failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_pipe(Model<T> &m, const ShapeId &profile, const ShapeId &path, std::type_identity_t<PipeOptions<T>> o = {})
        -> BuildResult<SweepResult<T>>
    {
        using namespace pipe_detail;
        if (!m.alive(path) || !m.alive(profile))
            return detail::fail(BuildErrc::InvalidId, "make_pipe: dead or invalid profile or path");
        // the path, in the direction of the wire
        std::vector<Piece<T>> pieces;
        bool closed = false;
        if (const auto *e = std::get_if<EdgeId>(&path))
        {
            const auto &ed = m.edge(*e);
            if (!ed.degenerate && ed.curve)
                pieces.push_back({*e, ed.curve, ed.u1, ed.u2, false});
            closed = ed.v1 == ed.v2;
        }
        else if (const auto *w = std::get_if<WireId>(&path))
        {
            for (const auto &ce : m.wire(*w).coedges)
            {
                const auto &ed = m.edge(ce.edge);
                if (!ed.degenerate && ed.curve)
                    pieces.push_back({ce.edge, ed.curve, ed.u1, ed.u2, ce.orient == Orientation::Reversed});
            }
            closed = m.wire(*w).closed;
        }
        else
            return detail::fail(BuildErrc::InvalidId, "make_pipe: the path must be an edge or a wire");
        if (pieces.empty())
            return detail::fail(BuildErrc::EmptyInput, "make_pipe: empty path");
        const auto curve = std::make_shared<PathCurve<T>>(pieces);
        const std::size_t n = pieces.size();
        const auto x0 = curve->value(T(0));

        // points of the profile: its size, the corners, the closure and the folds are measured on them
        std::vector<point<T, 3>> pts;
        for (auto eid : explore<EdgeId>(m, profile))
        {
            const auto &e = m.edge(eid);
            if (e.degenerate || !e.curve)
                continue;
            for (int i = 0; i <= 16; ++i)
                pts.push_back(e.curve->value(e.u1 + (e.u2 - e.u1) * T(i) / T(16)));
        }
        if (pts.empty())
            return detail::fail(BuildErrc::EmptyInput, "make_pipe: the profile sweeps no surface");
        T rho{};
        for (const auto &p : pts)
            rho = std::max(rho, distance(p, x0));
        rho = std::max(rho, o.tol);

        // G1: the corners at the junctions, the profile carried across them off by angle x size
        auto unit = [](const point<T, 3> &v) { return v / norm(v); };
        for (std::size_t k = 0; k < n; ++k)
        {
            if (k + 1 == n && !closed)
                break;
            const auto ta = unit(curve->at(k, T(k + 1), 1)), tb = unit(curve->at((k + 1) % n, T((k + 1) % n), 1));
            if (std::atan2(norm(cross(ta, tb)), ta * tb) * rho > o.tol)
                return detail::fail(BuildErrc::NonGeneric, "make_pipe: the path is not G1 (a corner at a junction of its edges)");
        }

        std::optional<MovingFrame<T>> frame;
        try
        {
            std::vector<T> breaks;
            for (std::size_t k = 1; k < n; ++k)
                breaks.push_back(T(k));
            frame.emplace(curve, T(0), T(n), o.frame, o.approx_tol / rho, breaks, closed, o.binormal);
        }
        catch (const std::domain_error &e)
        {
            return detail::fail(BuildErrc::NonGeneric, std::string("make_pipe: ") + e.what());
        }
        catch (const std::invalid_argument &e)
        {
            return detail::fail(BuildErrc::InvalidBounds, std::string("make_pipe: ") + e.what());
        }
        const auto &F = *frame;

        // closed path: the profile must come back onto itself
        if (closed)
        {
            const auto back = motion(F, T(0), T(n));
            for (const auto &p : pts)
                if (distance(back.apply(p), p) > o.tol)
                    return detail::fail(BuildErrc::NonGeneric, "make_pipe: closed path along which the profile does not come back onto itself");
        }
        // folds: the speed of a point q of the profile along the path is t (1 - k . (q - c)) across
        for (std::size_t k = 0; k < n; ++k)
            for (int i = 0; i <= 64; ++i)
            {
                const T s = T(k) + T(i) / T(64);
                const auto d1 = curve->at(k, s, 1), d2 = curve->at(k, s, 2);
                const T l2 = d1 * d1;
                const auto kv = (d2 - ((d2 * d1) / l2) * d1) / l2;
                const auto p = motion(F, T(0), s);
                const auto c = curve->at(k, s);
                for (const auto &q : pts)
                    if (T(1) - kv * (p.apply(q) - c) <= T(1e-9))
                        return detail::fail(BuildErrc::NonGeneric, "make_pipe: the profile is beyond the radius of curvature of the path (the swept surface folds)");
            }

        const SweepOptions<T> so{o.tol, o.approx_tol};
        auto surfaces_of = [&](const Edge<T> &e, std::size_t &approximated, T &deviation) {
            std::vector<std::shared_ptr<Surface<T, 3>>> out;
            for (std::size_t k = 0; k < n; ++k)
            {
                const auto &piece = pieces[k];
                const T s0 = T(k);
                const auto p0 = motion(F, T(0), s0);
                const T size = T(1) + rho + norm(curve->at(k, s0));
                // along a segment, the frame not turning: a prism
                if (sweep_detail::segment(m.edge(piece.edge), o.tol))
                {
                    const auto d = curve->at(k, s0 + T(1)) - curve->at(k, s0);
                    bool translation = true;
                    for (int i = 1; i <= 4 && translation; ++i)
                    {
                        const T s = s0 + T(i) / T(4);
                        translation = same(motion(F, T(0), s), Placement<T>{.t = (T(i) / T(4)) * d} * p0, T(1e-9), size);
                    }
                    if (translation)
                    {
                        out.push_back(sweep_detail::prism_surface(placed(e, p0), d, so, approximated, deviation));
                        continue;
                    }
                }
                // along an arc of circle, the frame turning with it: a revolution
                if (const auto *el = dynamic_cast<const EllipseCurve<T> *>(piece.curve.get()); el && el->is_circle())
                {
                    const auto O = el->frame().o;
                    const auto Z = unit(cross(curve->at(k, s0) - O, curve->at(k, s0, 1)));
                    const T angle = piece.u2 - piece.u1;
                    bool rotation = true;
                    for (int i = 1; i <= 4 && rotation; ++i)
                    {
                        const T s = s0 + T(i) / T(4);
                        rotation = same(motion(F, T(0), s), sweep_detail::rotation<T>(O, Z, angle * T(i) / T(4)) * p0, T(1e-9), size);
                    }
                    if (rotation)
                    {
                        out.push_back(sweep_detail::revolution_surface(placed(e, p0), O, Z, angle, so, approximated, deviation));
                        continue;
                    }
                }
                out.push_back(skinned(e, F, k, o, approximated, deviation));
            }
            return out;
        };
        return sweep_detail::sweep(m, profile, so, surfaces_of, closed ? std::optional<Placement<T>>{} : std::optional{motion(F, T(0), T(n))});
    }
}
