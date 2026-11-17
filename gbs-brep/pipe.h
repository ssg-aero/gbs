#pragma once
/**
 * @file pipe.h
 * @brief Sweep of a profile along a path (pipe): an edge gives a face per path edge, a wire a
 * shell, a face a solid.
 *
 * The profile is carried by a moving frame of the path (`gbs/movingframe.h`): rotation-minimizing
 * by default, Frenet or with a fixed binormal on request. It is taken as it is: the frame at the
 * start of the path is the reference, the profile need not lie in the normal plane nor touch the
 * path. Each edge of the profile and each edge of the path give a face, swept by the frame of that
 * path edge (at its end, the limit from the left):
 * - along a segment, where the frame does not turn: the prism of the placed profile edge (plane,
 *   cylinder, exact NURBS extrusion, as `make_prism`);
 * - along an arc of circle (recognized on its points, whatever the type of its curve), where the
 *   frame turns with the arc (any planar path for a rotation-minimizing or Frenet frame): the
 *   revolution of the placed profile edge (cylinder, cone, sphere, torus, exact NURBS revolution,
 *   as `make_revolution`);
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
 * profile, exceeds the tolerance of the junction, max of `tol` and of its vertex); a frame that
 * moves the profile at a junction beyond that tolerance (Frenet across a change of the principal
 * normal); a profile beyond the radius of curvature of the path, where the swept surface folds; a
 * sweep that intersects itself away from the folds (the path coming back within the size of the
 * profile); a frame undefined on the path (Frenet where the curvature vanishes or the binormal
 * reverses, a fixed binormal along the tangent). An accuracy that cannot be reached is a
 * `PCurveApproximation`. Design note: docs/sources/design/brep4_pr02_pipe.md.
 */
#include <gbs-brep/sweep.h>
#include <gbs/loft/loftGeneric.h>
#include <gbs/movingframe.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
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
            auto has_second_derivatives() const -> bool override
            {
                return std::ranges::all_of(p_, [](const Piece<T> &q) { return q.curve->has_second_derivatives(); });
            }
        };

        /// The end of the piece k, as the limit from the left: its frame is that of the piece k.
        template <std::floating_point T>
        T end_of(std::size_t k) { return std::nextafter(T(k + 1), T(k)); }

        /// The rigid motions carrying the frame at the start of the path to the frame at s.
        template <std::floating_point T>
        class Motion
        {
            const MovingFrame<T> &f_;
            Placement<T> inv0_;

        public:
            explicit Motion(const MovingFrame<T> &f) : f_{f}
            {
                const auto a = f(f.bounds()[0]);
                inv0_ = Placement<T>::frame(a.o, a.x, a.y, a.z).inverse();
            }
            auto operator()(T s) const -> Placement<T>
            {
                const auto b = f_(s);
                return Placement<T>::frame(b.o, b.x, b.y, b.z) * inv0_;
            }
        };

        /// Placements of a piece for its skinning at n intervals: at the stations (the last one at the
        /// end of the piece), and at the quarters between them, where the deviation is measured.
        template <std::floating_point T>
        struct Stations
        {
            std::vector<Placement<T>> at, between;
        };

        /// The stations of each piece and number of intervals, computed once for all the edges of the profile.
        template <std::floating_point T>
        class StationCache
        {
            const Motion<T> &mo_;
            std::vector<std::map<std::size_t, Stations<T>>> c_;

        public:
            StationCache(const Motion<T> &mo, std::size_t pieces) : mo_{mo}, c_(pieces) {}
            auto get(std::size_t k, std::size_t n) -> const Stations<T> &
            {
                auto &levels = c_[k];
                if (auto it = levels.find(n); it != levels.end())
                    return it->second;
                Stations<T> st;
                st.at.reserve(n + 1);
                st.between.reserve(3 * n);
                const T s0 = T(k);
                for (std::size_t i = 0; i <= n; ++i)
                    st.at.push_back(mo_(i == n ? end_of<T>(k) : s0 + T(i) / T(n)));
                for (std::size_t i = 0; i < n; ++i)
                {
                    const T a = s0 + T(i) / T(n), b = s0 + T(i + 1) / T(n);
                    for (int j = 1; j < 4; ++j)
                        st.between.push_back(mo_(a + (b - a) * T(j) / T(4)));
                }
                return levels.emplace(n, std::move(st)).first->second;
            }
        };

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

        /// An arc of circle: its center, its axis oriented with the travel, the angle travelled at the
        /// quarters of the piece (from its start; the last one is the angle of the arc).
        template <std::floating_point T>
        struct Arc
        {
            point<T, 3> o, z;
            std::array<T, 5> angle;
        };

        /**
         * The piece k as an arc of circle, whatever the type of its curve: the circle through three
         * of 49 points of the piece, all of them within tol of it, travelled one way.
         */
        template <std::floating_point T>
        auto arc_of(const PathCurve<T> &c, std::size_t k, T tol) -> std::optional<Arc<T>>
        {
            constexpr int N = 48;
            std::array<point<T, 3>, N + 1> p;
            for (int i = 0; i <= N; ++i)
                p[i] = c.at(k, T(k) + T(i) / T(N));
            // the circumcenter of p[0], p[N/3], p[2N/3] (three distinct points, also on a full circle)
            const auto ab = p[N / 3] - p[0], ac = p[2 * N / 3] - p[0], nn = cross(ab, ac);
            const T n2 = nn * nn;
            if (!(n2 > T(0)))
                return std::nullopt;
            const auto o = p[0] + ((ac * ac) * cross(nn, ab) + (ab * ab) * cross(ac, nn)) / (T(2) * n2);
            auto z = nn / std::sqrt(n2);
            const T r = norm(p[0] - o);
            for (const auto &q : p)
                if (std::abs(norm(q - o) - r) > tol || std::abs((q - o) * z) > tol)
                    return std::nullopt;
            if (cross(p[0] - o, p[1] - p[0]) * z < T(0))
                z = T(-1) * z;
            Arc<T> a{o, z, {}};
            T acc{};
            for (int i = 1; i <= N; ++i)
            {
                const auto r0 = p[i - 1] - o, r1 = p[i] - o;
                const T step = std::atan2(cross(r0, r1) * z, r0 * r1);
                if (!(step > T(0)))
                    return std::nullopt;
                acc += step;
                if (i % (N / 4) == 0)
                    a.angle[std::size_t(i / (N / 4))] = acc;
            }
            return a;
        }

        template <std::floating_point T>
        bool segment_crosses(const point<T, 3> &p, const point<T, 3> &q, const point<T, 3> &a, const point<T, 3> &b, const point<T, 3> &c)
        {
            const auto e1 = b - a, e2 = c - a, d = q - p;
            const auto h = cross(d, e2);
            const T det = e1 * h;
            if (!(std::abs(det) > T(1e-14) * norm(e1) * norm(e2) * norm(d)))
                return false; // parallel: left to the other edges
            const T f = T(1) / det;
            const auto s = p - a;
            const T u = f * (s * h);
            if (u < T(0) || u > T(1))
                return false;
            const auto sq = cross(s, e1);
            const T v = f * (d * sq);
            if (v < T(0) || u + v > T(1))
                return false;
            const T t = f * (e2 * sq);
            return t >= T(0) && t <= T(1);
        }

        template <std::floating_point T>
        bool triangles_cross(const std::array<point<T, 3>, 3> &a, const std::array<point<T, 3>, 3> &b)
        {
            for (int i = 0; i < 3; ++i)
                if (segment_crosses(a[i], a[(i + 1) % 3], b[0], b[1], b[2]) || segment_crosses(b[i], b[(i + 1) % 3], a[0], a[1], a[2]))
                    return true;
            return false;
        }

        /**
         * Does the sweep intersect itself away from the folds? The polylines of the profile placed by
         * P at stations along the path make quads; two quads whose stations are more than `gap` apart
         * along the path (lengths `ell`, around the path when closed) must not cross. A grid of cells
         * as large as the largest quad gathers the candidates.
         */
        template <std::floating_point T>
        bool self_intersects(const std::vector<std::vector<point<T, 3>>> &polylines, const std::vector<Placement<T>> &P,
                             const std::vector<T> &ell, T gap, bool closed)
        {
            struct Quad
            {
                std::array<point<T, 3>, 4> p;
                point<T, 3> lo, hi;
                std::size_t band;
            };
            auto place = [&](std::size_t j) {
                std::vector<std::vector<point<T, 3>>> r;
                for (const auto &pl : polylines)
                {
                    auto &q = r.emplace_back();
                    for (const auto &x : pl)
                        q.push_back(P[j].apply(x));
                }
                return r;
            };
            std::vector<Quad> quads;
            T cell{};
            auto a = place(0);
            for (std::size_t j = 0; j + 1 < P.size(); ++j)
            {
                auto b = place(j + 1);
                for (std::size_t l = 0; l < a.size(); ++l)
                    for (std::size_t i = 0; i + 1 < a[l].size(); ++i)
                    {
                        Quad q{{a[l][i], a[l][i + 1], b[l][i + 1], b[l][i]}, a[l][i], a[l][i], j};
                        for (const auto &x : q.p)
                            for (std::size_t d = 0; d < 3; ++d)
                                q.lo[d] = std::min(q.lo[d], x[d]), q.hi[d] = std::max(q.hi[d], x[d]);
                        for (std::size_t d = 0; d < 3; ++d)
                            cell = std::max(cell, q.hi[d] - q.lo[d]);
                        quads.push_back(q);
                    }
                a = std::move(b);
            }
            if (!(cell > T(0)))
                return false;
            const T L = ell.back();
            auto near = [&](std::size_t i, std::size_t j) {
                T d = std::abs(ell[i] - ell[j]);
                if (closed)
                    d = std::min(d, L - d);
                return d <= gap;
            };
            using Key = std::array<std::int64_t, 3>;
            auto key = [cell](const point<T, 3> &x) {
                return Key{std::int64_t(std::floor(x[0] / cell)), std::int64_t(std::floor(x[1] / cell)), std::int64_t(std::floor(x[2] / cell))};
            };
            std::map<Key, std::vector<std::size_t>> grid;
            for (std::size_t i = 0; i < quads.size(); ++i)
            {
                const auto lo = key(quads[i].lo), hi = key(quads[i].hi);
                for (auto x = lo[0]; x <= hi[0]; ++x)
                    for (auto y = lo[1]; y <= hi[1]; ++y)
                        for (auto z = lo[2]; z <= hi[2]; ++z)
                            grid[Key{x, y, z}].push_back(i);
            }
            for (const auto &[c, ids] : grid)
                for (std::size_t i = 0; i < ids.size(); ++i)
                    for (std::size_t j = i + 1; j < ids.size(); ++j)
                    {
                        const auto &A = quads[ids[i]], &B = quads[ids[j]];
                        if (near(A.band, B.band))
                            continue;
                        point<T, 3> lo, hi;
                        bool overlap = true;
                        for (std::size_t d = 0; d < 3; ++d)
                        {
                            lo[d] = std::max(A.lo[d], B.lo[d]);
                            hi[d] = std::min(A.hi[d], B.hi[d]);
                            overlap = overlap && lo[d] <= hi[d];
                        }
                        // each pair once: in the cell of the lowest corner of the common box
                        if (!overlap || key(lo) != c)
                            continue;
                        for (const auto &ta : {std::array{A.p[0], A.p[1], A.p[2]}, std::array{A.p[0], A.p[2], A.p[3]}})
                            for (const auto &tb : {std::array{B.p[0], B.p[1], B.p[2]}, std::array{B.p[0], B.p[2], B.p[3]}})
                                if (triangles_cross(ta, tb))
                                    return true;
                    }
            return false;
        }

        /**
         * B-spline sweep of the profile edge along the piece k: its NURBS placed at stations and
         * skinned, the stations doubled until the surface is within approx_tol of the exact sweep.
         */
        template <std::floating_point T>
        auto skinned(const Edge<T> &e, StationCache<T> &cache, std::size_t k, const PipeOptions<T> &o, std::size_t &approximated,
                     T &deviation) -> std::shared_ptr<Surface<T, 3>>
        {
            std::shared_ptr<Surface<T, 3>> r;
            sweep_detail::with_nurbs(e, o.approx_tol, approximated, deviation, [&](const auto &c) {
                using C = std::remove_cvref_t<decltype(c)>;
                // parameters of the profile where the gap is measured: the knots and the quarters of each span
                std::vector<T> us;
                const auto kn = c.knots();
                for (std::size_t i = 0; i + 1 < kn.size(); ++i)
                    for (int j = 0; j < 4; ++j)
                        us.push_back(kn[i] + (kn[i + 1] - kn[i]) * T(j) / T(4));
                us.push_back(kn.back());
                const T s0 = T(k);
                for (std::size_t n = 4; n <= 4096; n *= 2)
                {
                    const auto &st = cache.get(k, n);
                    std::vector<T> v(n + 1);
                    std::vector<C> sections;
                    sections.reserve(n + 1);
                    for (std::size_t i = 0; i <= n; ++i)
                    {
                        v[i] = s0 + T(i) / T(n);
                        auto ci = c;
                        const auto &p = st.at[i];
                        gbs::transform(ci, [&p](const point<T, 3> &x) { return p.apply(x); });
                        sections.push_back(std::move(ci));
                    }
                    auto srf = loft_generic<T, 3>(sections, v, std::min<std::size_t>(3, n));
                    T dev{};
                    for (std::size_t i = 0; i < n && dev <= o.approx_tol; ++i)
                        for (int j = 1; j < 4; ++j)
                        {
                            const T s = v[i] + (v[i + 1] - v[i]) * T(j) / T(4);
                            const auto &p = st.between[3 * i + std::size_t(j - 1)];
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
        using sweep_detail::unit;
        if (!m.alive(path) || !m.alive(profile))
            return detail::fail(BuildErrc::InvalidId, "make_pipe: dead or invalid profile or path");
        if (!std::holds_alternative<FaceId>(profile) && !std::holds_alternative<WireId>(profile) && !std::holds_alternative<EdgeId>(profile))
            return detail::fail(BuildErrc::InvalidId, "make_pipe: the profile must be an edge, a wire or a face");
        if (!detail::valid_tolerance(o.tol) || !detail::valid_tolerance(o.approx_tol))
            return detail::fail(BuildErrc::InvalidBounds, "make_pipe: tol and approx_tol must be positive");
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

        // points of the profile: its size, the corners, the closure and the folds are measured on them;
        // its polylines (one segment for a straight edge) for the self-intersection
        std::vector<point<T, 3>> pts;
        std::vector<std::vector<point<T, 3>>> polylines;
        for (auto eid : explore<EdgeId>(m, profile))
        {
            const auto &e = m.edge(eid);
            if (e.degenerate || !e.curve)
                continue;
            auto &pl = polylines.emplace_back();
            for (int i = 0; i <= 16; ++i)
                pl.push_back(pts.emplace_back(e.curve->value(e.u1 + (e.u2 - e.u1) * T(i) / T(16))));
            if (sweep_detail::segment(e, o.tol))
                pl = {pl.front(), pl.back()};
        }
        if (pts.empty())
            return detail::fail(BuildErrc::EmptyInput, "make_pipe: the profile sweeps no surface");
        T rho{};
        for (const auto &p : pts)
            rho = std::max(rho, distance(p, x0));
        rho = std::max(rho, o.tol);

        // tolerance of the junction at the start of the piece k (k = 0: the closure): of its vertex
        auto junction_tol = [&](std::size_t k) {
            const auto &prev = pieces[(k + n - 1) % n];
            const auto &ed = m.edge(prev.edge);
            return std::max(o.tol, m.vertex(prev.reversed ? ed.v1 : ed.v2).tol);
        };
        // G1: the corners at the junctions, the profile carried across them off by angle x size
        for (std::size_t k = 0; k < n; ++k)
        {
            if (k + 1 == n && !closed)
                break;
            const auto ta = unit<T>(curve->at(k, T(k + 1), 1)), tb = unit<T>(curve->at((k + 1) % n, T((k + 1) % n), 1));
            if (std::atan2(norm(cross(ta, tb)), ta * tb) * rho > junction_tol((k + 1) % n))
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
        catch (const std::runtime_error &e)
        {
            return detail::fail(BuildErrc::PCurveApproximation, std::string("make_pipe: ") + e.what());
        }
        const Motion<T> mo{*frame};

        // the junctions (and the closure): the profile at the end of a piece and at the start of the
        // next one, within the tolerance of the junction; the sewing covers these gaps
        T sew_tol{};
        for (std::size_t k = closed ? 0 : 1; k < n; ++k)
        {
            const auto left = mo(end_of<T>((k + n - 1) % n)), right = mo(T(k));
            T gap{};
            for (const auto &p : pts)
                gap = std::max(gap, distance(left.apply(p), right.apply(p)));
            if (gap > junction_tol(k))
                return detail::fail(BuildErrc::NonGeneric, k == 0 ? "make_pipe: closed path along which the profile does not come back onto itself"
                                                                  : "make_pipe: the frame turns the profile at a junction of the path (Frenet across a change of the principal normal)");
            sew_tol = std::max(sew_tol, gap);
        }

        // folds: the speed of a point q of the profile along the path is t (1 - k . (q - c)) across;
        // the curvature by differences of the tangent where the path has no second derivatives
        const bool second = curve->has_second_derivatives();
        std::vector<T> lengths(n);
        for (std::size_t k = 0; k < n; ++k)
        {
            auto prev = curve->at(k, T(k));
            for (int i = 0; i <= 64; ++i)
            {
                const T s = i == 64 ? end_of<T>(k) : T(k) + T(i) / T(64);
                const auto d1 = curve->at(k, s, 1);
                point<T, 3> d2;
                if (second)
                    d2 = curve->at(k, s, 2);
                else
                {
                    const T a = std::max(T(k), s - T(1e-4)), b = std::min(T(k + 1), s + T(1e-4));
                    d2 = (curve->at(k, b, 1) - curve->at(k, a, 1)) / (b - a);
                }
                const T l2 = d1 * d1;
                const auto kv = (d2 - ((d2 * d1) / l2) * d1) / l2;
                const auto p = mo(s);
                const auto c = curve->at(k, s);
                lengths[k] += norm(c - prev);
                prev = c;
                for (const auto &q : pts)
                    if (T(1) - kv * (p.apply(q) - c) <= T(1e-9))
                        return detail::fail(BuildErrc::NonGeneric, "make_pipe: the profile is beyond the radius of curvature of the path (the swept surface folds)");
            }
        }

        // self-intersection away from the folds: stations about half the size of the profile apart,
        // at most 200000 quads; the quads closer than pi times the size along the path are left to the folds
        {
            std::size_t segs{};
            for (const auto &pl : polylines)
                segs += pl.size() - 1;
            std::vector<std::size_t> per(n);
            std::size_t total{};
            for (std::size_t k = 0; k < n; ++k)
                total += per[k] = std::max<std::size_t>(2, std::size_t(std::ceil(lengths[k] / (rho / 2))));
            const std::size_t cap = std::max<std::size_t>(2 * n, 200000 / std::max<std::size_t>(segs, 1));
            if (total > cap)
                for (auto &p : per)
                    p = std::max<std::size_t>(1, p * cap / total);
            std::vector<Placement<T>> P;
            std::vector<T> ell;
            auto prev = x0;
            auto station = [&](T s) {
                const auto c = curve->value(s);
                ell.push_back(ell.empty() ? T(0) : ell.back() + norm(c - prev));
                prev = c;
                P.push_back(mo(s));
            };
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t i = 0; i < per[k]; ++i)
                    station(T(k) + T(i) / T(per[k]));
            station(end_of<T>(n - 1));
            T step{};
            for (std::size_t j = 1; j < ell.size(); ++j)
                step = std::max(step, ell[j] - ell[j - 1]);
            if (self_intersects(polylines, P, ell, std::numbers::pi_v<T> * rho + T(2) * step, closed))
                return detail::fail(BuildErrc::NonGeneric, "make_pipe: the swept surface intersects itself (the path comes back within the size of the profile)");
        }

        // along each path edge, the same for all the edges of the profile: a prism along a segment where
        // the frame translates, a revolution along an arc where it turns with it, else a skinning
        struct Along
        {
            Placement<T> p0;
            std::optional<point<T, 3>> translation;
            std::optional<Arc<T>> rotation;
        };
        std::vector<Along> along(n);
        for (std::size_t k = 0; k < n; ++k)
        {
            const T s0 = T(k);
            // the quarters of the piece, the last one at its end (the frame of this piece)
            auto quarter = [&](int i) { return i == 4 ? end_of<T>(k) : s0 + T(i) / T(4); };
            auto &a = along[k];
            a.p0 = mo(s0);
            const T size = T(1) + rho + norm(curve->at(k, s0));
            if (sweep_detail::segment(m.edge(pieces[k].edge), o.tol))
            {
                const auto d = curve->at(k, s0 + T(1)) - curve->at(k, s0);
                bool translation = true;
                for (int i = 1; i <= 4 && translation; ++i)
                    translation = same(mo(quarter(i)), Placement<T>{.t = (T(i) / T(4)) * d} * a.p0, T(1e-9), size);
                if (translation)
                    a.translation = d;
            }
            else if (const auto arc = arc_of(*curve, k, o.tol))
            {
                bool rotation = true;
                for (int i = 1; i <= 4 && rotation; ++i)
                    rotation = same(mo(quarter(i)), sweep_detail::rotation<T>(arc->o, arc->z, arc->angle[std::size_t(i)]) * a.p0, T(1e-9), size);
                if (rotation)
                    a.rotation = arc;
            }
        }

        const SweepOptions<T> so{o.tol, o.approx_tol};
        StationCache<T> cache{mo, n};
        auto surfaces_of = [&](const Edge<T> &e, std::size_t &approximated, T &deviation) {
            std::vector<std::shared_ptr<Surface<T, 3>>> out;
            for (std::size_t k = 0; k < n; ++k)
            {
                const auto &a = along[k];
                if (a.translation)
                    out.push_back(sweep_detail::prism_surface(placed(e, a.p0), *a.translation, so, approximated, deviation));
                else if (a.rotation)
                    out.push_back(sweep_detail::revolution_surface(placed(e, a.p0), a.rotation->o, a.rotation->z, a.rotation->angle[4], so,
                                                                   approximated, deviation));
                else
                    out.push_back(skinned(e, cache, k, o, approximated, deviation));
            }
            return out;
        };
        return sweep_detail::sweep(m, profile, so, surfaces_of, closed ? std::optional<Placement<T>>{} : std::optional{mo(end_of<T>(n - 1))},
                                   sew_tol);
    }
}
