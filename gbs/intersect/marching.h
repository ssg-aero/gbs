#pragma once
/**
 * @file marching.h
 * @brief Intersection of two surfaces by marching, for the pairs without a closed form.
 *
 * 1. Starting points: both surfaces are cut into patches (`SurfacePieces`); pairs of patches
 *    whose boxes overlap are halved until both are flat, then a Gauss-Newton on the four
 *    parameters (minimum norm) puts the centre of the pair on the intersection. Nearly
 *    parallel flat patches may hide a small closed loop or a tangency: they are halved further.
 * 2. Marching from each starting point that is not on a curve already traced, in both
 *    directions: prediction along n1 x n2, correction by Newton on (u1, v1, u2, v2) in the plane
 *    normal to the step, step adapted to the turn of the tangent. A branch stops on the boundary
 *    of a box (the last point is put on it), when it comes back to its start (closed curve) or on
 *    a singular point where the normals are parallel (tangency, branch point).
 * 3. Curves: a cubic B-spline through the 3D points and two through the parameters on each
 *    surface, on the same chord-length parameter, densified until they agree within the
 *    tolerance.
 *
 * Periodic parameters are not unwrapped: a curve crossing the seam of a box stops there, as the
 * curves of the closed forms are cut at the seams.
 */
#include <gbs/intersect/intersection_curve.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/bscinterp.h>
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    namespace intersect_detail
    {
        /// Solves A x = b (n x n, row major) by Gaussian elimination with partial pivoting.
        template <std::floating_point T, std::size_t n>
        auto solve(std::array<std::array<T, n>, n> A, std::array<T, n> b) -> std::optional<std::array<T, n>>
        {
            T scale{};
            for (const auto &r : A)
                for (T x : r)
                    scale = std::max(scale, std::abs(x));
            for (std::size_t c = 0; c < n; ++c)
            {
                std::size_t piv = c;
                for (std::size_t r = c + 1; r < n; ++r)
                    if (std::abs(A[r][c]) > std::abs(A[piv][c]))
                        piv = r;
                if (!(std::abs(A[piv][c]) > T(1e-14) * scale))
                    return std::nullopt;
                std::swap(A[c], A[piv]);
                std::swap(b[c], b[piv]);
                for (std::size_t r = c + 1; r < n; ++r)
                {
                    const T f = A[r][c] / A[c][c];
                    for (std::size_t k = c; k < n; ++k)
                        A[r][k] -= f * A[c][k];
                    b[r] -= f * b[c];
                }
            }
            std::array<T, n> x{};
            for (std::size_t c = n; c-- > 0;)
            {
                T s = b[c];
                for (std::size_t k = c + 1; k < n; ++k)
                    s -= A[c][k] * x[k];
                x[c] = s / A[c][c];
            }
            return x;
        }

        /// A point of the intersection: 3D point and parameters (u1, v1, u2, v2).
        template <std::floating_point T>
        struct SSPoint
        {
            point<T, 3> p;
            std::array<T, 4> q;
        };

        template <std::floating_point T>
        class Marcher
        {
            const Surface<T, 3> &s1_, &s2_;
            std::array<T, 4> b1_, b2_;
            T tol_, sing_, max_turn_, max_step_;

        public:
            Marcher(const Surface<T, 3> &s1, const std::array<T, 4> &b1, const Surface<T, 3> &s2, const std::array<T, 4> &b2,
                    const IntersectOptions<T> &o, T size)
                : s1_{s1}, s2_{s2}, b1_{b1}, b2_{b2}, tol_{o.tol}, sing_{o.angular_tol}, max_turn_{o.max_turn}, max_step_{size / T(10)} {}

            [[nodiscard]] T tol() const { return tol_; }

            /// Lower and upper bounds of the four parameters.
            [[nodiscard]] T lo(int k) const { return k < 2 ? b1_[2 * k] : b2_[2 * (k - 2)]; }
            [[nodiscard]] T hi(int k) const { return k < 2 ? b1_[2 * k + 1] : b2_[2 * (k - 2) + 1]; }
            [[nodiscard]] T slo(int k) const { return k < 2 ? s1_.bounds()[2 * k] : s2_.bounds()[2 * (k - 2)]; }
            [[nodiscard]] T shi(int k) const { return k < 2 ? s1_.bounds()[2 * k + 1] : s2_.bounds()[2 * (k - 2) + 1]; }
            [[nodiscard]] bool inside(const std::array<T, 4> &q) const { return in_box(b1_, q[0], q[1]) && in_box(b2_, q[2], q[3]); }

            /// Derivatives: columns of J = d(S1 - S2) / d(u1, v1, u2, v2), and the normals.
            struct Local
            {
                point<T, 3> p1, p2;
                std::array<point<T, 3>, 4> J;
                point<T, 3> n1, n2;
            };
            [[nodiscard]] Local local(const std::array<T, 4> &q) const
            {
                Local l;
                l.p1 = s1_.value(q[0], q[1]);
                l.p2 = s2_.value(q[2], q[3]);
                l.J[0] = s1_.value(q[0], q[1], 1, 0);
                l.J[1] = s1_.value(q[0], q[1], 0, 1);
                l.J[2] = T(-1) * s2_.value(q[2], q[3], 1, 0);
                l.J[3] = T(-1) * s2_.value(q[2], q[3], 0, 1);
                l.n1 = cross(l.J[0], l.J[1]);
                l.n2 = cross(l.J[2], l.J[3]);
                return l;
            }

            /// Sine of the angle between the normals (0 at a degenerate point).
            [[nodiscard]] T sine(const Local &l) const
            {
                const T a = norm(l.n1), b = norm(l.n2);
                if (!(a > T(0)) || !(b > T(0)))
                    return T(0);
                return norm(cross(l.n1, l.n2)) / (a * b);
            }

            /// Unit tangent n1 x n2.
            [[nodiscard]] point<T, 3> tangent(const Local &l) const
            {
                const auto t = cross(l.n1, l.n2);
                return t / norm(t);
            }

            /**
             * Newton onto the intersection. Constraint: none (minimum-norm step, parameters kept in
             * the boxes), the plane through `c` normal to `t`, or the parameter `fixed` set to `value`.
             */
            [[nodiscard]] auto correct(std::array<T, 4> q, const point<T, 3> *c, const point<T, 3> *t, int fixed = -1, T value = T(0)) const
                -> std::optional<SSPoint<T>>
            {
                if (fixed >= 0)
                    q[fixed] = value;
                // a point off a surface is not evaluated: the caller lands on the boundary of the box
                auto off = [&](const std::array<T, 4> &x) -> std::optional<SSPoint<T>> {
                    std::array<T, 4> qc;
                    bool out = false;
                    for (int k = 0; k < 4; ++k)
                    {
                        qc[k] = std::clamp(x[k], slo(k), shi(k));
                        out = out || qc[k] != x[k];
                    }
                    if (!out)
                        return std::nullopt;
                    const auto lc = local(qc);
                    return SSPoint<T>{(lc.p1 + lc.p2) / T(2), x};
                };
                if (!c && fixed < 0)
                    for (int k = 0; k < 4; ++k)
                        q[k] = std::clamp(q[k], lo(k), hi(k));
                else if (auto o = off(q))
                    return o;
                for (int it = 0; it < (c || fixed >= 0 ? 30 : 100); ++it)
                {
                    const auto l = local(q);
                    const auto F = l.p1 - l.p2;
                    std::array<T, 4> dx{};
                    if (fixed >= 0)
                    {
                        std::array<std::array<T, 3>, 3> A;
                        std::array<int, 3> col;
                        for (int k = 0, m = 0; k < 4; ++k)
                            if (k != fixed)
                                col[m++] = k;
                        for (int r = 0; r < 3; ++r)
                            for (int m = 0; m < 3; ++m)
                                A[r][m] = l.J[col[m]][r];
                        const auto x = solve<T, 3>(A, {-F[0], -F[1], -F[2]});
                        if (!x)
                            return std::nullopt;
                        for (int m = 0; m < 3; ++m)
                            dx[col[m]] = (*x)[m];
                    }
                    else if (c)
                    {
                        std::array<std::array<T, 4>, 4> A{};
                        for (int r = 0; r < 3; ++r)
                            for (int k = 0; k < 4; ++k)
                                A[r][k] = l.J[k][r];
                        A[3] = {l.J[0] * *t, l.J[1] * *t, T(0), T(0)};
                        const auto x = solve<T, 4>(A, {-F[0], -F[1], -F[2], -((l.p1 - *c) * *t)});
                        if (!x)
                            return std::nullopt;
                        dx = *x;
                    }
                    else
                    {
                        // minimum norm least squares by a truncated SVD: where the surfaces are tangent
                        // J has rank 2, and the step goes to the closest points
                        Eigen::Matrix<T, 3, 4> Jm;
                        for (int r = 0; r < 3; ++r)
                            for (int k = 0; k < 4; ++k)
                                Jm(r, k) = l.J[k][r];
                        Eigen::JacobiSVD<Eigen::Matrix<T, 3, 4>> svd(Jm, Eigen::ComputeFullU | Eigen::ComputeFullV);
                        svd.setThreshold(T(1e-8));
                        const Eigen::Matrix<T, 4, 1> x = svd.solve(Eigen::Matrix<T, 3, 1>{-F[0], -F[1], -F[2]});
                        for (int k = 0; k < 4; ++k)
                            dx[k] = x(k);
                    }
                    T moved{};
                    bool out = false;
                    for (int k = 0; k < 4; ++k)
                    {
                        q[k] += dx[k];
                        if (!std::isfinite(q[k]))
                            return std::nullopt;
                        if (!c && fixed < 0)
                            q[k] = std::clamp(q[k], lo(k), hi(k)); // a seed stays in the boxes
                        else if (q[k] < slo(k) || q[k] > shi(k))
                            out = true;
                        moved = std::max(moved, std::abs(dx[k]) * norm(l.J[k]));
                    }
                    if (out)
                        return off(q);
                    if (moved <= tol_ * T(1e-4))
                    {
                        const auto l2 = local(q);
                        if (norm(l2.p1 - l2.p2) > tol_ * T(1e-2))
                            return std::nullopt;
                        return SSPoint<T>{(l2.p1 + l2.p2) / T(2), q};
                    }
                }
                return std::nullopt;
            }

            enum class End
            {
                boundary,
                closed,
                singular,
                capped
            };

            /// Marching from `start` in the direction `sign` of n1 x n2. The first point is `start`.
            [[nodiscard]] auto march(const SSPoint<T> &start, T sign, End &end) const -> std::vector<SSPoint<T>>
            {
                std::vector<SSPoint<T>> pts{start};
                auto l = local(start.q);
                if (sine(l) < sing_)
                {
                    end = End::singular;
                    return pts;
                }
                point<T, 3> t = sign * tangent(l);
                T h = max_step_ / T(4);
                const T h_min = T(10) * tol_;
                for (std::size_t n = 0; n < 100000; ++n)
                {
                    const auto &cur = pts.back();
                    // prediction: the step on both surfaces, by least squares on their tangent planes
                    std::array<T, 4> q = cur.q;
                    for (int s = 0; s < 2; ++s)
                    {
                        const auto &a = l.J[2 * s], &b = l.J[2 * s + 1];
                        const T sg = s == 0 ? T(1) : T(-1); // J holds -dS2
                        const T A = a * a, B = a * b, C = b * b, det = A * C - B * B;
                        if (!(det > T(0)))
                            continue;
                        const auto r = sg * h * t;
                        q[2 * s] += (C * (a * r) - B * (b * r)) / det;
                        q[2 * s + 1] += (A * (b * r) - B * (a * r)) / det;
                    }
                    const auto target = cur.p + h * t;
                    auto next = correct(q, &target, &t);
                    bool accepted = false;
                    if (next && norm(next->p - cur.p) <= T(2) * h && !inside(next->q))
                    {
                        // leaving a box: land on its boundary, at the first parameter to cross it
                        int k_out = -1;
                        T ratio = T(2);
                        for (int k = 0; k < 4; ++k)
                        {
                            const T bound = next->q[k] < lo(k) ? lo(k) : (next->q[k] > hi(k) ? hi(k) : next->q[k]);
                            if (bound != next->q[k])
                            {
                                const T r = (bound - cur.q[k]) / (next->q[k] - cur.q[k]);
                                if (r < ratio)
                                {
                                    ratio = r;
                                    k_out = k;
                                }
                            }
                        }
                        std::array<T, 4> guess;
                        for (int k = 0; k < 4; ++k)
                            guess[k] = cur.q[k] + ratio * (next->q[k] - cur.q[k]);
                        const T bound = next->q[k_out] < lo(k_out) ? lo(k_out) : hi(k_out);
                        auto last = correct(guess, nullptr, nullptr, k_out, bound);
                        // the last step obeys the same turn limit as the others
                        auto turn_to = [&](const SSPoint<T> &x) {
                            const auto lx = local(x.q);
                            if (sine(lx) < sing_)
                                return T(0);
                            const auto tx = tangent(lx);
                            return std::acos(std::clamp(std::abs(tx * t), T(0), T(1)));
                        };
                        if (last && inside(last->q) && norm(last->p - cur.p) <= T(2) * h && turn_to(*last) <= max_turn_)
                        {
                            if (norm(last->p - cur.p) > tol_)
                                pts.push_back(*last);
                            end = End::boundary;
                            return pts;
                        }
                    }
                    else if (next && norm(next->p - cur.p) <= T(2) * h)
                    {
                        const auto ln = local(next->q);
                        if (sine(ln) < sing_)
                        {
                            end = End::singular;
                            return pts;
                        }
                        auto tn = tangent(ln);
                        if (tn * t < T(0))
                            tn = T(-1) * tn;
                        const T turn = std::acos(std::clamp(tn * t, T(-1), T(1)));
                        // a dip of the angle between the normals inside the step: a singular point
                        // may be stepped over, shorten the step until it is reached
                        std::array<T, 4> qm;
                        for (int k = 0; k < 4; ++k)
                            qm[k] = (cur.q[k] + next->q[k]) / T(2);
                        const T s_mid = sine(local(qm)), s_ends = std::min(sine(l), sine(ln));
                        const bool dip = s_mid < T(0.1) && s_mid < T(0.5) * s_ends;
                        if (turn <= max_turn_ && !dip)
                        {
                            accepted = true;
                            // back to the start: the loop is closed
                            // (the parameters must come back too: a curve may meet its start across a seam)
                            const auto to_start = start.p - cur.p;
                            bool params_back = true;
                            for (int k = 0; k < 4; ++k)
                                params_back = params_back && std::abs(start.q[k] - cur.q[k]) <= T(1.5) * std::abs(next->q[k] - cur.q[k]) + T(1e-12);
                            if (pts.size() > 3 && params_back && norm(to_start) <= T(1.05) * h && to_start * t > T(0))
                            {
                                pts.push_back(start);
                                end = End::closed;
                                return pts;
                            }
                            pts.push_back(*next);
                            l = ln;
                            t = tn;
                            h = std::min(max_step_, h * std::clamp(max_turn_ / std::max(turn, T(1e-12)) * T(0.8), T(0.5), T(2)));
                        }
                    }
                    if (!accepted)
                    {
                        h /= T(2);
                        if (h < h_min)
                        {
                            end = End::singular;
                            return pts;
                        }
                    }
                }
                end = End::capped;
                return pts;
            }

            /// A point on the intersection between two consecutive points (for densification).
            [[nodiscard]] auto between(const SSPoint<T> &a, const SSPoint<T> &b) const -> std::optional<SSPoint<T>>
            {
                std::array<T, 4> q;
                for (int k = 0; k < 4; ++k)
                    q[k] = (a.q[k] + b.q[k]) / T(2);
                const auto c = (a.p + b.p) / T(2);
                const auto d = b.p - a.p;
                const T nd = norm(d);
                if (!(nd > T(0)))
                    return std::nullopt;
                const auto t = d / nd;
                auto r = correct(q, &c, &t);
                if (r && !inside(r->q))
                    return std::nullopt;
                return r;
            }
        };

        /// Distance of p to the polyline of a branch.
        template <std::floating_point T>
        T polyline_distance(const std::vector<SSPoint<T>> &pts, const point<T, 3> &p)
        {
            T d = std::numeric_limits<T>::max();
            for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                d = std::min(d, distance_to_segment(p, pts[k].p, pts[k + 1].p) - T(0.01) * norm(pts[k + 1].p - pts[k].p));
            if (pts.size() == 1)
                d = norm(p - pts[0].p);
            return d;
        }

        /// Curves through the points of a branch, densified until they agree within the tolerance.
        template <std::floating_point T>
        auto branch_curve(const Marcher<T> &m, const Surface<T, 3> &s1, const Surface<T, 3> &s2, std::vector<SSPoint<T>> pts, bool closed)
            -> IntersectionCurve<T>
        {
            IntersectionCurve<T> ic;
            const T tol = m.tol();
            for (int round = 0; round < 8; ++round)
            {
                std::vector<T> t{T(0)};
                for (std::size_t k = 1; k < pts.size(); ++k)
                    t.push_back(t.back() + std::max(norm(pts[k].p - pts[k - 1].p), std::numeric_limits<T>::epsilon()));
                const std::size_t p = std::min<std::size_t>(3, pts.size() - 1);
                std::vector<constrType<T, 3, 1>> Q(pts.size());
                std::vector<constrType<T, 2, 1>> U1(pts.size()), U2(pts.size());
                for (std::size_t k = 0; k < pts.size(); ++k)
                {
                    Q[k] = {pts[k].p};
                    U1[k] = {point<T, 2>{pts[k].q[0], pts[k].q[1]}};
                    U2[k] = {point<T, 2>{pts[k].q[2], pts[k].q[3]}};
                }
                auto c = std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(Q, t, p));
                auto pc1 = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(U1, t, p));
                auto pc2 = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(U2, t, p));
                // deviation on each interval, at its quarter points
                T dev{};
                std::vector<bool> refine(pts.size(), false);
                for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                    for (T f : {T(0.25), T(0.5), T(0.75)})
                    {
                        const T x = t[k] + f * (t[k + 1] - t[k]);
                        const auto P = c->value(x);
                        const auto a = pc1->value(x), b = pc2->value(x);
                        const T d = std::max(norm(s1.value(a[0], a[1]) - P), norm(s2.value(b[0], b[1]) - P));
                        dev = std::max(dev, d);
                        if (d > tol / T(2))
                            refine[k] = true;
                    }
                ic.curve = c;
                ic.pc1 = pc1;
                ic.pc2 = pc2;
                ic.a = t.front();
                ic.b = t.back();
                ic.deviation = dev;
                if (dev <= tol / T(2) || round == 7)
                    break;
                std::vector<SSPoint<T>> denser{pts.front()};
                for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                {
                    if (refine[k])
                        if (auto mid = m.between(pts[k], pts[k + 1]))
                            denser.push_back(*mid);
                    denser.push_back(pts[k + 1]);
                }
                if (denser.size() == pts.size())
                    break;
                pts = std::move(denser);
            }
            ic.closed = closed || norm(pts.front().p - pts.back().p) <= tol;
            return ic;
        }

        /// Intersection of two surfaces by marching.
        template <std::floating_point T>
        auto marching(const Surface<T, 3> &s1, const std::array<T, 4> &b1, const Surface<T, 3> &s2, const std::array<T, 4> &b2,
                      const IntersectOptions<T> &o) -> SurfaceSurfaceResult<T>
        {
            SurfaceSurfaceResult<T> res;
            const T tol = o.tol, small = T(10) * tol;
            SurfacePieces<T> P1{s1}, P2{s2};
            const auto pieces1 = P1.pieces(b1), pieces2 = P2.pieces(b2);
            Box<T, 3> box1, box2;
            for (const auto &p : pieces1)
                box1.add(p.box);
            for (const auto &p : pieces2)
                box2.add(p.box);
            const Marcher<T> m{s1, b1, s2, b2, o, std::min(box1.diagonal(), box2.diagonal())};

            // ---- starting points
            using Piece = SurfacePiece<T>;
            auto normal = [](const Piece &pc) {
                const auto n = cross(pc.pole(pc.nu - 1, pc.nv - 1) - pc.pole(0, 0), pc.pole(0, pc.nv - 1) - pc.pole(pc.nu - 1, 0));
                const T nn = norm(n);
                return nn > T(0) ? n / nn : n;
            };
            auto cone = [](const Piece &pc) { return std::asin(std::min(T(1), T(4) * pc.flatness / std::max(pc.size(), std::numeric_limits<T>::min()))); };
            std::vector<SSPoint<T>> seeds;
            std::vector<std::pair<Piece, Piece>> stack;
            for (const auto &p : pieces1)
                for (const auto &q : pieces2)
                    if (p.box.overlaps(q.box, tol))
                        stack.emplace_back(p, q);
            std::size_t pairs{};
            while (!stack.empty())
            {
                auto [p, q] = std::move(stack.back());
                stack.pop_back();
                if (++pairs > o.max_pairs)
                {
                    res.complete = false;
                    break;
                }
                if (!p.box.overlaps(q.box, tol))
                    continue;
                const bool p_flat = p.flatness <= T(0.02) * p.size() || p.size() <= small;
                const bool q_flat = q.flatness <= T(0.02) * q.size() || q.size() <= small;
                if (p_flat && q_flat)
                {
                    // nearly parallel patches may hide a small closed loop or a tangency: when no
                    // starting point is found at their centre, they are halved further
                    const T angle = std::acos(std::clamp(std::abs(normal(p) * normal(q)), T(0), T(1)));
                    const bool parallel = angle <= cone(p) + cone(q) + T(0.05);
                    const auto cp = p.center(), cq = q.center();
                    auto sp = m.correct({cp[0], cp[1], cq[0], cq[1]}, nullptr, nullptr);
                    if (sp && m.inside(sp->q))
                    {
                        seeds.push_back(*sp);
                        continue;
                    }
                    if (!parallel || (p.size() <= T(100) * small && q.size() <= T(100) * small))
                        continue;
                }
                // halve the larger patch (both when parallel and flat)
                const bool split_p = q_flat && !p_flat ? true : (p_flat && !q_flat ? false : p.size() >= q.size());
                if (split_p)
                {
                    auto [l, r] = P1.split(p, SurfacePieces<T>::longer_along_u(p));
                    stack.emplace_back(std::move(l), q);
                    stack.emplace_back(std::move(r), std::move(q));
                }
                else
                {
                    auto [l, r] = P2.split(q, SurfacePieces<T>::longer_along_u(q));
                    stack.emplace_back(p, std::move(l));
                    stack.emplace_back(std::move(p), std::move(r));
                }
            }

            // ---- marching from the seeds not on a traced branch
            std::vector<std::vector<SSPoint<T>>> branches;
            std::vector<Box<T, 3>> boxes;
            auto traced = [&](const point<T, 3> &p) {
                for (std::size_t b = 0; b < branches.size(); ++b)
                {
                    auto bb = boxes[b];
                    bb.inflate(small);
                    if (bb.contains(p) && polyline_distance(branches[b], p) <= small)
                        return true;
                }
                return false;
            };
            for (const auto &seed : seeds)
            {
                if (traced(seed.p))
                    continue;
                typename Marcher<T>::End e1{}, e2{};
                auto fwd = m.march(seed, T(1), e1);
                std::vector<SSPoint<T>> pts;
                bool closed = e1 == Marcher<T>::End::closed;
                if (closed)
                    pts = std::move(fwd);
                else
                {
                    auto bwd = m.march(seed, T(-1), e2);
                    pts.assign(bwd.rbegin(), bwd.rend());
                    pts.insert(pts.end(), fwd.begin() + 1, fwd.end());
                }
                if (e1 == Marcher<T>::End::singular || e2 == Marcher<T>::End::singular)
                    res.singular = true;
                if (e1 == Marcher<T>::End::capped || e2 == Marcher<T>::End::capped)
                    res.complete = false;
                // box of the branch, grown by the sag of the curve between its points
                Box<T, 3> bb;
                T longest{};
                for (std::size_t k = 0; k < pts.size(); ++k)
                {
                    bb.add(pts[k].p);
                    if (k > 0)
                        longest = std::max(longest, norm(pts[k].p - pts[k - 1].p));
                }
                bb.inflate(T(0.01) * longest);
                if (pts.size() == 1 || bb.diagonal() <= small)
                {
                    // an isolated contact point
                    const auto &x = pts.front();
                    const bool known = std::ranges::any_of(res.points, [&](const auto &h) { return norm(h.p - x.p) <= small; });
                    if (!known)
                        res.points.push_back({x.q[0], x.q[1], x.q[2], x.q[3], x.p});
                    branches.push_back({x});
                    boxes.push_back(bb);
                    continue;
                }
                branches.push_back(pts);
                boxes.push_back(bb);
                res.curves.push_back(branch_curve(m, s1, s2, std::move(pts), closed));
            }
            return res;
        }
    }
}
