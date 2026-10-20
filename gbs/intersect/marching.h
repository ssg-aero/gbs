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
 *    a singular point where the normals are parallel. A branch point (curves crossing where the
 *    surfaces are tangent) is stepped over straight on: the step that lands beyond it, on the
 *    same branch, with the sine of the angle between the normals dipping in between, is kept and
 *    its interval reported (`IntersectionCurve::branch_spans`). A landing where the curve only
 *    touches a bound is put on the contact.
 * 2b. A start where the surfaces are tangent may lie on a curve of tangency (a cylinder lying on
 *    a face, a torus in a cylinder): it is traced by `TangentMarcher`, on the equations of the
 *    foot of S1 on S2 and of parallel normals, in the direction where the surfaces bend alike
 *    (`IntersectionCurve::tangent`). Where the tangency ends inside the boxes, the branch ends
 *    there if the surfaces part beyond (a free end), and is singular if they cross (two
 *    transversal curves start there).
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
#include <span>
#include <utility>
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
            std::array<T, 4> b1_, b2_, sb1_, sb2_;
            T tol_, sing_, max_turn_, max_step_;

        public:
            Marcher(const Surface<T, 3> &s1, const std::array<T, 4> &b1, const Surface<T, 3> &s2, const std::array<T, 4> &b2,
                    const IntersectOptions<T> &o, T size)
                : s1_{s1}, s2_{s2}, b1_{b1}, b2_{b2}, sb1_{s1.bounds()}, sb2_{s2.bounds()}, tol_{o.tol}, sing_{o.angular_tol}, max_turn_{o.max_turn}, max_step_{size / T(10)} {}

            [[nodiscard]] T tol() const { return tol_; }
            [[nodiscard]] T angular_tol() const { return sing_; }
            [[nodiscard]] T max_turn() const { return max_turn_; }
            [[nodiscard]] T max_step() const { return max_step_; }

            /// Lower and upper bounds of the four parameters.
            [[nodiscard]] T lo(int k) const { return k < 2 ? b1_[2 * k] : b2_[2 * (k - 2)]; }
            [[nodiscard]] T hi(int k) const { return k < 2 ? b1_[2 * k + 1] : b2_[2 * (k - 2) + 1]; }
            [[nodiscard]] T slo(int k) const { return k < 2 ? sb1_[2 * k] : sb2_[2 * (k - 2)]; }
            [[nodiscard]] T shi(int k) const { return k < 2 ? sb1_[2 * k + 1] : sb2_[2 * (k - 2) + 1]; }
            [[nodiscard]] bool inside(const std::array<T, 4> &q) const { return in_box(b1_, q[0], q[1]) && in_box(b2_, q[2], q[3]); }
            /// The rounding margin of the parameter k, that of `in_box`.
            [[nodiscard]] T margin(int k) const { return T(1e-10) * (T(1) + std::abs(lo(k)) + std::abs(hi(k))); }
            /// The parameter k is out of its box, beyond the rounding margin.
            [[nodiscard]] bool out_of(int k, T x) const { return x < lo(k) - margin(k) || x > hi(k) + margin(k); }
            /**
             * A curve along the border of a surface goes beyond it by rounding errors: the parameters
             * within the margin of a bound of the surface are put back on it, not taken as leaving
             * the surface. With the margin of `inside`, a point left beyond a bound is out of the box.
             */
            void snap(std::array<T, 4> &q) const
            {
                for (int k = 0; k < 4; ++k)
                    if (q[k] < slo(k) && q[k] >= slo(k) - margin(k))
                        q[k] = slo(k);
                    else if (q[k] > shi(k) && q[k] <= shi(k) + margin(k))
                        q[k] = shi(k);
            }

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
                snap(q);
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
                        // truncated SVD: the boundary may cross the curve at a branch point, where
                        // the surfaces are tangent and the system has rank 2
                        std::array<int, 3> col;
                        for (int k = 0, m = 0; k < 4; ++k)
                            if (k != fixed)
                                col[m++] = k;
                        Eigen::Matrix<T, 3, 3> A;
                        for (int r = 0; r < 3; ++r)
                            for (int m = 0; m < 3; ++m)
                                A(r, m) = l.J[col[m]][r];
                        Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
                        svd.setThreshold(T(1e-6));
                        const Eigen::Matrix<T, 3, 1> x = svd.solve(Eigen::Matrix<T, 3, 1>{-F[0], -F[1], -F[2]});
                        for (int m = 0; m < 3; ++m)
                            dx[col[m]] = x(m);
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
                        // J has rank 2, and the step goes to the closest points. A parameter on a bound
                        // of its box that the step pushes out is held there, the others solved again
                        // (the iterate slides along the border instead of stalling on it)
                        auto solve_free = [&](const std::array<bool, 4> &held) {
                            Eigen::Matrix<T, 3, 4> Jm;
                            for (int r = 0; r < 3; ++r)
                                for (int k = 0; k < 4; ++k)
                                    Jm(r, k) = held[k] ? T(0) : l.J[k][r];
                            Eigen::JacobiSVD<Eigen::Matrix<T, 3, 4>> svd(Jm, Eigen::ComputeFullU | Eigen::ComputeFullV);
                            svd.setThreshold(T(1e-8));
                            const Eigen::Matrix<T, 4, 1> x = svd.solve(Eigen::Matrix<T, 3, 1>{-F[0], -F[1], -F[2]});
                            for (int k = 0; k < 4; ++k)
                                dx[k] = held[k] ? T(0) : x(k);
                        };
                        std::array<bool, 4> held{};
                        solve_free(held);
                        for (bool more = true; more;) // at most four passes
                        {
                            more = false;
                            for (int k = 0; k < 4; ++k)
                                if (!held[k] && ((q[k] <= lo(k) && dx[k] < T(0)) || (q[k] >= hi(k) && dx[k] > T(0))))
                                    more = held[k] = true;
                            if (more)
                                solve_free(held);
                        }
                    }
                    T moved{};
                    bool out = false;
                    for (int k = 0; k < 4; ++k)
                    {
                        q[k] += dx[k];
                        if (!std::isfinite(q[k]))
                            return std::nullopt;
                    }
                    snap(q);
                    for (int k = 0; k < 4; ++k)
                    {
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
                capped,
                parted ///< a tangency ends where the surfaces part (`TangentMarcher` only)
            };

            /**
             * The point where the step from `cur` to `next` (out of the boxes) reaches their
             * boundary: at the first parameter to leave its box (beyond the rounding margin), with
             * that parameter fixed on the bound. `cur` may itself be out of the box by a rounding
             * error: the ratios are kept in [0, 1].
             *
             * A curve that only touches a bound, tangent to it, does not leave the box there: the
             * landing on such a bound is ill-conditioned and stops short, the next parameter to
             * leave is tried instead (the bound of a plane disk around a circle, which the circle
             * touches where it crosses the seam of the other surface). When every landing inside
             * the boxes is tangent, the first one ends the branch: a section has a vertex where a
             * curve touches the boundary of a face (the double point of Viviani's curve on the
             * border of a face). A landing off the boxes is no landing: the next exit is tried.
             */
            [[nodiscard]] auto land(const SSPoint<T> &cur, const SSPoint<T> &next) const -> std::optional<SSPoint<T>>
            {
                struct Exit
                {
                    T ratio, bound;
                    int k;
                };
                std::array<Exit, 4> buf;
                std::size_t n{};
                for (int k = 0; k < 4; ++k)
                {
                    if (!out_of(k, next.q[k]))
                        continue;
                    const T b = next.q[k] < lo(k) ? lo(k) : hi(k);
                    const T d = next.q[k] - cur.q[k];
                    buf[n++] = {d != T(0) ? std::clamp((b - cur.q[k]) / d, T(0), T(1)) : T(0), b, k};
                }
                const std::span exits{buf.data(), n};
                // (ties broken by the parameter, for the same order on every standard library)
                std::ranges::sort(exits, {}, [](const Exit &e) { return std::pair{e.ratio, e.k}; });
                std::optional<SSPoint<T>> tangent_landing;
                int k_tangent = -1;
                T b_tangent{};
                for (const auto &[ratio, bound, k_out] : exits)
                {
                    // (in the boxes: a guess off a surface is not corrected)
                    std::array<T, 4> guess;
                    for (int k = 0; k < 4; ++k)
                        guess[k] = std::clamp(cur.q[k] + ratio * (next.q[k] - cur.q[k]), lo(k), hi(k));
                    auto r = correct(guess, nullptr, nullptr, k_out, bound);
                    if (!r || !inside(r->q))
                        continue;
                    if (crosses(local(r->q), k_out))
                        return r;
                    if (!tangent_landing)
                    {
                        tangent_landing = r;
                        k_tangent = k_out;
                        b_tangent = bound;
                    }
                }
                // the landing on a touched bound stops short: the contact is where the parameter
                // is extremal along the curve
                if (tangent_landing)
                    if (auto c = touch(*tangent_landing, k_tangent, b_tangent); c && inside(c->q))
                        return c;
                return tangent_landing;
            }

            /**
             * The point near x where the curve touches the iso q_k = b: S1 = S2 and the tangent
             * along the iso, (n1 x n2) . (n x S_o) = 0 with n the normal of the surface of k and
             * S_o its derivative along the other parameter; regular where the curve bends away
             * from the iso, unlike the landing with q_k fixed. The contact is put on the bound
             * when it is within the tolerance of it.
             */
            [[nodiscard]] auto touch(const SSPoint<T> &x, int k, T b) const -> std::optional<SSPoint<T>>
            {
                auto g = [&](const Local &l) {
                    const auto &n = k < 2 ? l.n1 : l.n2;
                    const auto m = cross(n, l.J[k ^ 1]), t = cross(l.n1, l.n2);
                    const T s = norm(m) * norm(t);
                    return s > T(0) ? (t * m) / s : T(0);
                };
                auto q = x.q;
                for (int it = 0; it < 30; ++it)
                {
                    const auto l = local(q);
                    const auto F = l.p1 - l.p2;
                    const T g0 = g(l);
                    Eigen::Matrix<T, 4, 4> A;
                    Eigen::Matrix<T, 4, 1> r;
                    for (int c = 0; c < 4; ++c)
                    {
                        for (int i = 0; i < 3; ++i)
                            A(i, c) = l.J[c][i];
                        // (the derivatives of g by differences, which need first derivatives only:
                        // central, one-sided inwards against a bound of the surface)
                        const T h = T(1e-6) * (T(1) + std::abs(q[c]));
                        auto qp = q, qm = q;
                        qp[c] = std::min(q[c] + h, shi(c));
                        qm[c] = std::max(q[c] - h, slo(c));
                        const T gp = qp[c] != q[c] ? g(local(qp)) : g0, gm = qm[c] != q[c] ? g(local(qm)) : g0;
                        A(3, c) = (gp - gm) / (qp[c] - qm[c]);
                    }
                    // the last row in lengths: g is a sine, scaled by a step
                    A.row(3) *= max_step_;
                    r << -F[0], -F[1], -F[2], -g0 * max_step_;
                    const Eigen::Matrix<T, 4, 1> dx = A.fullPivLu().solve(r);
                    T moved{};
                    for (int c = 0; c < 4; ++c)
                    {
                        if (!std::isfinite(dx(c)))
                            return std::nullopt;
                        q[c] += dx(c);
                        moved = std::max(moved, std::abs(dx(c)) * norm(l.J[c]));
                    }
                    snap(q);
                    for (int c = 0; c < 4; ++c)
                        if (q[c] < slo(c) || q[c] > shi(c))
                            return std::nullopt;
                    if (moved <= tol_ * T(1e-4))
                    {
                        // on the bound, within the tolerance in 3D
                        if (std::abs(q[k] - b) * norm(l.J[k]) > tol_)
                            return std::nullopt;
                        q[k] = b;
                        const auto lt = local(q);
                        if (norm(lt.p1 - lt.p2) > tol_ * T(1e-2))
                            return std::nullopt;
                        return SSPoint<T>{(lt.p1 + lt.p2) / T(2), q};
                    }
                }
                return std::nullopt;
            }

            /**
             * Does the curve cross the iso of the parameter k at the point of `l`, rather than touch
             * it? The sine of the angle between the tangent and the iso, |t x S_o| / |S_o| with S_o
             * the derivative along the other parameter of the same surface, is above
             * `angular_tol`. At a singular point the tangent is unknown, and where the surface is
             * degenerate (a pole, an apex: S_k or S_o null next to the other, or both parallel)
             * the isos are not curves through the point: both taken as a crossing.
             */
            [[nodiscard]] bool crosses(const Local &l, int k) const
            {
                if (sine(l) < sing_)
                    return true;
                const auto &Sk = l.J[k], &So = l.J[k ^ 1];
                const T nk = norm(Sk), no = norm(So);
                if (!(nk > sing_ * no) || !(no > sing_ * nk) || !(norm(cross(Sk, So)) > sing_ * nk * no))
                    return true;
                return norm(cross(tangent(l), So)) > sing_ * no;
            }

            /// The parameters after a step h along t from q, by least squares on the tangent planes.
            [[nodiscard]] auto predict(std::array<T, 4> q, const Local &l, const point<T, 3> &t, T h) const -> std::array<T, 4>
            {
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
                return q;
            }

            /**
             * A step over a branch point ahead of `cur` (tangent t): the shortest step along t that
             * lands on the same branch beyond a dip of the angle between the normals; or, when the
             * branch point lies on the boundary of a box (a seam through it), the point where the
             * branch reaches that boundary (second member true).
             */
            [[nodiscard]] auto step_over(const SSPoint<T> &cur, const Local &l, const point<T, 3> &t) const
                -> std::optional<std::pair<SSPoint<T>, bool>>
            {
                const T s_cur = sine(l);
                for (T H = T(20) * tol_; H <= max_step_; H *= T(2))
                {
                    const auto target = cur.p + H * t;
                    const auto next = correct(predict(cur.q, l, t, H), &target, &t);
                    if (!next)
                        continue;
                    if (!inside(next->q))
                    {
                        // on the boundary at the branch point (the surfaces are tangent there: the
                        // correction takes minimum norm steps)
                        const auto last = land(cur, *next);
                        if (last && inside(last->q) && (last->p - cur.p) * t > T(0) && norm(last->p - cur.p) <= T(1.5) * H)
                            return std::pair{*last, true};
                        continue;
                    }
                    const auto d = next->p - cur.p;
                    if (norm(d) < H / T(2) || norm(d) > T(1.5) * H || d * t <= T(0))
                        continue;
                    const auto ln = local(next->q);
                    const T s_next = sine(ln);
                    if (s_next < sing_)
                        continue;
                    const auto tn = tangent(ln);
                    if (std::acos(std::clamp(std::abs(tn * t), T(0), T(1))) > T(4) * max_turn_)
                        continue;
                    std::array<T, 4> qm;
                    for (int k = 0; k < 4; ++k)
                        qm[k] = (cur.q[k] + next->q[k]) / T(2);
                    if (!(sine(local(qm)) < std::min(s_cur, s_next)))
                        continue; // no singular point in between
                    return std::pair{*next, false};
                }
                return std::nullopt;
            }

            /**
             * Marching from `start` in the direction `sign` of n1 x n2. The first point is `start`.
             * The steps over branch points are appended to `jumps` (from, to).
             */
            [[nodiscard]] auto march(const SSPoint<T> &start, T sign, End &end, std::vector<std::pair<SSPoint<T>, SSPoint<T>>> &jumps) const
                -> std::vector<SSPoint<T>>
            {
                std::vector<SSPoint<T>> pts{start};
                auto l = local(start.q);
                if (sine(l) < sing_)
                {
                    end = End::singular;
                    return pts;
                }
                point<T, 3> t = sign * tangent(l);
                T orient = sign; // t = orient n1 x n2 / |n1 x n2|
                T h = max_step_ / T(4);
                const T h_min = T(10) * tol_;
                for (std::size_t n = 0; n < 100000; ++n)
                {
                    const auto cur = pts.back();
                    // a singular point ahead: step over it if it is a branch point
                    // (0: not a branch point, 1: stepped over, 2: landed on a boundary at the branch point)
                    auto over = [&]() {
                        if (jumps.size() >= 64)
                            return 0;
                        const auto r = step_over(cur, l, t);
                        if (!r)
                            return 0;
                        const auto &next = r->first;
                        pts.push_back(next);
                        if (r->second)
                            return 2;
                        jumps.emplace_back(cur, next);
                        l = local(next.q);
                        auto tn = tangent(l);
                        orient = tn * t < T(0) ? T(-1) : T(1);
                        t = orient * tn;
                        h = std::max(h, norm(next.p - cur.p));
                        return 1;
                    };
                    // prediction: the step on both surfaces, by least squares on their tangent planes
                    const auto q = predict(cur.q, l, t, h);
                    const auto target = cur.p + h * t;
                    auto next = correct(q, &target, &t);
                    bool accepted = false;
                    if (next && norm(next->p - cur.p) <= T(2) * h && !inside(next->q))
                    {
                        // leaving a box: land on its boundary, at the first parameter to cross it
                        const auto last = land(cur, *next);
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
                            if (const int o = over(); o == 1)
                                continue;
                            else if (o == 2)
                            {
                                end = End::boundary;
                                return pts;
                            }
                            end = End::singular;
                            return pts;
                        }
                        auto tn = tangent(ln);
                        // n1 x n2 turning back: the step went through a branch point; it is kept
                        // short, and reported
                        const T o_next = tn * t < T(0) ? T(-1) : T(1);
                        const bool through = o_next != orient;
                        tn = o_next * tn;
                        const T turn = std::acos(std::clamp(tn * t, T(-1), T(1)));
                        // a dip of the angle between the normals inside the step: a singular point
                        // may be stepped over, shorten the step until it is reached
                        std::array<T, 4> qm;
                        for (int k = 0; k < 4; ++k)
                            qm[k] = (cur.q[k] + next->q[k]) / T(2);
                        const T s_mid = sine(local(qm)), s_ends = std::min(sine(l), sine(ln));
                        const bool dip = s_mid < T(0.1) && s_mid < T(0.5) * s_ends;
                        if (turn <= max_turn_ && !dip && (!through || h <= T(64) * tol_))
                        {
                            accepted = true;
                            if (through && jumps.size() < 64)
                                jumps.emplace_back(cur, *next);
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
                            orient = o_next;
                            h = std::min(max_step_, h * std::clamp(max_turn_ / std::max(turn, T(1e-12)) * T(0.8), T(0.5), T(2)));
                        }
                    }
                    if (!accepted)
                    {
                        h /= T(2);
                        if (h < h_min)
                        {
                            if (const int o = over(); o == 1)
                                continue;
                            else if (o == 2)
                            {
                                end = End::boundary;
                                return pts;
                            }
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

        /**
         * Tracing of a curve along which the surfaces are tangent. There n1 x n2 vanishes and the
         * system of `Marcher` is singular at every point. The curve is the solution set of
         *
         *     (S1 - S2) . S2_u = 0,  (S1 - S2) . S2_v = 0   (S2 is the foot of S1)
         *     n2 . S1_u = 0,         n2 . S1_v = 0           (parallel normals)
         *
         * which has rank 3 on it, with |S1 - S2| within the tolerance. Its direction is the null
         * direction of the difference of the second fundamental forms in the common tangent
         * plane: along the curve the surfaces bend alike. Each point is corrected on the plane
         * normal to the step (five equations, four unknowns, by least squares).
         */
        template <std::floating_point T>
        class TangentMarcher
        {
            const Marcher<T> &m_;
            const Surface<T, 3> &s1_, &s2_;

            struct Second
            {
                point<T, 3> p, u, v, uu, uv, vv;
            };
            static Second second(const Surface<T, 3> &s, T u, T v)
            {
                // (one pass on B-spline surfaces; SKL[a * 3 + b] is the derivative (a, b))
                std::array<point<T, 3>, 9> d;
                s.derivatives(u, v, 2, 2, d.data());
                return {d[0], d[3], d[1], d[6], d[4], d[2]};
            }

        public:
            TangentMarcher(const Marcher<T> &m, const Surface<T, 3> &s1, const Surface<T, 3> &s2) : m_{m}, s1_{s1}, s2_{s2} {}
            [[nodiscard]] bool inside(const std::array<T, 4> &q) const { return m_.inside(q); }

            /// Direction of the tangency at q (unit, in 3D), and the eigenvalues of the difference
            /// of the second fundamental forms, smallest first in absolute value.
            struct Direction
            {
                point<T, 3> t;
                T small, large;
                T bend; ///< the large eigenvalue, signed: (S1 - S2) . n2 grows as bend s^2 / 2 across the curve
            };
            [[nodiscard]] auto direction(const std::array<T, 4> &q) const -> std::optional<Direction>
            {
                const auto a = second(s1_, q[0], q[1]), b = second(s2_, q[2], q[3]);
                auto n = cross(b.u, b.v);
                const T nn = norm(n);
                if (!(nn > T(0)) || !(norm(b.u) > T(0)))
                    return std::nullopt;
                n = n / nn;
                const auto e1 = b.u / norm(b.u), e2 = cross(n, e1);
                // second fundamental form of a surface on the common basis (e1, e2)
                auto form = [&](const Second &s) -> std::optional<std::array<T, 3>> {
                    const T g11 = s.u * s.u, g12 = s.u * s.v, g22 = s.v * s.v, det = g11 * g22 - g12 * g12;
                    if (!(det > T(0)))
                        return std::nullopt;
                    // parameters of e1, e2 on the surface: (J^T J)^-1 J^T e
                    auto coef = [&](const point<T, 3> &e) {
                        const T r1 = s.u * e, r2 = s.v * e;
                        return std::array<T, 2>{(g22 * r1 - g12 * r2) / det, (g11 * r2 - g12 * r1) / det};
                    };
                    const auto c1 = coef(e1), c2 = coef(e2);
                    const T L = s.uu * n, M = s.uv * n, N = s.vv * n;
                    auto II = [&](const std::array<T, 2> &x, const std::array<T, 2> &y) { return L * x[0] * y[0] + M * (x[0] * y[1] + x[1] * y[0]) + N * x[1] * y[1]; };
                    return std::array<T, 3>{II(c1, c1), II(c1, c2), II(c2, c2)};
                };
                const auto f1 = form(a), f2 = form(b);
                if (!f1 || !f2)
                    return std::nullopt;
                const T A = (*f1)[0] - (*f2)[0], B = (*f1)[1] - (*f2)[1], C = (*f1)[2] - (*f2)[2];
                // eigenvalues of [[A, B], [B, C]]
                const T mean = (A + C) / T(2), rad = std::hypot((A - C) / T(2), B);
                T l1 = mean - rad, l2 = mean + rad;
                if (std::abs(l1) > std::abs(l2))
                    std::swap(l1, l2);
                // no direction where the surfaces bend alike every way (coincident, or umbilic contact)
                if (!(std::abs(l2) * m_.max_step() > T(1e-6)))
                    return std::nullopt;
                // eigenvector of l1: (B, l1 - A) or (l1 - C, B), the longer
                const T x1 = B, y1 = l1 - A, x2 = l1 - C, y2 = B;
                const auto t = std::hypot(x1, y1) >= std::hypot(x2, y2) ? x1 * e1 + y1 * e2 : x2 * e1 + y2 * e2;
                const T nt = norm(t);
                if (!(nt > T(0)))
                    return std::nullopt;
                return Direction{t / nt, std::abs(l1), std::abs(l2), l2};
            }

            /**
             * Gauss-Newton onto the tangency, on the plane through `c` normal to `t` (none: the
             * four equations only, minimum norm), or with the parameter `fixed` set to `value`.
             * Returns the point when |S1 - S2| is within the tolerance (at any distance with
             * `any_distance`: the extremum of the gap across the curve); or, when an iterate leaves
             * a surface, that iterate (its parameters beyond the bound, its point clamped), on
             * which the caller lands.
             */
            struct Correction
            {
                int fixed{-1};          ///< a parameter set to `value` (none: -1)
                T value{};
                bool any_distance{};    ///< the extremum of the gap, beyond the tolerance too
            };
            [[nodiscard]] auto correct(std::array<T, 4> q, const point<T, 3> *c, const point<T, 3> *t, Correction how = {}) const
                -> std::optional<SSPoint<T>>
            {
                const auto [fixed, value, any_distance] = how;
                if (fixed >= 0)
                    q[fixed] = value;
                const T tol = m_.tol();
                auto off = [&](const std::array<T, 4> &x) -> std::optional<SSPoint<T>> {
                    std::array<T, 4> qc;
                    bool out = false;
                    for (int k = 0; k < 4; ++k)
                    {
                        qc[k] = std::clamp(x[k], m_.slo(k), m_.shi(k));
                        out = out || qc[k] != x[k];
                    }
                    if (!out)
                        return std::nullopt;
                    return SSPoint<T>{(s1_.value(qc[0], qc[1]) + s2_.value(qc[2], qc[3])) / T(2), x};
                };
                for (int it = 0; it < 50; ++it)
                {
                    m_.snap(q);
                    if (auto o = off(q))
                        return fixed >= 0 ? std::nullopt : o;
                    const auto a = second(s1_, q[0], q[1]), b = second(s2_, q[2], q[3]);
                    const auto d = a.p - b.p;
                    auto N = cross(b.u, b.v);
                    const T nN = norm(N);
                    if (!(nN > T(0)))
                        return std::nullopt;
                    const auto n = N / nN;
                    const auto Nu = cross(b.uu, b.v) + cross(b.u, b.uv), Nv = cross(b.uv, b.v) + cross(b.u, b.vv);
                    const auto nu = (Nu - (n * Nu) * n) / nN, nv = (Nv - (n * Nv) * n) / nN;
                    // rows scaled to lengths: the foot equations by |S2_u|, |S2_v|, the parallel
                    // normals (sines) by the length of a step
                    const T su = norm(b.u), sv = norm(b.v), L = m_.max_step();
                    const T w1 = T(1) / su, w2 = T(1) / sv, w3 = L / norm(a.u), w4 = L / norm(a.v);
                    Eigen::Matrix<T, 5, 4> J;
                    Eigen::Matrix<T, 5, 1> F;
                    J.row(0) << w1 * (a.u * b.u), w1 * (a.v * b.u), w1 * (d * b.uu - b.u * b.u), w1 * (d * b.uv - b.v * b.u);
                    J.row(1) << w2 * (a.u * b.v), w2 * (a.v * b.v), w2 * (d * b.uv - b.u * b.v), w2 * (d * b.vv - b.v * b.v);
                    J.row(2) << w3 * (n * a.uu), w3 * (n * a.uv), w3 * (nu * a.u), w3 * (nv * a.u);
                    J.row(3) << w4 * (n * a.uv), w4 * (n * a.vv), w4 * (nu * a.v), w4 * (nv * a.v);
                    F << w1 * (d * b.u), w2 * (d * b.v), w3 * (n * a.u), w4 * (n * a.v), T(0);
                    int rows = 4;
                    if (c)
                    {
                        J.row(4) << a.u * *t, a.v * *t, T(0), T(0);
                        F(4) = (a.p - *c) * *t;
                        rows = 5;
                    }
                    if (fixed >= 0)
                        J.col(fixed).setZero();
                    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, 4>> svd(J.topRows(rows), Eigen::ComputeThinU | Eigen::ComputeThinV);
                    svd.setThreshold(T(1e-10));
                    const Eigen::Matrix<T, 4, 1> dx = svd.solve(-F.head(rows));
                    T moved{};
                    for (int k = 0; k < 4; ++k)
                    {
                        if (!std::isfinite(dx(k)))
                            return std::nullopt;
                        q[k] += dx(k);
                        moved = std::max(moved, std::abs(dx(k)) * (k < 2 ? norm(k == 0 ? a.u : a.v) : norm(k == 2 ? b.u : b.v)));
                    }
                    if (moved <= tol * T(1e-4))
                    {
                        m_.snap(q);
                        if (auto o = off(q))
                            return fixed >= 0 ? std::nullopt : o;
                        const auto p1 = s1_.value(q[0], q[1]), p2 = s2_.value(q[2], q[3]);
                        if (!any_distance && norm(p1 - p2) > tol)
                            return std::nullopt;
                        return SSPoint<T>{(p1 + p2) / T(2), q};
                    }
                }
                return std::nullopt;
            }

            /// The point where the step from `cur` to `next` (out of the boxes) reaches their boundary.
            [[nodiscard]] auto land(const SSPoint<T> &cur, const SSPoint<T> &next) const -> std::optional<SSPoint<T>>
            {
                std::optional<SSPoint<T>> best;
                T best_ratio = T(2);
                for (int k = 0; k < 4; ++k)
                {
                    if (!m_.out_of(k, next.q[k]))
                        continue;
                    const T b = next.q[k] < m_.lo(k) ? m_.lo(k) : m_.hi(k);
                    const T d = next.q[k] - cur.q[k];
                    const T ratio = d != T(0) ? std::clamp((b - cur.q[k]) / d, T(0), T(1)) : T(0);
                    if (ratio >= best_ratio)
                        continue;
                    std::array<T, 4> guess;
                    for (int j = 0; j < 4; ++j)
                        guess[j] = std::clamp(cur.q[j] + ratio * (next.q[j] - cur.q[j]), m_.lo(j), m_.hi(j));
                    if (auto r = correct(guess, nullptr, nullptr, {.fixed = k, .value = b}); r && m_.inside(r->q))
                    {
                        best = r;
                        best_ratio = ratio;
                    }
                }
                return best;
            }

            /**
             * Tracing from `start` in the direction `sign` of its tangency direction. The first point
             * is `start`. The end: `boundary` on a box, `closed` back at the start, `parted` where
             * the tangency ends and the surfaces part, `singular` where they cross beyond it (two
             * transversal curves start there) or the direction is undefined.
             */
            [[nodiscard]] auto march(const SSPoint<T> &start, T sign, typename Marcher<T>::End &end) const -> std::vector<SSPoint<T>>
            {
                using End = typename Marcher<T>::End;
                std::vector<SSPoint<T>> pts{start};
                auto dir = direction(start.q);
                if (!dir)
                {
                    end = End::singular;
                    return pts;
                }
                point<T, 3> t = sign * dir->t;
                const T h_max = m_.max_step(), h_min = T(10) * m_.tol();
                T h = h_max / T(4);
                for (std::size_t n = 0; n < 100000; ++n)
                {
                    const auto cur = pts.back();
                    const auto l = m_.local(cur.q);
                    const auto target = cur.p + h * t;
                    auto next = correct(m_.predict(cur.q, l, t, h), &target, &t);
                    bool accepted = false;
                    if (next && norm(next->p - cur.p) <= T(2) * h)
                    {
                        if (!m_.inside(next->q))
                        {
                            if (auto last = land(cur, *next); last && norm(last->p - cur.p) <= T(2) * h)
                            {
                                if (norm(last->p - cur.p) > m_.tol())
                                    pts.push_back(*last);
                                end = End::boundary;
                                return pts;
                            }
                        }
                        else if (auto dn = direction(next->q))
                        {
                            auto tn = dn->t * t < T(0) ? T(-1) * dn->t : dn->t;
                            const T turn = std::acos(std::clamp(tn * t, T(-1), T(1)));
                            // (still tangent along a curve: the surfaces bend alike along it; where
                            // they stop, the tangency ends, whatever the tolerance)
                            const bool alike = dn->small <= m_.angular_tol() * dn->large;
                            if (alike && turn <= m_.max_turn() && (next->p - cur.p) * t > T(0))
                            {
                                accepted = true;
                                // (the parameters must come back too: across a seam, the branch lands
                                // on it instead, with monotone parameters)
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
                                t = tn;
                                h = std::min(h_max, h * std::clamp(m_.max_turn() / std::max(turn, T(1e-12)) * T(0.8), T(0.5), T(2)));
                            }
                        }
                    }
                    if (!accepted)
                    {
                        h /= T(2);
                        if (h < h_min)
                        {
                            // the end of the tangency: the surfaces part, cross, or keep within the
                            // tolerance up to the boundary of a box, where the branch lands
                            const auto [beyond, off] = parting(pts.back(), t);
                            end = End::singular;
                            if (beyond == Beyond::parts)
                                end = End::parted;
                            else if (beyond == Beyond::boundary)
                                if (auto last = land(pts.back(), *off))
                                {
                                    if (norm(last->p - pts.back().p) > m_.tol())
                                        pts.push_back(*last);
                                    end = End::boundary;
                                }
                            return pts;
                        }
                    }
                }
                end = End::capped;
                return pts;
            }

            /**
             * What follows `cur`, the last point of a tangency, in the direction t? Across the curve
             * the gap (S1 - S2) . n2 grows as bend s^2 / 2; further on, its extremum across (the
             * foot with parallel normals, at any distance) keeps the sign of bend when the surfaces
             * part, and has the other sign when they cross (the gap then vanishes twice across:
             * two transversal curves). The probe goes further, from 10 tol doubling up to a step,
             * until the gap exceeds half the tolerance; leaving the boxes first, the surfaces keep
             * within the tolerance up to their boundary (the off iterate, to land on). Otherwise
             * (no solution, the gap too small up to a step): unknown.
             */
            enum class Beyond
            {
                unknown,
                parts,
                crosses,
                boundary
            };
            [[nodiscard]] auto parting(const SSPoint<T> &cur, const point<T, 3> &t) const -> std::pair<Beyond, std::optional<SSPoint<T>>>
            {
                const auto d = direction(cur.q);
                if (!d)
                    return {Beyond::unknown, std::nullopt};
                const auto l = m_.local(cur.q);
                for (T H = T(10) * m_.tol(); H <= m_.max_step(); H *= T(2))
                {
                    const auto target = cur.p + H * t;
                    const auto x = correct(m_.predict(cur.q, l, t, H), &target, &t, {.any_distance = true});
                    if (!x)
                        return {Beyond::unknown, std::nullopt};
                    if (!m_.inside(x->q))
                        return {Beyond::boundary, x};
                    const auto b = second(s2_, x->q[2], x->q[3]);
                    const auto n = cross(b.u, b.v);
                    const T nn = norm(n);
                    if (!(nn > T(0)))
                        return {Beyond::unknown, std::nullopt};
                    const T gap = (s1_.value(x->q[0], x->q[1]) - b.p) * n / nn;
                    if (std::abs(gap) > m_.tol() / T(2))
                        return {gap * d->bend > T(0) ? Beyond::parts : Beyond::crosses, std::nullopt};
                }
                return {Beyond::unknown, std::nullopt};
            }

            /// A point on the tangency between two consecutive points (for densification).
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
                if (r && !m_.inside(r->q))
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
        auto branch_curve(const Marcher<T> &m, const Surface<T, 3> &s1, const Surface<T, 3> &s2, std::vector<SSPoint<T>> pts, bool closed,
                          const std::vector<std::pair<SSPoint<T>, SSPoint<T>>> &jumps = {}, const TangentMarcher<T> *tangency = nullptr)
            -> IntersectionCurve<T>
        {
            // the steps over a branch point, which are not densified (no point between them)
            auto jump_at = [&](std::size_t k) {
                return std::ranges::any_of(jumps, [&](const auto &j) {
                    return (pts[k].p == j.first.p && pts[k + 1].p == j.second.p) || (pts[k].p == j.second.p && pts[k + 1].p == j.first.p);
                });
            };
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
                // pcurves kept in the parametric rectangles: a branch along a border is interpolated
                // with rounding errors beyond it
                auto pc1 = std::make_shared<BSCurve<T, 2>>(clamp_poles(interpolate<T, 2>(U1, t, p), s1.bounds()));
                auto pc2 = std::make_shared<BSCurve<T, 2>>(clamp_poles(interpolate<T, 2>(U2, t, p), s2.bounds()));
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
                        if (d > tol / T(2) && !jump_at(k))
                            refine[k] = true;
                    }
                ic.curve = c;
                ic.pc1 = pc1;
                ic.pc2 = pc2;
                ic.a = t.front();
                ic.b = t.back();
                ic.deviation = dev;
                ic.branch_spans.clear();
                for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                    if (jump_at(k))
                        ic.branch_spans.push_back({t[k], t[k + 1]});
                if (dev <= tol / T(2) || round == 7)
                    break;
                std::vector<SSPoint<T>> denser{pts.front()};
                for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                {
                    if (refine[k])
                        if (auto mid = tangency ? tangency->between(pts[k], pts[k + 1]) : m.between(pts[k], pts[k + 1]))
                            denser.push_back(*mid);
                    denser.push_back(pts[k + 1]);
                }
                if (denser.size() == pts.size())
                    break;
                pts = std::move(denser);
            }
            ic.closed = closed || norm(pts.front().p - pts.back().p) <= tol;
            ic.tangent = tangency != nullptr;
            return ic;
        }

        /// A branch of tangency through a start, both ways.
        template <std::floating_point T>
        struct TangentBranch
        {
            enum class Kind
            {
                not_a_curve, ///< no direction, an isolated contact: the regular path applies
                traced,      ///< the start, put on the tangency, is on a branch already traced
                branch
            } kind{Kind::not_a_curve};
            std::vector<SSPoint<T>> pts;
            bool closed{}, singular{}, capped{};
        };
        /**
         * The branch of tangency through `seed`. A curve of contact is certified by II1 - II2: its
         * smallest eigenvalue vanishes along the curve, both are of the same size at an isolated
         * contact (a sphere on a plane), which is not traced.
         */
        template <std::floating_point T, class Traced>
        auto tangent_branch(const TangentMarcher<T> &tm, const SSPoint<T> &seed, T small, T angular_tol, Traced &&traced) -> TangentBranch<T>
        {
            using End = typename Marcher<T>::End;
            using Kind = typename TangentBranch<T>::Kind;
            TangentBranch<T> b;
            const auto dir = tm.direction(seed.q);
            if (!dir || dir->small > angular_tol * dir->large)
                return b;
            // the seed is put on the tangency first (seeds near a tangency spread across it, by
            // the square root of the tolerance)
            const auto start = tm.correct(seed.q, &seed.p, &dir->t);
            if (!start || !tm.inside(start->q))
                return b;
            if (traced(start->p))
            {
                b.kind = Kind::traced;
                return b;
            }
            End e1{}, e2{};
            auto fwd = tm.march(*start, T(1), e1);
            if (e1 == End::closed)
                b.pts = std::move(fwd);
            else
            {
                auto bwd = tm.march(*start, T(-1), e2);
                b.pts.assign(bwd.rbegin(), bwd.rend());
                b.pts.insert(b.pts.end(), fwd.begin() + 1, fwd.end());
            }
            Box<T, 3> bb;
            for (const auto &x : b.pts)
                bb.add(x.p);
            if (b.pts.size() < 2 || bb.diagonal() <= small)
            {
                b.pts.clear();
                return b; // not a curve: an isolated contact
            }
            b.kind = Kind::branch;
            b.closed = e1 == End::closed;
            b.singular = e1 == End::singular || e2 == End::singular;
            b.capped = e1 == End::capped || e2 == End::capped;
            return b;
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
            const TangentMarcher<T> tm{m, s1, s2};
            const bool second = s1.has_second_derivatives() && s2.has_second_derivatives();

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

            // ---- marching from the seeds not on a traced branch, the best conditioned first (a seed
            //      near a branch point is then on a branch already traced through it)
            std::ranges::stable_sort(seeds, [&](const SSPoint<T> &a, const SSPoint<T> &b) { return m.sine(m.local(a.q)) > m.sine(m.local(b.q)); });
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
                // a start where the surfaces are tangent: on a curve of tangency, traced as such
                // (with second derivatives only)
                if (second && m.sine(m.local(seed.q)) < o.angular_tol)
                    if (auto tb = tangent_branch(tm, seed, small, o.angular_tol, traced); tb.kind != TangentBranch<T>::Kind::not_a_curve)
                    {
                        if (tb.kind == TangentBranch<T>::Kind::traced)
                            continue;
                        if (tb.singular)
                            res.singular = true;
                        if (tb.capped)
                            res.complete = false;
                        Box<T, 3> bb;
                        T longest{};
                        for (std::size_t k = 0; k < tb.pts.size(); ++k)
                        {
                            bb.add(tb.pts[k].p);
                            if (k > 0)
                                longest = std::max(longest, norm(tb.pts[k].p - tb.pts[k - 1].p));
                        }
                        bb.inflate(std::max(small, T(0.01) * longest));
                        branches.push_back(tb.pts);
                        boxes.push_back(bb);
                        res.curves.push_back(branch_curve(m, s1, s2, std::move(tb.pts), tb.closed, {}, &tm));
                        continue;
                    }
                std::vector<std::pair<SSPoint<T>, SSPoint<T>>> jumps;
                auto fwd = m.march(seed, T(1), e1, jumps);
                std::vector<SSPoint<T>> pts;
                bool closed = e1 == Marcher<T>::End::closed;
                if (closed)
                    pts = std::move(fwd);
                else
                {
                    auto bwd = m.march(seed, T(-1), e2, jumps);
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
                    // an isolated contact point where the surfaces are tangent; elsewhere, a branch
                    // reduced to a point is a start on a corner of the boxes (both seams at once)
                    const auto &x = pts.front();
                    const bool tangent = m.sine(m.local(x.q)) < o.angular_tol;
                    const bool known = std::ranges::any_of(res.points, [&](const auto &h) { return norm(h.p - x.p) <= small; });
                    if (tangent && !known)
                        res.points.push_back({x.q[0], x.q[1], x.q[2], x.q[3], x.p});
                    branches.push_back({x});
                    boxes.push_back(bb);
                    continue;
                }
                branches.push_back(pts);
                boxes.push_back(bb);
                res.curves.push_back(branch_curve(m, s1, s2, std::move(pts), closed, jumps));
            }
            return res;
        }
    }
}
