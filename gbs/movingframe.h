#pragma once
/**
 * @file movingframe.h
 * @brief Moving frames along a curve, for the sweeps: rotation-minimizing, Frenet, fixed binormal.
 *
 * A moving frame gives at each parameter u of a path c an orthonormal frame (o, x, y, z):
 * o = c(u), z the unit tangent, x and y in the normal plane. A profile carried by the frame (its
 * coordinates in the frame at the start kept along the path) sweeps a surface.
 *
 * - **Rotation-minimizing** (the default of the sweeps): the normal plane turns as little as it
 *   can, no rotation around the tangent. Computed by the double reflection of Wang, Jüttler, Zheng
 *   and Liu (*Computation of rotation minimizing frames*, ACM TOG 27(1), 2008): from a station to
 *   the next, a reflection in the bisector plane of the chord, then one that brings the reflected
 *   tangent onto the next tangent; fourth order. The stations are doubled until the frame at the
 *   common stations moves by less than `angle_tol`; between two stations, one double reflection
 *   from the previous one. On a closed path the frame does not come back to itself in general
 *   (the holonomy of the path): the residual angle is spread along the path, linearly in the
 *   length of the chords.
 * - **Frenet**: x the principal normal, y the binormal. Undefined where the curvature vanishes
 *   (a straight segment, an inflection): refused there, and where the binormal reverses between
 *   two samples or at a break (an inflection between the samples, two arcs turning opposite
 *   ways); refused on a curve without second derivatives.
 * - **Fixed binormal**: y the unit projection of a given direction on the normal plane, x = y × z;
 *   refused where the tangent is along that direction, also between the samples (y reverses).
 *
 * The continuity is checked on 32 samples per piece, an interval halved while the binormal (y for
 * the fixed binormal) turns across it by more than pi / 8: a reversal left after 40 halvings is
 * refused.
 *
 * The path may be made of pieces only G1 at given breaks (the junctions of the edges of a wire):
 * the stations include them. At a break the frame is that of the piece after it; the frame at the
 * end of the piece before is the limit from the left, at `std::nextafter(u, u1)`.
 */
#include <gbs/bscurve.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    enum class FrameLaw : std::uint8_t
    {
        RotationMinimizing,
        Frenet,
        FixedBinormal
    };

    template <std::floating_point T>
    class MovingFrame
    {
        std::shared_ptr<const Curve<T, 3>> c_;
        T u1_, u2_;
        FrameLaw law_;
        point<T, 3> binormal_{};
        // rotation-minimizing: stations, their points, tangents and normals; chord lengths for the closure
        std::vector<T> s_;
        std::vector<point<T, 3>> x_, t_, r_;
        std::vector<T> l_;
        T phi_{};

        static auto unit(const point<T, 3> &v) -> point<T, 3> { return v / norm(v); }

        /// Unit tangent: the first derivative; where it vanishes (a cusp of the parametrization), the
        /// second one, or the chord around u for a curve without second derivatives.
        auto tangent(T u) const -> point<T, 3>
        {
            const auto d1 = c_->value(u, 1);
            if (norm(d1) > T(0))
                return unit(d1);
            if (c_->has_second_derivatives())
                return unit(c_->value(u, 2));
            const T h = T(1e-6) * (u2_ - u1_);
            return unit(c_->value(std::min(u + h, u2_)) - c_->value(std::max(u - h, u1_)));
        }

        /// The vector a law turns with the path, which must not reverse: Frenet's binormal, the fixed
        /// binormal projected on the normal plane.
        auto steered(T u) const -> point<T, 3>
        {
            if (law_ == FrameLaw::Frenet)
                return unit(cross(c_->value(u, 1), c_->value(u, 2)));
            const auto t = tangent(u);
            return unit(binormal_ - (binormal_ * t) * t);
        }

        /// The steered vector turns continuously from a to b: the interval halved while it turns by more
        /// than pi / 8 across it, a reversal (or NaN) left after `depth` halvings refused.
        auto steady(T a, const point<T, 3> &va, T b, const point<T, 3> &vb, int depth) const -> bool
        {
            if (va * vb >= std::cos(std::numbers::pi_v<T> / 8))
                return true;
            if (depth == 0)
                return false;
            const T mid = (a + b) / 2;
            const auto vm = steered(mid);
            return steady(a, va, mid, vm, depth - 1) && steady(mid, vm, b, vb, depth - 1);
        }

        /// Frenet and the fixed binormal: the steered vector does not reverse along each piece nor at a break.
        void check_steady(std::span<const T> b, const char *what) const
        {
            point<T, 3> before{};
            for (std::size_t k = 0; k + 1 < b.size(); ++k)
            {
                // the piece from b[k] to the limit from the left at b[k + 1]
                const T end = std::nextafter(b[k + 1], b[k]);
                T a = b[k];
                auto va = steered(a);
                if (k > 0 && !(before * va > T(0)))
                    throw std::domain_error(std::string("MovingFrame: ") + what + " reversed at a break of the path");
                for (int i = 1; i <= 32; ++i)
                {
                    const T u = i == 32 ? end : b[k] + (b[k + 1] - b[k]) * T(i) / T(32);
                    const auto vu = steered(u);
                    if (!steady(a, va, u, vu, 40))
                        throw std::domain_error(std::string("MovingFrame: ") + what + " reversed along the path");
                    a = u;
                    va = vu;
                }
                before = va;
            }
        }

        /// Next normal by double reflection, from (x, t, r) to the point y of tangent ty.
        static auto reflect(const point<T, 3> &x, const point<T, 3> &t, const point<T, 3> &r, const point<T, 3> &y,
                            const point<T, 3> &ty) -> point<T, 3>
        {
            const auto v1 = y - x;
            const T c1 = v1 * v1;
            if (!(c1 > T(0)))
                return r;
            const auto rl = r - (T(2) / c1) * (v1 * r) * v1;
            const auto tl = t - (T(2) / c1) * (v1 * t) * v1;
            const auto v2 = ty - tl;
            const T c2 = v2 * v2;
            auto rn = c2 > T(0) ? rl - (T(2) / c2) * (v2 * rl) * v2 : rl;
            // back onto the normal plane, against the rounding
            return unit(rn - (rn * ty) * ty);
        }

        /// Normals at the stations, the reference r0 normal to the tangent at u1.
        auto march(const std::vector<T> &s, const point<T, 3> &r0, std::vector<point<T, 3>> &x, std::vector<point<T, 3>> &t) const
            -> std::vector<point<T, 3>>
        {
            x.resize(s.size());
            t.resize(s.size());
            std::vector<point<T, 3>> r(s.size());
            for (std::size_t i = 0; i < s.size(); ++i)
            {
                x[i] = c_->value(s[i]);
                t[i] = tangent(s[i]);
                r[i] = i == 0 ? r0 : reflect(x[i - 1], t[i - 1], r[i - 1], x[i], t[i]);
            }
            return r;
        }

        static auto stations(std::span<const T> breaks, std::size_t n) -> std::vector<T>
        {
            std::vector<T> s;
            for (std::size_t k = 0; k + 1 < breaks.size(); ++k)
                for (std::size_t i = 0; i < n; ++i)
                    s.push_back(breaks[k] + (breaks[k + 1] - breaks[k]) * T(i) / T(n));
            s.push_back(breaks.back());
            return s;
        }

        /// Angle of the rotation around t bringing a onto b (both normal to t).
        static auto angle(const point<T, 3> &a, const point<T, 3> &b, const point<T, 3> &t) -> T
        {
            return std::atan2(cross(a, b) * t, a * b);
        }

    public:
        /**
         * @param c the path, on [u1, u2]
         * @param law the frame law
         * @param angle_tol accuracy of the rotation-minimizing frame (radians)
         * @param breaks parameters where the path is only G1 (inside (u1, u2)), included in the stations
         * @param closed the path comes back to its start (rotation-minimizing: the residual angle spread)
         * @param binormal the direction kept by FixedBinormal
         * @throw std::invalid_argument an empty range, angle_tol not positive, a null binormal for FixedBinormal
         * @throw std::domain_error Frenet where the curvature vanishes, where the binormal reverses or
         * without second derivatives; FixedBinormal where the tangent is along the binormal
         * @throw std::runtime_error a rotation-minimizing frame not reaching angle_tol (within 65536
         * stations)
         */
        MovingFrame(std::shared_ptr<const Curve<T, 3>> c, T u1, T u2, FrameLaw law, T angle_tol = T(1e-8),
                    std::span<const T> breaks = {}, bool closed = false, const point<T, 3> &binormal = {})
            : c_{std::move(c)}, u1_{u1}, u2_{u2}, law_{law}
        {
            if (!c_ || !(u2 > u1))
                throw std::invalid_argument("MovingFrame: empty path");
            if (!(angle_tol > T(0)) || !std::isfinite(angle_tol))
                throw std::invalid_argument("MovingFrame: angle_tol must be positive");
            std::vector<T> b{u1};
            for (T x : breaks)
                if (x > u1 && x < u2)
                    b.push_back(x);
            b.push_back(u2);
            std::ranges::sort(b);
            const auto dup = std::ranges::unique(b);
            b.erase(dup.begin(), dup.end());
            const std::size_t pieces = b.size() - 1;
            // the checks of Frenet and of the fixed binormal, at points of every piece
            auto samples = [&](auto &&f) {
                for (T u : stations(b, 32))
                    f(u);
            };
            switch (law)
            {
            case FrameLaw::Frenet:
                if (!c_->has_second_derivatives())
                    throw std::domain_error("MovingFrame: Frenet frame of a curve without second derivatives");
                samples([&](T u) {
                    const auto d1 = c_->value(u, 1), d2 = c_->value(u, 2);
                    if (!(norm(cross(d1, d2)) > T(1e-10) * norm(d1) * norm(d2)) || !(norm(d2) > T(0)))
                        throw std::domain_error("MovingFrame: Frenet frame where the curvature vanishes");
                });
                check_steady(b, "Frenet binormal");
                return;
            case FrameLaw::FixedBinormal:
                if (!(norm(binormal) > T(0)))
                    throw std::invalid_argument("MovingFrame: null binormal");
                binormal_ = unit(binormal);
                samples([&](T u) {
                    if (!(norm(cross(tangent(u), binormal_)) > T(1e-9)))
                        throw std::domain_error("MovingFrame: tangent along the fixed binormal");
                });
                check_steady(b, "fixed binormal projected on the normal plane");
                return;
            case FrameLaw::RotationMinimizing:
                break;
            }
            // reference normal: the axis least along the tangent, projected
            const auto t0 = tangent(u1);
            point<T, 3> a{};
            std::size_t k = 0;
            for (std::size_t i = 1; i < 3; ++i)
                if (std::abs(t0[i]) < std::abs(t0[k]))
                    k = i;
            a[k] = T(1);
            const auto r0 = unit(a - (a * t0) * t0);
            std::vector<point<T, 3>> xs, ts;
            auto r = march(stations(b, 8), r0, xs, ts);
            for (std::size_t n = 16;; n *= 2)
            {
                auto s = stations(b, n);
                auto rf = march(s, r0, x_, t_);
                T gap{};
                for (std::size_t i = 0; i < r.size(); ++i)
                    gap = std::max(gap, norm(rf[2 * i] - r[i]));
                r = std::move(rf);
                s_ = std::move(s);
                if (gap <= angle_tol)
                    break;
                if (n * pieces > 65536)
                    throw std::runtime_error("MovingFrame: the rotation-minimizing frame does not reach angle_tol");
            }
            r_ = std::move(r);
            l_.assign(s_.size(), T(0));
            for (std::size_t i = 1; i < s_.size(); ++i)
                l_[i] = l_[i - 1] + norm(x_[i] - x_[i - 1]);
            if (closed)
                phi_ = angle(r_.front(), r_.back() - (r_.back() * t_.front()) * t_.front(), t_.front());
        }

        [[nodiscard]] auto law() const -> FrameLaw { return law_; }
        [[nodiscard]] auto bounds() const -> std::array<T, 2> { return {u1_, u2_}; }
        /// Residual angle of the rotation-minimizing frame around a closed path, spread along it.
        [[nodiscard]] auto closure_angle() const -> T { return phi_; }

        /// The frame at u: o on the path, z the unit tangent, (x, y, z) direct.
        [[nodiscard]] auto operator()(T u) const -> Frame3<T>
        {
            u = std::clamp(u, u1_, u2_);
            Frame3<T> f;
            f.o = c_->value(u);
            f.z = tangent(u);
            switch (law_)
            {
            case FrameLaw::Frenet:
            {
                const auto d2 = c_->value(u, 2);
                f.x = unit(d2 - (d2 * f.z) * f.z);
                break;
            }
            case FrameLaw::FixedBinormal:
                f.y = unit(binormal_ - (binormal_ * f.z) * f.z);
                f.x = cross(f.y, f.z);
                return f;
            case FrameLaw::RotationMinimizing:
            {
                // the station at or before u (upper_bound - 1 is at most the last one)
                const auto j = std::size_t(std::max<std::ptrdiff_t>(0, std::ranges::upper_bound(s_, u) - s_.begin() - 1));
                auto r = u == s_[j] ? r_[j] : reflect(x_[j], t_[j], r_[j], f.o, f.z);
                if (phi_ != T(0))
                {
                    const T w = (l_[j] + (u == s_[j] ? T(0) : norm(f.o - x_[j]))) / l_.back();
                    const T a = -phi_ * w, ca = std::cos(a), sa = std::sin(a);
                    r = ca * r + sa * cross(f.z, r);
                }
                f.x = r;
                break;
            }
            }
            f.y = cross(f.z, f.x);
            return f;
        }
    };
}
