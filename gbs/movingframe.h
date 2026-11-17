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
 *   (a straight segment, an inflection): refused there.
 * - **Fixed binormal**: y the unit projection of a given direction on the normal plane, x = y × z;
 *   refused where the tangent is along that direction.
 *
 * The path may be made of pieces only G1 at given breaks (the junctions of the edges of a wire):
 * the stations include them.
 */
#include <gbs/bscurve.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
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

        /// Unit tangent: the first derivative, the second one where the first vanishes (a cusp of the parametrization).
        auto tangent(T u) const -> point<T, 3>
        {
            const auto d1 = c_->value(u, 1);
            if (norm(d1) > T(0))
                return unit(d1);
            return unit(c_->value(u, 2));
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
         * @throw std::invalid_argument an empty range, a null binormal for FixedBinormal
         * @throw std::domain_error Frenet where the curvature vanishes, FixedBinormal where the tangent
         * is along the binormal, a rotation-minimizing frame not reaching angle_tol
         */
        MovingFrame(std::shared_ptr<const Curve<T, 3>> c, T u1, T u2, FrameLaw law, T angle_tol = T(1e-8),
                    std::span<const T> breaks = {}, bool closed = false, const point<T, 3> &binormal = {})
            : c_{std::move(c)}, u1_{u1}, u2_{u2}, law_{law}
        {
            if (!c_ || !(u2 > u1))
                throw std::invalid_argument("MovingFrame: empty path");
            std::vector<T> b{u1};
            for (T x : breaks)
                if (x > u1 && x < u2)
                    b.push_back(x);
            b.push_back(u2);
            std::ranges::sort(b);
            const std::size_t pieces = b.size() - 1;
            // the checks of Frenet and of the fixed binormal, at points of every piece
            auto samples = [&](auto &&f) {
                for (T u : stations(b, 32))
                    f(u);
            };
            switch (law)
            {
            case FrameLaw::Frenet:
                samples([&](T u) {
                    const auto d1 = c_->value(u, 1), d2 = c_->value(u, 2);
                    if (!(norm(cross(d1, d2)) > T(1e-10) * norm(d1) * norm(d2)) || !(norm(d2) > T(0)))
                        throw std::domain_error("MovingFrame: Frenet frame where the curvature vanishes");
                });
                return;
            case FrameLaw::FixedBinormal:
                if (!(norm(binormal) > T(0)))
                    throw std::invalid_argument("MovingFrame: null binormal");
                binormal_ = unit(binormal);
                samples([&](T u) {
                    if (!(norm(cross(tangent(u), binormal_)) > T(1e-9)))
                        throw std::domain_error("MovingFrame: tangent along the fixed binormal");
                });
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
                    throw std::domain_error("MovingFrame: the rotation-minimizing frame does not reach angle_tol");
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
                const auto i = std::size_t(std::max<std::ptrdiff_t>(0, std::ranges::upper_bound(s_, u) - s_.begin() - 1));
                const auto j = std::min(i, s_.size() - 1);
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
