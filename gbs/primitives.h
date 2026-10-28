#pragma once
/**
 * @file primitives.h
 * @brief Analytic primitives: line, ellipse (and circle), plane, cylinder, cone, sphere, torus.
 *
 * Design: docs/sources/design/step_writer.md § 12.1; architecture note
 * docs/sources/design/step_p1_primitives.md.
 *
 * The parametrizations are those of ISO 10303-42 (STEP), in a frame
 * ax2 = {origin, axis, reference direction} made orthonormal (Y = Z x X):
 *
 *   line      C(t) = O + t V
 *   ellipse   C(t) = O + r1 cos t X + r2 sin t Y                  (circle: r1 = r2)
 *   plane     S(u, v) = O + u X + v Y
 *   cylinder  S(u, v) = O + R (cos u X + sin u Y) + v Z
 *   cone      S(u, v) = O + (R + v tan a)(cos u X + sin u Y) + v Z
 *   sphere    S(u, v) = O + R cos v (cos u X + sin u Y) + R sin v Z
 *   torus     S(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
 *
 * Every primitive has a finite parameter domain given at construction (the
 * range of an edge, the box of a face), as the BREP requires; angles are not
 * restricted to [0, 2 pi]: an arc from 5.9 to 6.6 rad crosses 0. Derivatives
 * of any order are exact. parameters() inverts a point in closed form and
 * to_nurbs() gives the exact NURBS over the domain (gbs/bselementary.h), whose
 * parameters are not the analytic ones (nurbs_parameters() converts).
 */
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bselementary.h>

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <stdexcept>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include "vecop.ixx"
#endif

namespace gbs
{
    namespace primitives_detail
    {
        template <std::floating_point T>
        constexpr T two_pi = T(2) * std::numbers::pi_v<T>;

        /// n-th derivative of cos and sin.
        template <std::floating_point T>
        T dcos(T x, std::size_t n) { return std::cos(x + T(n) * T(0.5) * std::numbers::pi_v<T>); }
        template <std::floating_point T>
        T dsin(T x, std::size_t n) { return std::sin(x + T(n) * T(0.5) * std::numbers::pi_v<T>); }

        /// t + 2 pi k in [lo, lo + 2 pi)
        template <std::floating_point T>
        T wrap(T t, T lo) { return t - two_pi<T> * std::floor((t - lo) / two_pi<T>); }
    }

    template <std::floating_point T>
    struct Frame3;

    namespace primitives_detail
    {
        /// Frame of the meridian plane at angle u, centered at c: x = rho(u), y = Z.
        template <std::floating_point T>
        Frame3<T> meridian(const Frame3<T> &f, T u, const point<T, 3> &c);
    }

    /// Orthonormal frame of an ax2 {origin, axis, reference direction}.
    template <std::floating_point T>
    struct Frame3
    {
        point<T, 3> o, x, y, z;

        Frame3() : o{}, x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1} {}
        explicit Frame3(const ax2<T, 3> &ax)
        {
            const T nz = norm(ax[1]);
            if (!(nz > T(0)))
                throw std::invalid_argument("frame: null axis");
            z = ax[1] / nz;
            x = ax[2] - (ax[2] * z) * z;
            const T nx = norm(x);
            if (!(nx > T(1e-12) * norm(ax[2])))
                throw std::invalid_argument("frame: reference direction parallel to the axis");
            x = x / nx;
            y = cross(z, x);
            o = ax[0];
        }
        [[nodiscard]] ax2<T, 3> ax() const { return {o, z, x}; }
        /// Frame moved by a rigid map of points.
        [[nodiscard]] Frame3 moved(const std::function<point<T, 3>(const point<T, 3> &)> &f) const
        {
            Frame3 r;
            r.o = f(o);
            r.x = f(o + x) - r.o;
            r.y = f(o + y) - r.o;
            r.z = f(o + z) - r.o;
            return r;
        }
        [[nodiscard]] point<T, 3> radial(T u, std::size_t d = 0) const
        {
            return primitives_detail::dcos(u, d) * x + primitives_detail::dsin(u, d) * y;
        }
    };

    template <std::floating_point T>
    Frame3<T> primitives_detail::meridian(const Frame3<T> &f, T u, const point<T, 3> &c)
    {
        Frame3<T> m;
        m.o = c;
        m.x = f.radial(u);
        m.y = f.z;
        m.z = cross(m.x, m.y);
        return m;
    }

    // =========================================================================
    // Curves
    // =========================================================================

    /// Analytic curve with a finite parameter range.
    template <std::floating_point T>
    class ElementaryCurve : public Curve<T, 3>
    {
    protected:
        std::array<T, 2> range_;

    public:
        explicit ElementaryCurve(std::array<T, 2> range) : range_{range}
        {
            if (!(range[1] > range[0]) || !std::isfinite(range[0]) || !std::isfinite(range[1]))
                throw std::invalid_argument("primitive curve: invalid parameter range");
        }
        auto bounds() const -> std::array<T, 2> override { return range_; }
        /// Parameter of the point of the curve nearest to p, within a turn of the range for a conic.
        [[nodiscard]] virtual T parameter(const point<T, 3> &p) const = 0;
        /// Exact NURBS over the range.
        [[nodiscard]] virtual std::shared_ptr<Curve<T, 3>> to_nurbs() const = 0;
        /// NURBS parameter of the point at parameter t.
        [[nodiscard]] virtual T nurbs_parameter(T t) const = 0;
        /// The same curve moved by a rigid map of points.
        [[nodiscard]] virtual std::shared_ptr<ElementaryCurve<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &f) const = 0;
        /// The same point set run backwards: parameter t becomes -t.
        [[nodiscard]] virtual std::shared_ptr<ElementaryCurve<T>> reversed() const = 0;
    };

    /// C(t) = O + t V (STEP line), t in [t1, t2].
    template <std::floating_point T>
    class LineCurve : public ElementaryCurve<T>
    {
        point<T, 3> o_, v_;

    public:
        LineCurve(const point<T, 3> &origin, const point<T, 3> &vector, T t1, T t2) : ElementaryCurve<T>{{t1, t2}}, o_{origin}, v_{vector}
        {
            if (!(norm(vector) > T(0)))
                throw std::invalid_argument("line: null vector");
        }
        [[nodiscard]] const point<T, 3> &origin() const { return o_; }
        [[nodiscard]] const point<T, 3> &vector() const { return v_; }
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            if (d == 0)
                return o_ + t * v_;
            if (d == 1)
                return v_;
            return {};
        }
        T parameter(const point<T, 3> &p) const override { return ((p - o_) * v_) / (v_ * v_); }
        std::shared_ptr<Curve<T, 3>> to_nurbs() const override
        {
            const auto [a, b] = this->range_;
            return std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{value(a), value(b)}, std::vector<T>{a, a, b, b}, 1);
        }
        T nurbs_parameter(T t) const override { return t; }
        std::shared_ptr<ElementaryCurve<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &f) const override
        {
            const auto o = f(o_);
            return std::make_shared<LineCurve<T>>(o, f(o_ + v_) - o, this->range_[0], this->range_[1]);
        }
        std::shared_ptr<ElementaryCurve<T>> reversed() const override
        {
            return std::make_shared<LineCurve<T>>(o_, T(-1) * v_, -this->range_[1], -this->range_[0]);
        }
    };

    /// C(t) = O + r1 cos t X + r2 sin t Y (STEP ellipse; a circle when r1 = r2), t in [t1, t2], t2 - t1 <= 2 pi.
    template <std::floating_point T>
    class EllipseCurve : public ElementaryCurve<T>
    {
        Frame3<T> f_;
        T r1_, r2_;

    public:
        EllipseCurve(T r1, T r2, const Frame3<T> &frame, T t1 = T(0), T t2 = primitives_detail::two_pi<T>)
            : ElementaryCurve<T>{{t1, t2}}, f_{frame}, r1_{r1}, r2_{r2}
        {
            if (!(r1 > T(0)) || !(r2 > T(0)))
                throw std::invalid_argument("ellipse: radii must be positive");
            if (t2 - t1 > primitives_detail::two_pi<T> * (T(1) + T(1e-12)))
                throw std::invalid_argument("ellipse: range larger than a turn");
        }
        EllipseCurve(T r1, T r2, const ax2<T, 3> &ax, T t1 = T(0), T t2 = primitives_detail::two_pi<T>)
            : EllipseCurve{r1, r2, Frame3<T>{ax}, t1, t2} {}
        [[nodiscard]] const Frame3<T> &frame() const { return f_; }
        [[nodiscard]] T radius1() const { return r1_; }
        [[nodiscard]] T radius2() const { return r2_; }
        [[nodiscard]] bool is_circle() const { return r1_ == r2_; }
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            auto p = (r1_ * primitives_detail::dcos(t, d)) * f_.x + (r2_ * primitives_detail::dsin(t, d)) * f_.y;
            return d == 0 ? f_.o + p : p;
        }
        T parameter(const point<T, 3> &p) const override
        {
            const auto q = p - f_.o;
            const T t = std::atan2((q * f_.y) / r2_, (q * f_.x) / r1_);
            // the turn closest to the middle of the range
            return primitives_detail::wrap(t, T(0.5) * (this->range_[0] + this->range_[1]) - std::numbers::pi_v<T>);
        }
        std::shared_ptr<Curve<T, 3>> to_nurbs() const override
        {
            return std::make_shared<BSCurveRational<T, 3>>(
                build_ellipse_arc<T, 3>(r1_, r2_, this->range_[0], this->range_[1], f_.o, f_.x, f_.y));
        }
        T nurbs_parameter(T t) const override { return ellipse_arc_parameter(this->range_[0], this->range_[1], t); }
        std::shared_ptr<ElementaryCurve<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &f) const override
        {
            return std::make_shared<EllipseCurve<T>>(r1_, r2_, f_.moved(f), this->range_[0], this->range_[1]);
        }
        std::shared_ptr<ElementaryCurve<T>> reversed() const override
        {
            auto g = f_; // C(-t) with Y -> -Y
            g.y = T(-1) * g.y;
            g.z = T(-1) * g.z;
            return std::make_shared<EllipseCurve<T>>(r1_, r2_, g, -this->range_[1], -this->range_[0]);
        }
    };

    // =========================================================================
    // Surfaces
    // =========================================================================

    /// Analytic surface on a finite box of parameters.
    template <std::floating_point T>
    class ElementarySurface : public Surface<T, 3>
    {
    protected:
        Frame3<T> f_;
        std::array<T, 4> box_;

        /// u of the angle a in the turn of the box.
        T turn_u(T a) const { return primitives_detail::wrap(a, T(0.5) * (box_[0] + box_[1]) - std::numbers::pi_v<T>); }

    public:
        ElementarySurface(const Frame3<T> &f, std::array<T, 4> box) : f_{f}, box_{box}
        {
            for (T x : box)
                if (!std::isfinite(x))
                    throw std::invalid_argument("primitive surface: infinite parameter box");
            if (!(box[1] > box[0]) || !(box[3] > box[2]))
                throw std::invalid_argument("primitive surface: empty parameter box");
        }
        auto bounds() const -> std::array<T, 4> override { return box_; }
        [[nodiscard]] const Frame3<T> &frame() const { return f_; }
        /// Parameters of a point of the surface (closed form; u in the turn of the box for an angle).
        [[nodiscard]] virtual std::array<T, 2> parameters(const point<T, 3> &p) const = 0;
        /// Exact NURBS over the box.
        [[nodiscard]] virtual std::shared_ptr<Surface<T, 3>> to_nurbs() const = 0;
        /// NURBS parameters of the point at (u, v).
        [[nodiscard]] virtual std::array<T, 2> nurbs_parameters(T u, T v) const = 0;
        /// The same surface moved by a rigid map of points.
        [[nodiscard]] virtual std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &f) const = 0;
        /// The iso-u line between v1 and v2, an analytic curve parametrized by v.
        [[nodiscard]] virtual std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const = 0;
        /// The iso-v line between u1 and u2, an analytic curve parametrized by u (nullptr if it reduces to a point).
        [[nodiscard]] virtual std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const = 0;

    protected:
        /// Circle of radius r at height h on the axis, parametrized by the angle u like the surface.
        std::shared_ptr<ElementaryCurve<T>> parallel(T r, T h, T u1, T u2) const
        {
            if (!(r > T(1e-14)))
                return nullptr;
            Frame3<T> c = f_;
            c.o = f_.o + h * f_.z;
            return std::make_shared<EllipseCurve<T>>(r, r, c, u1, u2);
        }
    };

    /// S(u, v) = O + u X + v Y
    template <std::floating_point T>
    class Plane : public ElementarySurface<T>
    {
    public:
        Plane(const Frame3<T> &f, std::array<T, 4> box) : ElementarySurface<T>{f, box} {}
        Plane(const ax2<T, 3> &ax, std::array<T, 4> box) : Plane{Frame3<T>{ax}, box} {}
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const auto &f = this->f_;
            if (du == 0 && dv == 0)
                return f.o + u * f.x + v * f.y;
            if (du == 1 && dv == 0)
                return f.x;
            if (du == 0 && dv == 1)
                return f.y;
            return {};
        }
        std::array<T, 2> parameters(const point<T, 3> &p) const override
        {
            const auto q = p - this->f_.o;
            return {q * this->f_.x, q * this->f_.y};
        }
        std::shared_ptr<Surface<T, 3>> to_nurbs() const override
        {
            const auto [u1, u2, v1, v2] = this->box_;
            return std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{value(u1, v1), value(u2, v1), value(u1, v2), value(u2, v2)},
                                                     std::vector<T>{u1, u1, u2, u2}, std::vector<T>{v1, v1, v2, v2}, 1, 1);
        }
        std::array<T, 2> nurbs_parameters(T u, T v) const override { return {u, v}; }
        std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &m) const override
        {
            return std::make_shared<Plane<T>>(this->f_.moved(m), this->box_);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const override
        {
            return std::make_shared<LineCurve<T>>(this->f_.o + u * this->f_.x, this->f_.y, v1, v2);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const override
        {
            return std::make_shared<LineCurve<T>>(this->f_.o + v * this->f_.y, this->f_.x, u1, u2);
        }
    };

    /// S(u, v) = O + R (cos u X + sin u Y) + v Z
    template <std::floating_point T>
    class CylindricalSurface : public ElementarySurface<T>
    {
        T r_;

    public:
        CylindricalSurface(T radius, const Frame3<T> &f, std::array<T, 4> box) : ElementarySurface<T>{f, box}, r_{radius}
        {
            if (!(radius > T(0)))
                throw std::invalid_argument("cylinder: radius must be positive");
        }
        CylindricalSurface(T radius, const ax2<T, 3> &ax, std::array<T, 4> box) : CylindricalSurface{radius, Frame3<T>{ax}, box} {}
        [[nodiscard]] T radius() const { return r_; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const auto &f = this->f_;
            if (dv == 0)
                return du == 0 ? f.o + r_ * f.radial(u) + v * f.z : r_ * f.radial(u, du);
            if (dv == 1 && du == 0)
                return f.z;
            return {};
        }
        std::array<T, 2> parameters(const point<T, 3> &p) const override
        {
            const auto q = p - this->f_.o;
            return {this->turn_u(std::atan2(q * this->f_.y, q * this->f_.x)), q * this->f_.z};
        }
        std::shared_ptr<Surface<T, 3>> to_nurbs() const override
        {
            const auto [u1, u2, v1, v2] = this->box_;
            return std::make_shared<BSSurfaceRational<T, 3>>(build_cylinder<T>(r_, this->f_.ax(), v1, v2, u1, u2));
        }
        std::array<T, 2> nurbs_parameters(T u, T v) const override { return {ellipse_arc_parameter(this->box_[0], this->box_[1], u), v}; }
        std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &m) const override
        {
            return std::make_shared<CylindricalSurface<T>>(r_, this->f_.moved(m), this->box_);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const override
        {
            return std::make_shared<LineCurve<T>>(this->f_.o + r_ * this->f_.radial(u), this->f_.z, v1, v2);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const override { return this->parallel(r_, v, u1, u2); }
    };

    /// S(u, v) = O + (R + v tan a)(cos u X + sin u Y) + v Z
    template <std::floating_point T>
    class ConicalSurface : public ElementarySurface<T>
    {
        T r_, a_;

    public:
        ConicalSurface(T radius, T semi_angle, const Frame3<T> &f, std::array<T, 4> box)
            : ElementarySurface<T>{f, box}, r_{radius}, a_{semi_angle}
        {
            if (r_ < T(0) || !(std::abs(a_) < T(0.5) * std::numbers::pi_v<T>) || a_ == T(0))
                throw std::invalid_argument("cone: invalid radius or semi angle");
        }
        ConicalSurface(T radius, T semi_angle, const ax2<T, 3> &ax, std::array<T, 4> box)
            : ConicalSurface{radius, semi_angle, Frame3<T>{ax}, box} {}
        [[nodiscard]] T radius() const { return r_; }
        [[nodiscard]] T semi_angle() const { return a_; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const auto &f = this->f_;
            const T t = std::tan(a_);
            point<T, 3> p{};
            if (dv == 0)
                p = (r_ + v * t) * f.radial(u, du);
            else if (dv == 1)
                p = t * f.radial(u, du);
            if (du == 0 && dv == 0)
                p = p + f.o + v * f.z;
            else if (du == 0 && dv == 1)
                p = p + f.z;
            return p;
        }
        std::array<T, 2> parameters(const point<T, 3> &p) const override
        {
            const auto q = p - this->f_.o;
            const T v = q * this->f_.z, u = std::atan2(q * this->f_.y, q * this->f_.x);
            // beyond the apex (R + v tan a < 0) the point is on the other nappe, at u + pi
            return {this->turn_u(r_ + v * std::tan(a_) < T(0) ? u + std::numbers::pi_v<T> : u), v};
        }
        std::shared_ptr<Surface<T, 3>> to_nurbs() const override
        {
            const auto [u1, u2, v1, v2] = this->box_;
            return std::make_shared<BSSurfaceRational<T, 3>>(build_cone<T>(r_, a_, this->f_.ax(), v1, v2, u1, u2));
        }
        std::array<T, 2> nurbs_parameters(T u, T v) const override { return {ellipse_arc_parameter(this->box_[0], this->box_[1], u), v}; }
        std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &m) const override
        {
            return std::make_shared<ConicalSurface<T>>(r_, a_, this->f_.moved(m), this->box_);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const override
        {
            const auto rho = this->f_.radial(u);
            return std::make_shared<LineCurve<T>>(this->f_.o + r_ * rho, std::tan(a_) * rho + this->f_.z, v1, v2);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const override { return this->parallel(r_ + v * std::tan(a_), v, u1, u2); }
    };

    /// S(u, v) = O + R cos v (cos u X + sin u Y) + R sin v Z, v in [-pi/2, pi/2]
    template <std::floating_point T>
    class SphericalSurface : public ElementarySurface<T>
    {
        T r_;

    public:
        SphericalSurface(T radius, const Frame3<T> &f, std::array<T, 4> box) : ElementarySurface<T>{f, box}, r_{radius}
        {
            if (!(radius > T(0)))
                throw std::invalid_argument("sphere: radius must be positive");
            const T h = T(0.5) * std::numbers::pi_v<T> * (T(1) + T(1e-12));
            if (box[2] < -h || box[3] > h)
                throw std::invalid_argument("sphere: latitude out of [-pi/2, pi/2]");
        }
        SphericalSurface(T radius, const ax2<T, 3> &ax, std::array<T, 4> box) : SphericalSurface{radius, Frame3<T>{ax}, box} {}
        [[nodiscard]] T radius() const { return r_; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            using namespace primitives_detail;
            const auto &f = this->f_;
            point<T, 3> p = (r_ * dcos(v, dv)) * f.radial(u, du);
            if (du == 0)
                p = p + (r_ * dsin(v, dv)) * f.z;
            return du == 0 && dv == 0 ? f.o + p : p;
        }
        std::array<T, 2> parameters(const point<T, 3> &p) const override
        {
            const auto q = p - this->f_.o;
            const T h = q * this->f_.z;
            const T rho = norm(q - h * this->f_.z);
            return {this->turn_u(std::atan2(q * this->f_.y, q * this->f_.x)), std::atan2(h, rho)};
        }
        std::shared_ptr<Surface<T, 3>> to_nurbs() const override
        {
            const auto [u1, u2, v1, v2] = this->box_;
            return std::make_shared<BSSurfaceRational<T, 3>>(build_sphere<T>(r_, this->f_.ax(), u1, u2, v1, v2));
        }
        std::array<T, 2> nurbs_parameters(T u, T v) const override
        {
            return {ellipse_arc_parameter(this->box_[0], this->box_[1], u), ellipse_arc_parameter(this->box_[2], this->box_[3], v)};
        }
        std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &m) const override
        {
            return std::make_shared<SphericalSurface<T>>(r_, this->f_.moved(m), this->box_);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const override
        {
            // meridian: O + R cos v rho(u) + R sin v Z
            return std::make_shared<EllipseCurve<T>>(r_, r_, primitives_detail::meridian(this->f_, u, this->f_.o), v1, v2);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const override
        {
            return this->parallel(r_ * std::cos(v), r_ * std::sin(v), u1, u2);
        }
    };

    /// S(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
    template <std::floating_point T>
    class ToroidalSurface : public ElementarySurface<T>
    {
        T major_, minor_;

    public:
        ToroidalSurface(T major, T minor, const Frame3<T> &f, std::array<T, 4> box)
            : ElementarySurface<T>{f, box}, major_{major}, minor_{minor}
        {
            if (!(minor > T(0)))
                throw std::invalid_argument("torus: minor radius must be positive");
        }
        ToroidalSurface(T major, T minor, const ax2<T, 3> &ax, std::array<T, 4> box)
            : ToroidalSurface{major, minor, Frame3<T>{ax}, box} {}
        [[nodiscard]] T major_radius() const { return major_; }
        [[nodiscard]] T minor_radius() const { return minor_; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            using namespace primitives_detail;
            const auto &f = this->f_;
            const T radial = dv == 0 ? major_ + minor_ * std::cos(v) : minor_ * dcos(v, dv);
            point<T, 3> p = radial * f.radial(u, du);
            if (du == 0)
                p = p + (minor_ * dsin(v, dv)) * f.z;
            return du == 0 && dv == 0 ? f.o + p : p;
        }
        std::array<T, 2> parameters(const point<T, 3> &p) const override
        {
            const auto q = p - this->f_.o;
            const T h = q * this->f_.z;
            const T rho = norm(q - h * this->f_.z);
            T u = std::atan2(q * this->f_.y, q * this->f_.x), v = std::atan2(h, rho - major_);
            // a spindle torus (R < r) also passes across the axis, where R + r cos v < 0: the point
            // is there at u + pi when the tube circle centred at -R is the nearer
            if (major_ < minor_ && std::abs(std::hypot(rho + major_, h) - minor_) < std::abs(std::hypot(rho - major_, h) - minor_))
            {
                u += std::numbers::pi_v<T>;
                v = std::atan2(h, -rho - major_);
            }
            return {this->turn_u(u), primitives_detail::wrap(v, T(0.5) * (this->box_[2] + this->box_[3]) - std::numbers::pi_v<T>)};
        }
        std::shared_ptr<Surface<T, 3>> to_nurbs() const override
        {
            const auto [u1, u2, v1, v2] = this->box_;
            return std::make_shared<BSSurfaceRational<T, 3>>(build_torus<T>(major_, minor_, this->f_.ax(), u1, u2, v1, v2));
        }
        std::array<T, 2> nurbs_parameters(T u, T v) const override
        {
            return {ellipse_arc_parameter(this->box_[0], this->box_[1], u), ellipse_arc_parameter(this->box_[2], this->box_[3], v)};
        }
        std::shared_ptr<ElementarySurface<T>> moved(const std::function<point<T, 3>(const point<T, 3> &)> &m) const override
        {
            return std::make_shared<ToroidalSurface<T>>(major_, minor_, this->f_.moved(m), this->box_);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_u(T u, T v1, T v2) const override
        {
            // section: O + R rho(u) + r cos v rho(u) + r sin v Z
            return std::make_shared<EllipseCurve<T>>(minor_, minor_, primitives_detail::meridian(this->f_, u, this->f_.o + major_ * this->f_.radial(u)), v1, v2);
        }
        std::shared_ptr<ElementaryCurve<T>> iso_v(T v, T u1, T u2) const override
        {
            return this->parallel(major_ + minor_ * std::cos(v), minor_ * std::sin(v), u1, u2);
        }
    };
}
