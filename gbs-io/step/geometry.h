#pragma once

/**
 * @file geometry.h
 * @brief Translation of STEP geometry (ISO 10303-42) to gbs geometry.
 *
 * Design: docs/sources/design/step_reader.md, sections 3.2, 3.3 and 4;
 * architecture note docs/sources/design/step_pr04_geometry.md.
 *
 * Two stages:
 *  1. GeometryReader reads an instance into a *definition* (CurveDef, SurfaceDef)
 *     that keeps the STEP parametrization, already converted to the target units
 *     (lengths) and to radians (angles): lines and planes are infinite, circles
 *     are parametrized by the angle. A definition evaluates points, inverts
 *     points to parameters (analytically for lines, conics and elementary
 *     surfaces) and knows the kind of each parameter direction.
 *  2. to_curve() / to_surface() build the geometry over a parameter range: an
 *     analytic primitive (gbs/primitives.h) for a line, a conic or an
 *     elementary surface, an exact NURBS otherwise; to_nurbs() always gives
 *     the exact NURBS (gbs/bselementary.h). The range of an edge comes from
 *     its vertices, the box of a face from its boundary (parameter_box()).
 */

#include <gbs-io/step/p21.h>
#include <gbs-io/step/units.h>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bsctools.h>
#include <gbs/bselementary.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <optional>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace gbs::step
{
    using Real = double;
    template <std::size_t dim>
    using Point = point<Real, dim>;

    /// What a parameter measures, which decides its unit conversion and its periodicity.
    enum class ParamKind : std::uint8_t
    {
        None,   ///< dimensionless (B-spline, line with a vector of given magnitude)
        Length, ///< a length (plane, height of a cylinder or cone)
        Angle,  ///< an angle in radians, periodic on a full turn
    };

    // =========================================================================
    // Curves
    // =========================================================================

    template <std::size_t dim>
    struct CurveDef;

    /// C(t) = origin + t V
    template <std::size_t dim>
    struct LineDef
    {
        Point<dim> origin, vector;
    };

    /// C(t) = center + r1 cos t X + r2 sin t Y, X and Y orthonormal (circle: r1 = r2)
    template <std::size_t dim>
    struct ConicDef
    {
        Point<dim> center, x, y;
        Real r1, r2;
    };

    /// A bounded NURBS in its own parameter (B-spline, Bezier, polyline, composite curve).
    template <std::size_t dim>
    struct SplineDef
    {
        std::shared_ptr<const Curve<Real, dim>> curve;
    };

    /// Basis curve restricted to [t1, t2]; reversed when the curve runs from t2 to t1.
    template <std::size_t dim>
    struct TrimmedDef
    {
        std::shared_ptr<const CurveDef<dim>> basis;
        Real t1, t2;
        bool reversed;
    };

    /// NURBS of an edge with the NURBS parameters of the requested STEP range.
    template <std::size_t dim>
    struct NurbsCurve
    {
        std::shared_ptr<Curve<Real, dim>> curve;
        Real u1, u2;
    };

    namespace detail
    {
        constexpr Real two_pi = 2. * std::numbers::pi;
        constexpr Real inf = std::numeric_limits<Real>::infinity();

        /**
         * Relative precision of the angles of a file. A degree factor written with 10 digits
         * (0.01745329252, SolidWorks) puts 360 degrees 2e-11 rad above 2 pi, and up to 1.8e-9 rad
         * away for the worst rounding of the last digit: a span within angle_precision x max(2 pi, span)
         * of a whole number of turns is that number of turns (turn_end()).
         */
        constexpr Real angle_precision = 1e-9;

        inline Real wrap_to(Real t, Real lo) // t + 2 pi k in [lo, lo + 2 pi)
        {
            return t - two_pi * std::floor((t - lo) / two_pi);
        }

        /**
         * End of the angular range running from `from` to `to` in the increasing sense, in
         * (from, from + 2 pi]: a full turn, exactly, when `to` is `from` up to whole turns within the
         * precision of the file (angle_precision), never a span slightly above 2 pi nor a sliver.
         */
        inline Real turn_end(Real from, Real to)
        {
            const Real d = to - from;
            if (std::abs(d - two_pi * std::round(d / two_pi)) <= angle_precision * std::max(two_pi, std::abs(d)))
                return from + two_pi;
            return wrap_to(to, from);
        }

        /**
         * Parameters at which a curve is sampled over [a, b] before a local search: p + 2 per non-empty
         * knot span of a B-spline (64 in all at least), so that each span, hence each branch of a
         * hairpin narrower than the sampling step, has its own samples; 64 uniform ones otherwise.
         */
        template <std::size_t dim>
        std::vector<Real> sample_parameters(const Curve<Real, dim> &c, Real a, Real b)
        {
            const std::vector<Real> *k = nullptr;
            std::size_t p = 1;
            if (auto s = dynamic_cast<const BSCurveGeneral<Real, dim, false> *>(&c))
                k = &s->knotsFlats(), p = s->degree();
            else if (auto r = dynamic_cast<const BSCurveGeneral<Real, dim, true> *>(&c))
                k = &r->knotsFlats(), p = r->degree();
            std::vector<Real> breaks{a};
            if (k)
                for (Real x : *k)
                    if (x > breaks.back() && x < b)
                        breaks.push_back(x);
            breaks.push_back(b);
            const std::size_t spans = breaks.size() - 1;
            const std::size_t m = std::max<std::size_t>(p + 2, (64 + spans - 1) / spans);
            std::vector<Real> t;
            t.reserve(spans * m + 1);
            for (std::size_t i = 0; i < spans; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    t.push_back(breaks[i] + (breaks[i + 1] - breaks[i]) * Real(j) / Real(m));
            t.push_back(b);
            return t;
        }

        /**
         * Parameter minimising f over [s.front(), s.back()], from the samples s (increasing).
         *
         * The best local minima of the sampled values (8 at most) are refined: on each side of the
         * sample where f decreases (its one-sided derivative, taken one ulp inside, so that a knot
         * gives the derivative of the span searched), Newton on f' = 0 safeguarded by bisection on
         * the sign of f' up to the neighbouring sample (a kink of a degree 1 curve converges to its
         * knot). The least f wins; values whose square roots (the distances) are within `tie` of the
         * least are ties, won by the smallest parameter. f(t) gives f, fd(t) gives {f, f', f''}.
         */
        template <typename F, typename FD>
        Real minimize_sampled(F &&f, FD &&fd, const std::vector<Real> &s, Real tie)
        {
            const std::size_t n = s.size();
            std::vector<Real> fs(n);
            for (std::size_t i = 0; i < n; ++i)
                fs[i] = f(s[i]);
            std::vector<std::size_t> minima;
            for (std::size_t i = 0; i < n; ++i)
                if ((i == 0 || fs[i] < fs[i - 1]) && (i + 1 == n || fs[i] <= fs[i + 1]))
                    minima.push_back(i);
            std::ranges::stable_sort(minima, {}, [&](std::size_t i) { return fs[i]; });
            if (minima.size() > 8)
                minima.resize(8);

            auto df = [&](Real x) { return fd(x)[1]; };
            // zero of f' in [lo, hi], f'(lo) < 0 < f'(hi), from t
            auto newton = [&](Real lo, Real hi, Real t) {
                constexpr Real eps = 4. * std::numeric_limits<Real>::epsilon();
                for (int it = 0; it < 100; ++it)
                {
                    const auto e = fd(t);
                    if (e[1] == 0.)
                        break;
                    (e[1] < 0. ? lo : hi) = t;
                    Real next = t - e[1] / e[2];
                    if (!(e[2] > 0.) || !(next > lo && next < hi))
                        next = 0.5 * (lo + hi);
                    const Real step = eps * std::max({std::abs(lo), std::abs(hi), hi - lo});
                    const bool done = std::abs(next - t) <= step || hi - lo <= step;
                    t = next;
                    if (done)
                        break;
                }
                return t;
            };
            // minimum in [a, b] where f decreases from a (right) or from b (left)
            auto side = [&](Real a, Real b, bool right) -> std::optional<Real> {
                Real lo = a, hi = b;
                Real &outer = right ? hi : lo;
                const Real gf = df(right ? std::nextafter(b, a) : std::nextafter(a, b));
                if (right ? !(gf > 0.) : !(gf < 0.))
                {
                    // f' of one sign at both ends: a minimum then a maximum in between, if any
                    bool found = false;
                    for (int j = 1; j < 16 && !found; ++j)
                    {
                        const Real x = right ? a + (b - a) * j / 16. : b - (b - a) * j / 16.;
                        const Real gx = df(x);
                        if (right ? gx > 0. : gx < 0.)
                            outer = x, found = true;
                    }
                    if (!found)
                        return std::nullopt;
                }
                return newton(lo, hi, right ? std::nextafter(lo, hi) : std::nextafter(hi, lo));
            };

            std::vector<std::pair<Real, Real>> found; // (t, f)
            Real best = s.front(), fbest = inf;
            for (auto i : minima)
            {
                Real t = s[i], ft = fs[i];
                auto keep = [&](std::optional<Real> x) {
                    if (!x)
                        return;
                    if (const Real fx = f(*x); fx < ft)
                        t = *x, ft = fx;
                };
                if (i + 1 < n && df(std::nextafter(s[i], s[i + 1])) < 0.)
                    keep(side(s[i], s[i + 1], true));
                if (i > 0 && df(std::nextafter(s[i], s[i - 1])) > 0.)
                    keep(side(s[i - 1], s[i], false));
                found.emplace_back(t, ft);
                if (ft < fbest)
                    fbest = ft, best = t;
            }
            const Real dbest = std::sqrt(std::max<Real>(fbest, 0.));
            for (auto [t, ft] : found)
                if (t < best && std::sqrt(std::max<Real>(ft, 0.)) <= dbest + tie)
                    best = t;
            return best;
        }

        /**
         * Parameter of the point of c nearest to p, seen through the linear projection `proj`
         * (identity: in space; I - V V^T / |V|^2: in the plane normal to V). Ties (within 1e-12
         * relative) go to the smallest parameter: the closure point of a closed curve gives its start.
         */
        template <std::size_t dim, typename Proj = std::identity>
        Real closest_parameter(const Curve<Real, dim> &c, const Point<dim> &p, Proj proj = {})
        {
            const auto [a, b] = c.bounds();
            auto f = [&](Real t) {
                const Point<dim> r = proj(c.value(t) - p);
                return r * r;
            };
            auto fd = [&](Real t) -> std::array<Real, 3> {
                const Point<dim> r = proj(c.value(t) - p), d1 = proj(c.value(t, 1)), d2 = proj(c.value(t, 2));
                return {r * r, 2. * (r * d1), 2. * (d1 * d1 + r * d2)};
            };
            return minimize_sampled(f, fd, sample_parameters(c, a, b), 1e-12 * std::max<Real>(1., norm(p)));
        }
    }

    /**
     * @brief A STEP curve in its own parametrization (see the file comment).
     */
    template <std::size_t dim>
    struct CurveDef
    {
        std::variant<LineDef<dim>, ConicDef<dim>, SplineDef<dim>, TrimmedDef<dim>> def;
        std::uint64_t id{};

        /// Point at t, or its derivative of order d.
        [[nodiscard]] Point<dim> value(Real t, std::size_t d = 0) const
        {
            return std::visit([t, d](const auto &c) -> Point<dim> {
                using D = std::decay_t<decltype(c)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return d == 0 ? c.origin + t * c.vector : d == 1 ? c.vector : Point<dim>{};
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                {
                    // (cos, sin) turned by d quarter turns: their derivative of order d
                    Real co = std::cos(t), si = std::sin(t);
                    for (std::size_t k = 0; k < d % 4; ++k)
                        std::tie(co, si) = std::pair{-si, co};
                    const Point<dim> v = (c.r1 * co) * c.x + (c.r2 * si) * c.y;
                    return d == 0 ? c.center + v : v;
                }
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return c.curve->value(t, d);
                else
                    return c.basis->value(t, d);
            }, def);
        }

        /**
         * @brief Parameter of the curve point nearest to p (exact for a point on a line or conic).
         *
         * For a B-spline, ties (points equally near several parts of the curve, within 1e-12
         * relative) go to the smallest parameter: the closure point of a closed B-spline gives the
         * start of its domain. A conic gives an angle in [0, 2 pi), a trimmed conic the turn of the
         * angle closest to its range.
         */
        [[nodiscard]] Real parameter(const Point<dim> &p) const
        {
            return std::visit([&p](const auto &d) -> Real {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return ((p - d.origin) * d.vector) / sq_norm(d.vector);
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                    return detail::wrap_to(std::atan2(((p - d.center) * d.y) / d.r2, ((p - d.center) * d.x) / d.r1), 0.);
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return detail::closest_parameter(*d.curve, p);
                else
                {
                    const Real t = d.basis->parameter(p);
                    if (!d.basis->periodic())
                        return t;
                    // the turn of t closest to [t1, t2]
                    const Real m = 0.5 * (d.t1 + d.t2);
                    return detail::wrap_to(t, m - std::numbers::pi);
                }
            }, def);
        }

        /**
         * @brief Parameter of the curve point nearest to p seen along w, that is in the plane normal
         * to w: the u of a point of the extrusion C(u) + v w (exact for a point on it, and for a line
         * or a conic not parallel to w). Ties as parameter().
         */
        [[nodiscard]] Real parameter_along(const Point<dim> &p, const Point<dim> &w) const
            requires(dim == 3)
        {
            const Real ww = sq_norm(w);
            return std::visit([&](const auto &d) -> Real {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                {
                    // p - O = t L + s w, least squares
                    const auto q = p - d.origin;
                    const Real a = sq_norm(d.vector), b = d.vector * w, det = a * ww - b * b;
                    return (ww * (q * d.vector) - b * (q * w)) / det;
                }
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                {
                    // slide p along w onto the plane of the conic
                    const auto n = cross(d.x, d.y);
                    const Real wn = w * n;
                    return parameter(std::abs(wn) > 1e-12 * std::sqrt(ww) ? p - (((p - d.center) * n) / wn) * w : p);
                }
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return detail::closest_parameter(*d.curve, p, [&](const Point<dim> &x) { return x - ((x * w) / ww) * w; });
                else
                {
                    const Real t = d.basis->parameter_along(p, w);
                    return d.basis->periodic() ? detail::wrap_to(t, 0.5 * (d.t1 + d.t2) - std::numbers::pi) : t;
                }
            }, def);
        }

        /// Parameters at which to sample the curve over [t1, t2] (finite) before a local search.
        [[nodiscard]] std::vector<Real> samples(Real t1, Real t2) const
        {
            if (auto s = std::get_if<SplineDef<dim>>(&def))
                return detail::sample_parameters(*s->curve, t1, t2);
            if (auto tr = std::get_if<TrimmedDef<dim>>(&def))
                return tr->basis->samples(t1, t2);
            // a line, or a conic over a turn at most: 16 intervals
            constexpr std::size_t n = 16;
            std::vector<Real> t(n + 1);
            for (std::size_t i = 0; i <= n; ++i)
                t[i] = t1 + (t2 - t1) * Real(i) / Real(n);
            return t;
        }

        /// Parameter domain: infinite for a line, [0, 2 pi] for a full conic.
        [[nodiscard]] std::array<Real, 2> domain() const
        {
            return std::visit([](const auto &d) -> std::array<Real, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                    return {-detail::inf, detail::inf};
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                    return {0., detail::two_pi};
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return d.curve->bounds();
                else
                    return {d.t1, d.t2};
            }, def);
        }

        /// True for a full circle or ellipse: any range of at most 2 pi is valid.
        [[nodiscard]] bool periodic() const { return std::holds_alternative<ConicDef<dim>>(def); }

        /// True for a trimmed curve running against its basis parameter.
        [[nodiscard]] bool reversed() const
        {
            auto t = std::get_if<TrimmedDef<dim>>(&def);
            return t && t->reversed;
        }

        [[nodiscard]] ParamKind kind() const
        {
            if (std::holds_alternative<ConicDef<dim>>(def))
                return ParamKind::Angle;
            if (auto t = std::get_if<TrimmedDef<dim>>(&def))
                return t->basis->kind();
            return ParamKind::None;
        }

        /// NURBS parameter of the point at parameter t, for the NURBS built by to_nurbs(t1, t2).
        [[nodiscard]] Real nurbs_parameter(Real t1, Real t2, Real t) const
        {
            if (std::holds_alternative<ConicDef<dim>>(def))
                return ellipse_arc_parameter(t1, t2, t);
            if (auto tr = std::get_if<TrimmedDef<dim>>(&def))
                return tr->basis->nurbs_parameter(t1, t2, t);
            return t;
        }

        /**
         * @brief Exact NURBS of the curve between t1 < t2 (STEP parameters).
         *
         * A line gives a segment parametrized by t, a conic an arc whose ends are
         * at the parameters t1 and t2 (2 pi at most), a B-spline itself with u1 = t1, u2 = t2.
         */
        [[nodiscard]] NurbsCurve<dim> to_nurbs(Real t1, Real t2) const
        {
            if (!(t2 > t1) || !std::isfinite(t1) || !std::isfinite(t2))
                throw StepError(id, "invalid curve range");
            return std::visit([&](const auto &d) -> NurbsCurve<dim> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, LineDef<dim>>)
                {
                    points_vector<Real, dim> poles{value(t1), value(t2)};
                    return {std::make_shared<BSCurve<Real, dim>>(poles, std::vector<Real>{t1, t1, t2, t2}, 1), t1, t2};
                }
                else if constexpr (std::is_same_v<D, ConicDef<dim>>)
                {
                    if (t2 - t1 > detail::two_pi * (1. + 1e-12))
                        throw StepError(id, "conic range larger than a turn");
                    auto c = build_ellipse_arc<Real, dim>(d.r1, d.r2, t1, std::min(t2, t1 + detail::two_pi), d.center, d.x, d.y);
                    return {std::make_shared<BSCurveRational<Real, dim>>(std::move(c)), t1, t2};
                }
                else if constexpr (std::is_same_v<D, SplineDef<dim>>)
                    return {std::const_pointer_cast<Curve<Real, dim>>(d.curve), t1, t2};
                else
                    return d.basis->to_nurbs(t1, t2);
            }, def);
        }

        /**
         * @brief Curve of an edge between t1 < t2 (STEP parameters): a primitive (gbs/primitives.h)
         * for a line or a conic, parametrized like STEP (u1 = t1, u2 = t2); the NURBS otherwise.
         */
        [[nodiscard]] NurbsCurve<dim> to_curve(Real t1, Real t2) const
        {
            if constexpr (dim == 3)
            {
                if (!(t2 > t1) || !std::isfinite(t1) || !std::isfinite(t2))
                    throw StepError(id, "invalid curve range");
                if (auto l = std::get_if<LineDef<3>>(&def))
                    return {std::make_shared<LineCurve<Real>>(l->origin, l->vector, t1, t2), t1, t2};
                if (auto c = std::get_if<ConicDef<3>>(&def))
                {
                    if (t2 - t1 > detail::two_pi * (1. + 1e-12))
                        throw StepError(id, "conic range larger than a turn");
                    Frame3<Real> f;
                    f.o = c->center, f.x = c->x, f.y = c->y, f.z = cross(c->x, c->y);
                    return {std::make_shared<EllipseCurve<Real>>(c->r1, c->r2, f, t1, std::min(t2, t1 + detail::two_pi)), t1, t2};
                }
                if (auto tr = std::get_if<TrimmedDef<3>>(&def))
                    return tr->basis->to_curve(t1, t2);
            }
            return to_nurbs(t1, t2);
        }
    };

    // =========================================================================
    // Surfaces
    // =========================================================================

    /// Orthonormal frame of an axis2_placement_3d.
    struct Frame
    {
        Point<3> origin, x, y, z;
        [[nodiscard]] ax2<Real, 3> ax() const { return {origin, z, x}; }
    };

    struct PlaneDef { Frame f; };                          ///< S = O + u X + v Y
    struct CylinderDef { Frame f; Real r; };               ///< S = O + R (cos u X + sin u Y) + v Z
    struct ConeDef { Frame f; Real r, semi_angle; };       ///< S = O + (R + v tan a)(cos u X + sin u Y) + v Z
    struct SphereDef { Frame f; Real r; };                 ///< S = O + R cos v (cos u X + sin u Y) + R sin v Z
    /**
     * Part of a degenerate torus (minor radius r > major radius R, self-intersecting) selected by
     * select_outer (ISO 10303-42): the apple, outer part, v in [-t, t]; the lemon, inner part, v in
     * [t, 2 pi - t], where R + r cos v <= 0 (a point at u lies on the side of u + pi); cos t = -R / r.
     */
    enum class TorusPart : std::uint8_t
    {
        Whole, ///< a torus, v periodic
        Apple, ///< select_outer .T.
        Lemon, ///< select_outer .F.
    };
    /// S = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
    struct TorusDef { Frame f; Real major, minor; TorusPart part = TorusPart::Whole; };
    /// S(u, v) = generatrix(v) rotated by u around (origin, z)
    struct RevolutionDef { std::shared_ptr<const CurveDef<3>> generatrix; Point<3> origin, z; };
    /// S(u, v) = C(u) + v V
    struct ExtrusionDef { std::shared_ptr<const CurveDef<3>> curve; Point<3> v; };
    /// B-spline surface in its own parameters
    struct SplineSurfaceDef { std::shared_ptr<const Surface<Real, 3>> surface; };
    struct SurfaceDef;
    /// Basis surface restricted to [u1, u2] x [v1, v2]
    struct TrimmedSurfaceDef { std::shared_ptr<const SurfaceDef> basis; Real u1, u2, v1, v2; };

    /// Parameter box of a surface: [u1, u2] x [v1, v2].
    struct UVBox
    {
        Real u1, u2, v1, v2;
    };

    namespace detail
    {
        /**
         * (u, v) of a point p of a surface of revolution, whatever its generatrix (skew to the axis,
         * across it). v minimises the distance in the meridian half-plane between the generatrix
         * point and p, (a(v) - a_p)^2 + (rho(v) - rho_p)^2, a being the axial coordinate and rho the
         * distance to the axis; u then turns the generatrix point onto p (0 when either is on the axis).
         */
        inline std::array<Real, 2> revolution_parameters(const RevolutionDef &d, const Point<3> &p)
        {
            const CurveDef<3> &g = *d.generatrix;
            // axial coordinate and radial part of a vector
            auto split = [&](const Point<3> &x) {
                const Real a = x * d.z;
                return std::pair{a, x - a * d.z};
            };
            const auto mp = split(p - d.origin);
            const Real ap = mp.first;
            const Point<3> rp = mp.second;
            const Real rhop = norm(rp);
            auto f = [&](Real v) {
                const auto [a, r] = split(g.value(v) - d.origin);
                const Real da = a - ap, dr = norm(r) - rhop;
                return da * da + dr * dr;
            };
            auto fd = [&](Real v) -> std::array<Real, 3> {
                const auto [a, r] = split(g.value(v) - d.origin);
                const auto [a1, r1] = split(g.value(v, 1));
                const auto [a2, r2] = split(g.value(v, 2));
                const Real rho = norm(r), da = a - ap, dr = rho - rhop;
                // rho' = r.r' / rho, rho'' = (r'.r' + r.r'' - rho'^2) / rho; a kink on the axis
                Real rho1 = 0., rho2 = 0.;
                if (rho > 0.)
                {
                    rho1 = (r * r1) / rho;
                    rho2 = (r1 * r1 + r * r2 - rho1 * rho1) / rho;
                }
                return {da * da + dr * dr, 2. * (da * a1 + dr * rho1), 2. * (a1 * a1 + da * a2 + rho1 * rho1 + dr * rho2)};
            };
            Real v;
            const auto [g1, g2] = g.domain();
            if (std::isfinite(g1) && std::isfinite(g2))
                v = minimize_sampled(f, fd, g.samples(g1, g2), 1e-12 * std::max<Real>(1., norm(p - d.origin)));
            else if (auto l = std::get_if<LineDef<3>>(&g.def))
            {
                // an unbounded line O + t L: a(t) = a0 + t a1 = a_p, exact on the surface
                const auto [a0, r0] = split(l->origin - d.origin);
                const auto [a1, r1] = split(l->vector);
                const Real rr = sq_norm(r1);
                if (std::abs(a1) > 1e-12 * norm(l->vector))
                    v = (ap - a0) / a1;
                else // normal to the axis: a plane annulus covered twice, the root of rho(t) = rho_p after the foot
                {
                    const Real foot = -(r0 * r1) / rr;
                    v = foot + std::sqrt(std::max<Real>(0., rhop * rhop - sq_norm(r0 + foot * r1)) / rr);
                }
            }
            else
                v = g.parameter(p); // not reached: only a line is unbounded
            // u turns the radial direction of the generatrix point onto that of p
            const auto rg = split(g.value(v) - d.origin).second;
            if (!(norm(rg) > 0.) || !(rhop > 0.))
                return {0., v};
            return {wrap_to(std::atan2(d.z * cross(rg, rp), rg * rp), 0.), v};
        }
    }

    /**
     * @brief A STEP surface in its own parametrization (see the file comment).
     */
    struct SurfaceDef
    {
        std::variant<PlaneDef, CylinderDef, ConeDef, SphereDef, TorusDef, RevolutionDef, ExtrusionDef, SplineSurfaceDef, TrimmedSurfaceDef> def;
        std::uint64_t id{};

        [[nodiscard]] Point<3> value(Real u, Real v) const
        {
            return std::visit([&](const auto &d) -> Point<3> {
                using D = std::decay_t<decltype(d)>;
                auto radial = [u](const Frame &f) { return std::cos(u) * f.x + std::sin(u) * f.y; };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return d.f.origin + u * d.f.x + v * d.f.y;
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return d.f.origin + d.r * radial(d.f) + v * d.f.z;
                else if constexpr (std::is_same_v<D, ConeDef>)
                    return d.f.origin + (d.r + v * std::tan(d.semi_angle)) * radial(d.f) + v * d.f.z;
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return d.f.origin + (d.r * std::cos(v)) * radial(d.f) + (d.r * std::sin(v)) * d.f.z;
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return d.f.origin + (d.major + d.minor * std::cos(v)) * radial(d.f) + (d.minor * std::sin(v)) * d.f.z;
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    const auto p = d.generatrix->value(v) - d.origin;
                    const auto h = (p * d.z) * d.z;
                    const auto r = p - h;
                    return d.origin + h + std::cos(u) * r + std::sin(u) * cross(d.z, r);
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return d.curve->value(u) + v * d.v;
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return d.surface->value(u, v);
                else
                    return d.basis->value(u, v);
            }, def);
        }

        /// Kinds of the u and v parameters.
        [[nodiscard]] std::array<ParamKind, 2> kinds() const
        {
            using K = ParamKind;
            return std::visit([](const auto &d) -> std::array<K, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {K::Length, K::Length};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return {K::Angle, K::Length};
                else if constexpr (std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                    return {K::Angle, K::Angle};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {K::Angle, d.generatrix->kind()};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->kind(), K::None};
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return {K::None, K::None};
                else
                    return d.basis->kinds();
            }, def);
        }

        /// Natural parameter domain (infinite in the open directions, [-pi/2, pi/2] for a sphere latitude…).
        [[nodiscard]] UVBox domain() const
        {
            using detail::inf;
            constexpr Real tp = detail::two_pi, hp = 0.5 * std::numbers::pi;
            return std::visit([](const auto &d) -> UVBox {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {-inf, inf, -inf, inf};
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return {0., tp, -inf, inf};
                else if constexpr (std::is_same_v<D, ConeDef>)
                {
                    // the half cone of the reference radius, up to the apex
                    const Real t = std::tan(d.semi_angle);
                    if (t == 0.)
                        return {0., tp, -inf, inf};
                    const Real apex = -d.r / t;
                    return t > 0. ? UVBox{0., tp, apex, inf} : UVBox{0., tp, -inf, apex};
                }
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return {0., tp, -hp, hp};
                else if constexpr (std::is_same_v<D, TorusDef>)
                {
                    if (d.part == TorusPart::Whole)
                        return {0., tp, 0., tp};
                    const Real t = std::acos(-d.major / d.minor);
                    return d.part == TorusPart::Apple ? UVBox{0., tp, -t, t} : UVBox{0., tp, t, tp - t};
                }
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    auto [a, b] = d.generatrix->domain();
                    return {0., tp, a, b};
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    auto [a, b] = d.curve->domain();
                    return {a, b, -inf, inf};
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                {
                    auto b = d.surface->bounds();
                    return {b[0], b[1], b[2], b[3]};
                }
                else
                    return {d.u1, d.u2, d.v1, d.v2};
            }, def);
        }

        /// Periodicity of u and v (full turn, not trimmed).
        [[nodiscard]] std::array<bool, 2> periodic() const
        {
            return std::visit([](const auto &d) -> std::array<bool, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef> || std::is_same_v<D, SphereDef>)
                    return {true, false};
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return {true, d.part == TorusPart::Whole};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {true, d.generatrix->periodic()};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->periodic(), false};
                else
                    return {false, false};
            }, def);
        }

        /**
         * @brief Parameters (u, v) of a point on the surface, analytic except for
         * B-spline surfaces (for which this returns nullopt: use the natural domain).
         *
         * u is in [0, 2 pi) in a periodic direction; a trimmed surface gives the turn closest to its
         * range. A point on the axis or at a pole has u = 0. A surface of revolution is inverted in
         * its meridian half-plane, numerically unless its generatrix is an unbounded line; an
         * extrusion of a B-spline by the B-spline seen along the extrusion vector.
         */
        [[nodiscard]] std::optional<std::array<Real, 2>> parameters(const Point<3> &p) const
        {
            return std::visit([&p](const auto &d) -> std::optional<std::array<Real, 2>> {
                using D = std::decay_t<decltype(d)>;
                auto angle = [&p](const Frame &f) {
                    const auto q = p - f.origin;
                    return detail::wrap_to(std::atan2(q * f.y, q * f.x), 0.);
                };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return std::array{(p - d.f.origin) * d.f.x, (p - d.f.origin) * d.f.y};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return std::array{angle(d.f), (p - d.f.origin) * d.f.z};
                else if constexpr (std::is_same_v<D, SphereDef>)
                {
                    const auto q = p - d.f.origin;
                    const Real h = q * d.f.z;
                    const Real rho = norm(q - h * d.f.z);
                    return std::array{angle(d.f), std::atan2(h, rho)};
                }
                else if constexpr (std::is_same_v<D, TorusDef>)
                {
                    const auto q = p - d.f.origin;
                    const Real h = q * d.f.z;
                    const Real rho = norm(q - h * d.f.z);
                    if (d.part == TorusPart::Lemon) // R + r cos v = -rho: the point is on the side of u + pi
                        return std::array{rho > 0. ? detail::wrap_to(angle(d.f) + std::numbers::pi, 0.) : 0.,
                                          detail::wrap_to(std::atan2(h, -rho - d.major), 0.)};
                    if (d.part == TorusPart::Apple) // v in [-t, t]
                        return std::array{angle(d.f), std::atan2(h, rho - d.major)};
                    return std::array{angle(d.f), detail::wrap_to(std::atan2(h, rho - d.major), 0.)};
                }
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return detail::revolution_parameters(d, p);
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    const Real u = d.curve->parameter_along(p, d.v);
                    return std::array{u, ((p - d.curve->value(u)) * d.v) / sq_norm(d.v)};
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return std::nullopt;
                else
                {
                    auto uv = d.basis->parameters(p);
                    // in a periodic direction of the basis, the turn closest to the trimmed range
                    const auto per = d.basis->periodic();
                    if (uv && per[0])
                        (*uv)[0] = detail::wrap_to((*uv)[0], 0.5 * (d.u1 + d.u2) - std::numbers::pi);
                    if (uv && per[1])
                        (*uv)[1] = detail::wrap_to((*uv)[1], 0.5 * (d.v1 + d.v2) - std::numbers::pi);
                    return uv;
                }
            }, def);
        }

        /// NURBS parameters of the point at (u, v), for the NURBS built by to_nurbs(box).
        [[nodiscard]] std::array<Real, 2> nurbs_parameters(const UVBox &b, Real u, Real v) const
        {
            return std::visit([&](const auto &d) -> std::array<Real, 2> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return {u, v};
                else if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), v};
                else if constexpr (std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), ellipse_arc_parameter(b.v1, b.v2, v)};
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                    return {ellipse_arc_parameter(b.u1, b.u2, u), d.generatrix->nurbs_parameter(b.v1, b.v2, v)};
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                    return {d.curve->nurbs_parameter(b.u1, b.u2, u), v};
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return {u, v};
                else
                    return d.basis->nurbs_parameters(b, u, v);
            }, def);
        }

        /// Exact NURBS of the surface over the box (STEP parameters; a B-spline surface is returned whole).
        [[nodiscard]] std::shared_ptr<Surface<Real, 3>> to_nurbs(const UVBox &b) const
        {
            if (!(b.u2 > b.u1) || !(b.v2 > b.v1) || !std::isfinite(b.u1 + b.u2 + b.v1 + b.v2))
                throw StepError(id, "invalid surface parameter box");
            return std::visit([&](const auto &d) -> std::shared_ptr<Surface<Real, 3>> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                {
                    points_vector<Real, 3> poles{value(b.u1, b.v1), value(b.u2, b.v1), value(b.u1, b.v2), value(b.u2, b.v2)};
                    return std::make_shared<BSSurface<Real, 3>>(poles, std::vector<Real>{b.u1, b.u1, b.u2, b.u2},
                                                                std::vector<Real>{b.v1, b.v1, b.v2, b.v2}, 1, 1);
                }
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_cylinder<Real>(d.r, d.f.ax(), b.v1, b.v2, b.u1, b.u2));
                else if constexpr (std::is_same_v<D, ConeDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_cone<Real>(d.r, d.semi_angle, d.f.ax(), b.v1, b.v2, b.u1, b.u2));
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_sphere<Real>(d.r, d.f.ax(), b.u1, b.u2, b.v1, b.v2));
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return std::make_shared<BSSurfaceRational<Real, 3>>(build_torus<Real>(d.major, d.minor, d.f.ax(), b.u1, b.u2, b.v1, b.v2));
                else if constexpr (std::is_same_v<D, RevolutionDef>)
                {
                    auto g = d.generatrix->to_nurbs(b.v1, b.v2).curve;
                    const ax1<Real, 3> axis{d.origin, d.z};
                    if (auto r = std::dynamic_pointer_cast<BSCurveRational<Real, 3>>(g))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_revolution(*r, axis, b.u1, b.u2));
                    if (auto p = std::dynamic_pointer_cast<BSCurve<Real, 3>>(g))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_revolution(*p, axis, b.u1, b.u2));
                    throw StepError(id, "generatrix is not a NURBS");
                }
                else if constexpr (std::is_same_v<D, ExtrusionDef>)
                {
                    auto c = d.curve->to_nurbs(b.u1, b.u2).curve;
                    if (auto r = std::dynamic_pointer_cast<BSCurveRational<Real, 3>>(c))
                        return std::make_shared<BSSurfaceRational<Real, 3>>(build_extrusion(*r, d.v, b.v1, b.v2));
                    if (auto p = std::dynamic_pointer_cast<BSCurve<Real, 3>>(c))
                        return std::make_shared<BSSurface<Real, 3>>(build_extrusion(*p, d.v, b.v1, b.v2));
                    throw StepError(id, "swept curve is not a NURBS");
                }
                else if constexpr (std::is_same_v<D, SplineSurfaceDef>)
                    return std::const_pointer_cast<Surface<Real, 3>>(d.surface);
                else
                    return d.basis->to_nurbs(b);
            }, def);
        }

        /// True for a plane, cylinder, cone, sphere or torus (possibly trimmed): kept as a primitive.
        [[nodiscard]] bool elementary() const
        {
            if (auto t = std::get_if<TrimmedSurfaceDef>(&def))
                return t->basis->elementary();
            return std::holds_alternative<PlaneDef>(def) || std::holds_alternative<CylinderDef>(def) || std::holds_alternative<ConeDef>(def) ||
                   std::holds_alternative<SphereDef>(def) || std::holds_alternative<TorusDef>(def);
        }

        /**
         * @brief Surface of a face over the box (STEP parameters): a primitive (gbs/primitives.h),
         * parametrized like STEP, for an elementary surface; the exact NURBS otherwise.
         */
        [[nodiscard]] std::shared_ptr<Surface<Real, 3>> to_surface(const UVBox &b) const
        {
            if (!(b.u2 > b.u1) || !(b.v2 > b.v1) || !std::isfinite(b.u1 + b.u2 + b.v1 + b.v2))
                throw StepError(id, "invalid surface parameter box");
            const std::array<Real, 4> box{b.u1, b.u2, b.v1, b.v2};
            auto frame = [](const Frame &f) {
                Frame3<Real> g;
                g.o = f.origin, g.x = f.x, g.y = f.y, g.z = f.z;
                return g;
            };
            return std::visit([&](const auto &d) -> std::shared_ptr<Surface<Real, 3>> {
                using D = std::decay_t<decltype(d)>;
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return std::make_shared<Plane<Real>>(frame(d.f), box);
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return std::make_shared<CylindricalSurface<Real>>(d.r, frame(d.f), box);
                else if constexpr (std::is_same_v<D, ConeDef>)
                    return std::make_shared<ConicalSurface<Real>>(d.r, d.semi_angle, frame(d.f), box);
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return std::make_shared<SphericalSurface<Real>>(d.r, frame(d.f), box);
                else if constexpr (std::is_same_v<D, TorusDef>)
                    return std::make_shared<ToroidalSurface<Real>>(d.major, d.minor, frame(d.f), box);
                else if constexpr (std::is_same_v<D, TrimmedSurfaceDef>)
                    return d.basis->to_surface(b);
                else
                    return to_nurbs(b);
            }, def);
        }

        /// Parameters, on the surface built by to_surface(box), of the point at STEP parameters (u, v).
        [[nodiscard]] std::array<Real, 2> surface_parameters(const UVBox &b, Real u, Real v) const
        {
            return elementary() ? std::array<Real, 2>{u, v} : nurbs_parameters(b, u, v);
        }
    };

    namespace detail
    {
        /**
         * Step of an angle between two points along a loop, d = a2 - a1 modulo 2 pi: the smallest one,
         * except across a pole (`pole`: +1 at the top of v, -1 at the bottom, 0 none), where the loop
         * runs along the pole in the sense that keeps the face on its side. With the face on the left
         * of its loops in (u, v) (`sense` > 0), the face is below a pole at the top of v: the step is
         * in (-2 pi, 0] there, in [0, 2 pi) at the bottom; the opposite for `sense` < 0. A lune between
         * two meridians is wider or narrower than pi by its sense only.
         */
        inline Real angle_step(Real d, int pole, int sense)
        {
            if (pole == 0 || sense == 0)
                return d - two_pi * std::round(d / two_pi);
            return (pole > 0) == (sense > 0) ? d - two_pi * std::ceil(d / two_pi) : d - two_pi * std::floor(d / two_pi);
        }
    }

    /**
     * @brief Parameter box of the face of a surface, from points sampled on its boundary loops
     * (in the target units), enlarged by `margin` of its size (design question 4).
     *
     * - A B-spline or trimmed surface gives its own domain.
     * - In a periodic direction, a loop that winds around the axis makes the face cover a full
     *   turn: the box is [0, 2 pi], the seam of the STEP surface. Otherwise the box is the angular
     *   extent of the loops, which may cross 0 (e.g. [-0.3, 0.4]): such a face does not meet the seam.
     * - In other directions, the extent of the points, clamped to the natural domain (sphere
     *   latitude, cone apex).
     * Points on the axis or at a pole carry no angle and are ignored in u. Across a pole, the angle
     * of a loop steps by `detail::angle_step`: `sense` is +1 for a face on the left of its loops in
     * (u, v) (STEP: `same_sense`, the normal of the surface being S_u x S_v), -1 on the right, 0
     * unknown (the smallest step).
     */
    inline UVBox parameter_box(const SurfaceDef &s, std::span<const std::vector<Point<3>>> loops, Real margin = 0.1, int sense = 0)
    {
        const auto dom = s.domain();
        if (std::holds_alternative<SplineSurfaceDef>(s.def) || std::holds_alternative<TrimmedSurfaceDef>(s.def))
            return dom;
        const auto per = s.periodic();
        // axis of rotation, to drop points without angle
        std::optional<std::pair<Point<3>, Point<3>>> axis;
        std::visit([&](const auto &d) {
            using D = std::decay_t<decltype(d)>;
            if constexpr (std::is_same_v<D, CylinderDef> || std::is_same_v<D, ConeDef> || std::is_same_v<D, SphereDef> || std::is_same_v<D, TorusDef>)
                axis = std::pair{d.f.origin, d.f.z};
            else if constexpr (std::is_same_v<D, RevolutionDef>)
                axis = std::pair{d.origin, d.z};
        }, s.def);

        Real scale = 0.;
        for (const auto &l : loops)
            for (const auto &p : l)
                scale = std::max(scale, norm(p));
        const Real eps = 1e-9 * std::max<Real>(1., scale);

        std::array<Real, 2> lo{detail::inf, detail::inf}, hi{-detail::inf, -detail::inf};
        std::array<bool, 2> full{false, false};
        std::array<std::optional<Real>, 2> ref; // centre of the first loop in a periodic direction
        for (const auto &l : loops)
        {
            std::array<Real, 2> llo{detail::inf, detail::inf}, lhi{-detail::inf, -detail::inf};
            std::array<std::optional<Real>, 2> first, prev;
            std::array<Real, 2> unwrapped{};
            int pole = 0, first_pole = 0; // pole passed since the last angle, and before the first one
            for (const auto &p : l)
            {
                auto uv = s.parameters(p);
                if (!uv)
                    return dom;
                for (int k = 0; k < 2; ++k)
                {
                    if (k == 0 && per[0] && axis)
                    {
                        const auto q = p - axis->first;
                        if (norm(q - (q * axis->second) * axis->second) < eps)
                        {
                            pole = std::abs((*uv)[1] - dom.v2) < std::abs((*uv)[1] - dom.v1) ? 1 : -1;
                            continue; // on the axis: no angle
                        }
                    }
                    Real t = (*uv)[k];
                    if (per[k])
                    {
                        if (prev[k])
                            unwrapped[k] += detail::angle_step(t - *prev[k], k == 0 ? pole : 0, sense);
                        else
                        {
                            unwrapped[k] = t, first[k] = t;
                            if (k == 0)
                                first_pole = pole;
                        }
                        if (k == 0)
                            pole = 0;
                        prev[k] = t;
                        t = unwrapped[k];
                    }
                    llo[k] = std::min(llo[k], t);
                    lhi[k] = std::max(lhi[k], t);
                }
            }
            for (int k = 0; k < 2; ++k)
            {
                if (!(lhi[k] >= llo[k]))
                    continue;
                if (per[k])
                {
                    // closing step back to the first point: a winding loop has turned by 2 pi
                    const Real d = detail::angle_step(*first[k] - *prev[k], k == 0 ? (pole != 0 ? pole : first_pole) : 0, sense);
                    if (std::abs(unwrapped[k] + d - *first[k]) > std::numbers::pi)
                        full[k] = true;
                    // bring the loop next to the first one
                    const Real c = 0.5 * (llo[k] + lhi[k]);
                    if (ref[k])
                    {
                        const Real shift = detail::two_pi * std::round((*ref[k] - c) / detail::two_pi);
                        llo[k] += shift, lhi[k] += shift;
                    }
                    else
                        ref[k] = c;
                }
                lo[k] = std::min(lo[k], llo[k]);
                hi[k] = std::max(hi[k], lhi[k]);
            }
        }

        std::array<Real, 4> box{dom.u1, dom.u2, dom.v1, dom.v2};
        const Real size = std::max({hi[0] - lo[0], hi[1] - lo[1], Real(0)});
        for (int k = 0; k < 2; ++k)
        {
            Real &b1 = box[2 * k], &b2 = box[2 * k + 1];
            if (!(hi[k] >= lo[k])) // no usable point (a face reduced to a pole): keep the domain if finite
            {
                if (!std::isfinite(b1) || !std::isfinite(b2))
                    throw StepError(s.id, "cannot bound the surface: no boundary point");
                continue;
            }
            if (per[k])
            {
                // the turn whose centre is in [0, 2 pi): [-20, 30] degrees rather than [340, 390]
                const Real shift = detail::two_pi * std::floor(0.5 * (lo[k] + hi[k]) / detail::two_pi);
                lo[k] -= shift, hi[k] -= shift;
                const Real span = hi[k] - lo[k];
                const Real m = std::min(margin * span, 0.5 * (detail::two_pi - span));
                if (full[k] || m <= 0.)
                    b1 = 0., b2 = detail::two_pi;
                else
                    b1 = lo[k] - m, b2 = hi[k] + m;
                continue;
            }
            Real m = margin * (hi[k] - lo[k]);
            if (!(m > 0.))
                m = margin * size > 0. ? margin * size : Real(1);
            b1 = std::max(b1, lo[k] - m);
            b2 = std::min(b2, hi[k] + m);
        }
        return {box[0], box[1], box[2], box[3]};
    }

    // =========================================================================
    // Reader
    // =========================================================================

    namespace detail
    {
        // Insert the knot u once into a clamped or unclamped knot vector, for every row of poles (Boehm).
        template <typename P>
        void insert_knot(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p, Real u)
        {
            const auto span = static_cast<std::size_t>(std::upper_bound(k.begin(), k.end(), u) - k.begin()) - 1;
            const auto s = static_cast<std::size_t>(std::count(k.begin(), k.end(), u));
            for (auto &poles : rows)
            {
                std::vector<P> q;
                q.reserve(poles.size() + 1);
                for (std::size_t i = 0; i <= poles.size(); ++i)
                {
                    if (i + p <= span)
                        q.push_back(poles[i]);
                    else if (i + s > span)
                        q.push_back(poles[i - 1]);
                    else
                    {
                        const Real a = (u - k[i]) / (k[i + p] - k[i]);
                        q.push_back(a * poles[i] + (1. - a) * poles[i - 1]);
                    }
                }
                poles = std::move(q);
            }
            k.insert(k.begin() + static_cast<std::ptrdiff_t>(span) + 1, u);
        }

        // Clamp the start of a knot vector at its domain start k[p] (multiplicity p + 1).
        template <typename P>
        void clamp_start(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p)
        {
            const Real a = k[p];
            while (static_cast<std::size_t>(std::count(k.begin(), k.end(), a)) < p)
                insert_knot(k, rows, p, a);
            const auto j = static_cast<std::size_t>(std::lower_bound(k.begin(), k.end(), a) - k.begin());
            if (j == 0)
                return; // already clamped
            k.erase(k.begin(), k.begin() + static_cast<std::ptrdiff_t>(j));
            k.insert(k.begin(), a);
            for (auto &poles : rows)
                poles.erase(poles.begin(), poles.begin() + static_cast<std::ptrdiff_t>(j - 1));
        }

        // Clamp both ends of a (possibly unclamped, e.g. periodic) knot vector.
        template <typename P>
        void clamp(std::vector<Real> &k, std::vector<std::vector<P>> &rows, std::size_t p)
        {
            clamp_start(k, rows, p);
            std::ranges::reverse(k);
            for (auto &x : k)
                x = -x;
            for (auto &poles : rows)
                std::ranges::reverse(poles);
            clamp_start(k, rows, p);
            std::ranges::reverse(k);
            for (auto &x : k)
                x = -x;
            for (auto &poles : rows)
                std::ranges::reverse(poles);
        }

        // Degree of a B-spline with n poles: 1 at least, and n >= degree + 1 (checked before any
        // knot is generated or allocated).
        inline std::size_t checked_degree(std::int64_t degree, std::size_t n, std::uint64_t id)
        {
            if (degree < 1 || degree >= static_cast<std::int64_t>(n))
                throw StepError(id, "B-spline of degree " + std::to_string(degree) + " with " + std::to_string(n) + " poles");
            return static_cast<std::size_t>(degree);
        }

        // Flat knots of the B-spline subtypes defined by their degree p >= 1 and number of poles n >= p + 1.
        inline std::vector<Real> implicit_knots(std::string_view kind, std::size_t p, std::size_t n, std::uint64_t id)
        {
            if (p < 1 || n < p + 1)
                throw StepError(id, "B-spline with too few poles for its degree");
            std::vector<Real> k;
            if (kind == "BEZIER")
            {
                // n = p * segments + 1 ; knots 0..segments, ends p + 1 times, inner p times
                if (p == 0 || (n - 1) % p != 0)
                    throw StepError(id, "Bezier poles do not match the degree");
                const std::size_t segs = (n - 1) / p;
                for (std::size_t i = 0; i <= segs; ++i)
                    k.insert(k.end(), i == 0 || i == segs ? p + 1 : p, Real(i));
            }
            else if (kind == "QUASI_UNIFORM")
            {
                const std::size_t segs = n - p;
                for (std::size_t i = 0; i <= segs; ++i)
                    k.insert(k.end(), i == 0 || i == segs ? p + 1 : 1, Real(i));
            }
            else // UNIFORM: unclamped, k_i = i - p, clamped afterwards
                for (std::size_t i = 0; i < n + p + 1; ++i)
                    k.push_back(Real(i) - Real(p));
            return k;
        }

        // Flat knots of a B-spline of degree p with n poles from multiplicities (1 to p + 1 each,
        // n + p + 1 in all, checked before allocating) and distinct knots.
        inline std::vector<Real> flat_knots_of(ListView mults, ListView knots, std::size_t p, std::size_t n, std::uint64_t id)
        {
            if (mults.size() != knots.size())
                throw StepError(id, "knot multiplicities and knots differ in size");
            std::size_t total = 0;
            for (std::size_t i = 0; i < mults.size(); ++i)
            {
                const auto m = mults[i].as_int();
                if (m < 1 || m > static_cast<std::int64_t>(p) + 1)
                    throw StepError(id, "knot multiplicity " + std::to_string(m) + " out of [1, degree + 1]");
                total += static_cast<std::size_t>(m);
            }
            if (total != n + p + 1)
                throw StepError(id, "B-spline poles, knots and degree do not match");
            std::vector<Real> k;
            k.reserve(total);
            for (std::size_t i = 0; i < knots.size(); ++i)
                k.insert(k.end(), static_cast<std::size_t>(mults[i].as_int()), knots[i].as_real());
            return k;
        }

        inline bool clamped(const std::vector<Real> &k, std::size_t p)
        {
            return k.size() > 2 * p + 1 && k[0] == k[p] && k[k.size() - 1] == k[k.size() - 1 - p];
        }
    }

    /**
     * @brief Reads STEP geometric entities of a file into definitions, in target units.
     *
     * Results are cached by #id. Errors: StepUnsupported for a valid but unsupported entity,
     * StepError for invalid content, P21AccessError for arguments of the wrong kind.
     * 2D entities (pcurves of definitional representations) are read in the parameter space of
     * their STEP surface, without unit conversion.
     */
    class GeometryReader
    {
        const P21File &f_;
        Units u_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const CurveDef<2>>> curves2_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const CurveDef<3>>> curves3_;
        std::unordered_map<std::uint64_t, std::shared_ptr<const SurfaceDef>> surfaces_;
        std::unordered_set<std::uint64_t> reading_; // curves and surfaces being read

        // Marks an instance as being read while it lives: reading it again before the end is a cycle.
        class Reading
        {
            std::unordered_set<std::uint64_t> &set_;
            std::uint64_t id_;

        public:
            Reading(std::unordered_set<std::uint64_t> &set, std::uint64_t id) : set_{set}, id_{id}
            {
                if (!set_.insert(id).second)
                    throw StepError(id, "cyclic reference");
            }
            ~Reading() { set_.erase(id_); }
            Reading(const Reading &) = delete;
            Reading &operator=(const Reading &) = delete;
        };

        template <std::size_t dim>
        auto &cache()
        {
            if constexpr (dim == 2)
                return curves2_;
            else
                return curves3_;
        }

        // Own attributes of entity `type` in a simple instance of a subtype (after `skip`
        // inherited ones) or in the partial entity of a complex instance.
        static std::vector<ValueView> own(const InstanceView &in, std::string_view type, std::size_t skip, std::size_t count)
        {
            ListView a = in.is_complex() ? in.args(type) : in.args();
            const std::size_t first = in.is_complex() ? 0 : skip;
            if (a.size() < first + count)
                throw StepError(in.id(), "missing arguments of " + std::string(type));
            std::vector<ValueView> r;
            for (std::size_t i = 0; i < count; ++i)
                r.push_back(a[first + i]);
            return r;
        }

        // A simple instance of one of the types, or a complex one with one of them as partial type
        // (a simple instance does not list its supertypes).
        static bool is_one_of(const InstanceView &in, std::initializer_list<std::string_view> types)
        {
            return std::ranges::any_of(types, [&](std::string_view t) { return in.has_type(t); });
        }

        template <std::size_t dim>
        Point<dim> coords(std::uint64_t id, Real scale)
        {
            const auto in = f_.instance(id);
            const auto l = in.args()[1].as_list();
            if (l.size() < dim)
                throw StepError(id, "expected " + std::to_string(dim) + " coordinates");
            Point<dim> p{};
            for (std::size_t i = 0; i < dim; ++i)
                p[i] = l[i].as_real() * scale;
            return p;
        }

    public:
        GeometryReader(const P21File &f, Units u) : f_{f}, u_{std::move(u)} {}

        [[nodiscard]] const Units &units() const noexcept { return u_; }
        [[nodiscard]] const P21File &file() const noexcept { return f_; }

        /// Scale of a length in 3D, none in 2D (parameter space).
        template <std::size_t dim>
        [[nodiscard]] Real length_scale() const { return dim == 3 ? u_.length : 1.; }

        template <std::size_t dim>
        Point<dim> point(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            if (!in.has_type("CARTESIAN_POINT"))
                throw StepUnsupported(id, in.type());
            return coords<dim>(id, length_scale<dim>());
        }

        /// Unit direction (DIRECTION).
        template <std::size_t dim>
        Point<dim> direction(std::uint64_t id)
        {
            auto d = coords<dim>(id, 1.);
            const Real n = norm(d);
            if (!(n > 0.))
                throw StepError(id, "null direction");
            return d / n;
        }

        /// VECTOR(orientation, magnitude) as a vector in target units.
        template <std::size_t dim>
        Point<dim> vector(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            return (a[2].as_real() * length_scale<dim>()) * direction<dim>(a[1].as_ref());
        }

        /// AXIS2_PLACEMENT_3D: orthonormal frame, defaults as ISO 10303-42 build_axes (axis z,
        /// reference x, or y when the axis is along x).
        Frame frame(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            if (!in.has_type("AXIS2_PLACEMENT_3D"))
                throw StepUnsupported(id, in.type());
            const auto a = in.args();
            const auto o = point<3>(a[1].as_ref());
            const Point<3> z = a[2].is_null() ? Point<3>{0., 0., 1.} : direction<3>(a[2].as_ref());
            Point<3> r = a[3].is_null() ? Point<3>{1., 0., 0.} : direction<3>(a[3].as_ref());
            if (a[3].is_null() && norm(cross(z, r)) < 1e-12)
                r = {0., 1., 0.}; // first_proj_axis: (0, 1, 0) for an axis (1, 0, 0) or (-1, 0, 0)
            auto x = r - (r * z) * z;
            const Real n = norm(x);
            if (!(n > 1e-12))
                throw StepError(id, "reference direction parallel to the axis");
            x = x / n;
            return {o, x, cross(z, x), z};
        }

        /// AXIS2_PLACEMENT_2D: origin and orthonormal (x, y).
        std::array<Point<2>, 3> frame2d(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            const auto o = point<2>(a[1].as_ref());
            const Point<2> x = a[2].is_null() ? Point<2>{1., 0.} : direction<2>(a[2].as_ref());
            return {o, x, Point<2>{-x[1], x[0]}};
        }

        /// AXIS1_PLACEMENT: origin and unit axis.
        std::pair<Point<3>, Point<3>> axis1(std::uint64_t id)
        {
            const auto a = f_.instance(id).args();
            return {point<3>(a[1].as_ref()), a[2].is_null() ? Point<3>{0., 0., 1.} : direction<3>(a[2].as_ref())};
        }

        /// Curve of any supported type, 2D or 3D. A curve referring to itself, directly or not,
        /// throws StepError (cyclic reference).
        template <std::size_t dim>
        std::shared_ptr<const CurveDef<dim>> curve(std::uint64_t id)
        {
            auto &c = cache<dim>();
            if (auto it = c.find(id); it != c.end())
                return it->second;
            const Reading guard{reading_, id};
            auto r = std::make_shared<const CurveDef<dim>>(read_curve<dim>(id));
            c.emplace(id, r);
            return r;
        }

        /// Surface of any supported type. A surface referring to itself, directly or not, throws
        /// StepError (cyclic reference).
        std::shared_ptr<const SurfaceDef> surface(std::uint64_t id)
        {
            if (auto it = surfaces_.find(id); it != surfaces_.end())
                return it->second;
            const Reading guard{reading_, id};
            auto r = std::make_shared<const SurfaceDef>(read_surface(id));
            surfaces_.emplace(id, r);
            return r;
        }

    private:
        template <std::size_t dim>
        CurveDef<dim> read_curve(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            const auto a = in.args();
            const auto t = in.type();
            CurveDef<dim> c{.id = id};
            if (t == "LINE")
                c.def = LineDef<dim>{point<dim>(a[1].as_ref()), vector<dim>(a[2].as_ref())};
            else if (t == "CIRCLE" || t == "ELLIPSE")
            {
                const Real s = length_scale<dim>();
                const Real r1 = a[2].as_real() * s, r2 = (t == "CIRCLE" ? a[2] : a[3]).as_real() * s;
                if (!(r1 > 0.) || !(r2 > 0.))
                    throw StepError(id, "radius must be positive");
                if constexpr (dim == 3)
                {
                    const auto fr = frame(a[1].as_ref());
                    c.def = ConicDef<3>{fr.origin, fr.x, fr.y, r1, r2};
                }
                else
                {
                    const auto [o, x, y] = frame2d(a[1].as_ref());
                    c.def = ConicDef<2>{o, x, y, r1, r2};
                }
            }
            else if (is_one_of(in, {"B_SPLINE_CURVE", "B_SPLINE_CURVE_WITH_KNOTS", "BEZIER_CURVE", "QUASI_UNIFORM_CURVE", "UNIFORM_CURVE"}))
                c.def = SplineDef<dim>{read_bspline_curve<dim>(in)};
            else if (t == "POLYLINE")
            {
                points_vector<Real, dim> poles;
                for (auto p : a[1].as_list())
                    poles.push_back(point<dim>(p.as_ref()));
                if (poles.size() < 2)
                    throw StepError(id, "polyline with less than two points");
                std::vector<Real> k{0.};
                for (std::size_t i = 0; i < poles.size(); ++i)
                    k.push_back(Real(i));
                k.push_back(Real(poles.size() - 1));
                c.def = SplineDef<dim>{std::make_shared<BSCurve<Real, dim>>(poles, k, 1)};
            }
            else if (t == "TRIMMED_CURVE")
                c.def = read_trimmed<dim>(in);
            else if (t == "COMPOSITE_CURVE")
                c.def = SplineDef<dim>{read_composite<dim>(in)};
            else if (t == "SURFACE_CURVE" || t == "SEAM_CURVE" || t == "INTERSECTION_CURVE")
                return CurveDef<dim>{curve<dim>(a[1].as_ref())->def, id};
            else if (t == "PCURVE" && dim == 2)
                throw StepUnsupported(id, t); // read through its definitional representation
            else
                throw StepUnsupported(id, t);
            return c;
        }

        template <std::size_t dim>
        std::shared_ptr<const Curve<Real, dim>> read_bspline_curve(const InstanceView &in)
        {
            const auto id = in.id();
            const auto b = own(in, "B_SPLINE_CURVE", 1, 5);
            std::vector<Point<dim>> cps;
            for (auto r : b[1].as_list())
                cps.push_back(point<dim>(r.as_ref()));
            const auto n = cps.size();
            const auto p = detail::checked_degree(b[0].as_int(), n, id);
            std::vector<Real> k;
            if (in.has_type("B_SPLINE_CURVE_WITH_KNOTS"))
            {
                const auto kn = own(in, "B_SPLINE_CURVE_WITH_KNOTS", 6, 3);
                k = detail::flat_knots_of(kn[0].as_list(), kn[1].as_list(), p, n, id);
            }
            else if (in.has_type("BEZIER_CURVE"))
                k = detail::implicit_knots("BEZIER", p, n, id);
            else if (in.has_type("QUASI_UNIFORM_CURVE"))
                k = detail::implicit_knots("QUASI_UNIFORM", p, n, id);
            else if (in.has_type("UNIFORM_CURVE"))
                k = detail::implicit_knots("UNIFORM", p, n, id);
            else
                throw StepUnsupported(id, in.type());
            if (k.size() != n + p + 1)
                throw StepError(id, "B-spline poles, knots and degree do not match");

            std::vector<Real> w;
            if (auto rp = in.part("RATIONAL_B_SPLINE_CURVE"))
                for (auto x : rp->args()[0].as_list())
                    w.push_back(x.as_real());
            if (!w.empty() && w.size() != n)
                throw StepError(id, "weights and poles differ in size");
            if (!w.empty())
            {
                std::vector<std::vector<std::array<Real, dim + 1>>> rows(1);
                for (std::size_t i = 0; i < n; ++i)
                {
                    if (!(w[i] > 0.))
                        throw StepError(id, "weights must be positive");
                    std::array<Real, dim + 1> h;
                    for (std::size_t j = 0; j < dim; ++j)
                        h[j] = cps[i][j] * w[i]; // homogeneous coordinates
                    h[dim] = w[i];
                    rows[0].push_back(h);
                }
                if (!detail::clamped(k, p))
                    detail::clamp(k, rows, p);
                return std::make_shared<BSCurveRational<Real, dim>>(rows[0], k, p);
            }
            std::vector<std::vector<Point<dim>>> rows{cps};
            if (!detail::clamped(k, p))
                detail::clamp(k, rows, p);
            return std::make_shared<BSCurve<Real, dim>>(rows[0], k, p);
        }

        template <std::size_t dim>
        TrimmedDef<dim> read_trimmed(const InstanceView &in)
        {
            const auto id = in.id();
            const auto a = in.args();
            auto basis = curve<dim>(a[1].as_ref());
            if (std::holds_alternative<TrimmedDef<dim>>(basis->def))
                throw StepUnsupported(id, "TRIMMED_CURVE of a trimmed curve");
            const bool sense = a[4].as_bool();
            const bool prefer_point = a.size() > 5 && a[5].is(ValueKind::Enum) && a[5].as_enum() == "CARTESIAN";
            auto trim = [&](ValueView v) -> Real {
                std::optional<Real> param, by_point;
                for (auto s : v.as_list())
                {
                    if (s.is(ValueKind::Reference))
                        by_point = basis->parameter(point<dim>(s.as_ref()));
                    else if (s.is(ValueKind::Typed))
                    {
                        Real t = s.as_measure();
                        if (basis->kind() == ParamKind::Angle)
                            t *= u_.angle;
                        param = t;
                    }
                }
                if (by_point && (prefer_point || !param))
                    return *by_point;
                if (param)
                    return *param;
                throw StepError(id, "trimming select without parameter nor point");
            };
            Real t1 = trim(a[2]), t2 = trim(a[3]);
            if (basis->periodic())
            {
                // from t1 to t2 in the sense of the curve, at most a turn; trims equal up to whole
                // turns (within the angle precision of the file) give a full turn
                if (sense)
                    t2 = detail::turn_end(t1, t2);
                else
                    t1 = detail::turn_end(t2, t1);
            }
            if (!sense)
                std::swap(t1, t2);
            if (!(t2 > t1))
                throw StepError(id, "empty trimmed curve");
            return {basis, t1, t2, !sense};
        }

        template <std::size_t dim>
        std::shared_ptr<const Curve<Real, dim>> read_composite(const InstanceView &in)
        {
            std::optional<BSCurveRational<Real, dim>> joined;
            for (auto sref : in.args()[1].as_list())
            {
                const auto seg = f_.instance(sref.as_ref());
                if (!seg.has_type("COMPOSITE_CURVE_SEGMENT"))
                    throw StepUnsupported(seg.id(), seg.type());
                const auto sa = seg.args("COMPOSITE_CURVE_SEGMENT");
                const bool same = sa[1].as_bool();
                auto c = curve<dim>(sa[2].as_ref());
                const auto [t1, t2] = c->domain();
                if (!std::isfinite(t1) || !std::isfinite(t2))
                    throw StepError(seg.id(), "unbounded composite curve segment");
                auto nc = c->to_nurbs(t1, t2);
                BSCurveRational<Real, dim> r = [&] {
                    if (auto rr = std::dynamic_pointer_cast<BSCurveRational<Real, dim>>(nc.curve))
                        return *rr;
                    return BSCurveRational<Real, dim>{*std::dynamic_pointer_cast<BSCurve<Real, dim>>(nc.curve)};
                }();
                if (nc.u1 > r.bounds()[0] || nc.u2 < r.bounds()[1])
                    r.trim(nc.u1, nc.u2);
                if (same == c->reversed())
                    r.reverse();
                if (!joined)
                    joined = r;
                else
                {
                    if (norm(joined->end() - r.begin()) > 1e-6 * std::max<Real>(1., norm(r.begin())))
                        throw StepError(in.id(), "composite curve segments are not connected");
                    joined = BSCurveRational<Real, dim>{join(*joined, r)};
                }
            }
            if (!joined)
                throw StepError(in.id(), "empty composite curve");
            return std::make_shared<BSCurveRational<Real, dim>>(*joined);
        }

        SurfaceDef read_surface(std::uint64_t id)
        {
            const auto in = f_.instance(id);
            const auto a = in.args();
            const auto t = in.type();
            const Real L = u_.length;
            SurfaceDef s{.id = id};
            auto positive = [&](Real r) {
                if (!(r > 0.))
                    throw StepError(id, "radius must be positive");
                return r;
            };
            if (t == "PLANE")
                s.def = PlaneDef{frame(a[1].as_ref())};
            else if (t == "CYLINDRICAL_SURFACE")
                s.def = CylinderDef{frame(a[1].as_ref()), positive(a[2].as_real() * L)};
            else if (t == "CONICAL_SURFACE")
            {
                const Real r = a[2].as_real() * L, alpha = a[3].as_measure() * u_.angle;
                if (r < 0. || !(std::abs(alpha) < 0.5 * std::numbers::pi) || alpha == 0.)
                    throw StepError(id, "invalid cone radius or semi angle");
                s.def = ConeDef{frame(a[1].as_ref()), r, alpha};
            }
            else if (t == "SPHERICAL_SURFACE")
                s.def = SphereDef{frame(a[1].as_ref()), positive(a[2].as_real() * L)};
            else if (t == "TOROIDAL_SURFACE" || t == "DEGENERATE_TOROIDAL_SURFACE")
            {
                TorusDef tor{frame(a[1].as_ref()), a[2].as_real() * L, positive(a[3].as_real() * L)};
                if (!(tor.major >= 0.))
                    throw StepError(id, "negative major radius");
                // minor > major: self-intersecting, the part chosen by select_outer
                if (t == "DEGENERATE_TOROIDAL_SURFACE" && tor.minor > tor.major)
                    tor.part = a[4].as_bool() ? TorusPart::Apple : TorusPart::Lemon;
                s.def = tor;
            }
            else if (is_one_of(in, {"B_SPLINE_SURFACE", "B_SPLINE_SURFACE_WITH_KNOTS", "BEZIER_SURFACE", "QUASI_UNIFORM_SURFACE", "UNIFORM_SURFACE"}))
                s.def = SplineSurfaceDef{read_bspline_surface(in)};
            else if (t == "SURFACE_OF_REVOLUTION")
            {
                auto g = curve<3>(a[1].as_ref());
                const auto [o, z] = axis1(a[2].as_ref());
                // a generatrix on the axis sweeps no surface
                const auto [g1, g2] = g->domain();
                const Real lo = std::isfinite(g1) ? g1 : (std::isfinite(g2) ? g2 - 1. : 0.);
                const Real hi = std::isfinite(g2) ? g2 : lo + 1.;
                bool off = false;
                for (int i = 0; i <= 16 && !off; ++i)
                {
                    const auto q = g->value(lo + (hi - lo) * i / 16.) - o;
                    off = norm(q - (q * z) * z) > 1e-9 * std::max<Real>(1., norm(q));
                }
                if (!off)
                    throw StepError(id, "generatrix on the axis of revolution");
                s.def = RevolutionDef{g, o, z};
            }
            else if (t == "SURFACE_OF_LINEAR_EXTRUSION")
            {
                const auto v = vector<3>(a[2].as_ref());
                if (!(norm(v) > 0.))
                    throw StepError(id, "null extrusion vector");
                s.def = ExtrusionDef{curve<3>(a[1].as_ref()), v};
            }
            else if (t == "RECTANGULAR_TRIMMED_SURFACE")
            {
                auto basis = surface(a[1].as_ref());
                const auto k = basis->kinds();
                const auto per = basis->periodic();
                auto conv = [&](ValueView v, ParamKind kind) {
                    const Real x = v.as_measure();
                    return kind == ParamKind::Angle ? x * u_.angle : kind == ParamKind::Length ? x * L : x;
                };
                // ISO 10303-42: in a periodic direction the range runs from the first bound to the
                // second in the sense given (usense, vsense), across the seam if need be, a turn at
                // most (bounds equal up to whole turns: a full turn); in another one they are sorted
                auto range = [](Real &b1, Real &b2, bool sense, bool periodic) {
                    if (!periodic)
                    {
                        if (b1 > b2)
                            std::swap(b1, b2);
                    }
                    else if (sense)
                        b2 = detail::turn_end(b1, b2);
                    else
                    {
                        const Real end = detail::turn_end(b2, b1); // decreasing from b1 to b2
                        b1 = b2, b2 = end;
                    }
                };
                Real u1 = conv(a[2], k[0]), u2 = conv(a[3], k[0]), v1 = conv(a[4], k[1]), v2 = conv(a[5], k[1]);
                range(u1, u2, a[6].as_bool(), per[0]);
                range(v1, v2, a[7].as_bool(), per[1]);
                s.def = TrimmedSurfaceDef{basis, u1, u2, v1, v2};
            }
            else if (t == "OFFSET_SURFACE")
                s.def = read_offset(in);
            else
                throw StepUnsupported(id, t);
            return s;
        }

        // Offset of an elementary surface, S + d N (N = Su x Sv / |Su x Sv|): the elementary surface of
        // offset radius, with the same parameters (u, v).
        decltype(SurfaceDef::def) read_offset(const InstanceView &in)
        {
            const auto a = in.args();
            const auto basis = surface(a[1].as_ref());
            const Real d = a[2].as_measure() * u_.length;
            return std::visit([&](const auto &b) -> decltype(SurfaceDef::def) {
                using D = std::decay_t<decltype(b)>;
                auto positive = [&](Real r) {
                    if (!(r > 0.))
                        throw StepError(in.id(), "offset gives a null radius");
                    return r;
                };
                if constexpr (std::is_same_v<D, PlaneDef>)
                    return PlaneDef{{b.f.origin + d * b.f.z, b.f.x, b.f.y, b.f.z}};
                else if constexpr (std::is_same_v<D, CylinderDef>)
                    return CylinderDef{b.f, positive(b.r + d)};
                else if constexpr (std::is_same_v<D, SphereDef>)
                    return SphereDef{b.f, positive(b.r + d)};
                else if constexpr (std::is_same_v<D, TorusDef>)
                {
                    // the normal of a degenerate torus flips on the lemon (R + r cos v < 0), and the
                    // range of v of a part is that of the basis: not a part of the offset torus
                    if (b.part != TorusPart::Whole)
                        throw StepUnsupported(in.id(), "OFFSET_SURFACE of a degenerate torus");
                    return TorusDef{b.f, b.major, positive(b.minor + d)};
                }
                else if constexpr (std::is_same_v<D, ConeDef>)
                {
                    // N = cos a X(u) - sin a Z: S + d N = (O - d sin a Z) + (R + d cos a + v tan a) X(u) + v Z,
                    // the cone of the same half-angle, origin moved by -d sin a along the axis,
                    // radius R + d cos a
                    const Real r = b.r + d * std::cos(b.semi_angle);
                    if (r < 0.) // a negative reference radius: no such cone with the same parameters
                        throw StepUnsupported(in.id(), "OFFSET_SURFACE of a cone across its apex");
                    Frame f = b.f;
                    f.origin = f.origin - (d * std::sin(b.semi_angle)) * f.z;
                    return ConeDef{f, r, b.semi_angle};
                }
                else
                    throw StepUnsupported(in.id(), "OFFSET_SURFACE of a non elementary surface");
            }, basis->def);
        }

        std::shared_ptr<const Surface<Real, 3>> read_bspline_surface(const InstanceView &in)
        {
            const auto id = in.id();
            const auto b = own(in, "B_SPLINE_SURFACE", 1, 7);
            const auto grid = b[2].as_list();
            const std::size_t nu = grid.size();
            const std::size_t nv = nu ? grid[0].as_list().size() : 0;
            if (nu == 0 || nv == 0)
                throw StepError(id, "empty control point grid");
            const auto pu = detail::checked_degree(b[0].as_int(), nu, id), pv = detail::checked_degree(b[1].as_int(), nv, id);
            std::vector<Real> ku, kv;
            if (in.has_type("B_SPLINE_SURFACE_WITH_KNOTS"))
            {
                const auto kn = own(in, "B_SPLINE_SURFACE_WITH_KNOTS", 8, 5);
                ku = detail::flat_knots_of(kn[0].as_list(), kn[2].as_list(), pu, nu, id);
                kv = detail::flat_knots_of(kn[1].as_list(), kn[3].as_list(), pv, nv, id);
            }
            else
            {
                std::string_view kind = in.has_type("BEZIER_SURFACE")          ? "BEZIER"
                                        : in.has_type("QUASI_UNIFORM_SURFACE") ? "QUASI_UNIFORM"
                                        : in.has_type("UNIFORM_SURFACE")       ? "UNIFORM"
                                                                               : "";
                if (kind.empty())
                    throw StepUnsupported(id, in.type());
                ku = detail::implicit_knots(kind, pu, nu, id);
                kv = detail::implicit_knots(kind, pv, nv, id);
            }
            if (ku.size() != nu + pu + 1 || kv.size() != nv + pv + 1)
                throw StepError(id, "B-spline surface poles, knots and degrees do not match");
            std::optional<ListView> weights;
            if (auto rp = in.part("RATIONAL_B_SPLINE_SURFACE"))
                weights = rp->args()[0].as_list();

            // rows along u (one per v), homogeneous 4D poles; STEP stores [u][v]
            std::vector<std::vector<std::array<Real, 4>>> rows(nv, std::vector<std::array<Real, 4>>(nu));
            for (std::size_t i = 0; i < nu; ++i)
            {
                const auto col = grid[i].as_list();
                if (col.size() != nv)
                    throw StepError(id, "ragged control point grid");
                for (std::size_t j = 0; j < nv; ++j)
                {
                    const auto p = point<3>(col[j].as_ref());
                    const Real w = weights ? (*weights)[i].as_list()[j].as_real() : 1.;
                    if (!(w > 0.))
                        throw StepError(id, "weights must be positive");
                    rows[j][i] = {p[0] * w, p[1] * w, p[2] * w, w};
                }
            }
            if (!detail::clamped(ku, pu))
                detail::clamp(ku, rows, pu);
            if (!detail::clamped(kv, pv))
            {
                // transpose: rows along v (one per u)
                std::vector<std::vector<std::array<Real, 4>>> cols(rows[0].size(), std::vector<std::array<Real, 4>>(rows.size()));
                for (std::size_t j = 0; j < rows.size(); ++j)
                    for (std::size_t i = 0; i < rows[j].size(); ++i)
                        cols[i][j] = rows[j][i];
                detail::clamp(kv, cols, pv);
                rows.assign(cols[0].size(), std::vector<std::array<Real, 4>>(cols.size()));
                for (std::size_t i = 0; i < cols.size(); ++i)
                    for (std::size_t j = 0; j < cols[i].size(); ++j)
                        rows[j][i] = cols[i][j];
            }
            std::vector<std::array<Real, 4>> hp; // u fastest
            for (const auto &r : rows)
                hp.insert(hp.end(), r.begin(), r.end());
            if (weights)
                return std::make_shared<BSSurfaceRational<Real, 3>>(hp, ku, kv, pu, pv);
            points_vector<Real, 3> poles;
            poles.reserve(hp.size());
            for (const auto &h : hp)
                poles.push_back({h[0], h[1], h[2]});
            return std::make_shared<BSSurface<Real, 3>>(poles, ku, kv, pu, pv);
        }
    };
}
