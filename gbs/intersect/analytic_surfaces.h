#pragma once
/**
 * @file analytic_surfaces.h
 * @brief Closed-form intersection of two primitive surfaces: the carrier curves (lines, circles,
 * ellipses, unbounded or full turns) and the isolated contact points, before their restriction
 * to the parameter boxes of the surfaces (`surface_surface.h`).
 *
 * Coaxial pairs (surfaces of revolution sharing an axis, a sphere centred on the axis, a plane
 * normal to it, two spheres) intersect along circles: their meridians are lines and circles of
 * the half-plane (rho >= 0, h), intersected in 2D. The other handled pairs: two planes (a line),
 * a plane and a cylinder or a cone (conic of the plane: ellipse, or lines), a plane containing
 * the axis of a torus (two meridian circles), two cylinders of parallel axes (lines). Anything
 * else (parabolas, hyperbolas, skew cylinders, a torus cut obliquely…) has no closed form here.
 */
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    namespace intersect_detail
    {
        /// A curve carrying a part of the intersection, before restriction to the boxes.
        template <std::floating_point T>
        struct Carrier
        {
            std::shared_ptr<ElementaryCurve<T>> curve; ///< line on a placeholder range, or a full turn
            bool tangent{};                            ///< the surfaces are tangent along it
        };

        template <std::floating_point T>
        struct AnalyticIntersection
        {
            std::vector<Carrier<T>> carriers;
            std::vector<point<T, 3>> points;
            bool coincident{};
        };

        template <std::floating_point T>
        auto unit(const point<T, 3> &v) -> point<T, 3> { return v / norm(v); }

        /// A unit vector orthogonal to z.
        template <std::floating_point T>
        auto orthogonal(const point<T, 3> &z) -> point<T, 3>
        {
            const point<T, 3> a = std::abs(z[0]) < T(0.9) ? point<T, 3>{1, 0, 0} : point<T, 3>{0, 1, 0};
            return unit<T>(a - (a * z) * z);
        }

        template <std::floating_point T>
        auto frame_of(const point<T, 3> &o, const point<T, 3> &z, const point<T, 3> &x) -> Frame3<T>
        {
            Frame3<T> f;
            f.o = o;
            f.z = z;
            f.x = unit<T>(x - (x * z) * z);
            f.y = cross(f.z, f.x);
            return f;
        }

        // ---- meridians: lines a rho + b h + c = 0 and circles of centre (rho, h), in the half-plane rho >= 0

        template <std::floating_point T>
        struct Meridian
        {
            bool circle{};
            T a{}, b{}, c{}; ///< line
            T cr{}, ch{}, r{}; ///< circle
        };

        /// 2D intersection points of two meridian elements (rho, h, tangent); `same` when they coincide.
        template <std::floating_point T>
        auto meet(const Meridian<T> &m1, const Meridian<T> &m2, T tol, bool &same) -> std::vector<std::array<T, 3>>
        {
            std::vector<std::array<T, 3>> pts;
            auto line_circle = [&](const Meridian<T> &l, const Meridian<T> &c) {
                const T n = std::hypot(l.a, l.b);
                const T a = l.a / n, b = l.b / n, cc = l.c / n;
                const T d = a * c.cr + b * c.ch + cc; // signed distance of the centre
                const T fr = c.cr - d * a, fh = c.ch - d * b; // foot of the centre on the line
                if (std::abs(std::abs(d) - c.r) <= tol)
                    pts.push_back({fr, fh, T(1)});
                else if (std::abs(d) < c.r)
                {
                    const T w = std::sqrt(c.r * c.r - d * d);
                    pts.push_back({fr - w * b, fh + w * a, T(0)});
                    pts.push_back({fr + w * b, fh - w * a, T(0)});
                }
            };
            if (!m1.circle && !m2.circle)
            {
                const T det = m1.a * m2.b - m1.b * m2.a;
                const T n1 = std::hypot(m1.a, m1.b), n2 = std::hypot(m2.a, m2.b);
                if (std::abs(det) <= T(1e-12) * n1 * n2)
                {
                    // parallel: the same line when the distance between them vanishes
                    const T d = std::abs(m1.c / n1 - (m1.a * m2.a + m1.b * m2.b > 0 ? T(1) : T(-1)) * m2.c / n2);
                    same = same || d <= tol;
                    return pts;
                }
                pts.push_back({(m1.b * m2.c - m2.b * m1.c) / det, (m2.a * m1.c - m1.a * m2.c) / det, T(0)});
            }
            else if (m1.circle != m2.circle)
                line_circle(m1.circle ? m2 : m1, m1.circle ? m1 : m2);
            else
            {
                const T dr = m2.cr - m1.cr, dh = m2.ch - m1.ch, d = std::hypot(dr, dh);
                if (d <= tol)
                {
                    same = same || std::abs(m1.r - m2.r) <= tol;
                    return pts;
                }
                const T x = (d * d + m1.r * m1.r - m2.r * m2.r) / (T(2) * d); // along the centres
                const T y2 = m1.r * m1.r - x * x;
                const T er = dr / d, eh = dh / d;
                const bool tangent = std::abs(d - (m1.r + m2.r)) <= tol || std::abs(d - std::abs(m1.r - m2.r)) <= tol;
                if (tangent)
                    pts.push_back({m1.cr + x * er, m1.ch + x * eh, T(1)});
                else if (y2 > T(0))
                {
                    const T y = std::sqrt(y2);
                    pts.push_back({m1.cr + x * er - y * eh, m1.ch + x * eh + y * er, T(0)});
                    pts.push_back({m1.cr + x * er + y * eh, m1.ch + x * eh - y * er, T(0)});
                }
            }
            return pts;
        }

        /// Is the surface a surface of revolution with its own axis (cylinder, cone, torus)?
        template <std::floating_point T>
        bool has_axis(const ElementarySurface<T> &s)
        {
            return dynamic_cast<const CylindricalSurface<T> *>(&s) || dynamic_cast<const ConicalSurface<T> *>(&s) ||
                   dynamic_cast<const ToroidalSurface<T> *>(&s);
        }

        /// Meridian of a surface around the axis (o, z); empty when the surface is not of revolution around it.
        template <std::floating_point T>
        auto meridian(const ElementarySurface<T> &s, const point<T, 3> &o, const point<T, 3> &z) -> std::vector<Meridian<T>>
        {
            const auto &f = s.frame();
            const T hc = (f.o - o) * z;
            if (dynamic_cast<const Plane<T> *>(&s))
                return {Meridian<T>{.a = T(0), .b = T(1), .c = -hc}};
            if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(&s))
                return {Meridian<T>{.circle = true, .cr = T(0), .ch = hc, .r = sp->radius()}};
            if (const auto *cy = dynamic_cast<const CylindricalSurface<T> *>(&s))
                return {Meridian<T>{.a = T(1), .b = T(0), .c = -cy->radius()}};
            if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(&s))
                return {Meridian<T>{.circle = true, .cr = to->major_radius(), .ch = hc, .r = to->minor_radius()}};
            if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(&s))
            {
                // rho = +-(R + k s (h - hc)), s = +-1 the orientation of the cone axis on z
                const T k = std::tan(co->semi_angle()) * (f.z * z > 0 ? T(1) : T(-1)), R = co->radius();
                return {Meridian<T>{.a = T(1), .b = -k, .c = -R + k * hc}, Meridian<T>{.a = T(1), .b = k, .c = R - k * hc}};
            }
            return {};
        }

        /// The common axis of two surfaces, when they are coaxial: origin, direction, reference x.
        template <std::floating_point T>
        auto common_axis(const ElementarySurface<T> &s1, const ElementarySurface<T> &s2, T tol)
            -> std::optional<std::array<point<T, 3>, 3>>
        {
            const bool p1 = dynamic_cast<const Plane<T> *>(&s1) != nullptr, p2 = dynamic_cast<const Plane<T> *>(&s2) != nullptr;
            const bool sp1 = dynamic_cast<const SphericalSurface<T> *>(&s1) != nullptr, sp2 = dynamic_cast<const SphericalSurface<T> *>(&s2) != nullptr;
            const auto &f1 = s1.frame(), &f2 = s2.frame();
            const T ang = T(1e-10);
            auto on_axis = [&](const point<T, 3> &p, const point<T, 3> &o, const point<T, 3> &z) {
                const auto q = p - o;
                return norm(q - (q * z) * z) <= tol;
            };
            if (has_axis(s1) || has_axis(s2))
            {
                const auto &fa = has_axis(s1) ? f1 : f2;
                const auto &other = has_axis(s1) ? s2 : s1;
                const auto &fo = other.frame();
                const bool o_plane = has_axis(s1) ? p2 : p1, o_sphere = has_axis(s1) ? sp2 : sp1;
                if (o_plane)
                {
                    if (norm(cross(fa.z, fo.z)) <= ang)
                        return std::array<point<T, 3>, 3>{fa.o, fa.z, fa.x};
                    return std::nullopt;
                }
                if (o_sphere)
                {
                    if (on_axis(fo.o, fa.o, fa.z))
                        return std::array<point<T, 3>, 3>{fa.o, fa.z, fa.x};
                    return std::nullopt;
                }
                if (norm(cross(fa.z, fo.z)) <= ang && on_axis(fo.o, fa.o, fa.z))
                    return std::array<point<T, 3>, 3>{fa.o, fa.z, fa.x};
                return std::nullopt;
            }
            if (sp1 && sp2)
            {
                const auto d = f2.o - f1.o;
                if (norm(d) <= tol)
                    return std::array<point<T, 3>, 3>{f1.o, f1.z, f1.x}; // concentric
                const auto z = unit<T>(d);
                return std::array<point<T, 3>, 3>{f1.o, z, orthogonal<T>(z)};
            }
            if (sp1 != sp2 && (p1 || p2)) // sphere and plane: the axis through the centre, normal to the plane
            {
                const auto &fs = sp1 ? f1 : f2, &fp = sp1 ? f2 : f1;
                return std::array<point<T, 3>, 3>{fs.o, fp.z, fp.x};
            }
            return std::nullopt;
        }

        template <std::floating_point T>
        auto full_circle(T r, const Frame3<T> &f, T t0) -> std::shared_ptr<ElementaryCurve<T>>
        {
            return std::make_shared<EllipseCurve<T>>(r, r, f, t0, t0 + T(2) * std::numbers::pi_v<T>);
        }

        template <std::floating_point T>
        auto line(const point<T, 3> &o, const point<T, 3> &d) -> std::shared_ptr<ElementaryCurve<T>>
        {
            return std::make_shared<LineCurve<T>>(o, d, T(-1), T(1));
        }

        /// Plane and cylinder or cone, not coaxial: the conic of the plane (ellipse or lines).
        template <std::floating_point T>
        auto plane_quadric(const Plane<T> &pl, const ElementarySurface<T> &q, T tol) -> std::optional<AnalyticIntersection<T>>
        {
            const auto &fq = q.frame(), &fp = pl.frame();
            T R{}, k{};
            if (const auto *cy = dynamic_cast<const CylindricalSurface<T> *>(&q))
                R = cy->radius();
            else if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(&q))
            {
                R = co->radius();
                k = std::tan(co->semi_angle());
            }
            else
                return std::nullopt;
            AnalyticIntersection<T> r;
            const auto n = fp.z;
            auto in_q = [&](const point<T, 3> &v) { return point<T, 3>{v * fq.x, v * fq.y, v * fq.z}; };
            if (k == T(0) && std::abs(n * fq.z) <= T(1e-10))
            {
                // plane parallel to the axis: generatrices at the distance of the plane to the axis
                const T d = (fq.o - fp.o) * n; // axis point above the plane
                const auto foot = fq.o - d * n;
                const auto w = unit<T>(cross(n, fq.z));
                if (std::abs(std::abs(d) - R) <= tol)
                    r.carriers.push_back({line<T>(foot, fq.z), true});
                else if (std::abs(d) < R)
                {
                    const T h = std::sqrt(R * R - d * d);
                    r.carriers.push_back({line<T>(foot + h * w, fq.z), false});
                    r.carriers.push_back({line<T>(foot - h * w, fq.z), false});
                }
                return r;
            }
            // the plane P(x, y) = c + x E1 + y E2 in the frame of the quadric; g = X² + Y² - (R + k Z)²
            const auto c = in_q(fp.o - fq.o), E1 = in_q(fp.x), E2 = in_q(fp.y);
            auto quad = [&](const point<T, 3> &a, const point<T, 3> &b) { return a[0] * b[0] + a[1] * b[1] - k * k * a[2] * b[2]; };
            const T A = quad(E1, E1), B = T(2) * quad(E1, E2), C = quad(E2, E2);
            const T D = T(2) * quad(c, E1) - T(2) * R * k * E1[2], E = T(2) * quad(c, E2) - T(2) * R * k * E2[2];
            const T F = quad(c, c) - T(2) * R * k * c[2] - R * R;
            const T disc = B * B - T(4) * A * C;
            const T scale = std::max({std::abs(A), std::abs(B), std::abs(C)});
            auto at = [&](T x, T y) { return fp.o + x * fp.x + y * fp.y; };
            if (k != T(0))
            {
                // plane through the apex: generatrices (or the apex alone)
                const auto apex = fq.o - (R / k) * fq.z;
                if (std::abs((apex - fp.o) * n) <= tol)
                {
                    if (disc < -T(1e-12) * scale * scale)
                        r.points.push_back(apex);
                    else if (disc <= T(1e-12) * scale * scale)
                    {
                        const auto d = std::abs(A) >= std::abs(C) ? point<T, 3>{-B / (T(2) * A), T(1), T(0)} : point<T, 3>{T(1), -B / (T(2) * C), T(0)};
                        r.carriers.push_back({line<T>(apex, unit<T>(d[0] * fp.x + d[1] * fp.y)), true});
                    }
                    else
                    {
                        // directions (x, y) with A x² + B x y + C y² = 0
                        const T s = std::sqrt(disc);
                        std::vector<point<T, 3>> dirs;
                        if (std::abs(A) >= std::abs(C))
                            for (T sg : {T(-1), T(1)})
                                dirs.push_back(unit<T>(((-B + sg * s) / (T(2) * A)) * fp.x + fp.y));
                        else
                            for (T sg : {T(-1), T(1)})
                                dirs.push_back(unit<T>(fp.x + ((-B + sg * s) / (T(2) * C)) * fp.y));
                        for (const auto &d : dirs)
                            r.carriers.push_back({line<T>(apex, d), false});
                    }
                    return r;
                }
            }
            if (disc >= -T(1e-12) * scale * scale)
                return std::nullopt; // parabola or hyperbola
            // ellipse: centre, axes from the eigenvectors of [[A, B/2], [B/2, C]]
            const T det = T(4) * A * C - B * B;
            const T x0 = (B * E - T(2) * C * D) / det, y0 = (B * D - T(2) * A * E) / det;
            const T F0 = F + (D * x0 + E * y0) / T(2);
            const T th = T(0.5) * std::atan2(B, A - C);
            const T cs = std::cos(th), sn = std::sin(th);
            const T l1 = A * cs * cs + B * cs * sn + C * sn * sn, l2 = A * sn * sn - B * cs * sn + C * cs * cs;
            const T r1s = -F0 / l1, r2s = -F0 / l2;
            const auto centre = at(x0, y0);
            if (r1s <= tol * tol || r2s <= tol * tol)
            {
                if (std::max(std::abs(r1s), std::abs(r2s)) <= tol * tol)
                    r.points.push_back(centre); // the plane touches the quadric at one point
                return r;
            }
            Frame3<T> fe;
            fe.o = centre;
            fe.z = n;
            fe.x = cs * fp.x + sn * fp.y;
            fe.y = cross(fe.z, fe.x);
            r.carriers.push_back({std::make_shared<EllipseCurve<T>>(std::sqrt(r1s), std::sqrt(r2s), fe, T(0),
                                                                    T(2) * std::numbers::pi_v<T>),
                                  false});
            return r;
        }

        /**
         * @brief Closed-form intersection of two primitives, or nullopt when the pair has none here.
         */
        template <std::floating_point T>
        auto analytic_intersection(const ElementarySurface<T> &s1, const ElementarySurface<T> &s2, T tol) -> std::optional<AnalyticIntersection<T>>
        {
            const auto &f1 = s1.frame(), &f2 = s2.frame();
            const auto *pl1 = dynamic_cast<const Plane<T> *>(&s1), *pl2 = dynamic_cast<const Plane<T> *>(&s2);
            AnalyticIntersection<T> r;
            if (pl1 && pl2)
            {
                const auto n1 = f1.z, n2 = f2.z;
                const auto d = cross(n1, n2);
                if (norm(d) <= T(1e-10))
                {
                    r.coincident = std::abs((f2.o - f1.o) * n1) <= tol;
                    return r;
                }
                const T h1 = n1 * f1.o, h2 = n2 * f2.o, c12 = n1 * n2, det = T(1) - c12 * c12;
                const auto p = ((h1 - h2 * c12) / det) * n1 + ((h2 - h1 * c12) / det) * n2;
                r.carriers.push_back({line<T>(p, unit<T>(d)), false});
                return r;
            }
            if (const auto ax = common_axis(s1, s2, tol))
            {
                const auto &[o, z, x] = *ax;
                const auto fa = frame_of(o, z, x);
                const auto m1 = meridian(s1, o, z), m2 = meridian(s2, o, z);
                bool same = false;
                // the parameter of the circles is the angle u of the surface giving the axis
                const auto &fs = has_axis(s1) ? s1 : s2;
                const T t0 = has_axis(s1) || has_axis(s2) ? fs.bounds()[0] : T(0);
                std::vector<std::array<T, 3>> pts;
                for (const auto &a : m1)
                    for (const auto &b : m2)
                        for (const auto &p : meet(a, b, tol, same))
                            if (p[0] >= -tol)
                                pts.push_back(p);
                if (same)
                {
                    r.coincident = true;
                    return r;
                }
                for (const auto &[rho, h, tg] : pts)
                {
                    Frame3<T> fc = fa;
                    fc.o = o + h * z;
                    if (rho <= tol)
                        r.points.push_back(fc.o);
                    else
                        r.carriers.push_back({full_circle<T>(rho, fc, t0), tg != T(0)});
                }
                // a cone meets itself on its axis: duplicates of the apex
                std::ranges::sort(r.points, [](const auto &a, const auto &b) { return a < b; });
                r.points.erase(std::unique(r.points.begin(), r.points.end(), [&](const auto &a, const auto &b) { return norm(a - b) <= tol; }),
                               r.points.end());
                return r;
            }
            if (pl1 || pl2)
            {
                const auto &pl = pl1 ? *pl1 : *pl2;
                const auto &q = pl1 ? s2 : s1;
                if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(&q))
                {
                    // plane containing the axis: the meridian circles at u and u + pi
                    const auto &ft = to->frame(), &fp = pl.frame();
                    if (std::abs(fp.z * ft.z) > T(1e-10) || std::abs((ft.o - fp.o) * fp.z) > tol)
                        return std::nullopt;
                    const auto w = unit<T>(cross(fp.z, ft.z)); // in the plane, normal to the axis
                    const T u = std::atan2(w * ft.y, w * ft.x);
                    for (T a : {u, u + std::numbers::pi_v<T>})
                    {
                        auto c = to->iso_u(a, T(0), T(2) * std::numbers::pi_v<T>);
                        r.carriers.push_back({c, false});
                    }
                    return r;
                }
                return plane_quadric(pl, q, tol);
            }
            const auto *c1 = dynamic_cast<const CylindricalSurface<T> *>(&s1), *c2 = dynamic_cast<const CylindricalSurface<T> *>(&s2);
            if (c1 && c2 && norm(cross(f1.z, f2.z)) <= T(1e-10))
            {
                // parallel axes: the circles of the cross-section of the first cylinder
                const auto q = f2.o - f1.o;
                Meridian<T> a{.circle = true, .cr = T(0), .ch = T(0), .r = c1->radius()};
                Meridian<T> b{.circle = true, .cr = q * f1.x, .ch = q * f1.y, .r = c2->radius()};
                bool same = false;
                for (const auto &[x, y, tg] : meet(a, b, tol, same))
                    r.carriers.push_back({line<T>(f1.o + x * f1.x + y * f1.y, f1.z), tg != T(0)});
                r.coincident = same;
                return r;
            }
            return std::nullopt;
        }
    }
}
