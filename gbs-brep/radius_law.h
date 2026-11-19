#pragma once
/**
 * @file radius_law.h
 * @brief The law of the radius of a fillet along an edge or a chain of tangent edges (stage 5,
 * `docs/sources/design/brep_palier5_variable_fillets.md`, § 2): control points (an edge, a
 * fraction of its length from its first vertex, a radius), and between them a linear or a smooth
 * variation along the arc length.
 *
 * The law is a B-spline function of the arc length (`BSCfunction`):
 * - **linear**: degree 1 through the points (its slope jumps at them);
 * - **smooth**: the cubic Hermite through the points, the slopes of Fritsch and Carlson (C1,
 *   monotone between two points: never beyond the radii on either side), written as a cubic
 *   B-spline with double knots; zero slopes at the ends of an open law, so that it goes on into
 *   the constant radius beyond its first and last points; a closed law with its slope across its
 *   closure.
 *
 * Before the first point and after the last one, the law goes on with its slope there: constant
 * for a smooth law (its slope 0 there), straight on for a linear one (no kink there).
 *
 * Source: Fritsch and Carlson, "Monotone piecewise cubic interpolation" (SIAM Journal on
 * Numerical Analysis 17(2), 1980), with the weighted harmonic mean of Fritsch and Butland for the
 * slopes.
 */
#include <gbs-brep/model.h>
#include <gbs/bscurve.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gbs::brep
{
    /// A radius at a point of an edge: `at` the fraction of its length from its first vertex (0 or 1: at a vertex).
    template <std::floating_point T>
    struct RadiusPoint
    {
        EdgeId edge;
        T at{};
        T radius{};
    };

    /// How the radius varies between two control points.
    enum class RadiusVariation
    {
        Linear, ///< linearly (its slope jumps at the points)
        Smooth  ///< monotone cubic Hermite, C1 (Fritsch and Carlson)
    };

    /// The radius of a fillet along its edges: control points and the variation between them.
    template <std::floating_point T>
    struct RadiusLaw
    {
        std::vector<RadiusPoint<T>> points;
        RadiusVariation variation{RadiusVariation::Smooth};

        /// From r1 at the first vertex of the edge to r2 at its second, linearly.
        static auto linear(EdgeId e, T r1, T r2) -> RadiusLaw
        {
            return {{{e, T(0), r1}, {e, T(1), r2}}, RadiusVariation::Linear};
        }
    };

    namespace radius_detail
    {
        /**
         * The radius as a function of the arc length through the points (s, r), s increasing
         * (`closed`: the law of a loop of length `period`, the points in [s0, s0 + period), its
         * last point followed by the first one a period further). Throws std::invalid_argument
         * for no point, abscissae not increasing, a radius not positive and finite.
         */
        template <std::floating_point T>
        auto radius_function(std::vector<std::pair<T, T>> pts, RadiusVariation variation, bool closed = false, T period = T(0)) -> BSCfunction<T>
        {
            if (pts.empty())
                throw std::invalid_argument("radius law: no point");
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                if (!(pts[i].second > T(0)) || !std::isfinite(pts[i].second) || !std::isfinite(pts[i].first))
                    throw std::invalid_argument("radius law: a radius not positive and finite");
                if (i > 0 && !(pts[i].first > pts[i - 1].first))
                    throw std::invalid_argument("radius law: abscissae not increasing");
            }
            if (closed)
            {
                if (!(period > pts.back().first - pts.front().first))
                    throw std::invalid_argument("radius law: points beyond a period");
                pts.emplace_back(pts.front().first + period, pts.front().second);
            }
            if (pts.size() == 1)
                return BSCfunction<T>{std::vector<T>{pts[0].second, pts[0].second}, std::vector<T>{pts[0].first, pts[0].first, pts[0].first + T(1), pts[0].first + T(1)}, 1};
            const std::size_t n = pts.size() - 1; // spans
            std::vector<T> h(n), d(n);
            for (std::size_t i = 0; i < n; ++i)
            {
                h[i] = pts[i + 1].first - pts[i].first;
                d[i] = (pts[i + 1].second - pts[i].second) / h[i];
            }
            if (variation == RadiusVariation::Linear)
            {
                std::vector<T> poles, knots{pts.front().first};
                for (const auto &p : pts)
                {
                    poles.push_back(p.second);
                    knots.push_back(p.first);
                }
                knots.push_back(pts.back().first);
                return BSCfunction<T>{poles, knots, 1};
            }
            // the slopes: weighted harmonic mean of the slopes of the spans on either side where
            // they agree in sign, 0 where they do not (an extremum) and at the ends of an open law
            auto slope = [&](std::size_t a, std::size_t b) { // spans a (before) and b (after)
                if (!(d[a] * d[b] > T(0)))
                    return T(0);
                return T(3) * (h[a] + h[b]) / ((T(2) * h[b] + h[a]) / d[a] + (h[b] + T(2) * h[a]) / d[b]);
            };
            std::vector<T> m(n + 1, T(0));
            for (std::size_t i = 1; i < n; ++i)
                m[i] = slope(i - 1, i);
            if (closed)
                m[0] = m[n] = slope(n - 1, 0);
            std::vector<T> poles{pts[0].second}, knots(4, pts[0].first);
            for (std::size_t i = 0; i < n; ++i)
            {
                poles.push_back(pts[i].second + m[i] * h[i] / T(3));
                poles.push_back(pts[i + 1].second - m[i + 1] * h[i] / T(3));
                const std::size_t mult = i + 1 == n ? 4 : 2;
                knots.insert(knots.end(), mult, pts[i + 1].first);
            }
            poles.push_back(pts.back().second);
            return BSCfunction<T>{poles, knots, 3};
        }
    }
}
