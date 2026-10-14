#pragma once
/**
 * @file intersection_curve.h
 * @brief Results of the intersection of two surfaces.
 */
#include <gbs/bscurve.h>

#include <memory>
#include <vector>

namespace gbs
{
    /// A curve of the intersection of two surfaces, on [a, b].
    template <std::floating_point T>
    struct IntersectionCurve
    {
        std::shared_ptr<Curve<T, 3>> curve;
        T a{}, b{};
        std::shared_ptr<Curve<T, 2>> pc1, pc2; ///< pcurves on each surface, same parameter as the curve
        T deviation{};                         ///< largest distance between the curve and its pcurves on the surfaces
        bool exact{};                          ///< the three representations are exact
        bool closed{};                         ///< the curve is closed
        bool tangent{};                        ///< the surfaces are tangent along the curve
    };

    /// An isolated contact point of two surfaces.
    template <std::floating_point T>
    struct SurfaceSurfaceHit
    {
        T u1{}, v1{}, u2{}, v2{};
        point<T, 3> p{};
    };

    template <std::floating_point T>
    struct SurfaceSurfaceResult
    {
        std::vector<IntersectionCurve<T>> curves;
        std::vector<SurfaceSurfaceHit<T>> points;
        bool coincident{};    ///< the surfaces coincide on their boxes
        bool singular{};      ///< a curve stopped on a singular point (tangency, branch point), or the surfaces touch along a curve
        bool complete{true};  ///< false when the search was capped
    };

}
