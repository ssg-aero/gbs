#pragma once
#include "baseGeom.h"

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include <gbs/vecop.ixx>
#endif
namespace gbs
{
/**
 * @brief Determines whether a given point D lies inside a triangle defined by points A, B, and C.
 *
 * @tparam T A floating-point type.
 * @param a The coordinates of the triangle's first vertex, A.
 * @param b The coordinates of the triangle's second vertex, B.
 * @param c The coordinates of the triangle's third vertex, C.
 * @param d The coordinates of the point D to be tested.
 * @return `true` if the point D lies inside the triangle, `false` otherwise.
 */
    template<std::floating_point T>
    [[nodiscard]] bool in_triangle(
        const std::array<T, 2>& a,
        const std::array<T, 2>& b,
        const std::array<T, 2>& c,
        const std::array<T, 2>& d)
    {
        // Lambda function to compute the signed area of a triangle defined by three points
        auto sign = [](const auto& p1, const auto& p2, const auto& p3) {
            return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);
        };

        // Calculate the signed areas of the triangles formed by point D and each edge of the input triangle
        auto d1 = sign(d, a, b);
        auto d2 = sign(d, b, c);
        auto d3 = sign(d, c, a);

        // Determine if the point D has a consistent orientation with respect to the input triangle's vertices
        auto has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
        auto has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

        return !(has_neg && has_pos);
    }
/**
 * @brief Computes the barycentric coordinates of a point in a triangle in 2D of non degenerated triangle.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @param a The first vertex of the triangle.
 * @param b The second vertex of the triangle.
 * @param c The third vertex of the triangle.
 * @param d The point for which to compute the barycentric coordinates.
 * @return A pair of barycentric coordinates (u, v) such that d = u*a + v*b + (1-u-v)*c,
 *         or false if the denominator of the calculation is zero or very close to zero.
 */
    template<std::floating_point T>
    [[nodiscard]] auto triangle_barycentric_coordinates(
            const std::array<T, 2>& a,
            const std::array<T, 2>& b,
            const std::array<T, 2>& c,
            const std::array<T, 2>& d)
    {
        // Compute vectors
        std::array<T, 2> v0{c[0] - a[0], c[1] - a[1]};
        std::array<T, 2> v1{b[0] - a[0], b[1] - a[1]};
        std::array<T, 2> v2{d[0] - a[0], d[1] - a[1]};

        // Compute dot products
        T dot00{v0[0] * v0[0] + v0[1] * v0[1]};
        T dot01{v0[0] * v1[0] + v0[1] * v1[1]};
        T dot02{v0[0] * v2[0] + v0[1] * v2[1]};
        T dot11{v1[0] * v1[0] + v1[1] * v1[1]};
        T dot12{v1[0] * v2[0] + v1[1] * v2[1]};

        // Compute the denominator
        T denom{dot00 * dot11 - dot01 * dot01};

        // Compute barycentric coordinates
        T invDenom{ 1.0 / denom };
        T u{(dot11 * dot02 - dot01 * dot12) * invDenom};
        T v{(dot00 * dot12 - dot01 * dot02) * invDenom};
        return std::make_pair(u, v);

    }
/**
 * @brief Determines whether a given point D lies inside a triangle defined by points A, B, and C using barycentric coordinates.
 *
 * @tparam T A floating-point type.
 * @param a The coordinates of the triangle's first vertex, A.
 * @param b The coordinates of the triangle's second vertex, B.
 * @param c The coordinates of the triangle's third vertex, C.
 * @param d The coordinates of the point D to be tested.
 * @return `true` if the point D lies inside the triangle, `false` otherwise.
 */
    template<std::floating_point T>
    [[nodiscard]] bool in_triangle_barycentric(
        const std::array<T, 2>& a,
        const std::array<T, 2>& b,
        const std::array<T, 2>& c,
        const std::array<T, 2>& d)
    {
        // Check if d is equal to one of the triangle's vertices
        if (d == a || d == b || d == c)
        {
            return true;
        }

        auto [u,v] = triangle_barycentric_coordinates(a, b, c, d);

        return (u >= static_cast<T>(0)) && (v >= static_cast<T>(0)) && (u + v <= static_cast<T>(1));
    }
/**
 * @brief Determines the position of the point D (dx, dy) with respect to the circle passing through points A (ax, ay), B (bx, by), and C (cx, cy).
 * 
 * @tparam T A floating-point type.
 * @param ax The x-coordinate of point A.
 * @param ay The y-coordinate of point A.
 * @param bx The x-coordinate of point B.
 * @param by The y-coordinate of point B.
 * @param cx The x-coordinate of point C.
 * @param cy The y-coordinate of point C.
 * @param dx The x-coordinate of point D.
 * @param dy The y-coordinate of point D.
 * @return A positive value if the point D is inside the circle, a negative value if the point D is outside the circle, and zero if the point D is on the circle.
 */
    template <std::floating_point T>
    [[nodiscard]] T in_circle(T ax, T ay, T bx, T by, T cx, T cy, T dx, T dy)
    {
        T A{ax - dx};
        T B{ay - dy};
        T C{A * A + B * B};

        T D{bx - dx};
        T E{by - dy};
        T F{D * D + E * E};

        T G{cx - dx};
        T H{cy - dy};
        T I{G * G + H * H};

        return A * E * I + D * H * C + B * F * G - (G * E * C + D * B * I + A * H * F);
    }

    /**
     * @brief  returns > 0 if d inside circle formed by (a, b, c) return 0. if on
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @param d 
     * @return T 
     */
    template <std::floating_point T>
    [[nodiscard]] T in_circle(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &c,
        const std::array<T, 2> &d)
    {
        auto [ax, ay] = a;
        auto [bx, by] = b;
        auto [cx, cy] = c;
        auto [dx, dy] = d;

        return in_circle( ax, ay, bx, by, cx, cy, dx, dy);
    }

/**
 * @brief Checks if two line segments strictly intersect (excluding endpoints).
 * @tparam T Numeric type (e.g. float, double)
 * @param a First point of the first line segment
 * @param b Second point of the first line segment
 * @param c First point of the second line segment
 * @param d Second point of the second line segment
 * @param tol Tolerance for floating-point comparisons (default: 1e-10)
 * @return true if the line segments strictly intersect, false otherwise
 */
    template <std::floating_point T>
    [[nodiscard]] bool seg_seg_strict_intersection(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &c,
        const std::array<T, 2> &d,
        T tol = T{1e-10}
    )
    {
        T s1_x = b[0] - a[0];
        T s1_y = b[1] - a[1];
        T s2_x = d[0] - c[0];
        T s2_y = d[1] - c[1];

        T delta = -s2_x * s1_y + s1_x * s2_y;

        // Check if the segments are parallel or if they are degenerate (zero length)
        if (std::abs(delta) < tol || (std::abs(s1_x) < tol && std::abs(s1_y) < tol) || (std::abs(s2_x) < tol && std::abs(s2_y) < tol))
        {
            return false;
        }

        T s = (-s1_y * (a[0] - c[0]) + s1_x * (a[1] - c[1])) / delta;
        T t = (s2_x * (a[1] - c[1]) - s2_y * (a[0] - c[0])) / delta;

        // Check if s and t are within the strict intersection range
        const T max = T{1} - tol;
        return (s > tol && s < max && t > tol && t < max);
    }

    template <std::floating_point T>
    bool on_segment(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &p,
        T tol = T(1e-10))
    {

        T cross_product = (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1]);

        // Check if the point is on the line (tolerance for floating-point errors)
        if (std::abs(cross_product) > tol)
        {
            return false;
        }

        // Compute the parameter t to check if the point is inside the segment
        T t;
        if (std::fabs(b[0] - a[0]) >= std::fabs(b[1] - a[1]))
        {
            t = (p[0] - a[0]) / (b[0] - a[0]);
        }
        else
        {
            t = (p[1] - a[1]) / (b[1] - a[1]);
        }

        // Check if t is within the range [0, 1] to ensure the point is inside the segment
        return t >= T(0) && t <= T(1);
    }

    template <std::floating_point T>
    bool seg_D_strict_end_intersection(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c,
        const std::array<T,2> &D,
        T tol = 1e-10
    )
    {
        T s1_x, s1_y;
        s1_x = b[0] - a[0];
        s1_y = b[1] - a[1];
        auto [s2_x, s2_y] = D;

        T s, t, delta{-s2_x * s1_y + s1_x * s2_y};

        const T max = 1-tol;
        if(
            std::fabs(delta)<tol 
            || ( std::fabs(s1_x)<tol && std::fabs(s1_y)<tol))
        {
            return false;
        }
        s = (-s1_y * (a[0] - c[0]) + s1_x * (a[1] - c[1])) / delta;
        t = ( s2_x * (a[1] - c[1]) - s2_y * (a[0] - c[0])) / delta;

        if (
            s > 0. 
                && 
            t >= 0. 
                && 
            t < 1.
            )
        {
            return true;
        }

        return false; 
    }

    template <std::floating_point T>
    bool seg_H_strict_end_intersection(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &c,
        T tol = 1e-10)
    {
        return seg_D_strict_end_intersection(a, b, c, {1, 0}, tol);
    }
}