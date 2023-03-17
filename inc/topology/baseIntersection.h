#pragma once
#include "baseGeom.h"
#include <gbs/vecop.h>
namespace gbs
{

    template<typename T>
    bool in_triangle(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &c,
        const std::array<T, 2> &d)
    {
        auto sign = [](const auto &p1, const auto &p2, const auto &p3){
            return (p1[0]-p3[0])*(p2[1]-p3[1]) - (p2[0]-p3[0])*(p1[1]-p3[1]);
        };

        auto d1 = sign(d, a, b);
        auto d2 = sign(d, b, c);
        auto d3 = sign(d, c, a);

        auto has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
        auto has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

        return !(has_neg && has_pos);
    }
    /**
     * @brief returns > 0 if d inside circle formed by (a, b, c) return 0. if on
     * 
     * @tparam T 
     * @param ax 
     * @param ay 
     * @param bx 
     * @param by 
     * @param cx 
     * @param cy 
     * @param dx 
     * @param dy 
     * @return T 
     */
    template <typename T>
    T in_circle( T ax, T ay, T bx, T by, T cx, T cy, T dx, T dy)
    {
        auto A = ax - dx;
        auto B = ay - dy;
        auto C = A * A + B * B;

        auto D = bx - dx;
        auto E = by - dy;
        auto F = D * D + E * E;

        auto G = cx - dx;
        auto H = cy - dy;
        auto I = G * G + H * H;

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
    template <typename T>
    T in_circle(
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
     * @brief Intersection of segments ]a,b[ and ]c,d[
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @param d 
     * @param tol 
     * @return true 
     * @return false 
     */
    template <typename T>
    bool seg_seg_strict_intersection(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c,
        const std::array<T,2> &d,
        T tol = 1e-10
    )
    {
        T s1_x, s1_y, s2_x, s2_y;
        s1_x = b[0] - a[0];
        s1_y = b[1] - a[1];
        s2_x = d[0] - c[0];
        s2_y = d[1] - c[1];

        T s, t, delta{-s2_x * s1_y + s1_x * s2_y};

        const T max = 1-tol;
        if(
            std::fabs(delta)<tol 
            || ( std::fabs(s1_x)<tol && std::fabs(s1_y)<tol) 
            || ( std::fabs(s2_x)<tol && std::fabs(s2_y)<tol))
        {
            return false;
        }
        s = (-s1_y * (a[0] - c[0]) + s1_x * (a[1] - c[1])) / delta;
        t = ( s2_x * (a[1] - c[1]) - s2_y * (a[0] - c[0])) / delta;

        if (s > tol && s < max && t > tol && t < max)
        {
            return true;
        }

        return false; 
    }

    template <typename T>
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

    template <typename T>
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

    template <typename T>
    bool seg_H_strict_end_intersection(
        const std::array<T, 2> &a,
        const std::array<T, 2> &b,
        const std::array<T, 2> &c,
        T tol = 1e-10)
    {
        return seg_D_strict_end_intersection(a, b, c, {1, 0}, tol);
    }
}