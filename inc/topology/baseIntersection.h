#pragma once

namespace gbs
{
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

    template <typename T>
    T orient_2d( T ax, T ay, T bx, T by, T cx, T cy)
    {
        return (ax-cx)*(by-cy) - (ay-cy)*(bx-cx);
    }
    /**
     * @brief Returns > 0 if the triangle (a,b,c) is counter clockwise 0. if degenerated
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @return T 
     */
    template <typename T>
    T orient_2d(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c)
    {
        auto [ax,ay] = a;
        auto [bx,by] = b;
        auto [cx,cy] = c;
        return (ax-cx)*(by-cy) - (ay-cy)*(bx-cx);
    }
    /**
     * @brief Check if the triangle (a,b,c) is counter clockwise
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @return T 
     */
    template <typename T>
    T are_ccw(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c)
    {
        return orient_2d(a,b,c) >= 0;
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

}