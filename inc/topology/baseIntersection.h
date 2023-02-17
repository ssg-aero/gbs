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
}