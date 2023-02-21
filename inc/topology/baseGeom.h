#pragma once
#include <array>

namespace gbs
{
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
     * @brief Compute (a, b, c) triangle area
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @return auto 
     */
    template <typename T>
    inline auto tri_area(
        const std::array<T,3> &a,
        const std::array<T,3> &b,
        const std::array<T,3> &c )
        {
            return T(0.5)*norm((b-a)^(c-a));
        }
    /**
     * @brief Compute (a, b, c) triangle area
     * 
     * @tparam T 
     * @param a 
     * @param b 
     * @param c 
     * @return auto 
     */
    template <typename T>
    inline auto tri_area(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c )
        {
            return T(0.5) * std::abs(
                (b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0])
            );
        }
}