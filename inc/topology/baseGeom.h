#pragma once
#include <array>
#include <vector>
#include <utility>
#include <limits>

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

    template < typename T, size_t dim> 
    auto getCoordsMinMax(const std::vector< std::array<T, dim> > &X_lst)
    {
        auto Xmin{X_lst.front()};
        auto Xmax{X_lst.front()};

        for(const auto &X : X_lst)
        {
            for(size_t i{}; i < dim; i++)
            {
                Xmin[i] = std::min(Xmin[i], X[i]);
                Xmax[i] = std::max(Xmax[i], X[i]);
            }
        }

        return std::make_pair(Xmin,Xmax);
    }

    template <typename T>
    std::array<T, 2> circle_center(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c)
    {
        
        T A = b[0] - a[0], B = b[1] - a[1];
        T C = c[0] - a[0], D = c[1] - a[1];
        T E = A * (a[0] + b[0]) + B * (a[1] + b[1]);
        T F = C * (a[0] + c[0]) + D * (a[1] + c[1]);
        T G = static_cast<T>(2) * (A * (c[1] - b[1]) - B * (c[0] - b[0]));

        if (std::abs(G) < std::numeric_limits<T>::epsilon()) {
            // Points are colinear, no unique solution exists
            return {{std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()}};
        }

        T center_x = (D * E - B * F) / G;
        T center_y = (A * F - C * E) / G;

        return {{center_x, center_y}};
    }

}