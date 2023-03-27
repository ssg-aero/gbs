#pragma once
#include <array>
#include <vector>
#include <utility>
#include <limits>
#include <gbs/surfaces>
namespace gbs
{
    template <std::floating_point T>
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
    template <std::floating_point T>
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
    template <std::floating_point T>
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
    template <std::floating_point T>
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
    template <std::floating_point T>
    inline auto tri_area(
        const std::array<T,2> &a,
        const std::array<T,2> &b,
        const std::array<T,2> &c )
        {
            return T(0.5) * std::abs(
                (b[0]-a[0]) * (c[1]-a[1]) - (b[1]-a[1]) * (c[0]-a[0])
            );
        }

    template < std::floating_point T, size_t dim> 
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
/**
 * @brief Calculates the center of a circle that passes through three points.
 * 
 * Given three points a, b, and c in 2D space, this function calculates the center of
 * the circle that passes through all three points. The result is returned as a 2D
 * array of type T with x and y coordinates.
 * 
 * @tparam T The type of the elements in the input and output arrays.
 * @param a The first point.
 * @param b The second point.
 * @param c The third point.
 * @return A 2D array containing the x and y coordinates of the circle's center.
 * If the points are collinear (i.e. no unique solution exists), the result
 * will contain NaN values.
 */
    template <std::floating_point T>
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

/**
 * @brief Calculates the center of a sphere that passes through three points.
 * 
 * Given three points a, b, and c in 3D space, this function calculates the center of
 * the sphere that passes through all three points. The result is returned as a 3D
 * array of type T with x, y, and z coordinates.
 * 
 * @tparam T The type of the elements in the input and output arrays.
 * @param a The first point.
 * @param b The second point.
 * @param c The third point.
 * @return A 3D array containing the x, y, and z coordinates of the sphere's center.
 * If the points are collinear (i.e. no unique solution exists), the result
 * will contain NaN values.
 */
    template <std::floating_point T>
    std::array<T, 3> sphere_center(const std::array<T, 3>& a, const std::array<T, 3>& b, const std::array<T, 3>& c)
    {
        // Calculate some intermediate values
        T A = b[0] - a[0], B = b[1] - a[1], C = b[2] - a[2];
        T D = c[0] - a[0], E = c[1] - a[1], F = c[2] - a[2];
        T G = A * (E * F - B * F) + B * (D * F - C * F) + C * (B * D - E * D);
        T H = (A * A + B * B + C * C) * (E * F - D * F) + (D * D + E * E + F * F) * (B * F - C * E) + (A * A + C * C + B * B) * (D * E - B * C);

        if (std::abs(G) < std::numeric_limits<T>::epsilon()) {
            // Points are collinear, no unique solution exists
            return {{std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()}};
        }

        // Calculate the x, y, and z coordinates of the sphere's center
        T center_x = H / (2 * G);
        T center_y = ((A * A + B * B + C * C) * (D * F - E * F) + (D * D + E * E + F * F) * (C * A - D * B) + (A * A + C * C + B * B) * (E * D - C * B)) / (2 * G);
        T center_z = ((A * A + B * B + C * C) * (B * E - C * D) + (D * D + E * E + F * F) * (A * C - B * B) + (A * A + C * C + B * B) * (B * D - A * E)) / (2 * G);

        // Return the result as a 3D array
        return {{center_x, center_y, center_z}};
    }



    template <std::floating_point T, size_t dim>
    T surface_triangle_deviation(const Surface<T, dim> &srf, const std::array<T, 2> &uv1, const std::array<T, 2> &uv2, const std::array<T, 2> &uv3)
    {
        // Calculate the centroid of the triangle
        auto v1 = srf(uv1);
        auto v2 = srf(uv2);
        auto v3 = srf(uv3);
        auto centroid   = (v1 + v2 + v3) / static_cast<T>(3);
        auto centroid_uv= ( uv1 + uv2 + uv3 ) / static_cast<T>(3);

        // Calculate the normal vectors for the triangle and the surface
        auto tri_normal = normalized(cross(v2 - v1, v3 - v1));
        auto srf_normal = normalized(cross(srf(centroid_uv,1,0), srf(centroid_uv,0,1) ) );

        // Calculate the deviation between the triangle and the surface using the difference in normal vectors
        T dev = norm(tri_normal - srf_normal);

        return dev;
    }

    /**
     * @brief Computes the barycentric coordinates (u, v) for a given point `p` with respect to a triangle defined by vertices `a`, `b`, and `c`.
     * 
     * @tparam T The floating-point type used for coordinates.
     * @tparam dim The dimension of the input points (2 or 3).
     * @tparam TCoord The desired floating-point type for the output barycentric coordinates.
     * @param a The first vertex of the triangle.
     * @param b The second vertex of the triangle.
     * @param c The third vertex of the triangle.
     * @param p The point for which the barycentric coordinates are to be computed.
     * @return std::array<TCoord, 2> The (u, v) barycentric coordinates of the point `p`.
     * @throws std::runtime_error If the triangle vertices are collinear.
     */
    template <std::floating_point T, size_t dim>
    std::array<T, 2> computeUV(const std::array<T, dim> &a,
                                    const std::array<T, dim> &b,
                                    const std::array<T, dim> &c,
                                    const std::array<T, dim> &p)
    {
        // Compute the vector from a to b and a to c
        std::array<T, dim> ab, ac;
        std::transform(a.begin(), a.end(), b.begin(), ab.begin(), std::minus<T>());
        std::transform(a.begin(), a.end(), c.begin(), ac.begin(), std::minus<T>());

        // Compute the vector from a to p
        std::array<T, dim> ap;
        std::transform(a.begin(), a.end(), p.begin(), ap.begin(), std::minus<T>());

        // Compute the dot products
        T ab_ab = std::inner_product(ab.begin(), ab.end(), ab.begin(), T(0));
        T ab_ac = std::inner_product(ab.begin(), ab.end(), ac.begin(), T(0));
        T ac_ac = std::inner_product(ac.begin(), ac.end(), ac.begin(), T(0));
        T ab_ap = std::inner_product(ab.begin(), ab.end(), ap.begin(), T(0));
        T ac_ap = std::inner_product(ac.begin(), ac.end(), ap.begin(), T(0));

        // Compute the determinant
        T det = ab_ab * ac_ac - ab_ac * ab_ac;
        if (std::abs(det) < std::numeric_limits<T>::epsilon())
        {
            throw std::runtime_error("The triangle vertices are collinear.");
        }

        // Compute the barycentric coordinates (u, v)
        T u = (ac_ac * ab_ap - ab_ac * ac_ap) / det;
        T v = (ab_ab * ac_ap - ab_ac * ab_ap) / det;

        return {u, v};
    }

}