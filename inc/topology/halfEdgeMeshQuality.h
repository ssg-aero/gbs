#pragma once

#include "halfEdgeMeshData.h"
#include <execution>
#include <algorithm>
namespace gbs
{
/**
 * @brief Utility function to sum the area of a given triangle face with a given area value.
 * 
 * @tparam T Floating point type.
 * @param Area Existing area value.
 * @param h_face Shared pointer to a half-edge face with 3 vertices.
 * @return T The sum of the existing area value and the area of the given triangle face.
 */
    template <typename T>
    inline T sum_area(T Area, const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face)
    {
        auto coords = getFaceCoords(h_face);
        assert(coords.size()==3);
        return Area + tri_area(
                    *std::next(coords.begin(), 0),
                    *std::next(coords.begin(), 1),
                    *std::next(coords.begin(), 2)
                );
    }

    template <typename T>
    inline T sum_area(const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face, T Area)
    {
        return sum_area(Area, h_face);
    }

    template <typename T>
    inline T sum_area( T Area1, T Area2)
    {
        return Area1+Area2;
    }

    template <typename T>
    inline T sum_area(const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face1, const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face2)
    {
        return sum_area(T{},h_face1)+sum_area(T{},h_face2);
    }

/**
 * @brief Calculates the total area of a triangle mesh.
 * 
 * @tparam T Floating point type.
 * @param faces_lst List of shared pointers to half-edge faces.
 * @return T The total area of the triangle mesh.
 */
    auto getTriangle2dMeshArea(const auto &faces_lst)
    {
        return std::reduce(
            faces_lst.begin(),faces_lst.end(),
            0.,
            [](auto a1, auto a2)
            {
                return sum_area(a1, a2);
            }
        );
    }

/**
 * @brief Calculates the total area of a triangle mesh using parallel execution.
 * 
 * @tparam T Floating point type.
 * @param faces_lst List of shared pointers to half-edge faces.
 * @return T The total area of the triangle mesh.
 */
    auto getTriangle2dMeshAreaPar(const auto &faces_lst)
    {
        return std::reduce(
            std::execution::par,
            faces_lst.begin(),faces_lst.end(),
            0.,
            [](auto a1, auto a2)
            {
                return sum_area(a1, a2);
            }
        );
    }
    template <typename T, size_t dim>
    struct DistanceMeshSurface
    {
        // Reference to the input surface
        const Surface<T, dim> &srf;

        // Constructor taking the surface as input
        DistanceMeshSurface(const auto &srf_) : srf{srf_} {}

        // Function call operator overload to process a HalfEdgeFace
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            // Get the face coordinates in UV space
            auto coordsUV = getFaceCoords(hf);

            // Get the bounds of the surface in UV space
            auto [u1, u2, v1, v2] = srf.bounds();

            // Check if the face coordinates are within the surface bounds
            for (auto [u, v] : coordsUV)
            {
                if (u1 > u || u2 < u || v1 > v || v2 < v)
                {
                    // Return a zero distance and a dummy UV point for out-of-bounds coordinates
                    return std::make_pair(0., std::array<T, 2>{0., 0.});
                }
            }

            // Transform the face coordinates from UV space to 3D space
            std::vector<std::array<T, 3>> coords(coordsUV.size());
            std::transform(
                coordsUV.begin(), coordsUV.end(),
                coords.begin(),
                [&srf = this->srf](const auto &UV)
                {
                    return srf(UV[0], UV[1]);
                });

            // Compute the centroid (G) and midpoints (A, B, C) in both UV and 3D space
            auto G_UV = average(coordsUV.begin(), coordsUV.end());
            auto G = average(coords.begin(), coords.end());
            auto A_UV = midpoint(coordsUV[0], coordsUV[1]);
            auto A = midpoint(coords[0], coords[1]);
            auto B_UV = midpoint(coordsUV[1], coordsUV[2]);
            auto B = midpoint(coords[1], coords[2]);
            auto C_UV = midpoint(coordsUV[2], coordsUV[0]);
            auto C = midpoint(coords[2], coords[0]);

            // Compute the distances between the surface and the computed points
            auto test_lst = {
                std::make_pair(distance(G, srf(G_UV[0], G_UV[1])), G_UV),
                std::make_pair(distance(A, srf(A_UV[0], A_UV[1])), A_UV),
                std::make_pair(distance(B, srf(B_UV[0], B_UV[1])), B_UV),
                std::make_pair(distance(C, srf(C_UV[0], C_UV[1])), C_UV)};

            // Return the pair with the largest distance and its corresponding UV point
            return *std::max_element(
                test_lst.begin(), test_lst.end(),
                [](const auto &p1, const auto &p2)
                { return p1.first < p2.first; });
        }
    };
    // template <typename T, size_t dim>
    // struct DistanceMeshSurface
    // {
    //     const Surface<T, dim> &srf;
    //     DistanceMeshSurface(const auto &srf_) : srf{srf_} {}
    //     auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    //     {
    //         auto coordsUV = getFaceCoords(hf);
    //         auto [u1, u2, v1, v2] = srf.bounds();
    //         for (auto [u, v] : coordsUV)
    //         {
    //             if (u1 > u || u2 < u || v1 > v || v2 < v)
    //             // if (u1 >= u || u2 <= u || v1 >= v || v2 <= v)
    //             {
    //                 return std::make_pair(0., std::array<T, 2>{0., 0.}); // ext dummy are perfect
    //             }
    //         }
    //         std::vector<std::array<T, 3>> coords(coordsUV.size());
    //         std::transform(
    //             coordsUV.begin(), coordsUV.end(),
    //             coords.begin(),
    //             [&srf = this->srf](const auto &UV)
    //             {
    //                 return srf(UV[0], UV[1]);
    //             });
    //         auto G_UV = std::reduce(
    //                         coordsUV.begin(), coordsUV.end(),
    //                         std::array<T, 2>{0., 0.},
    //                         gbs::operator+<T, 2>) /
    //                     3.;
    //         auto G = std::reduce(
    //                     coords.begin(), coords.end(),
    //                     std::array<T, 3>{0., 0., 0.},
    //                     gbs::operator+<T, 3>) /
    //                 3.;

    //         auto A_UV = 0.5 * (*std::next(coordsUV.begin(), 0) + *std::next(coordsUV.begin(), 1));
    //         auto A = 0.5 * (*std::next(coords.begin(), 0) + *std::next(coords.begin(), 1));
    //         auto B_UV = 0.5 * (*std::next(coordsUV.begin(), 1) + *std::next(coordsUV.begin(), 2));
    //         auto B = 0.5 * (*std::next(coords.begin(), 1) + *std::next(coords.begin(), 2));
    //         auto C_UV = 0.5 * (*std::next(coordsUV.begin(), 2) + *std::next(coordsUV.begin(), 0));
    //         auto C = 0.5 * (*std::next(coords.begin(), 2) + *std::next(coords.begin(), 0));

    //         auto test_lst = {
    //             std::make_pair(distance(G, srf(G_UV[0], G_UV[1])), G_UV),
    //             std::make_pair(distance(A, srf(A_UV[0], A_UV[1])), A_UV),
    //             std::make_pair(distance(B, srf(B_UV[0], B_UV[1])), B_UV),
    //             std::make_pair(distance(C, srf(C_UV[0], C_UV[1])), C_UV)};

    //         return *std::max_element(
    //             test_lst.begin(), test_lst.end(),
    //             [](const auto &p1, const auto &p2)
    //             { return p1.first < p2.first; });
    //     }
    // };

    /**
     * @brief Struct to calculate the distance between a mesh surface and a parametric surface.
     * 
     * @tparam T Floating point type.
     * @tparam dim Dimension of the space.
     */
    template <typename T, size_t dim>
    struct DistanceMeshSurface2
    {
        const Surface<T, dim> &srf; // Reference to the parametric surface.

        /**
         * @brief Constructor for the DistanceMeshSurface2 struct.
         * 
         * @param srf_ Reference to the parametric surface.
         */
        DistanceMeshSurface2(const Surface<T, dim> &srf_) : srf(srf_) {}

        /**
         * @brief Calculate the distance between the parametric surface and a given half-edge face.
         * 
         * @param hf Shared pointer to a half-edge face.
         * @return std::pair<T, std::array<T, 2>> Pair containing the maximum distance and the corresponding UV coordinates.
         */
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};
            std::vector< std::pair<T, std::array<T,2> > > test_lst(std::distance(begin_, end_));

            std::transform(
                begin_, end_,
                std::back_inserter(test_lst),
                [&srf = this->srf](const auto &he){
                    assert(he && he->vertex && he->previous && he->previous->vertex);
                    if (!he->opposite)
                    {
                        return std::make_pair(static_cast<T>(0), std::array<T, 2>{0, 0}); // boundary edges are supposed to respect criteria
                    }
                    const auto &p1_uv = he->previous->vertex->coords;
                    const auto &p2_uv = he->vertex->coords;
                    auto mid_uv = p1_uv + static_cast<T>(0.5) * (p2_uv - p1_uv);
                    auto mid = srf(mid_uv[0], mid_uv[1]);
                    auto p1  = srf(p1_uv[0], p1_uv[1]);
                    auto p2  = srf(p2_uv[0], p2_uv[1]);
                    auto l1m = distance( p1, mid );
                    auto lm2 = distance( p2, mid );
                    auto l12 = distance( p2, p1 );
                    // auto dev = distance(mid, static_cast<T>(0.5) * (srf(p1_uv[0], p1_uv[1]) + srf(p2_uv[0], p2_uv[1])));
                    auto dev = ( l1m + lm2 - l12) / l12 ;
                    return std::make_pair( dev, mid_uv );
                }
            );

            const auto max_it = std::max_element(
                test_lst.cbegin(), test_lst.cend(),
                [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

            return *max_it;
        }
    };

    /**
     * @brief Struct to calculate the normal deviation between a mesh surface and a parametric surface.
     * 
     * @tparam T Floating point type.
     * @tparam dim Dimension of the space.
     */
    template <typename T, size_t dim>
    struct DeviationMeshSurface
    {
        const Surface<T, dim> &srf; // Reference to the parametric surface.

        /**
         * @brief Constructor for the DistanceMeshSurface2 struct.
         * 
         * @param srf_ Reference to the parametric surface.
         */
        DeviationMeshSurface(const Surface<T, dim> &srf_) : srf(srf_) {}

        /**
         * @brief Calculate the deviation between the parametric surface and a given half-edge triangular face.
         * 
         * @param hf Shared pointer to a half-edge face.
         * @return std::pair<T, std::array<T, 2>> Pair containing the maximum distance and the corresponding UV coordinates.
         */
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            auto [uv1, uv2, uv3] = getTriangleCoords(hf);

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

            return std::make_pair(dev, centroid_uv);
            // return std::make_pair(dev, circle_center(uv1, uv2, uv3));
        }
    };

    template <typename T, size_t dim>
    struct DeviationMeshSurface2
    {
        const Surface<T, dim> &srf; // Reference to the parametric surface.


        DeviationMeshSurface2(const Surface<T, dim> &srf_) : srf(srf_) {}

        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};
            std::vector< std::pair<T, std::array<T,2> > > test_lst();

            std::transform(
                begin_, end_,
                std::back_inserter(test_lst),
                [&srf = this->srf](const auto &he){
                    assert(he && he->vertex && he->previous && he->previous->vertex);
                    if (!he->opposite)
                    {
                        return std::make_pair(static_cast<T>(0), std::array<T, 2>{0, 0}); // boundary edges are supposed to respect criteria
                    }
                    const auto &p1_uv = he->previous->vertex->coords;
                    const auto &p2_uv = he->vertex->coords;
                    auto [delta_u, delta_v] = (p2_uv - p1_uv);
                    auto mid_uv = p1_uv + static_cast<T>(0.5) * (p2_uv - p1_uv);
                    auto v1 = srf.tangent(p1_uv[0], p1_uv[1], delta_u, delta_v);
                    auto v2 = srf.tangent(p2_uv[0], p2_uv[1], delta_u, delta_v);
                    auto dev_ = norm(cross(v1, v2)) / (norm(v1) * norm(v2));
                    return std::make_pair(dev_, mid_uv);
                }
            );

            const auto max_it = std::max_element(
                test_lst.cbegin(), test_lst.cend(),
                [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

            return *max_it;
        }
    };

    template <typename T, size_t dim>
    struct MaxEdgeSizeSurface
    {
        const Surface<T, dim> &srf; // Reference to the parametric surface.


        MaxEdgeSizeSurface(const Surface<T, dim> &srf_) : srf(srf_) {}

        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};

            // compute edges lenght
            std::vector<std::pair<T, std::array<T, 2>>> l_lst;
            std::transform(
                begin_, end_,
                std::back_inserter(l_lst),
                [&srf = this->srf](const auto &he)
                {
                    return he->opposite ? std::make_pair(sq_distance(srf(he->vertex->coords), srf(he->previous->vertex->coords)), edge_midpoint(he)) : std::make_pair(static_cast<T>(0.), std::array<T, 2>{});
                });
            auto max = std::max_element(l_lst.cbegin(), l_lst.cend(), [](const auto &p1, const auto &p2)
                                        { return p1.first < p2.first; });
            return *max;
        }
    };

    template <typename T, size_t dim>
    struct MaxEdgeSize
    {
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, dim>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};

            // compute edges length
            std::vector<std::pair<T, std::array<T, dim>>> l_lst(std::distance(begin_, end_));
            std::transform(
                begin_, end_,
                std::back_inserter(l_lst),
                [](const auto &he)
                {
                    return he->opposite ? std::make_pair(edge_sq_length(he), edge_midpoint(he)) : std::make_pair(static_cast<T>(0.), std::array<T, dim>{});
                });
            auto max = std::max_element(l_lst.cbegin(), l_lst.cend(), [](const auto &p1, const auto &p2)
                                        { return p1.first < p2.first; });
            return *max;
        }
    };

    template <typename T, size_t dim>
    struct TriangleShapeEdgeSplit
    {
        T min_edge_length{0.02};
        TriangleShapeEdgeSplit(T min_ed) : min_edge_length{min_ed} {}
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, dim>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};

            // if( std::find_if(begin_, end_,[](const auto &he){return he->opposite==nullptr;}) != end_){
            //     return std::make_pair(std::numeric_limits<T>::lowest(), std::array<T,dim>{});
            // }
            
            // compute edges length
            std::vector < T > l_lst;
            std::transform(
                begin_, end_,
                std::back_inserter(l_lst),
                [](const auto &he)
                {
                    return edge_sq_length(he);
                });

            const auto [min, max] = std::minmax_element(l_lst.cbegin(), l_lst.cend(), [](const auto &l1, const auto &l2)
                                        { return l1 < l2; });
            // find largest edge outside boundary
            auto mid_max = edge_midpoint( *(std::max_element(
                begin_, end_,
                [](const auto &he1, const auto &he2){
                    auto l1 = he1->opposite ? edge_sq_length(he1) : std::numeric_limits<T>::lowest();
                    auto l2 = he2->opposite ? edge_sq_length(he2) : std::numeric_limits<T>::lowest();
                    // auto l1 = edge_sq_length(he1) ;
                    // auto l2 = edge_sq_length(he2) ;
                    return  l1 < l2;}
            )));

        // const auto [a, b, c] = getTriangleCoords(hf);
        // auto pt = 0.5 * ( (( a+ b +c) / static_cast<T>(3)) + mid_max );

        // Check if the minimum edge length is below the given threshold
        if (std::sqrt(*max) < min_edge_length)
        {
            // Set the shape quality to a small value, so it's less likely to be picked
            // return std::make_pair(std::numeric_limits<T>::min(), pt);
            return std::make_pair(std::numeric_limits<T>::lowest(), mid_max);
        }

            return std::make_pair(
                std::sqrt((*max)/(*min)), 
                // pt
                mid_max
            );
        }
    };

    template <typename T, size_t dim>
    struct TriangleShape
    {
        T min_edge_length{0.02};
        TriangleShape(T min_ed) : min_edge_length{min_ed} {}
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, dim>> &hf) const
        {
            auto begin_{begin(*hf)};
            auto end_{end(*hf)};

            // compute edges length
            std::vector<T> l_lst;
            std::transform(
                begin_, end_,
                std::back_inserter(l_lst),
                [](const auto &he)
                {
                    return edge_sq_length(he);
                });

            const auto [min, max] = std::minmax_element(l_lst.cbegin(), l_lst.cend());

            const auto [a, b, c] = getTriangleCoords(hf);


        // Check if the minimum edge length is below the given threshold
        if (std::sqrt(*max) < min_edge_length)
        {
            // Set the shape quality to a high value, so it's less likely to be picked
            return std::make_pair(std::numeric_limits<T>::lowest(), ((a + b + c) / static_cast<T>(3)));
        }

            return std::make_pair(
                std::sqrt((*max)/(*min)), 
                // 0.3 * (
                // circle_center(a, b, c)
                // ) + 0.7 * (
                (( a+ b +c) / static_cast<T>(3))
                
                // )
                // incenter(a, b, c)
            );
        }
    };

    template <typename T, size_t dim>
    struct TriangleShape2
    {
        T min_edge_length{0.02};
        // TriangleShape(T min_ed) : min_edge_length{min_ed}
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, dim>> &hf) const
        {
            const auto [a, b, c] = getTriangleCoords(hf);
            T sin60{std::sqrt(T(3))*T(0.5)}; 
            auto a_ = std::abs(sin60-std::abs((b-a) * (c-a)));
            auto b_ = std::abs(sin60-std::abs((a-b) * (c-b)));
            auto c_ = std::abs(sin60-std::abs((b-c) * (a-c)));
            auto max = std::max(a_,std::max(b_, c_));

        // Check if the minimum edge length is below the given threshold
        // if (std::sqrt(*max) < min_edge_length)
        // {
        //     // Set the shape quality to a high value, so it's less likely to be picked
        //     return std::make_pair(std::numeric_limits<T>::min(), ((a + b + c) / static_cast<T>(3)));
        // }

            return std::make_pair(
                // (*max)/(*min), 
                max,
                // 0.3 * (
                // circle_center(a, b, c)
                // ) + 0.7 * (
                (( a+ b +c) / static_cast<T>(3))
                // )
                // incenter(a, b, c)
            );
        }
    };

    // template <typename T, size_t dim>
    // struct TriangleShape3
    // {
    //     T min_edge_length{0.02};
    //     // TriangleShape(T min_ed) : min_edge_length{min_ed}
    //     auto operator()(const std::shared_ptr<HalfEdgeFace<T, dim>> &hf)
    //     {
    //         auto begin_{begin(*hf)};
    //         auto end_{end(*hf)};

    //         // compute edges length
    //         std::vector<T> l_lst;
    //         std::transform(
    //             begin_, end_,
    //             std::back_inserter(l_lst),
    //             [](const auto &he)
    //             {
    //                 return edge_length(he);
    //             });
    //         const auto [a, b, c] = getTriangleCoords(hf);
    //         auto area = tri_area(a,b,c);
    //         auto s = std::reduce(l_lst.begin(), l_lst.end()) * static_cast<T>(0.5); // half perimeter

    //         return std::make_pair(
    //             // (*max)/(*min), 
    //             max,
    //             // 0.3 * (
    //             // circle_center(a, b, c)
    //             // ) + 0.7 * (
    //             (( a+ b +c) / static_cast<T>(3))
    //             // )
    //             // incenter(a, b, c)
    //         );
    //     }
    // };

    template <typename T, size_t dim>
    struct DistanceMeshSurface3
    {
        const Surface<T, dim> &srf;
        DistanceMeshSurface3(const auto &srf_) : srf{srf_} {}
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf) const
        {
            auto coords_lst = getFaceCoords(hf);
            auto ed_lst = getFaceEdges(hf);
            auto n = coords_lst.size();

            auto mid_uv = static_cast<T>(1./n) * std::reduce(
                coords_lst.begin(), coords_lst.end(),
                std::array<T,2>{},
                [](const auto &X1, const auto &X2){return X1+X2;}
            );

            // auto mid = static_cast<T>(1./n) * std::reduce(
            //     coords_lst.begin(), coords_lst.end(),
            //     std::array<T,dim>{},
            //     [&srf = this->srf](const auto &X1, const auto &X2)
            //     {
            //         return srf(X1[0],X1[1]) + srf(X2[0],X2[1]) );
            //     }
            // );
            std::array<T,dim> mid{};
            for(const auto &X : coords_lst)
            {
                mid = mid + srf(X[0],X[1]);
            }
            mid = static_cast<T>(1./n) * mid;

            std::vector< std::pair<T, std::array<T,2> > > test_lst(n+1);

            std::transform(
                ed_lst.begin(), ed_lst.end(),
                test_lst.begin(),
                [&srf = this->srf](const auto &he){
                    assert(he && he->vertex && he->previous && he->previous->vertex);
                    if(!he->opposite)
                    {
                        return std::make_pair(0., std::array<T, 2>{0., 0.}); // boundary edges are supposed to respect criteria
                    }
                    const auto &p1_uv = he->previous->vertex->coords;
                    const auto &p2_uv = he->vertex->coords;
                    auto mid_uv = p1_uv + static_cast<T>(0.5) *(p2_uv-p1_uv);
                    auto p1 = srf(p1_uv[0], p1_uv[1]);
                    auto p2 = srf(p2_uv[0], p2_uv[1]);
                    auto mid = p1 + static_cast<T>(0.5) *(p2-p1);
                    return std::make_pair(distance(mid, srf(mid_uv[0], mid_uv[1])), mid_uv);
                }
            );

            test_lst.back() = std::make_pair(distance(mid, srf(mid_uv[0], mid_uv[1])), mid_uv);
            
            return *std::max_element(
                test_lst.begin(), test_lst.end(),
                [](const auto &p1, const auto &p2)
                { return p1.first < p2.first; });

        }
    };
} // namespace gbs
