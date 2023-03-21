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
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
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
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
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
                    return std::make_pair(distance(mid, static_cast<T>(0.5) * (srf(p1_uv[0], p1_uv[1]) + srf(p2_uv[0], p2_uv[1]))), mid_uv);
                }
            );

            const auto max_it = std::max_element(
                test_lst.cbegin(), test_lst.cend(),
                [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

            return *max_it;
        }
    };

    template <typename T, size_t dim>
    struct DistanceMeshSurface3
    {
        const Surface<T, dim> &srf;
        DistanceMeshSurface3(const auto &srf_) : srf{srf_} {}
        auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
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
