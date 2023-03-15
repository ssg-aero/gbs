#pragma once

#include "halfEdgeMeshData.h"

namespace gbs
{
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
} // namespace gbs
