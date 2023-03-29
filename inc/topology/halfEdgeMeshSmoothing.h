#pragma once
#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
namespace gbs
{
    template <std::floating_point T, size_t dim>
    void laplacian_smoothing(auto &faces_lst, HalfEdgeVertex<T, dim> & h_v)
    {
        auto neighbors = getNeighboringVertices(h_v);

        std::array<T,dim> centroid{};
        for( const auto &vtx : neighbors)
        {
            centroid = centroid + vtx->coords;
        }
        centroid = centroid / static_cast<T>(neighbors.size());
        h_v.coords = centroid;
    }
}