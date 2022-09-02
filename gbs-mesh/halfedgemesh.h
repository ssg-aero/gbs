#pragma once
#include <gbs-mesh/types.h>

namespace gbs
{
    struct HalfEdge;
    struct HalfEdgeVertex;
    struct HalfEdgeFace;
    using ItHalfEdge       = vector<HalfEdge>::iterator;
    using ItHalfEdgeFace   = vector<HalfEdgeFace>::iterator;
    using ItHalfEdgeVertex = vector<HalfEdgeVertex>::iterator;

    struct HalfEdgeVertex
    {
        Real3 coords;
        HalfEdge &refHalfEdge;
    };

    struct HalfEdgeFace
    {
        HalfEdge &refHalfEdge;
    };

    struct HalfEdge
    {
        HalfEdgeVertex &vertex;
        HalfEdgeFace &face;
        HalfEdge      &next;
        HalfEdge      &previous;
        HalfEdge      &opposite;
    };

    struct HalfEdgeMesh
    {
        vector<HalfEdge> edges_;
        vector<HalfEdgeVertex> vertices_;
        vector<HalfEdgeFace> faces_;
    };

    // _CPU_GPU_
    // auto addEdge(const Real3 &p1)
    // {
    //     // auto n_vtx = vertices_.size();
    //     HalfEdge ed;
    //     HalfEdgeVertex vtx1{p1}; 
    // }
    // void make_half_edge()
    // _CPU_GPU_
    // void add_face(const Real3 &p1, const Real3 &p2, const Real3 &p3, HalfEdgeMesh &msh)
    // {
    //     // use c+20 asap
    //     // HalfEdge e1 {
    //     //     .vertex = HalfEdgeVertex{
    //     //         .coord{p1}
    //     //     }
    //     // };
    //     HalfEdge e1;
    //     HalfEdge e2;
    //     HalfEdge e3;
    //     HalfEdgeFace f;
    //     HalfEdgeVertex vtx1{p1,} 
    // }
    
}