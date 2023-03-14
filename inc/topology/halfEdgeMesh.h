#pragma once

namespace gbs
{
    template <typename T, size_t dim>
    struct HalfEdgeMesh
    {
        std::vector<std::shared_ptr<HalfEdge<T, dim>>> edges;
        std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> vertices;
        std::vector<std::shared_ptr<HalfEdgeFace<T, dim>>> faces;
        // public:
        void addFreeVertices(const auto &free_vertices)
        {
            vertices.insert(
                vertices.end(),
                free_vertices.begin(), free_vertices.end()
            );
        }
        
        /**
         * @brief Must be cw
         * 
         * @param i 
         * @param j 
         * @param k 
         * @return true 
         * @return false 
         */
        bool makeFace( size_t i, size_t j, size_t k)
        {
            // auto vtx1 = vertices[i];
            // auto vtx2 = vertices[j];
            // auto vtx3 = vertices[k];

            // auto p_face = std::make_shared<HalfEdgeFace<T, dim>>()

            // HalfEdge<T, dim> ed1{
            //     .vertex = vtx1,
            //     .face   = p_face,
            // };
            // HalfEdge<T, dim> ed2{
            //     .vertex = vtx2,
            //     .face   = p_face,
            // };
            // HalfEdge<T, dim> ed3{
            //     .vertex = vtx3,
            //     .face   = p_face,
            // };

            // // auto p_ed1 = std::make_shared< HalfEdge<T, dim> > (  HalfEdge<T, dim>{, .} 


            return true;
        }
    };
}