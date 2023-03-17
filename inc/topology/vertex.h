#pragma once

#include <array>

#include "basetopo.h"
#include "halfEdgeMeshData.h"

namespace gbs
{
    /**
     * @brief Vertex class representing a vertex in a mesh.
     *
     * @tparam T Floating-point type used for vertex coordinates.
     * @tparam dim Dimension of the vertex.
     */
    template <std::floating_point T, size_t dim>
    class Vertex : public BaseTopo<T, dim>
    {
        std::shared_ptr<HalfEdgeVertex<T, dim>> he_vertex{};

    public:
        Vertex() = delete;

        /// @brief Constructor with point coordinates.
        Vertex(const std::array<T, dim> &pnt) : he_vertex{std::make_shared<HalfEdgeVertex<T, dim>>(pnt)} {}

        /// @brief Constructor with point coordinates and tolerance reference.
        Vertex(const std::array<T, dim> &pnt, const BaseTopo<T, dim> &tolRef) : he_vertex{std::make_shared<HalfEdgeVertex<T, dim>>(pnt)}, BaseTopo<T, dim>{tolRef} {}

        /**
         * @brief Get the point coordinates.
         *
         * @return const auto& A reference to the point coordinates.
         */
        const auto &point() const noexcept
        {
            return he_vertex->coords;
        }

        /**
         * @brief Set the point coordinates.
         *
         * @param pnt The new point coordinates.
         */
        void setPoint(const std::array<T, dim> &pnt)
        {
            he_vertex->coords = pnt;
        }

        /// @brief Tessellate the vertex. Currently empty.
        void tessellate() override
        {
        }
    };
}