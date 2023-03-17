#pragma once

#include <array>
#include <memory>
#include <vector>
#include <algorithm>
#include <execution>

#include <gbs/curves>
#include <gbs/bscbuild.h>

#include "vertex.h"
#include "halfEdgeMeshData.h"

namespace gbs
{
    /**
     * @brief Edge class representing an edge in a mesh.
     *
     * @tparam T Floating-point type used for edge coordinates.
     * @tparam dim Dimension of the edge.
     */
    template <typename T, size_t dim>
    class Edge : public BaseTopo<T, dim>
    {
        std::shared_ptr<Curve<T, dim>> m_curve;
        std::shared_ptr<Vertex<T, dim>> m_vtx1;
        std::shared_ptr<Vertex<T, dim>> m_vtx2;

        std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> he_vertices{};

        T m_u1;
        T m_u2;

    public:
        /// @brief Deleted default constructor.
        Edge() = delete;

        /// @brief Constructor with two point coordinates.
        Edge(const std::array<T, dim>& pt1, const std::array<T, dim>& pt2);

        /// @brief Constructor with two vertices.
        Edge(const std::shared_ptr<Vertex<T, dim>>& vtx1, const std::shared_ptr<Vertex<T, dim>>& vtx2);

        /// @brief Constructor with a curve.
        Edge(const std::shared_ptr<Curve<T, dim>>& crv);

        /// @brief Constructor with a curve and tolerance reference.
        Edge(const std::shared_ptr<Curve<T, dim>>& crv, const BaseTopo<T, dim>& tolRef);

        /// @brief Get the first vertex.
        [[nodiscard]] const auto& vertex1() const;

        /// @brief Get the second vertex.
        [[nodiscard]] const auto& vertex2() const;

        /// @brief Get the curve.
        [[nodiscard]] const auto& curve() const;

        /// @brief Get the bounds.
        [[nodiscard]] auto bounds() const;

        /// @brief Set the first vertex.
        bool setVertex1(const std::shared_ptr<Vertex<T, dim>>& vtx);

        /// @brief Set the second vertex.
        bool setVertex2(const std::shared_ptr<Vertex<T, dim>>& vtx);

        /// @brief Tessellate the edge.
        void tessellate() override;
    };

    template <typename T, size_t dim>
    Edge<T, dim>::Edge(const std::array<T, dim> &pt1, const std::array<T, dim> &pt2) :
        m_curve{std::make_shared<BSCurve<T, dim>>(build_segment(pt1, pt2))},
        m_vtx1{std::make_shared<Vertex<T, dim>>(pt1)},
        m_vtx2{std::make_shared<Vertex<T, dim>>(pt2)},
        m_u1{m_curve->bounds().front()},
        m_u2{m_curve->bounds().back()}
    {
    }

    template <typename T, size_t dim>
    Edge<T, dim>::Edge(const std::shared_ptr<Vertex<T, dim>> &vtx1, const std::shared_ptr<Vertex<T, dim>> &vtx2) :
        m_curve{std::make_shared<BSCurve<T, dim>>(build_segment(vtx1->point(), vtx2->point()))},
        m_vtx1{vtx1},
        m_vtx2{vtx2},
        m_u1{m_curve->bounds().front()},
        m_u2{m_curve->bounds().back()}
    {
    }

    template <typename T, size_t dim>
    Edge<T, dim>::Edge(const std::shared_ptr<Curve<T, dim>> &crv) :
        m_curve{crv},
        m_vtx1{std::make_shared<Vertex<T, dim>>(crv->begin())},
        m_vtx2{std::make_shared<Vertex<T, dim>>(crv->end())},
        m_u1{m_curve->bounds().front()},
        m_u2{m_curve->bounds().back()}
    {
    }

    template <typename T, size_t dim>
    Edge<T, dim>::Edge(const std::shared_ptr<Curve<T, dim>> &crv, const BaseTopo<T, dim> &tolRef) :
        BaseTopo<T, dim>{tolRef},
        m_curve{crv},
        m_vtx1{std::make_shared<Vertex<T, dim>>(crv->begin())},
        m_vtx2{std::make_shared<Vertex<T, dim>>(crv->end())},
        m_u1{m_curve->bounds().front()},
        m_u2{m_curve->bounds().back()}
    {
    }

    template <typename T, size_t dim>
    const auto &Edge<T, dim>::vertex1() const
    {
        return m_vtx1;
    }

    template <typename T, size_t dim>
    const auto &Edge<T, dim>::vertex2() const
    {
        return m_vtx2;
    }

    template <typename T, size_t dim>
    const auto &Edge<T, dim>::curve() const
    {
        return m_curve;
    }

    template <typename T, size_t dim>
    auto Edge<T, dim>::bounds() const
    {
        return std::make_pair(m_u1, m_u2);
    }

    template <typename T, size_t dim>
    bool Edge<T, dim>::setVertex1(const std::shared_ptr<Vertex<T, dim>> &vtx)
    {
        auto [u, d] = extrema_curve_point(*m_curve, vtx->point(), knot_eps);

        if (d > this->approximationTolerance())
        {
            return false;
        }

        m_u1 = std::abs(u - m_curve->bounds()[1]) < knot_eps ? m_curve->bounds()[0] : u;
        m_vtx1 = vtx;

        return true;
    }

    template <typename T, size_t dim>
    bool Edge<T, dim>::setVertex2(const std::shared_ptr<Vertex<T, dim>> &vtx)
    {
        auto [u, d] = extrema_curve_point(*m_curve, vtx->point(), knot_eps);

        if (d > this->approximationTolerance())
        {
            return false;
        }

        m_u2 = std::abs(u - m_curve->bounds()[0]) < knot_eps ? m_curve->bounds()[1] : u;
        m_vtx2 = vtx;

        return true;
    }

    template <typename T, size_t dim>
    void Edge<T, dim>::tessellate()
    {
        auto points = discretize(*m_curve, 30, 0.05);

        points.front() = m_vtx1->point();
        points.back() = m_vtx2->point();

        // Clear and reserve space for the new he_vertices to minimize reallocations
        this->he_vertices.clear();
        this->he_vertices.reserve(points.size());

        std::transform(
            // std::execution::par,
            points.begin(), points.end(),
            std::back_inserter(this->he_vertices),
            [](const auto &pnt)
            {
                return std::make_shared<HalfEdgeVertex<T, dim>>(HalfEdgeVertex<T, dim>{pnt, nullptr});
            }
        );
    }
}
