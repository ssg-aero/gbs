#pragma once


#include <gbs/curves>
#include <gbs/bscbuild.h>
#include "vertex.h"

namespace gbs
{
    template <typename T, size_t dim>
    class Edge : public BaseTopo<T,dim>
    {
        std::shared_ptr<Curve<T,dim>> m_curve;
        std::shared_ptr<Vertex<T,dim>> m_vtx1;
        std::shared_ptr<Vertex<T,dim>> m_vtx2;

        T m_u1;
        T m_u2;

        public:

        Edge() = delete;

        Edge(const std::array<T,dim> &pt1, const std::array<T,dim> &pt2) :
            m_curve{std::make_shared<BSCurve<T,dim>>( build_segment(pt1, pt2) ) },
            m_vtx1{std::make_shared<Vertex<T,dim>>( pt1 )},
            m_vtx2{std::make_shared<Vertex<T,dim>>( pt2 )},
            m_u1{ m_curve->bounds().front() },
            m_u2{ m_curve->bounds().back() }
        {}

        Edge(const std::shared_ptr<Vertex<T,dim>> &vtx1, const std::shared_ptr<Vertex<T,dim>> &vtx2) :
            m_curve{std::make_shared<BSCurve<T,dim>>( build_segment(vtx1->point(), vtx2->point()) ) },
            m_vtx1{vtx1},
            m_vtx2{vtx2},
            m_u1{ m_curve->bounds().front() },
            m_u2{ m_curve->bounds().back() }
        {}

        Edge(const std::shared_ptr<Curve<T,dim>> &crv) : 
            m_curve{crv},
            m_vtx1{std::make_shared<Vertex<T,dim>>( crv->begin()) },
            m_vtx2{std::make_shared<Vertex<T,dim>>( crv->end()) },
            m_u1{ m_curve->bounds().front() },
            m_u2{ m_curve->bounds().back() }
        {}

        Edge(const std::shared_ptr<Curve<T,dim>> &crv, const BaseTopo<T,dim> &tolRef) : 
            m_curve{crv},
            m_vtx1{std::make_shared<Vertex<T,dim>>( crv->begin()) },
            m_vtx2{std::make_shared<Vertex<T,dim>>( crv->end()) },
            m_u1{ m_curve->bounds().front() },
            m_u2{ m_curve->bounds().back() },
            BaseTopo<T,dim>{tolRef} 
        {}

        const auto & vertex1() const 
        {
            return m_vtx1;
        }

        const auto & vertex2() const 
        {
            return m_vtx2;
        }

        const auto & curve() const 
        {
            return m_curve;
        }

        auto bounds() const
        {
            return std::make_pair(m_u1, m_u2);
        }

        bool setVertex1(const std::shared_ptr<Vertex<T,dim>> &vtx )
        {
            auto [u,d] = extrema_curve_point(*m_curve, vtx->point(), knot_eps );

            if(d > this->approximationTolerance())
            {
                return false;
            }
            
            m_u1   = std::abs(u-m_curve->bounds()[1]) < knot_eps ? m_curve->bounds()[0] : u;
            m_vtx1 = vtx;

            return true;
        }

        bool setVertex2(const std::shared_ptr<Vertex<T,dim>> &vtx )
        {
            auto [u,d] = extrema_curve_point(*m_curve, vtx->point(), knot_eps );

            if(d > this->approximationTolerance())
            {
                return false;
            }
            
            m_u2   = std::abs(u-m_curve->bounds()[0]) < knot_eps ? m_curve->bounds()[1] : u;
            m_vtx2 = vtx;

            return true;

        }

        void tessellate() override
        {
            this->m_msh.vertices.clear();

            auto points = discretize(*m_curve, 30, 0.05);

            points.front() = m_vtx1->point();
            points.back() = m_vtx2->point();

           this-> m_msh.vertices.resize(points.size());
            std::transform(
                std::execution::par,
                points.begin(), points.end(),
                this->m_msh.vertices.begin(),
                [](const auto &pnt)
                {
                    return  std::make_shared<HalfEdgeVertex<T,dim>>(HalfEdgeVertex<T,dim>{pnt,nullptr});
                }
            );

        }

    };
}