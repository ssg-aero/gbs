#pragma once

#include <array>

#include "basetopo.h"

namespace gbs
{
    template <typename T, size_t dim>
    class Vertex : public BaseTopo<T,dim>
    {
        std::array<T,dim> m_point;
        public:
        Vertex() = default;
        Vertex(T tolPrecision, T tolApproximation) : BaseTopo<T,dim>{tolPrecision,tolApproximation} {};
        Vertex(const std::array<T,dim> &pnt) : m_point{pnt} {}
        Vertex(const std::array<T,dim> &pnt, const BaseTopo<T,dim> &tolRef) : m_point{pnt}, BaseTopo<T,dim>{tolRef} {}
        const auto & point() const noexcept 
        {
            return m_point;
        }
        void setPoint(const std::array<T,dim> &pnt)
        {
            m_point = pnt;
        }

        void tessellate() override
        {
            this->m_msh.vertices.clear();
            this->m_msh.vertices.push_back( 
                std::make_shared<HalfEdgeVertex<T,dim>>( 
                    HalfEdgeVertex<T,dim>{m_point,nullptr} 
                )
            );
        }
    };
}