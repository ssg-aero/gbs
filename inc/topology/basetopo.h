#pragma once

#include "tessellations.h"

namespace gbs
{
    template <typename T, size_t dim>
    class BaseTopo
    {
        T m_precision;
        T m_approximation;
        
        protected:

        HalfEdgeMesh<T,dim> m_msh;

        public:

        BaseTopo() : m_precision{1e-6} , m_approximation{1e-5} {}

        BaseTopo(T tolPrecision, T tolApproximation) : m_precision{tolPrecision} , m_approximation{tolApproximation} {}

        BaseTopo(const BaseTopo<T,dim> &tolRef ) = default;

        void setPrecisionTolerance(T tol)
        {
            m_precision = tol;
        }

        void setApproximationTolerance(T tol)
        {
            m_approximation = tol;
        }

        T precisionTolerance() const noexcept
        {
            return m_precision;
        }

        T approximationTolerance() const noexcept
        {
            return m_approximation;
        }

        virtual void tessellate() = 0;

        const auto & mesh() const 
        {
            return m_msh;
        }

        // virtual bool near(const BaseTopo<T,d> &other) = 0;
    };
}