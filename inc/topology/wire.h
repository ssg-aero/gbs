#pragma once
#include <list>
#include "edge.h"

namespace gbs
{

    /**
     * @brief Wire class representing a sequence of connected edges.
     * 
     * @tparam T Floating-point type for coordinates.
     * @tparam dim Dimension of the space.
     */
    template <std::floating_point T, size_t dim>
    class Wire : public BaseTopo<T, dim>
    {
        std::list< std::shared_ptr< Edge<T,dim> > > m_edges;
        bool fuseVertex1(std::shared_ptr<Edge<T, dim>> &ed);
        bool fuseVertex2(std::shared_ptr<Edge<T, dim>> &ed);
    public:
        /// @brief Default constructor deleted.
        Wire() = delete;

        /// @brief Constructor with tolerance reference deleted.
        Wire(const BaseTopo<T, dim> &tolRef) = delete;

        /**
         * @brief Constructor with a single edge.
         * 
         * @param ed An edge to initialize the wire.
         */
        Wire(const Edge<T, dim> &ed);

        /**
         * @brief Constructor with a shared pointer to an edge.
         * 
         * @param ed Shared pointer to an edge to initialize the wire.
         */
        Wire(std::shared_ptr<Edge<T, dim>> &ed);

        /**
         * @brief Add an edge to the wire using a shared pointer.
         * 
         * @param ed Shared pointer to the edge to add.
         * @return true If the edge was added successfully.
         * @return false If the edge was not added.
         */
        bool addEdge(std::shared_ptr<Edge<T, dim>> &ed);

        /**
         * @brief Add an edge to the wire using an edge reference.
         * 
         * @param ed Reference to the edge to add.
         * @return true If the edge was added successfully.
         * @return false If the edge was not added.
         */
        bool addEdge(const Edge<T, dim> &ed);

        // Accessor methods
        [[nodiscard]] const auto begin() const;
        [[nodiscard]] const auto end() const;
        [[nodiscard]] const auto &front() const;
        [[nodiscard]] const auto &back() const;
        [[nodiscard]] const auto &vertex1() const;
        [[nodiscard]] const auto &vertex2() const;

        /**
         * @brief Check if the wire is closed.
         * 
         * @return true If the wire is closed.
         * @return false If the wire is open.
         */
        [[nodiscard]] bool isClosed() const;

        /**
         * @brief Close the wire by connecting the last edge to the first one.
         */
        void close();

        /**
         * @brief Tessellate the wire to generate a mesh representation.
         */
        void tessellate() override;
    };

    template <std::floating_point T, size_t dim>
    bool Wire<T, dim>::fuseVertex1(std::shared_ptr<Edge<T, dim>> &ed)
    {
        if (distance(
                m_edges.back()->vertex2()->point(),
                ed->vertex1()->point()) < this->approximationTolerance())
        {
            ed->setVertex1(m_edges.back()->vertex2());
            return true;
        }
        return false;
    }

    template <std::floating_point T, size_t dim>
    bool Wire<T, dim>::fuseVertex2(std::shared_ptr<Edge<T, dim>> &ed)
    {
        if (distance(
                m_edges.front()->vertex1()->point(),
                ed->vertex2()->point()) < this->approximationTolerance())
        {
            ed->setVertex2(m_edges.front()->vertex1());
            return true;
        }
        return false;
    }

    template <std::floating_point T, size_t dim>
    Wire<T, dim>::Wire(const Edge<T, dim> &ed) : m_edges{std::make_shared<Edge<T, dim>>(ed)}
    {
        fuseVertex2(m_edges.back()); // try to close
    }

    template <std::floating_point T, size_t dim>
    Wire<T, dim>::Wire(std::shared_ptr<Edge<T, dim>> &ed) : m_edges{ed}
    {
        fuseVertex2(m_edges.back()); // try to close
    }

    template <std::floating_point T, size_t dim>
    bool Wire<T, dim>::addEdge(std::shared_ptr<Edge<T, dim>> &ed)
    {
        if (m_edges.size() == 0)
        {
            m_edges.push_back(ed);
        }
        else
        {
            if (isClosed())
            {
                return false;
            }
            if (m_edges.front()->vertex1() == ed->vertex2())
            {
                m_edges.push_front(ed);
            }
            else if (m_edges.back()->vertex2() == ed->vertex1())
            {
                m_edges.push_back(ed);
            }
            else
            {
                bool added = false;
                if (distance(
                        m_edges.back()->vertex2()->point(),
                        ed->vertex1()->point()) < this->approximationTolerance())
                {
                    ed->setVertex1(m_edges.back()->vertex2());
                    m_edges.push_back(ed);
                    added = true;
                }
                if (distance(
                        m_edges.front()->vertex1()->point(),
                        ed->vertex2()->point()) < this->approximationTolerance())
                {
                    ed->setVertex2(m_edges.front()->vertex1());
                    if (!added)
                    {
                        m_edges.push_front(ed);
                        added = true;
                    }
                }
                return added;
            }
        }

        return true;
    }

    template <std::floating_point T, size_t dim>
    bool Wire<T, dim>::addEdge(const Edge<T, dim> &ed)
    {
        auto p_ed = std::make_shared<Edge<T, dim>>(ed);
        return addEdge(p_ed);
    }

    template <std::floating_point T, size_t dim>
    const auto Wire<T, dim>::begin() const
    {
        return m_edges.begin();
    }

    template <std::floating_point T, size_t dim>
    const auto Wire<T, dim>::end() const
    {
        return m_edges.end();
    }

    template <std::floating_point T, size_t dim>
    const auto &Wire<T, dim>::front() const
    {
        return m_edges.front();
    }

    template <std::floating_point T, size_t dim>
    const auto &Wire<T, dim>::back() const
    {
        return m_edges.back();
    }

    template <std::floating_point T, size_t dim>
    const auto &Wire<T, dim>::vertex1() const
    {
        return front()->vertex1();
    }

    template <std::floating_point T, size_t dim>
    const auto &Wire<T, dim>::vertex2() const
    {
        return back()->vertex2();
    }

    template <std::floating_point T, size_t dim>
    bool Wire<T, dim>::isClosed() const
    {
        if (m_edges.size())
        {
            return vertex1() == vertex2();
        }
        else
        {
            return false;
        }
    }

    template <std::floating_point T, size_t dim>
    void Wire<T, dim>::close()
    {
        addEdge(Edge<T, dim>{
            end()->vertex2(),
            begin()->vertex1()});
    }

    template <std::floating_point T, size_t dim>
    void Wire<T, dim>::tessellate()
    {
        std::ranges::for_each(m_edges, [](const auto &ed)
                              { ed->tessellate(); });
    }

    /*
    template <std::floating_point T, size_t dim>
    class Wire : public BaseTopo<T,dim>
    {

        std::list< std::shared_ptr< Edge<T,dim> > > m_edges;

        bool fuseVertex1(std::shared_ptr<Edge<T, dim>> &ed)
        {
            if (distance(
                    m_edges.back()->vertex2()->point(),
                    ed->vertex1()->point()) < this->approximationTolerance())
            {
                ed->setVertex1(m_edges.back()->vertex2());
                return true;
            }
            return false;
        }

        bool fuseVertex2(std::shared_ptr<Edge<T, dim>> &ed)
        {
            if (distance(
                    m_edges.front()->vertex1()->point(),
                    ed->vertex2()->point()) < this->approximationTolerance())
            {
                ed->setVertex2(m_edges.front()->vertex1());
                return true;
            }
            return false;
        }

        public:

        Wire() = delete;
        Wire(const BaseTopo<T,dim> &tolRef)  = delete;
        Wire(const Edge<T,dim> &ed) : m_edges{std::make_shared< Edge<T,dim> >( ed )} 
        {
            fuseVertex2(m_edges.back()); // try to close
        }
        Wire(std::shared_ptr< Edge<T,dim> > &ed) : m_edges{ed} 
        {
            fuseVertex2(m_edges.back()); // try to close
        }

        bool addEdge( std::shared_ptr< Edge<T,dim> > &ed)
        {
            if(m_edges.size() == 0)
            {
                m_edges.push_back( ed );
            }
            else
            {
                if(isClosed())
                {
                    return false;
                }
                if( m_edges.front()->vertex1() == ed->vertex2())
                {
                    m_edges.push_front( ed );
                }
                else if( m_edges.back()->vertex2() == ed->vertex1())
                {
                    m_edges.push_back( ed );
                }
                else
                {
                    bool added = false;
                    if( distance(
                        m_edges.back()->vertex2()->point(), 
                        ed->vertex1()->point() ) < this->approximationTolerance() )
                    {
                        ed->setVertex1( m_edges.back()->vertex2() );
                        m_edges.push_back( ed );
                        added = true;
                    }
                    if( distance(
                        m_edges.front()->vertex1()->point(), 
                        ed->vertex2()->point() ) < this->approximationTolerance() )
                    {
                        ed->setVertex2( m_edges.front()->vertex1() );
                        if(!added)
                        {
                            m_edges.push_front( ed );
                            added = true;
                        }
                    }
                    return added;
                }
            }

            return true;
        }

        bool addEdge(const Edge<T,dim> &ed )
        {
            auto p_ed = std::make_shared<Edge<T,dim>>(ed);
            return addEdge( p_ed );
        }

        const auto begin() const
        {
            return m_edges.begin();
        }

        const auto end() const
        {
            return m_edges.end();
        }

        const auto & front() const
        {
            return m_edges.front();
        }

        const auto & back() const
        {
            return m_edges.back();
        }

        const auto & vertex1() const 
        {
            return front()->vertex1();
        }

        const auto & vertex2() const 
        {
            return back()->vertex2();
        }

        bool isClosed() const
        {
            if(m_edges.size())
            {
                return  vertex1() == vertex2();
            }
            else
            {
                return false;
            }
        }

        void close()
        {
            addEdge( Edge<T,dim>{
                end()->vertex2(),
                begin()->vertex1() });
        }

        void tessellate() override
        {
            std::ranges::for_each(m_edges,[](const auto &ed){ed->tessellate();});
        }
    };

    template <typename T>
    class Wire2d : public Wire<T,2>
    {
        using Wire<T,2>::Wire;
        public:
        // Wire();
        // Wire(const BaseTopo<T,2> &tolRef);
        // Wire(const Edge<T,2> &ed);
        // Wire(const std::shared_ptr< Edge<T,2> > &ed);
        // bool addEdge( std::shared_ptr< Edge<T,2> > &ed);
        // bool addEdge(const Edge<T,2> &ed );
        // const auto begin() const;
        // const auto end() const;
        // const auto & front() const;
        // const auto & back() const;
        // const auto & vertex1() const;
        // const auto & vertex2() const;
        // bool isClosed() const;
        // void close();
        bool isInside(const std::array<T,2> &pnt) const
        {
            size_t intersections_count{};
            
            return intersections_count % 2 == 0;

        }
    };
*/
}