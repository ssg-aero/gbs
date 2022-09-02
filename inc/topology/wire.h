#pragma once

#include "edge.h"

namespace gbs
{
    template <typename T, size_t dim>
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

        Wire() = default;
        Wire(const BaseTopo<T,dim> &tolRef)  = delete;
        Wire(const Edge<T,dim> &ed) : m_edges{std::make_shared< Edge<T,dim> >( ed )} 
        {
            fuseVertex2(m_edges.back()); // try to close
        }
        Wire(std::shared_ptr< Edge<T,dim> > &ed) : m_edges{ed} 
        {
            fuseVertex2(m_edges.back()); // try to close
        }

        /**
         * @brief happend/prepend edge if possible. ed can be modified to update merged vertex
         * 
         * @param ed 
         * @return true 
         * @return false 
         */
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
            this->m_msh.vertices.clear();
            this->m_msh.edges.clear();
            this->m_msh.faces.clear();

            // Store each edge vertices but tail
            std::for_each(
                m_edges.begin(), m_edges.end(),
                [&msh=this->m_msh](const auto &ed)
                {
                    ed->tessellate();
                    msh.vertices.insert(
                        msh.vertices.end(),
                        std::next(ed->mesh().vertices.begin()),
                        std::next(ed->mesh().vertices.end(),-1)
                    );
                }
            );

            if(isClosed())// Build face
            {
                
                auto face = std::make_shared<HalfEdgeFace<T,dim>>();

                auto &msh = this->m_msh;

                // fill half edges loop
                auto n =  msh.vertices.size(); // initial vertices number
                msh.edges.resize(n);

                msh.vertices.push_back( msh.vertices.front() );// use head vertex for tail half edge

                std::transform(
                    std::next(msh.vertices.begin()), 
                    msh.vertices.end(),
                    msh.edges.begin(),
                    [&face](const auto &vtx)
                    {
                        HalfEdge<T,dim> hEd{
                                .vertex = vtx,
                                .face   = face,
                        };
                        return std::make_shared<HalfEdge<T,dim>>(hEd);
                    }
                );

                // build half edges connections
                for(size_t i{}; i < n; i++)
                {
                    auto next     = i == n-1 ? 0 : i+1;
                    auto previous = i == 0 ? n-1 : i-1;
                    msh.edges[i]->next = msh.edges[next];
                    msh.edges[i]->previous = msh.edges[previous];
                }

                // Store fce and cleanup data
                face->edge = msh.edges.front();
                msh.faces.push_back( face );
                msh.vertices.clear(); // no longer required

            }
            else
            {

            }




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

}