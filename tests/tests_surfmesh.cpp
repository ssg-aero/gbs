#include <gtest/gtest.h>
#include <limits>
#include <gbs/gbslib.h>
#include <gbs-render/vtkcurvesrender.h>

using namespace gbs;

namespace {
    const double tol = 1e-6;
    const double PI = acos(-1.);
    using index_t = size_t;
    // using index_t = unsigned long;
    // using index_t = unsigned int;
    // using index_t = unsigned short;
    const index_t invalid_index = std::numeric_limits<index_t>::max();

}

struct HalfEdge
{
    index_t v = invalid_index;
    index_t f = invalid_index;
    index_t next = invalid_index;
    index_t prev = invalid_index;
    index_t opp  = invalid_index;
};

template <typename T, size_t dim>
struct Vertex
{
    gbs::point<T, dim> p;
    index_t e = invalid_index;
};

struct Face
{
    index_t e;
};

// template <typename T, size_t dim>
// using vertex_it = typename std::unordered_map<index_t, Vertex<T, dim>>::iterator;

template <typename T, size_t dim>
class mesh
{
    std::unordered_map<index_t, Vertex<T, dim>> vt_;
    std::unordered_map<index_t, Face> fa_;
    std::unordered_map<index_t, HalfEdge> he_;

    index_t last_vt_free_id_ = 0;
    index_t last_fa_free_id_ = 0;
    index_t last_he_free_id_ = 0;

public:
    auto vertex_begin() const noexcept
    {
        return vt_.begin();
    }
    auto vertex_end() const noexcept
    {
        return vt_.end();
    }
    auto n_vertex() const -> index_t
    {
        return vt_.size();
    }
    auto vertex(index_t i) -> Vertex<T, dim> &
    {
        return vt_.find(i)->second;
    }
    auto vertex(index_t i) const -> Vertex<T, dim> 
    {
        return vt_.find(i)->second;
    }
    auto n_face() const -> index_t
    {
        return fa_.size();
    }
    auto face(index_t i) -> Face &
    {
        return fa_.find(i)->second;
    }
    auto face(index_t i) const -> Face 
    {
        return fa_.find(i)->second;
    }
    auto face_dim(index_t i) const -> index_t
    {
        return face_vertices_id(i).size();
    }
    auto face_begin() const noexcept
    {
        return fa_.begin();
    }
    auto face_end() const noexcept
    {
        return fa_.end();
    }
    auto n_half_edge() const -> index_t
    {
        return he_.size();
    }
    auto half_edge(index_t ihe) -> HalfEdge &
    {
        return he_.find(ihe)->second;
    }
    auto half_edge(index_t ihe) const -> HalfEdge
    {
        return he_.find(ihe)->second;
    }
    auto is_half_edge_bound(index_t ihe) const -> bool
    {
        return he_.find(ihe)->second.opp == invalid_index;
    }
    auto half_edge_start(index_t ihe) const -> index_t
    {
        return he_.find(he_.find(ihe)->second.prev)->second.v;
    }
    auto half_edge_end(index_t ihe) const -> index_t
    {
        return he_.find(ihe)->second.v;
    }
    auto are_half_edges_opposite(index_t ihe1,index_t ihe2) const -> bool
    {
        return (half_edge_start(ihe1) == half_edge_end(ihe2)) &&
               (half_edge_start(ihe2) == half_edge_end(ihe1));
    }
    auto vertex_secondary_half_edges(index_t ivt) const -> std::list<index_t>
    {

        std::list<index_t> he_lst;

        auto i_first_he = vertex(ivt).e;
        if (i_first_he == invalid_index)
        {
            return he_lst;
        }

        auto i_p_ed = half_edge(half_edge(i_first_he).next).opp;
        while ((i_p_ed != i_first_he)&&(i_p_ed!=invalid_index))
        {
            i_p_ed = half_edge(half_edge(i_p_ed).next).opp;
            if (i_p_ed != invalid_index)
            {
                he_lst.push_back(i_p_ed);
            }
            else
            {
                break;
            }
        }
        return he_lst;
    }

    auto add_vertex(const Vertex<T, dim> &vtx) -> index_t
    {
        vt_.insert({last_vt_free_id_,vtx});
        return last_vt_free_id_++;
    }
    auto add_face(const std::vector<index_t> &vtx_lst)
    {
        if (vtx_lst.size() < 2)
            std::invalid_argument("at least 3 vertex needed for a face!");

        std::vector<index_t> vtx_lst_{vtx_lst};
        vtx_lst_.push_back(vtx_lst.front());

        auto n_he = vtx_lst.size();
        auto ihe_    { last_he_free_id_};
        auto iF      {last_fa_free_id_};
        auto he_begin{last_he_free_id_};
        for (auto v1 = 0; v1 < n_he; v1++)// create and insert edges
        {
            auto v2 = (v1 + 1) % n_he; // for cyclic loop on provides vtx indices

            index_t i_next = (n_he + ihe_ + 1) % n_he + he_begin; // for cyclic loop on created edges
            index_t i_prev = (n_he + ihe_ - 1) % n_he + he_begin;
            index_t i_curr = ihe_ + he_begin;

            index_t i_opp = invalid_index;

            /*
            auto e_opp = vertex(vtx_lst[v1]).e;
            
            if (e_opp != invalid_index)
            {
                auto iF_opp = half_edge(e_opp).f;
                if (iF_opp != iF) // check if has a previous face sharing edge's vertices
                {
                    i_opp = vertex(vtx_lst[v1]).e;
                    if (half_edge(i_opp).opp == invalid_index)
                    {
                        half_edge(i_opp).opp = last_he_free_id_;
                    }
                }
            }
            */

            if (vertex(vtx_lst[v2]).e == invalid_index) // associate if not allready done
            {
                vertex(vtx_lst[v2]).e = ihe_;
            }

            he_.insert({last_he_free_id_++, {.v = vtx_lst[v2], .f = iF, .next = i_next, .prev = i_prev, .opp = i_opp}});

            ihe_++;
        }

        fa_.insert({last_fa_free_id_++,Face{.e =he_begin}}); // insert face

        auto fa_he = face_half_edges_id(last_fa_free_id_-1);
        for (auto i_curr : fa_he)
        {
            for (auto v : vtx_lst)
            {
                auto vt_edges = vertex_secondary_half_edges(v);
                for (auto he_ : vt_edges)
                {
                    if (are_half_edges_opposite(i_curr, he_))
                    {
                        if (!is_half_edge_bound(he_))
                        {
                            std::cerr << "non mainfoil!";
                        }
                        half_edge(he_).opp = i_curr;
                        half_edge(i_curr).opp = he_;
                    }
                }
            }
        }

        return last_fa_free_id_;
    }

    auto face_vertices_id(index_t iF) const 
    {
        auto fa = face(iF);
        auto first_he = fa.e;
        auto i_he = first_he;
        std::vector<index_t> vertices_;
        do
        {
            vertices_.push_back(  half_edge(i_he).v  );
            i_he = half_edge(i_he).next;
        } while (i_he != first_he);
        return vertices_;
    }

    auto face_vertices(index_t iF) const 
    {
        auto fa = face(iF);
        auto first_he = fa.e;
        auto i_he = first_he;
        std::vector<Vertex<T,dim>> vertices_;
        do
        {
            vertices_.push_back( vertex( half_edge(i_he).v ) );
            i_he = half_edge(i_he).next;
        } while (i_he != first_he);
        return vertices_;
    }

    auto face_half_edges(index_t iF) const
    {
        auto fa = face(iF);
        auto first_he = fa.e;
        auto i_he = first_he;
        std::vector<HalfEdge> half_edges_;
        do
        {
            auto he__ = half_edge(i_he);
            half_edges_.push_back(he__);
            i_he = he__.next;
        } while (i_he != first_he);
        return half_edges_;
    }

    auto face_half_edges_id(index_t iF) const
    {
        auto fa = face(iF);
        auto first_he = fa.e;
        auto i_he = first_he;
        std::vector<index_t> half_edges_;
        do
        {
            half_edges_.push_back(i_he);
            i_he = half_edge(i_he).next;
        } while (i_he != first_he);
        return half_edges_;
    }

    auto split_face(index_t iF)
    {
        auto vertices_ = face_vertices_id(iF);
        point<T, dim> c{T(0.)};
        auto w = T(vertices_.size());
        std::for_each(vertices_.begin(), vertices_.end(),
                      [&](const auto &vtx_id) {auto vtx_ = vertex( vtx_id ) ;c = c + vtx_.p; });
        // c = std::reduce(vertices_.begin(),vertices_.end());
        c = c / w;

        auto c_id = add_vertex(Vertex{.p=c});
        auto half_edges_ = face_half_edges(iF);
        std::for_each(
            half_edges_.begin(),half_edges_.end(),
            [&](const auto & he__)
            {
                add_face({
                    he__.v,
                    c_id,
                    half_edge(he__.prev).v
                });
            }
        );

    }

    auto tag_loop(const std::vector<index_t> &he_loop,index_t iF)
    {
        auto n_he = he_loop.size();
        for (auto he = 0; he < n_he; he++)// create and insert edges
        {
            auto he_next = he_loop[(n_he + he + 1) % n_he]; // for cyclic loop on provides indices
            auto he_prev = he_loop[(n_he + he - 1) % n_he]; 
            auto he_curr = he_loop[he];
            half_edge(he_curr).prev = he_prev;
            half_edge(he_curr).next = he_next;
            half_edge(he_curr).f = iF;
        }
    }

    auto flip_edge(index_t iE) // works only for triangle
    {
        auto he__ = half_edge(iE);
        if(he__.opp!=invalid_index)
        {
            auto he__opp = half_edge(he__.opp);
            // edges edition
            half_edge(iE).v = half_edge(he__.next).v;
            half_edge(he__.opp).v=half_edge(he__opp.next).v;
            // faces tag
            face(he__.f).e = iE;
            face(he__opp.f).e = he__.opp;
            // rebuild loop
            tag_loop({iE,half_edge(iE).prev,half_edge(he__.opp).next},he__.f);
            tag_loop({he__.opp,half_edge(he__.opp).prev,half_edge(iE).next},he__.f);
        }
    }
};

TEST(tests_surfmesh, limits)
{
    std::cerr << "index type:" << typeid(invalid_index).name() << std::endl;
    std::cerr << "index max:" << invalid_index-1 << std::endl;
}

TEST(tests_surfmesh, one_tri)
{
    mesh<double,2> msh;
    std::vector<index_t> vth_h {
        msh.add_vertex({0.,0.}),
        msh.add_vertex({1.,0.}),
        msh.add_vertex({1.,1.}),
    };
    msh.add_face(vth_h);
    ASSERT_EQ(msh.n_vertex(),3);
    ASSERT_EQ(msh.n_half_edge(),3);
    ASSERT_EQ(msh.n_face(),1);
    ASSERT_EQ(msh.half_edge(1).prev,0);
    ASSERT_EQ(msh.half_edge(1).next,2);
}

TEST(tests_surfmesh,  DISABLED_tow_tri)
{
    mesh<double, 2> msh;

    msh.add_vertex({0., 0.});
    msh.add_vertex({1., 0.});
    msh.add_vertex({1., 1.});
    msh.add_vertex({0., 1.});

    msh.add_face({0, 1, 2});
    msh.add_face({0, 2, 3});

    ASSERT_EQ(msh.half_edge(3).opp, 2);
}

template <typename T, size_t dim>
auto make_vtkPoint(const std::pair<const index_t,Vertex<T,dim>> &vtx) -> std::array<double, 3>
{
    std::array<double, 3> x = {0., 0., 0.};
    for (auto i = 0; i < fmin(dim, 3); i++)
        x[i] = vtx.second.p[i];
    return x;
}
#include <vtkXMLPolyDataWriter.h>
TEST(tests_surfmesh, tow_tri_vtk)
{
    mesh<double, 2> msh;

    msh.add_vertex({0., 0.});
    msh.add_vertex({1., 0.});
    msh.add_vertex({1., 1.});
    msh.add_vertex({0., 1.});
    // msh.add_vertex({2., 0.});
    // msh.add_vertex({2., 1.});

    msh.add_face({0, 1, 2});
    msh.add_face({0, 2, 3});
    // msh.add_face({1, 4, 5});
    // msh.add_face({1, 5, 2});

    // msh.split_face(1);
    // msh.split_face(4);

    // msh.vertex(2).p={1.,1.3};

    msh.flip_edge(2);

    auto points = gbs::make_vtkPoints(msh.vertex_begin(),msh.vertex_end());
    std::vector<std::array<vtkIdType, 3>> pts_tri;
    std::for_each(
        msh.face_begin(), msh.face_end(),
        [&](const auto &fa) {
            auto index = 0;
            auto first_he = fa.second.e;
            auto i_he = first_he;
            std::array<vtkIdType, 3> tri;
            do
            {
                tri[index] = msh.half_edge(i_he).v; // pas optimal
                index++;
                i_he = msh.half_edge(i_he).next;
            }while (i_he != first_he);
            pts_tri.push_back(tri);
        });

    vtkSmartPointer<vtkPolyData> ugrid =
        vtkSmartPointer<vtkPolyData>::New();
    ugrid->Allocate(pts_tri.size());

    std::for_each(pts_tri.begin(), pts_tri.end(),
                //   [&ugrid](const auto &tri) { ugrid->InsertNextCell(VTK_POLYGON, 3, tri.data()); });
                [&ugrid](const auto &tri) { ugrid->InsertNextCell(VTK_TRIANGLE, 3, tri.data()); });

    ugrid->SetPoints(points);

    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("tow_tri_vtk.vtp");

    writer->SetInputData(ugrid);

    writer->Write();
}