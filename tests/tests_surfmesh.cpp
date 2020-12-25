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

template <typename T, size_t dim>
class mesh
{
    std::list<Vertex<T, dim>> vt_;
    std::list<Face> fa_;
    std::list<HalfEdge> he_;

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
    auto vertex(index_t i) -> Vertex<T, dim>
    {
        return *std::next(vt_.begin(), i);
    }
    auto vertex(index_t i) const -> Vertex<T, dim>
    {
        return *std::next(vt_.begin(), i);
    }
    auto n_face() const -> index_t
    {
        return fa_.size();
    }
    auto face(index_t i) -> Face
    {
        return *std::next(fa_.begin(), i);
    }
    auto face(index_t i) const -> Face
    {
        return *std::next(fa_.begin(), i);
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
    auto half_edge(index_t i) -> HalfEdge
    {
        return *std::next(he_.begin(), i);
    }
    auto half_edge(index_t i) const -> HalfEdge
    {
        return *std::next(he_.begin(), i);
    }
    auto add_vertex(const Vertex<T, dim> &vtx) -> index_t
    {
        vt_.push_back(vtx);
        return index_t(vt_.size() - 1);
    }
    auto add_face(const std::vector<index_t> &vtx_lst)
    {
        if (vtx_lst.size() < 2)
            std::exception("at least 3 vertex needed for a face!");

        std::vector<index_t> vtx_lst_{vtx_lst};
        vtx_lst_.push_back(vtx_lst.front());

        auto n_he = vtx_lst.size();
        std::vector<HalfEdge> he_vec(n_he);

        index_t iF{fa_.size()};
        index_t he_begin{he_.size()};

        std::transform(
            vtx_lst.begin(),
            vtx_lst.end(),
            std::next(vtx_lst_.begin(), 1),
            he_vec.begin(),
            [&, ihe_ = he_.size()](index_t v1, index_t v2) mutable {

                index_t i_next = (n_he + ihe_ + 1) % n_he + he_begin;
                index_t i_prev = (n_he + ihe_ - 1) % n_he + he_begin;

                auto e_opp = (*std::next(vt_.begin(), v1)).e;
                index_t i_opp = invalid_index;
                if (e_opp != invalid_index)
                {
                    auto iF_opp = (*std::next(he_.begin(), e_opp)).f;
                    if (iF_opp != iF)
                    {
                        i_opp = (*std::next(vt_.begin(), v2)).e;
                    }
                }

                (*std::next(vt_.begin(), v1)).e = ihe_;
                ihe_++;
                return HalfEdge{.v = v2, .f = iF, .next = i_next, .prev = i_prev, .opp = i_opp};
            });

        fa_.push_back(Face{.e = he_.size()});
        he_vec.front().f = iF;

        he_.insert(he_.end(), he_vec.begin(), he_vec.end());
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

TEST(tests_surfmesh, tow_tri)
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
auto make_vtkPoint(const Vertex<T, dim> &vtx) -> std::array<double, 3>
{
    std::array<double, 3> x = {0., 0., 0.};
    for (auto i = 0; i < fmin(dim, 3); i++)
        x[i] = vtx.p[i];
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

    msh.add_face({0, 1, 2});
    msh.add_face({0, 2, 3});

    msh.split_face(1);
    msh.split_face(2);

    auto points = gbs::make_vtkPoints(msh.vertex_begin(),msh.vertex_end());
    std::vector<std::array<vtkIdType, 3>> pts_tri;
    std::for_each(
        msh.face_begin(), msh.face_end(),
        [&](const auto &fa) {
            auto index = 0;
            auto first_he = fa.e;
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