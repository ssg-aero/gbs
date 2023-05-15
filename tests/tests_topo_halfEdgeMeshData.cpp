#include <gtest/gtest.h>
#include <topology/halfEdgeMeshData.h>

using namespace gbs;
template <std::floating_point T, size_t dim>
using Vertex = HalfEdgeVertex<T, dim>;
template <std::floating_point T, size_t dim>
using Face = HalfEdgeFace<T, dim>;
template <std::floating_point T, size_t dim>
using Edge = HalfEdge<T, dim>;

TEST(HalfEdgeMeshDataTest, MakeSharedHEdgeVertices) {
    // Test case 1: Creating a half-edge with a vertex only
    {
        auto vertex = std::make_shared<Vertex<double, 3>>(Vertex<double, 3>{{0.0, 0.0, 0.0}});
        auto hedge = make_shared_h_edge(vertex);

        // Test if the half-edge has the correct vertex
        ASSERT_EQ(hedge->vertex, vertex);

        // Test if the vertex has the correct half-edge assigned
        ASSERT_EQ(vertex->edge, hedge);

        // Test if the half-edge does not have a face assigned
        ASSERT_EQ(hedge->face, nullptr);
    }

    // Test case 2: Creating a half-edge with a vertex and a face
    {
        auto vertex = std::make_shared<Vertex<double, 3>>(Vertex<double, 3>{{1.0, 1.0, 1.0}});
        auto face = std::make_shared<Face<double, 3>>();
        auto hedge = make_shared_h_edge(vertex, face);

        // Test if the half-edge has the correct vertex
        ASSERT_EQ(hedge->vertex, vertex);

        // Test if the vertex has the correct half-edge assigned
        ASSERT_EQ(vertex->edge, hedge);

        // Test if the half-edge has the correct face assigned
        ASSERT_EQ(hedge->face, face);
    }
    // Test case 3: Creating a vertex with given coordinates
    {
        std::array<double, 3> coords = {2.0, 2.0, 2.0};
        auto vertex = make_shared_h_vertex(coords);

        // Test if the vertex has the correct coordinates
        ASSERT_EQ(vertex->coords, coords);

        // Test if the vertex does not have an edge assigned
        ASSERT_EQ(vertex->edge, nullptr);
    }

    // Test case 4: Creating a half-edge with given vertex coordinates
    {
        std::array<double, 3> coords = {3.0, 3.0, 3.0};
        auto hedge = make_shared_h_edge(coords);

        // Test if the half-edge's vertex has the correct coordinates
        ASSERT_EQ(hedge->vertex->coords, coords);

        // Test if the vertex has the correct half-edge assigned
        ASSERT_EQ(hedge->vertex->edge, hedge);

        // Test if the half-edge does not have a face assigned
        ASSERT_EQ(hedge->face, nullptr);
    }
    // Test case 5: Creating a vector of vertices from a container of coordinate arrays
    {
        std::vector<std::array<double, 3>> coords = {
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0},
            {2.0, 2.0, 2.0}};
        auto vertices = make_shared_h_vertices<double, 3>(coords);

        // Test if the correct number of vertices is created
        ASSERT_EQ(vertices.size(), coords.size());

        // Test if each vertex has the correct coordinates
        for (size_t i = 0; i < coords.size(); i++) {
            ASSERT_EQ(vertices[i]->coords, coords[i]);
        }
    }

    // Test case 6: Creating a vector of half-edges from a container of coordinate arrays
    {
        std::vector<std::array<double, 3>> coords = {
            {0.0, 0.0, 0.0},
            {1.0, 1.0, 1.0},
            {2.0, 2.0, 2.0}};
        auto edges = make_shared_h_edges<double, 3>(coords);

        // Test if the correct number of half-edges is created
        ASSERT_EQ(edges.size(), coords.size());

        // Test if each half-edge's vertex has the correct coordinates
        for (size_t i = 0; i < coords.size(); i++) {
            ASSERT_EQ(edges[i]->vertex->coords, coords[i]);
        }
    }

    // Test case 7: Creating a half-edge loop from a vector of coordinates
    {
        std::vector<std::array<double, 3>> coords = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0},
            {0.0, 1.0, 0.0}
        };

        auto h_edges = make_HalfEdges(coords);
        size_t n = coords.size();

        // Check if the number of half-edges is equal to the number of coordinates
        ASSERT_EQ(h_edges.size(), n);

        // Check if each half-edge has the correct vertex coordinates and connections
        for (size_t i = 0; i < n; ++i) {
            ASSERT_EQ(h_edges[i]->vertex->coords, coords[i]);
            ASSERT_EQ(h_edges[i]->previous, h_edges[(i + n - 1) % n]);
            ASSERT_EQ(h_edges[i]->next, h_edges[(i + 1) % n]);
        }
    }
}

TEST(HalfEdgeMeshDataFaceTest, MakeLoopAndSharedHFace) {
    // Test case 1: Creating a half-edge loop and a face from a range of half-edges
    {
        std::vector<std::shared_ptr<Edge<double, 3>>> edges = {
            make_shared_h_edge<double, 3>(std::array<double, 3>{0.0, 0.0, 0.0}),
            make_shared_h_edge<double, 3>(std::array<double, 3>{1.0, 0.0, 0.0}),
            make_shared_h_edge<double, 3>(std::array<double, 3>{0.0, 1.0, 0.0})};

        auto face = make_shared_h_face<double, 3>(edges.begin(), edges.end());

        // Test if the face has the correct half-edge assigned
        ASSERT_EQ(face->edge, edges[0]);

        // Test if the half-edges have the correct face and next/previous half-edge relationships
        for (size_t i = 0; i < edges.size(); i++) {
            ASSERT_EQ(edges[i]->face, face);
            ASSERT_EQ(edges[i]->next, edges[(i + 1) % edges.size()]);
            ASSERT_EQ(edges[i]->previous, edges[(i + edges.size() - 1) % edges.size()]);
        }
    }

    // Test case 2: Creating a half-edge loop and a face from a container of half-edges
    {
        std::vector<std::shared_ptr<Edge<double, 3>>> edges = {
            make_shared_h_edge<double, 3>(std::array<double, 3>{0.0, 0.0, 0.0}),
            make_shared_h_edge<double, 3>(std::array<double, 3>{1.0, 0.0, 0.0}),
            make_shared_h_edge<double, 3>(std::array<double, 3>{0.0, 1.0, 0.0})};

        auto face = make_shared_h_face<double, 3>(edges);

        // Test if the face has the correct half-edge assigned
        ASSERT_EQ(face->edge, edges[0]);

        // Test if the half-edges have the correct face and next/previous half-edge relationships
        for (size_t i = 0; i < edges.size(); i++) {
            ASSERT_EQ(edges[i]->face, face);
            ASSERT_EQ(edges[i]->next, edges[(i + 1) % edges.size()]);
            ASSERT_EQ(edges[i]->previous, edges[(i + edges.size() - 1) % edges.size()]);
        }
    }


}

struct HalfEdgeFaceEdgeIteratorTest : public ::testing::Test {
    using Coord3D = std::array<double, 3>;
    std::vector<Coord3D> coords = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0}
    };
    std::vector<std::shared_ptr<HalfEdge<double, 3>>> h_edges;

    HalfEdgeFaceEdgeIteratorTest() {
        h_edges = make_HalfEdges(coords);
    }
};

TEST_F(HalfEdgeFaceEdgeIteratorTest, Circulation) {
    auto face = make_shared_h_face<double, 3>(h_edges);
    std::vector<Coord3D> collected_coords;
    for (const auto &edge : *face) {
        collected_coords.push_back(edge->vertex->coords);
    }

    ASSERT_EQ(coords, collected_coords);
}

TEST_F(HalfEdgeFaceEdgeIteratorTest, EmptyFace) {
    HalfEdgeFace<double, 3> empty_face;
    auto begin_it = begin(empty_face);
    auto end_it = end(empty_face);
    ASSERT_EQ(begin_it, end_it);
}

TEST_F(HalfEdgeFaceEdgeIteratorTest, InequalityOperator)
{
    // Create a face with 3 edges
    std::vector<std::array<double, 3>> coords = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}};
    auto h_edges = make_HalfEdges(coords);
    auto face = make_shared_h_face<double, 3>(h_edges);

    // Iterate through the face's half-edges and test the inequality operator
    auto begin_it = begin(*face);
    auto end_it = end(*face);

    ASSERT_NE(begin_it, end_it);

    ++begin_it;
    ASSERT_NE(begin_it, end_it);

    ++begin_it;
    ASSERT_NE(begin_it, end_it);

    ++begin_it;
    ASSERT_EQ(begin_it, end_it);
}

class CyclicHalfEdgeFaceEdgeRangeTest : public ::testing::Test {
protected:
    // CyclicHalfEdgeFaceEdgeRangeTest() : range_(HalfEdgeFaceEdgeIterator<double, 3>{}) {}
    void SetUp() override {
        // Define a square face in 2D
        auto coords = std::vector<std::array<double, 2>>{{{-1, -1}}, {{-1, 1}}, {{1, 1}}, {{1, -1}}};
        auto h_edges = make_HalfEdges(coords);
        auto face = make_shared_h_face<double, 2>(h_edges);
        range_ = CyclicHalfEdgeFaceEdgeRange<double, 2>(*face);
    }

    void TearDown() override {}

    CyclicHalfEdgeFaceEdgeRange<double, 2> range_;
};

// TEST_F(CyclicHalfEdgeFaceEdgeRangeTest, IterateOverEdges) {
//     std::vector<std::shared_ptr<HalfEdge<double, 2>>> expected_edges{range_.begin(), range_.end()};
//     std::vector<std::shared_ptr<HalfEdge<double, 2>>> actual_edges{
//         range_.begin(), std::prev(range_.end())};

//     ASSERT_EQ(expected_edges.size(), actual_edges.size());

//     for (size_t i = 0; i < expected_edges.size(); ++i) {
//         EXPECT_EQ(expected_edges[i], actual_edges[i]);
//     }
// }

// TEST_F(CyclicHalfEdgeFaceEdgeRangeTest, IterateOverEdgesReverse) {
//     std::vector<std::shared_ptr<HalfEdge<double, 2>>> expected_edges{range_.rbegin(), range_.rend()};
//     std::vector<std::shared_ptr<HalfEdge<double, 2>>> actual_edges{
//         range_.rbegin(), std::prev(range_.rend())};

//     ASSERT_EQ(expected_edges.size(), actual_edges.size());

//     for (size_t i = 0; i < expected_edges.size(); ++i) {
//         EXPECT_EQ(expected_edges[i], actual_edges[i]);
//     }
// }