#pragma once
#include <memory>
#include <algorithm>
#include <vector>
#include <ranges>
#include <iterator>


namespace gbs
{
    template <std::floating_point T, size_t dim>
    struct HalfEdge;
    template <std::floating_point T, size_t dim>
    struct HalfEdgeVertex;
    template <std::floating_point T, size_t dim>
    struct HalfEdgeFace;
/**
 * @brief A struct representing a vertex in a half-edge data structure
 *
 * @tparam T Numeric type of the vertex coordinates (must be a floating-point type)
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 */
    template <std::floating_point T, size_t dim>
    struct HalfEdgeVertex
    {
        std::array<T, dim> coords; ///< Coordinates of the vertex
        std::shared_ptr<HalfEdge<T, dim>> edge; ///< An outgoing half-edge from the vertex
    };

    /**
     * @brief Struct representing a face in the half-edge data structure
     *
     * @tparam T Numeric type of the vertex coordinates
     * @tparam dim Dimension of the vertex coordinates (2 or 3)
     */
    template <std::floating_point T, size_t dim>
    struct HalfEdgeFace
    {
        std::shared_ptr<HalfEdge<T, dim>> edge; ///< A half-edge on the boundary of the face
    };

    /**
     * @brief Struct representing a half-edge in the half-edge data structure
     *
     * @tparam T Numeric type of the vertex coordinates
     * @tparam dim Dimension of the vertex coordinates (2 or 3)
     */
    template <std::floating_point T, size_t dim>
    struct HalfEdge
    {
        std::shared_ptr<HalfEdgeVertex<T, dim>> vertex; ///< Starting vertex of the half-edge
        std::shared_ptr<HalfEdgeFace<T, dim>> face; ///< Face to the left of the half-edge
        std::shared_ptr<HalfEdge<T, dim>> next; ///< Next half-edge along the boundary of the face
        std::shared_ptr<HalfEdge<T, dim>> previous; ///< Previous half-edge along the boundary of the face
        std::shared_ptr<HalfEdge<T, dim>> opposite; ///< Opposite half-edge (pointing in the opposite direction)

        // HalfEdge(const std::array<T, dim> &coords) : vertex{std::make_shared<HalfEdgeVertex<T, dim>>({coords, this})} {}
        ///< Constructor for creating a HalfEdge object from vertex coordinates
        ///< @param coords The coordinates of the starting vertex of the half-edge
    };

/**
 * @brief Creates a shared pointer to a HalfEdge object with a given vertex and optional face
 *
 * @tparam T Numeric type of the vertex coordinates
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 * @param vertex The vertex at which the half-edge starts
 * @param face (optional) The face to the left of the half-edge
 * @return std::shared_ptr<HalfEdge<T, dim>> A shared pointer to the created HalfEdge object
 */
    template <std::floating_point T, size_t dim>
    auto make_shared_h_edge(const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex, const std::shared_ptr<HalfEdgeFace<T, dim>> &face = nullptr) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto hedge = std::make_shared<HalfEdge<T, dim>>(
            HalfEdge<T, dim>{
                .vertex = vertex,
                .face = face});

        if (!vertex->edge)
        {
            vertex->edge = hedge;
        }

        return hedge;
    }
/**
 * @brief Creates a shared pointer to a HalfEdgeVertex object with a given set of coordinates
 *
 * @tparam T Numeric type of the vertex coordinates
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 * @param coords The coordinates of the vertex
 * @return std::shared_ptr<HalfEdgeVertex<T, dim>> A shared pointer to the created HalfEdgeVertex object
 */
    template <std::floating_point T, size_t dim>
    auto make_shared_h_vertex(const std::array<T, dim> &coords)
    {
        return std::make_shared<HalfEdgeVertex<T, dim>>(HalfEdgeVertex<T, dim>{coords, nullptr});
    }

/**
 * @brief Creates a shared pointer to a HalfEdge object with a given set of coordinates for its vertex
 *
 * @tparam T Numeric type of the vertex coordinates
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 * @param coords The coordinates of the vertex
 * @return std::shared_ptr<HalfEdge<T, dim>> A shared pointer to the created HalfEdge object
 */
    template <std::floating_point T, size_t dim>
    auto make_shared_h_edge(const std::array<T, dim> &coords) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto vertex = std::make_shared<HalfEdgeVertex<T, dim>>(HalfEdgeVertex<T, dim>{coords, nullptr});
        return make_shared_h_edge<T, dim>(vertex);
    }
/**
 * @brief Creates a vector of shared pointers to HalfEdgeVertex objects from a container of coordinate arrays
 *
 * @tparam T Numeric type of the vertex coordinates (must be a floating-point type)
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 * @tparam _Container Type of the container holding the coordinate arrays
 * @param coords The container holding the coordinate arrays
 * @return std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> A vector of shared pointers to the created HalfEdgeVertex objects
 */
    template <std::floating_point T, size_t dim, typename _Container>
    auto make_shared_h_vertices(const _Container &coords)
    {
        auto n = coords.size();
        std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> vertices_cloud(n);
        std::transform(
            coords.begin(),
            coords.end(),
            vertices_cloud.begin(),
            [](const auto &xyz)
            { return make_shared_h_vertex<T, dim>(xyz); });
        return vertices_cloud;
    }
/**
 * @brief Creates a vector of shared pointers to HalfEdge objects from a container of coordinate arrays
 *
 * @tparam T Numeric type of the vertex coordinates (must be a floating-point type)
 * @tparam dim Dimension of the vertex coordinates (2 or 3)
 * @tparam _Container Type of the container holding the coordinate arrays
 * @param coords The container holding the coordinate arrays
 * @return std::vector<std::shared_ptr<HalfEdge<T, dim>>> A vector of shared pointers to the created HalfEdge objects
 */
    template <std::floating_point T, size_t dim, typename _Container>
    auto make_shared_h_edges(const _Container &coords)
    {
        auto n = coords.size();
        std::vector<std::shared_ptr<HalfEdge<T, dim>>> vertices_cloud(n);
        std::transform(
            coords.begin(),
            coords.end(),
            vertices_cloud.begin(),
            [](const auto &xyz)
            { return make_shared_h_edge<T, dim>(xyz); });
        return vertices_cloud;
    }
/**
 * @brief Constructs a loop of half-edges for a given face.
 *
 * This function takes a range of half-edge iterators and a face shared pointer.
 * It sets the face and next/previous half-edge relationships for each half-edge in the loop.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @tparam It The iterator type for the half-edge container.
 * @param begin An iterator pointing to the first half-edge in the range.
 * @param end An iterator pointing one past the last half-edge in the range.
 * @param p_face A shared pointer to the face for which the loop is constructed.
 */
    template <std::floating_point T, size_t dim, std::input_or_output_iterator It>
    void make_loop(const It &begin, const It &end, std::shared_ptr<HalfEdgeFace<T, dim>> &p_face)
    {
        auto n = std::ranges::distance(begin, end);
        assert(p_face);
        assert(n > 1);

        p_face->edge = (*begin);

        for (auto it = begin; it != end; ++it)
        {
            (*it)->face = p_face;
            (*it)->next = std::next(it) != end ? *std::next(it) : *begin;
            (*it)->previous = it != begin ? *std::prev(it) : *std::next(begin, n - 1);
        }
    }

/**
 * @brief Creates a shared pointer to a HalfEdgeFace and constructs a loop of half-edges for it.
 *
 * This function takes a range of half-edge iterators, creates a new HalfEdgeFace instance,
 * and sets up the half-edge loop for the face.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @tparam It The iterator type for the half-edge container.
 * @param begin An iterator pointing to the first half-edge in the range.
 * @param end An iterator pointing one past the last half-edge in the range.
 * @return A shared pointer to the created HalfEdgeFace with the half-edge loop constructed,
 *         or nullptr if there are less than 2 half-edges in the range.
 */
    template <std::floating_point T, size_t dim, std::input_or_output_iterator It>
    auto make_shared_h_face(const It &begin, const It &end) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {

        auto n = std::distance(begin, end);
        if (n < 2)
        {
            return nullptr;
        }

        auto p_face = std::make_shared<HalfEdgeFace<T, dim>>();
        make_loop(begin, end, p_face);

        return p_face;
    }

/**
 * @brief Creates a shared pointer to a HalfEdgeFace and constructs a loop of half-edges for it.
 *
 * This function takes a container of half-edges, creates a new HalfEdgeFace instance,
 * and sets up the half-edge loop for the face using the range defined by the container.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param he_lst A container of half-edges.
 * @return A shared pointer to the created HalfEdgeFace with the half-edge loop constructed,
 *         or nullptr if there are less than 2 half-edges in the container.
 */

    template <std::floating_point T, size_t dim>
    auto make_shared_h_face(const auto &he_lst) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        return make_shared_h_face<T, dim>(std::begin(he_lst), std::end(he_lst));
    }
/**
 * @brief Creates a new half-edge opposite to the given half-edge, and sets up the opposite relationship between them.
 *
 * This function takes a shared pointer to a HalfEdgeVertex and a shared pointer to a HalfEdge,
 * creates a new HalfEdge instance with the given vertex, and sets up the opposite relationship between the given edge and the new edge.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param vertex A shared pointer to the HalfEdgeVertex for the new opposite half-edge.
 * @param edge A shared pointer to the HalfEdge for which the opposite half-edge will be created.
 * @return A shared pointer to the created HalfEdge that is opposite to the given edge.
 */
    template <std::floating_point T, size_t dim>
    auto make_opposite(const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex, const std::shared_ptr<HalfEdge<T, dim>> &edge)
    {
        auto opposite = make_shared_h_edge(edge->previous->vertex);
        opposite->opposite = edge;
        edge->opposite = opposite;
        return opposite;
    }
/**
 * @brief Creates a vector of shared pointers to HalfEdge instances based on a vector of coordinates.
 *
 * This function takes a vector of coordinates and creates a loop of half-edges by connecting
 * consecutive coordinates. The resulting half-edges are returned as a vector of shared pointers.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param coords A vector of coordinates, represented as std::array<T, dim>.
 * @return A vector of shared pointers to HalfEdge instances, where each half-edge corresponds to a coordinate in the input vector.
 */
    template <std::floating_point T, size_t dim>
    auto make_HalfEdges(const std::vector<std::array<T, dim>> &coords)
    {
        size_t n = coords.size();

        // Preallocate memory for the half-edges vector
        std::vector<std::shared_ptr<HalfEdge<T, dim>>> h_edges;
        h_edges.reserve(n);

        // Create half-edges from coordinates
        for (const auto &coord : coords)
        {
            h_edges.push_back(std::make_shared<HalfEdge<T, dim>>(HalfEdge<T, dim>{.vertex = make_shared_h_vertex(coord)}));
        }

        // Connect half-edges in a loop
        for (size_t i = 0; i < n; ++i)
        {
            h_edges[i]->previous = h_edges[(i + n - 1) % n];
            h_edges[i]->next = h_edges[(i + 1) % n];
        }

        return h_edges;
    }
/**
 * @brief Custom iterator for circulating through half-edges of a face.
 *
 * @tparam T Floating-point type
 * @tparam dim Dimension of the half-edge data structure
 */
    template <std::floating_point T, size_t dim>
    class HalfEdgeFaceEdgeIterator
    {
    public:
        /// Iterator traits
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::shared_ptr<HalfEdge<T, dim>>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type *;
        using reference = value_type &;

        /**
         * @brief Constructor for creating a HalfEdgeFaceEdgeIterator object.
         *
         * @param start_edge Starting edge for the iterator
         * @param is_end Flag indicating if the iterator is at the end
         */
        explicit HalfEdgeFaceEdgeIterator(const value_type &start_edge, bool is_end = false)
            : current_edge(start_edge), initial_edge(start_edge), is_end(is_end) {}

        /// Dereference operator
        reference operator*() { return current_edge; }

        /// Arrow operator
        pointer operator->() { return &current_edge; }

        /// Prefix increment operator
        HalfEdgeFaceEdgeIterator &operator++()
        {
            if (current_edge && !is_end)
            {
                current_edge = current_edge->next;
                if (current_edge == initial_edge)
                {
                    is_end = true;
                }
            }
            return *this;
        }

        /// Postfix increment operator
        HalfEdgeFaceEdgeIterator operator++(int)
        {
            HalfEdgeFaceEdgeIterator tmp(*this);
            operator++();
            return tmp;
        }

        /// Equality operator
        bool operator==(const HalfEdgeFaceEdgeIterator &other) const
        {
            return (current_edge == other.current_edge) && (is_end == other.is_end);
        }

        /// Inequality operator
        bool operator!=(const HalfEdgeFaceEdgeIterator &other) const
        {
            return !(*this == other);
        }

    private:
        value_type current_edge;  ///< Current half-edge
        value_type initial_edge;  ///< Initial half-edge
        bool is_end;              ///< End flag
    };

    /**
     * @brief Begin function for HalfEdgeFaceEdgeIterator.
     *
     * @tparam T Floating-point type
     * @tparam dim Dimension of the half-edge data structure
     * @param face A constant reference to the face
     * @return HalfEdgeFaceEdgeIterator pointing to the beginning of the face
     */
    template <std::floating_point T, size_t dim>
    HalfEdgeFaceEdgeIterator<T, dim> begin(const HalfEdgeFace<T, dim> &face)
    {
        return HalfEdgeFaceEdgeIterator<T, dim>(face.edge);
    }

    /**
     * @brief End function for HalfEdgeFaceEdgeIterator.
     *
     * @tparam T Floating-point type
     * @tparam dim Dimension of the half-edge data structure
     * @param face A constant reference to the face
     * @return HalfEdgeFaceEdgeIterator pointing to the end of the face
     */
    template <std::floating_point T, size_t dim>
    HalfEdgeFaceEdgeIterator<T, dim> end(const HalfEdgeFace<T, dim> &face)
    {
        return HalfEdgeFaceEdgeIterator<T, dim>(face.edge, true);
    }

    /**
     * @brief A custom range representing a cycle of half-edges of a face.
     *
     * @tparam T Floating-point type
     * @tparam dim Dimension of the half-edge data structure
     */
    template <std::floating_point T, size_t dim>
    class CyclicHalfEdgeFaceEdgeRange
    {
    public:
        using iterator = HalfEdgeFaceEdgeIterator<T, dim>;
        using value_type = typename iterator::value_type;

        CyclicHalfEdgeFaceEdgeRange(const HalfEdgeFace<T, dim> &h_face) : start_edge_(h_face.edge) {}
        /**
         * @brief Constructor for creating a CyclicHalfEdgeFaceEdgeRange object.
         *
         * @param start_edge Starting edge for the range
         */
        CyclicHalfEdgeFaceEdgeRange(const value_type &start_edge) : start_edge_(start_edge) {}

        /// Returns an iterator to the beginning of the range
        iterator begin() const { return iterator{start_edge_}; }

        /// Returns an iterator to the end of the range
        iterator end() const { return iterator{start_edge_, true}; }

    private:
        value_type start_edge_;  ///< Starting half-edge
    };

    /**
     * @brief A custom view representing a cycle of half-edges of a face.
     *
     * @tparam T Floating-point type
     * @tparam dim Dimension of the half-edge data structure
     */
    template <std::floating_point T, size_t dim>
    class CyclicHalfEdgeFaceEdgeView : public std::ranges::view_interface<CyclicHalfEdgeFaceEdgeView<T, dim>>
    {
    public:
        using iterator = HalfEdgeFaceEdgeIterator<T, dim>;
        using value_type = typename iterator::value_type;

        CyclicHalfEdgeFaceEdgeView(const HalfEdgeFace<T, dim> &h_face) : start_edge_(h_face.edge) {}
        /**
         * @brief Constructor for creating a CyclicHalfEdgeFaceEdgeView object.
         *
         * @param start_edge Starting edge for the view
         */
        CyclicHalfEdgeFaceEdgeView(const value_type &start_edge) : start_edge_(start_edge) {}

        /// Returns an iterator to the beginning of the view
        iterator begin() const { return iterator{start_edge_}; }

        /// Returns an iterator to the end of the view
        iterator end() const { return iterator{start_edge_, true}; }

    private:
        value_type start_edge_;  ///< Starting half-edge
    };
}