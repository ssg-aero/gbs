#pragma once
/**
 * @file cage.h
 * @brief Control cage of a Catmull-Clark subdivision surface and its checked topology.
 *
 * A cage is a polygonal mesh: vertices, faces given as lists of vertex indices (any number of
 * sides from 3, counterclockwise seen from outside), infinitely sharp edges (creases) and
 * vertices tagged as corners. `Topology` builds its half-edges once and refuses, with a
 * `CageError` naming the culprit, what the subdivision cannot take: an edge shared by three
 * faces or more, a vertex whose faces make several fans, neighbouring faces of opposite
 * orientations, a face with a repeated vertex, an isolated vertex, an interior vertex of two
 * edges. Each vertex gets the rule of OpenSubdiv's Catmark scheme (boundary interpolation
 * "edge and corner"): corner when tagged, on three sharp edges or more, or on the boundary of a
 * single face; crease on two sharp edges (a boundary edge is sharp); dart on one; smooth
 * otherwise. Design: docs/sources/design/subd_catmull_clark.md, subd_pr01_cage.md.
 */
#include <gbs/gbslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace gbs::subd
{
    /// No half-edge, edge or face: the twin of a boundary half-edge, past the end of a fan.
    inline constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

    /// Control cage of a Catmull-Clark surface.
    template <std::floating_point T>
    struct Cage
    {
        points_vector<T, 3> vertices;
        std::vector<std::vector<std::size_t>> faces;     ///< counterclockwise seen from outside, 3 vertices or more
        std::vector<std::array<std::size_t, 2>> creases; ///< infinitely sharp edges, either way
        std::vector<std::size_t> corners;                ///< vertices interpolated by the limit surface
    };

    /// What makes a cage unfit for the subdivision; `CageError::index` names the element at fault.
    enum class CageIssue : std::uint8_t
    {
        NoFace,                  ///< no face at all (index: none)
        FaceTooSmall,            ///< a face of fewer than 3 vertices (face)
        VertexOutOfRange,        ///< a face refers to a vertex that does not exist (face)
        RepeatedVertex,          ///< a face goes through a vertex twice (face)
        NonFiniteVertex,         ///< a coordinate is infinite or NaN (vertex)
        NonManifoldEdge,         ///< an edge shared by three faces or more (face using it last)
        InconsistentOrientation, ///< two faces run along their common edge the same way (second face)
        NonManifoldVertex,       ///< the faces of a vertex make several fans (vertex)
        IsolatedVertex,          ///< a vertex of no face (vertex)
        InteriorValenceTwo,      ///< an interior vertex of two edges: two faces folded on each other (vertex)
        CreaseNotAnEdge,         ///< a crease joins two vertices that share no edge (crease)
        CornerOutOfRange,        ///< a corner refers to a vertex that does not exist (corner)
    };

    /// Refusal of a cage (`ValueError` in Python).
    class CageError : public std::invalid_argument
    {
    public:
        CageError(CageIssue issue, std::size_t index, const std::string &what)
            : std::invalid_argument{what}, issue_{issue}, index_{index} {}

        auto issue() const noexcept -> CageIssue { return issue_; }
        auto index() const noexcept -> std::size_t { return index_; } ///< face, vertex, crease or corner at fault (see CageIssue)

    private:
        CageIssue issue_;
        std::size_t index_;
    };

    /// Subdivision rule of a vertex, from its sharp edges (OpenSubdiv's Sdc::Crease::Rule).
    enum class VertexRule : std::uint8_t
    {
        Smooth, ///< no sharp edge
        Dart,   ///< one sharp edge: smooth vertex rule, the edge is still split at its middle
        Crease, ///< two sharp edges: (P₋ + 6 V + P₊) / 8, the limit follows the cubic B-spline of the crease
        Corner, ///< tagged, three sharp edges or more, or on the boundary of a single face: fixed
    };

    /**
     * @brief Half-edges of a cage, its edges, and the rule of each vertex.
     *
     * Half-edge h of face f, at position i, runs from faces[f][i] to faces[f][i + 1]; the
     * half-edges of a face are contiguous, in its order. Edges are numbered as they are first met
     * along the faces. Around a vertex, `outgoing` and `ccw` walk its faces counterclockwise
     * (seen from outside), from the boundary edge on a boundary.
     *
     * A strict topology checks a cage (`CageError`). A lenient one takes a piece of a cage, the
     * faces around a face or a vertex: a vertex of its border may then be pinched or look like
     * a boundary, its rule is meaningless; only the vertices whose faces are all in the piece
     * are used.
     */
    class Topology
    {
    public:
        Topology(std::size_t n_vertices, const std::vector<std::vector<std::size_t>> &faces,
                 std::span<const std::array<std::size_t, 2>> creases, std::span<const std::size_t> corners, bool strict = true)
            : n_vertices_{n_vertices}
        {
            if (faces.empty())
                throw CageError{CageIssue::NoFace, none, "subdivision cage: no face"};
            build_edges(build_half_edges(faces), creases);
            build_vertices(corners, strict);
        }

        template <std::floating_point T>
        explicit Topology(const Cage<T> &cage, bool strict = true)
            : Topology{cage.vertices.size(), cage.faces, cage.creases, cage.corners, strict}
        {
            if (strict)
                for (std::size_t v{}; v < cage.vertices.size(); ++v)
                    if (!std::ranges::all_of(cage.vertices[v], [](T x) { return std::isfinite(x); }))
                        throw CageError{CageIssue::NonFiniteVertex, v, "subdivision cage: vertex " + std::to_string(v) + " is not finite"};
        }

        auto n_vertices() const noexcept -> std::size_t { return n_vertices_; }
        auto n_faces() const noexcept -> std::size_t { return face_begin_.size() - 1; }
        auto n_edges() const noexcept -> std::size_t { return edge_vertices_.size(); }
        auto n_half_edges() const noexcept -> std::size_t { return origin_.size(); }

        auto face_size(std::size_t f) const -> std::size_t { return face_begin_[f + 1] - face_begin_[f]; }
        auto face_begin(std::size_t f) const -> std::size_t { return face_begin_[f]; } ///< half-edge from the first vertex of f

        auto origin(std::size_t h) const -> std::size_t { return origin_[h]; }
        auto target(std::size_t h) const -> std::size_t { return origin_[next(h)]; }
        auto face(std::size_t h) const -> std::size_t { return face_[h]; }
        auto next(std::size_t h) const -> std::size_t
        {
            const auto f = face_[h];
            return h + 1 == face_begin_[f + 1] ? face_begin_[f] : h + 1;
        }
        auto prev(std::size_t h) const -> std::size_t
        {
            const auto f = face_[h];
            return h == face_begin_[f] ? face_begin_[f + 1] - 1 : h - 1;
        }
        auto twin(std::size_t h) const -> std::size_t { return twin_[h]; } ///< none on a boundary
        auto edge(std::size_t h) const -> std::size_t { return edge_[h]; }

        auto edge_vertices(std::size_t e) const -> const std::array<std::size_t, 2> & { return edge_vertices_[e]; }
        auto edge_half_edge(std::size_t e) const -> std::size_t { return edge_half_edge_[e]; } ///< runs from edge_vertices(e)[0]
        auto is_boundary_edge(std::size_t e) const -> bool { return twin_[edge_half_edge_[e]] == none; }
        auto is_crease(std::size_t e) const -> bool { return crease_[e]; }
        auto is_sharp(std::size_t e) const -> bool { return crease_[e] || is_boundary_edge(e); } ///< crease or boundary

        /// First half-edge from v, counterclockwise: on a boundary, the one along the boundary.
        auto outgoing(std::size_t v) const -> std::size_t { return outgoing_[v]; }
        /// Next half-edge from origin(h), counterclockwise; none past a boundary.
        auto ccw(std::size_t h) const -> std::size_t { return twin_[prev(h)]; }
        /// Previous half-edge from origin(h), clockwise; none past a boundary.
        auto cw(std::size_t h) const -> std::size_t { return twin_[h] == none ? none : next(twin_[h]); }

        auto vertex_edges(std::size_t v) const -> std::span<const std::size_t>
        {
            return {vertex_edges_.data() + vertex_edges_begin_[v], vertex_edges_begin_[v + 1] - vertex_edges_begin_[v]};
        }
        auto vertex_faces(std::size_t v) const -> std::span<const std::size_t>
        {
            return {vertex_faces_.data() + vertex_faces_begin_[v], vertex_faces_begin_[v + 1] - vertex_faces_begin_[v]};
        }
        auto valence(std::size_t v) const -> std::size_t { return vertex_edges(v).size(); }
        auto face_count(std::size_t v) const -> std::size_t { return vertex_faces(v).size(); }
        auto is_boundary(std::size_t v) const -> bool { return twin_[outgoing_[v]] == none; }
        auto is_corner_tagged(std::size_t v) const -> bool { return corner_tag_[v]; }
        auto rule(std::size_t v) const -> VertexRule { return rule_[v]; }

        /// Half-edges from v counterclockwise, from outgoing(v): its faces in order.
        auto fan(std::size_t v) const -> std::vector<std::size_t>
        {
            std::vector<std::size_t> hs;
            const auto h0 = outgoing_[v];
            for (auto h = h0; h != none && hs.size() <= n_half_edges(); h = ccw(h))
            {
                hs.push_back(h);
                if (ccw(h) == h0)
                    break;
            }
            return hs;
        }

        /**
         * @brief The faces around origin(h0) whose limit is one smooth piece with the face of h0,
         * counterclockwise: the whole fan for a smooth vertex or a dart, the faces between two
         * sharp edges for a crease or a corner. Second: true when they go all around.
         */
        auto sector(std::size_t h0) const -> std::pair<std::vector<std::size_t>, bool>
        {
            const auto v = origin_[h0];
            const bool split = rule_[v] == VertexRule::Crease || rule_[v] == VertexRule::Corner;
            // back clockwise to the first face of the sector; all around: from h0
            auto first = h0;
            for (std::size_t n{}; n < n_half_edges(); ++n)
            {
                if (split && crease_[edge_[first]])
                    break;
                const auto h = cw(first);
                if (h == none)
                    break;
                if (h == h0)
                {
                    first = h0;
                    break;
                }
                first = h;
            }
            std::vector<std::size_t> hs;
            for (auto h = first; hs.size() <= n_half_edges();)
            {
                hs.push_back(h);
                const auto p = prev(h);
                if (split && crease_[edge_[p]])
                    return {hs, false};
                h = twin_[p];
                if (h == none)
                    return {hs, false};
                if (h == first)
                    return {hs, true};
            }
            return {hs, false};
        }

    private:
        std::size_t n_vertices_;
        std::vector<std::size_t> face_begin_{0}, origin_, face_, twin_, edge_;
        std::vector<std::array<std::size_t, 2>> edge_vertices_;
        std::vector<std::size_t> edge_half_edge_;
        std::vector<bool> crease_;
        std::vector<std::size_t> outgoing_, vertex_edges_begin_, vertex_edges_, vertex_faces_begin_, vertex_faces_;
        std::vector<bool> corner_tag_;
        std::vector<VertexRule> rule_;

        /// Origin, target and index of a half-edge; sorted, a pair of vertices is found by binary search.
        using Directed = std::array<std::size_t, 3>;

        static auto find(const std::vector<Directed> &directed, std::size_t a, std::size_t b) -> std::size_t
        {
            const auto it = std::ranges::lower_bound(directed, Directed{a, b, 0});
            return it != directed.end() && (*it)[0] == a && (*it)[1] == b ? (*it)[2] : none;
        }

        auto build_half_edges(const std::vector<std::vector<std::size_t>> &faces) -> std::vector<Directed>
        {
            for (std::size_t f{}; f < faces.size(); ++f)
            {
                const auto &fv = faces[f];
                if (fv.size() < 3)
                    throw CageError{CageIssue::FaceTooSmall, f, "subdivision cage: face " + std::to_string(f) + " has fewer than 3 vertices"};
                for (std::size_t i{}; i < fv.size(); ++i)
                {
                    if (fv[i] >= n_vertices_)
                        throw CageError{CageIssue::VertexOutOfRange, f, "subdivision cage: face " + std::to_string(f) + " refers to vertex " + std::to_string(fv[i]) + ", out of range"};
                    if (std::find(fv.begin(), fv.begin() + i, fv[i]) != fv.begin() + i)
                        throw CageError{CageIssue::RepeatedVertex, f, "subdivision cage: face " + std::to_string(f) + " goes through vertex " + std::to_string(fv[i]) + " twice"};
                    origin_.push_back(fv[i]);
                    face_.push_back(f);
                }
                face_begin_.push_back(origin_.size());
            }
            std::vector<Directed> directed(origin_.size());
            for (std::size_t h{}; h < origin_.size(); ++h)
                directed[h] = {origin_[h], target(h), h};
            std::ranges::sort(directed);
            // a pair of vertices run twice the same way: the first half-edge found doing so
            auto twice = none;
            for (std::size_t i = 1; i < directed.size(); ++i)
                if (directed[i][0] == directed[i - 1][0] && directed[i][1] == directed[i - 1][1])
                    twice = std::min(twice, directed[i][2]);
            if (twice != none)
            {
                const auto a = origin_[twice], b = target(twice);
                const bool third = find(directed, b, a) != none;
                throw CageError{third ? CageIssue::NonManifoldEdge : CageIssue::InconsistentOrientation, face_[twice],
                                "subdivision cage: edge (" + std::to_string(a) + ", " + std::to_string(b) + ") of face " + std::to_string(face_[twice]) +
                                    (third ? " is shared by three faces or more" : " runs the same way as in another face: orientations differ")};
            }
            // twins: the half-edge running the other way along the same edge
            twin_.resize(origin_.size());
            for (std::size_t h{}; h < origin_.size(); ++h)
                twin_[h] = find(directed, target(h), origin_[h]);
            return directed;
        }

        void build_edges(const std::vector<Directed> &directed, std::span<const std::array<std::size_t, 2>> creases)
        {
            edge_.assign(origin_.size(), none);
            for (std::size_t h{}; h < origin_.size(); ++h)
                if (edge_[h] == none)
                {
                    edge_[h] = edge_vertices_.size();
                    if (twin_[h] != none)
                        edge_[twin_[h]] = edge_vertices_.size();
                    edge_vertices_.push_back({origin_[h], target(h)});
                    edge_half_edge_.push_back(h);
                }
            crease_.assign(edge_vertices_.size(), false);
            for (std::size_t c{}; c < creases.size(); ++c)
            {
                const auto [a, b] = creases[c];
                auto h = find(directed, a, b);
                if (h == none)
                    h = find(directed, b, a);
                if (h == none)
                    throw CageError{CageIssue::CreaseNotAnEdge, c, "subdivision cage: crease " + std::to_string(c) + " (" + std::to_string(a) + ", " + std::to_string(b) + ") is not an edge"};
                crease_[edge_[h]] = true;
            }
        }

        void build_vertices(std::span<const std::size_t> corners, bool strict)
        {
            corner_tag_.assign(n_vertices_, false);
            for (std::size_t c{}; c < corners.size(); ++c)
            {
                if (corners[c] >= n_vertices_)
                    throw CageError{CageIssue::CornerOutOfRange, c, "subdivision cage: corner " + std::to_string(c) + " refers to vertex " + std::to_string(corners[c]) + ", out of range"};
                corner_tag_[corners[c]] = true;
            }
            // vertex -> edges and vertex -> faces, in increasing order
            auto csr = [&](std::vector<std::size_t> &begin, std::vector<std::size_t> &items, auto &&for_each_pair) {
                begin.assign(n_vertices_ + 1, 0);
                for_each_pair([&](std::size_t v, std::size_t) { ++begin[v + 1]; });
                for (std::size_t v{}; v < n_vertices_; ++v)
                    begin[v + 1] += begin[v];
                items.resize(begin[n_vertices_]);
                auto fill = begin;
                for_each_pair([&](std::size_t v, std::size_t x) { items[fill[v]++] = x; });
            };
            csr(vertex_edges_begin_, vertex_edges_, [&](auto &&add) {
                for (std::size_t e{}; e < edge_vertices_.size(); ++e)
                {
                    add(edge_vertices_[e][0], e);
                    add(edge_vertices_[e][1], e);
                }
            });
            csr(vertex_faces_begin_, vertex_faces_, [&](auto &&add) {
                for (std::size_t h{}; h < origin_.size(); ++h)
                    add(origin_[h], face_[h]);
            });

            outgoing_.assign(n_vertices_, none);
            for (std::size_t h{}; h < origin_.size(); ++h)
            {
                auto &o = outgoing_[origin_[h]];
                if (o == none || (twin_[h] == none && twin_[o] != none))
                    o = h;
            }
            rule_.assign(n_vertices_, VertexRule::Smooth);
            for (std::size_t v{}; v < n_vertices_; ++v)
            {
                if (outgoing_[v] == none)
                {
                    if (strict)
                        throw CageError{CageIssue::IsolatedVertex, v, "subdivision cage: vertex " + std::to_string(v) + " belongs to no face"};
                    continue;
                }
                if (strict && fan(v).size() != face_count(v))
                    throw CageError{CageIssue::NonManifoldVertex, v, "subdivision cage: the faces of vertex " + std::to_string(v) + " make several fans"};
                if (strict && !is_boundary(v) && valence(v) == 2)
                    throw CageError{CageIssue::InteriorValenceTwo, v, "subdivision cage: interior vertex " + std::to_string(v) + " has only two edges"};
                const auto sharp = std::ranges::count_if(vertex_edges(v), [&](std::size_t e) { return is_sharp(e); });
                if (corner_tag_[v] || sharp >= 3 || (is_boundary(v) && face_count(v) == 1))
                    rule_[v] = VertexRule::Corner;
                else if (sharp == 2)
                    rule_[v] = VertexRule::Crease;
                else if (sharp == 1)
                    rule_[v] = VertexRule::Dart;
            }
        }
    };
}
