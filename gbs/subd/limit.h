#pragma once
/**
 * @file limit.h
 * @brief Exact evaluation of the limit surface of a Catmull-Clark cage.
 *
 * The limit is a set of patches: one per quadrangle of the cage, n per face of n sides (one per
 * quadrangle of its first subdivision, like OpenSubdiv's ptex faces), each on [0, 1]². On a
 * quadrangle, (0, 0) is its first vertex, u runs towards the second, v towards the last; on
 * the quadrangle of corner i of an n-gon, (0, 0) is vertex i, u runs towards the middle of side
 * i, v towards the middle of side i - 1.
 *
 * - A **regular** patch is the uniform bicubic B-spline of its 16 neighbours: a quadrangle whose
 *   four corners are regular (an interior vertex of four quadrangles and no sharp edge; a crease
 *   or boundary vertex with two quadrangles on the side of the patch; a corner with the patch
 *   alone between its sharp edges). Beyond a sharp side, the missing points are the phantom
 *   points 2 B - I, reflected in both directions at a corner: the crease rule keeps that
 *   relation at every step.
 * - Elsewhere, the faces around the patch are subdivided and the evaluation goes down into the
 *   quarter that holds (u, v), until it is regular: about log2(1 / max(u, v)) steps from an
 *   extraordinary vertex, as in Stam (1998), who does the same with the eigenbasis.
 * - At an irregular vertex itself, its limit masks. Position: V at a corner, (4 V + P₋ + P₊) / 6
 *   on a crease, (n² V + 4 Σ eᵢ + Σ dᵢ) / (n (n + 5)) at a smooth vertex (Halstead et al. 1993,
 *   on quadrangles), and at a dart the left eigenvector of eigenvalue 1 of the subdivision
 *   matrix of its ring. Normal: the tangent masks of Halstead et al. at a smooth vertex; at a
 *   dart, a crease or a corner, the two dominant eigenvectors of the subdivision matrix of the
 *   sector of the patch, when its spectrum is clean (|λ₁| ≥ |λ₂| > |λ₃|, two independent
 *   eigenvectors); its sign is that of the normals of the patch next to the vertex. OpenSubdiv's
 *   masks are not exact there: the smooth mask at a dart, the masks of Biermann et al. across a
 *   crease of three faces or more on a side, which are not eigenvectors of the standard crease
 *   rules. Where λ₂ is tied to λ₃ (a crease of one or four faces on a side, a corner of two to
 *   four faces, a corner tagged without sharp edges), the limit has no tangent plane or reaches
 *   one like 1 / depth: `normal` throws std::domain_error there. So does a derivative at an
 *   irregular vertex.
 * - OpenSubdiv measures the derivatives on the patches of an n-gon in the parameters of the
 *   n-gon, where a patch spans half of them: twice ours per order. Ours are those of `value`.
 *
 * Design: docs/sources/design/subd_catmull_clark.md § 3, subd_pr01_cage.md.
 */
#include <gbs/subd/refine.h>
#include <gbs/exceptions.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gbs::subd
{
    namespace limit_detail
    {
        template <std::floating_point T>
        using Grid = std::array<std::array<point<T, 3>, 4>, 4>; ///< grid[i][j], i along u, j along v

        /// The faces `faces` of c as a cage of their own, faces[0] first, starting at its vertex `start`.
        template <std::floating_point T>
        auto piece(const Cage<T> &c, const Topology &t, std::span<const std::size_t> faces, std::size_t start) -> Cage<T>
        {
            Cage<T> p;
            std::vector<std::pair<std::size_t, std::size_t>> map; // cage vertex -> piece vertex
            auto local = [&](std::size_t v) {
                for (const auto &[from, to] : map)
                    if (from == v)
                        return to;
                map.emplace_back(v, p.vertices.size());
                p.vertices.push_back(c.vertices[v]);
                return p.vertices.size() - 1;
            };
            std::vector<std::size_t> creased;
            for (std::size_t i{}; i < faces.size(); ++i)
            {
                const auto f = faces[i], n = t.face_size(f), s = i == 0 ? start : 0;
                std::vector<std::size_t> fv(n);
                for (std::size_t j{}; j < n; ++j)
                    fv[j] = local(c.faces[f][(s + j) % n]);
                p.faces.push_back(std::move(fv));
                for (auto h = t.face_begin(f); h < t.face_begin(f) + n; ++h)
                    if (const auto e = t.edge(h); t.is_crease(e) && !std::ranges::contains(creased, e))
                        creased.push_back(e);
            }
            for (const auto e : creased)
            {
                const auto [a, b] = t.edge_vertices(e);
                p.creases.push_back({local(a), local(b)});
            }
            for (const auto &[from, to] : map)
                if (t.is_corner_tagged(from))
                    p.corners.push_back(to);
            return p;
        }

        /// Face f and the faces touching its vertices, f first from its vertex `start`.
        template <std::floating_point T>
        auto face_ring(const Cage<T> &c, const Topology &t, std::size_t f, std::size_t start) -> Cage<T>
        {
            std::vector<std::size_t> faces{f};
            for (const auto v : c.faces[f])
                for (const auto g : t.vertex_faces(v))
                    if (!std::ranges::contains(faces, g))
                        faces.push_back(g);
            return piece(c, t, faces, start);
        }

        /// The faces around vertex v, v first (vertex 0 of the piece).
        template <std::floating_point T>
        auto vertex_ring(const Cage<T> &c, const Topology &t, std::size_t v) -> Cage<T>
        {
            const auto h = t.outgoing(v), f = t.face(h);
            std::vector<std::size_t> faces{f};
            for (const auto g : t.vertex_faces(v))
                if (g != f)
                    faces.push_back(g);
            return piece(c, t, faces, h - t.face_begin(f));
        }

        /// Points of a sector of quadrangles around V: its neighbours along the edges, then the opposite corners.
        template <std::floating_point T>
        struct Ring
        {
            point<T, 3> V;
            bool closed;
            points_vector<T, 3> e{}, d{};
            std::vector<bool> sharp{}; ///< of each edge
        };

        template <std::floating_point T>
        auto ring(const Cage<T> &c, const Topology &t, const std::vector<std::size_t> &hs, bool closed) -> Ring<T>
        {
            Ring<T> r{c.vertices[t.origin(hs.front())], closed};
            for (const auto h : hs)
            {
                r.e.push_back(c.vertices[t.target(h)]);
                r.d.push_back(c.vertices[t.target(t.next(h))]);
                r.sharp.push_back(t.is_sharp(t.edge(h)));
            }
            if (!closed)
            {
                const auto p = t.prev(hs.back());
                r.e.push_back(c.vertices[t.origin(p)]);
                r.sharp.push_back(t.is_sharp(t.edge(p)));
            }
            return r;
        }

        /// True when the faces around v (all of them, or those of the sector of hs) are quadrangles.
        inline auto quads(const Topology &t, std::span<const std::size_t> faces) -> bool
        {
            return std::ranges::all_of(faces, [&](std::size_t f) { return t.face_size(f) == 4; });
        }

        /// Subdivision matrix of a ring, points [V, e₀…, d₀…]: row i gives the new point i from the old ones.
        template <std::floating_point T>
        auto ring_matrix(VertexRule rule, const Ring<T> &r) -> Eigen::MatrixX<T>
        {
            const auto k = r.d.size(), m = r.e.size(), N = 1 + m + k;
            const auto E = [](std::size_t i) { return 1 + i; };
            const auto D = [m](std::size_t i) { return 1 + m + i; };
            Eigen::MatrixX<T> S = Eigen::MatrixX<T>::Zero(N, N);
            const auto add_face = [&](std::size_t row, std::size_t i, T w) { // w times the centroid of quad i
                S(row, 0) += w / 4;
                S(row, E(i)) += w / 4;
                S(row, D(i)) += w / 4;
                S(row, E((i + 1) % m)) += w / 4;
            };
            switch (rule)
            {
            case VertexRule::Corner:
                S(0, 0) = 1;
                break;
            case VertexRule::Crease:
                S(0, 0) = T(6) / 8;
                S(0, E(0)) += T(1) / 8;
                S(0, E(m - 1)) += T(1) / 8;
                break;
            case VertexRule::Smooth:
            case VertexRule::Dart:
            {
                const auto n = T(k);
                S(0, 0) += (n - 2) / n;
                for (std::size_t i{}; i < k; ++i)
                {
                    S(0, E(i)) += 1 / (n * n);
                    add_face(0, i, 1 / (n * n));
                }
                break;
            }
            }
            for (std::size_t i{}; i < m; ++i)
            {
                if (r.sharp[i])
                {
                    S(E(i), 0) += T(0.5);
                    S(E(i), E(i)) += T(0.5);
                    continue;
                }
                S(E(i), 0) += T(0.25);
                S(E(i), E(i)) += T(0.25);
                add_face(E(i), (i + k - 1) % k, T(0.25));
                add_face(E(i), i, T(0.25));
            }
            for (std::size_t i{}; i < k; ++i)
                add_face(D(i), i, 1);
            return S;
        }

        template <std::floating_point T>
        auto ring_points(const Ring<T> &r) -> points_vector<T, 3>
        {
            points_vector<T, 3> P{r.V};
            P.insert(P.end(), r.e.begin(), r.e.end());
            P.insert(P.end(), r.d.begin(), r.d.end());
            return P;
        }

        /// Left eigenvector of eigenvalue 1, the weights of the limit position; std::runtime_error if the solver fails.
        template <std::floating_point T>
        auto limit_position(const Eigen::MatrixX<T> &S, const Ring<T> &r) -> point<T, 3>
        {
            Eigen::EigenSolver<Eigen::MatrixX<T>> es(S.transpose());
            if (es.info() != Eigen::Success)
                throw std::runtime_error{"subdivision: no eigenvalues for the limit mask of a dart"};
            const auto &w = es.eigenvalues();
            Eigen::Index best{};
            for (Eigen::Index i{1}; i < w.size(); ++i)
                if (std::abs(w[i] - T(1)) < std::abs(w[best] - T(1)))
                    best = i;
            const Eigen::VectorX<T> l = es.eigenvectors().col(best).real();
            if (std::abs(w[best] - T(1)) > T(1e-8) || std::abs(l.sum()) <= T(1e-12) * l.template lpNorm<1>())
                throw std::runtime_error{"subdivision: no eigenvalue 1 for the limit mask of a dart"};
            const auto P = ring_points(r);
            point<T, 3> p{};
            for (std::size_t j{}; j < P.size(); ++j)
                p += l[j] / l.sum() * P[j];
            return p;
        }

        /**
         * @brief Normal of the tangent plane of the two dominant eigenvectors of S (eigenvalue 1
         * aside), up to its sign; nullopt unless the spectrum is clean: |λ₁| ≥ |λ₂| > |λ₃|, two
         * independent eigenvectors for λ₁ and λ₂, whose directions span a plane.
         *
         * With λ₂ tied to λ₃ (a crease of one or four faces on a side, a corner of two to four
         * faces, a corner tagged without sharp edges), the limit has either no tangent plane or one
         * that the surface reaches like 1 / depth (a defective eigenvalue): not given. From five
         * faces in the sector of a crease or a corner, two inner modes dominate those of its sharp
         * edges: the plane they span, reached from inside the sector, is given.
         */
        template <std::floating_point T>
        auto eigen_normal(const Eigen::MatrixX<T> &S, const Ring<T> &r) -> std::optional<point<T, 3>>
        {
            Eigen::EigenSolver<Eigen::MatrixX<T>> es(S.transpose());
            if (es.info() != Eigen::Success)
                return std::nullopt;
            const auto &w = es.eigenvalues();
            const auto &X = es.eigenvectors();
            std::vector<Eigen::Index> order;
            for (Eigen::Index i{}; i < w.size(); ++i)
                if (std::abs(w[i] - T(1)) > T(1e-8))
                    order.push_back(i);
            std::ranges::sort(order, [&](auto a, auto b) { return std::abs(w[a]) > std::abs(w[b]); });
            if (order.size() < 2)
                return std::nullopt;
            const T m2 = std::abs(w[order[1]]), m3 = order.size() > 2 ? std::abs(w[order[2]]) : T(0);
            if (m2 - m3 <= T(1e-6) * m2)
                return std::nullopt;
            Eigen::VectorX<T> l1, l2;
            const auto x1 = X.col(order[0]), x2 = X.col(order[1]);
            if (std::abs(w[order[0]].imag()) > T(1e-9) * m2) // a complex pair: its real and imaginary parts
                l1 = x1.real(), l2 = x1.imag();
            else
            {
                if (std::abs(x1.dot(x2)) >= (1 - T(1e-6)) * x1.norm() * x2.norm())
                    return std::nullopt; // a defective eigenvalue
                l1 = x1.real(), l2 = x2.real();
            }
            const auto P = ring_points(r);
            point<T, 3> t1{}, t2{};
            for (std::size_t j{}; j < P.size(); ++j)
            {
                t1 += l1[j] * (P[j] - r.V);
                t2 += l2[j] * (P[j] - r.V);
            }
            const auto n = t1 ^ t2;
            if (norm(n) <= T(1e-9) * norm(t1) * norm(t2))
                return std::nullopt;
            return n / norm(n);
        }

        /// Normal of a smooth vertex from the tangent masks of Halstead et al., ring of quadrangles all around.
        template <std::floating_point T>
        auto smooth_normal(const Ring<T> &r) -> point<T, 3>
        {
            const auto n = r.d.size();
            const T theta = 2 * std::numbers::pi_v<T> / T(n);
            const T lambda = T(5) / 16 + (std::cos(theta) + std::cos(theta / 2) * std::sqrt(2 * (9 + std::cos(theta)))) / 16;
            const T c = 1 / (4 * lambda - 1);
            point<T, 3> t1{}, t2{};
            for (std::size_t i{}; i < n; ++i)
            {
                const T ci = std::cos(T(i) * theta), cn = std::cos(T(i + 1) * theta), cp = std::cos((T(i) - 1) * theta);
                const auto e = r.e[i] - r.V, d = r.d[i] - r.V;
                t1 += 4 * ci * e + c * (ci + cn) * d;
                t2 += 4 * cp * e + c * (cp + ci) * d;
            }
            return normalized(t1 ^ t2);
        }

        /// The cage with quadrangles all around vertex v: itself, or refined once (v keeps its index).
        template <std::floating_point T>
        auto quads_around(const Cage<T> &c, const Topology &t, std::size_t v) -> std::optional<Cage<T>>
        {
            if (quads(t, t.vertex_faces(v)))
                return std::nullopt;
            return refine(c, t);
        }

        /// Limit position of vertex v of a piece whose faces around v are all there.
        template <std::floating_point T>
        auto vertex_limit(const Cage<T> &c, const Topology &t, std::size_t v) -> point<T, 3>
        {
            const auto &V = c.vertices[v];
            switch (t.rule(v))
            {
            case VertexRule::Corner:
                return V;
            case VertexRule::Crease:
            {
                point<T, 3> s{};
                for (const auto e : t.vertex_edges(v))
                    if (t.is_sharp(e))
                    {
                        const auto &ev = t.edge_vertices(e);
                        s += c.vertices[ev[0] == v ? ev[1] : ev[0]];
                    }
                return (4 * V + s) / T(6);
            }
            case VertexRule::Smooth:
            case VertexRule::Dart:
                break;
            }
            const auto refined = quads_around(c, t, v);
            std::optional<Topology> tr;
            if (refined)
                tr.emplace(*refined, false);
            const auto &cq = refined ? *refined : c;
            const auto &tq = refined ? *tr : t;
            const auto r = ring(cq, tq, tq.fan(v), true);
            if (tq.rule(v) == VertexRule::Dart)
                return limit_position(ring_matrix(VertexRule::Dart, r), r);
            const auto n = T(r.d.size());
            point<T, 3> s = n * n * r.V;
            for (std::size_t i{}; i < r.d.size(); ++i)
                s += 4 * r.e[i] + r.d[i];
            return s / (n * (n + 5));
        }

        /// Limit normal at corner k of face 0 of a piece, on the side of that face, up to its sign; nullopt where not given (eigen_normal).
        template <std::floating_point T>
        auto vertex_normal(const Cage<T> &c, const Topology &t, std::size_t k) -> std::optional<point<T, 3>>
        {
            const auto v = c.faces[0][k];
            const auto refined = quads_around(c, t, v);
            std::optional<Topology> tr;
            if (refined)
                tr.emplace(*refined, false);
            const auto &cq = refined ? *refined : c;
            const auto &tq = refined ? *tr : t;
            // refined: the child of face 0 at corner k is face k, the child of v at its corner k
            const auto h0 = (refined ? tq.face_begin(k) : tq.face_begin(0)) + k;
            const auto [hs, closed] = tq.sector(h0);
            const auto r = ring(cq, tq, hs, closed);
            if (tq.rule(v) == VertexRule::Smooth)
                return smooth_normal(r);
            return eigen_normal(ring_matrix(tq.rule(v), r), r);
        }

        /// True when corner k of face 0 of a piece lets the face be a B-spline patch.
        inline auto regular_corner(const Topology &t, std::size_t k) -> bool
        {
            const auto h0 = t.face_begin(0) + k;
            const auto v = t.origin(h0);
            switch (t.rule(v))
            {
            case VertexRule::Smooth:
                return t.valence(v) == 4 && t.face_count(v) == 4 && quads(t, t.vertex_faces(v));
            case VertexRule::Dart:
                return false;
            case VertexRule::Crease:
            {
                const auto [hs, closed] = t.sector(h0);
                return !closed && hs.size() == 2 && t.face_size(t.face(hs[0])) == 4 && t.face_size(t.face(hs[1])) == 4;
            }
            case VertexRule::Corner:
            {
                const auto [hs, closed] = t.sector(h0);
                return !closed && hs.size() == 1;
            }
            }
            std::unreachable();
        }

        /// The 16 points of the B-spline patch of face 0 of a piece, all of whose corners are regular.
        template <std::floating_point T>
        auto regular_grid(const Cage<T> &c, const Topology &t) -> Grid<T>
        {
            Grid<T> g{};
            std::array<std::array<bool, 4>, 4> set{};
            constexpr std::array<std::array<int, 2>, 4> corner{{{1, 1}, {2, 1}, {2, 2}, {1, 2}}};
            constexpr std::array<std::array<int, 2>, 4> dir{{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}}; // counterclockwise
            const auto put = [&](int i, int j, std::size_t v) {
                g[i][j] = c.vertices[v];
                set[i][j] = true;
            };
            // face of half-edge h from corner k, m quarter turns counterclockwise from face 0
            const auto fill = [&](std::size_t k, int m, std::size_t h) {
                const auto r = dir[std::size_t((int(k) + m + 8) % 4)], s = dir[std::size_t((int(k) + m + 9) % 4)];
                const auto [i, j] = corner[k];
                put(i + r[0], j + r[1], t.target(h));
                put(i + r[0] + s[0], j + r[1] + s[1], t.target(t.next(h)));
                put(i + s[0], j + s[1], t.origin(t.prev(h)));
            };
            for (std::size_t k{}; k < 4; ++k)
            {
                const auto h0 = t.face_begin(0) + k;
                auto h = h0;
                fill(k, 0, h);
                for (int m = 1; m < 4; ++m)
                {
                    const auto p = t.prev(h);
                    if (t.is_sharp(t.edge(p)))
                        break;
                    h = t.twin(p);
                    fill(k, m, h);
                }
                h = h0;
                for (int m = -1; m > -4; --m)
                {
                    if (t.is_sharp(t.edge(h)))
                        break;
                    h = t.next(t.twin(h));
                    fill(k, m, h);
                }
            }
            // phantom points beyond the sharp sides: 0 bottom (v = 0), 1 right (u = 1), 2 top, 3 left
            std::array<bool, 4> sharp;
            for (std::size_t s{}; s < 4; ++s)
                sharp[s] = t.is_sharp(t.edge(t.face_begin(0) + s));
            const auto phantom = [&](int i, int j, int i1, int j1, int i2, int j2) {
                g[i][j] = 2 * g[i1][j1] - g[i2][j2];
                set[i][j] = true;
            };
            for (int j = sharp[0] ? 1 : 0; j <= (sharp[2] ? 2 : 3); ++j)
            {
                if (sharp[3])
                    phantom(0, j, 1, j, 2, j);
                if (sharp[1])
                    phantom(3, j, 2, j, 1, j);
            }
            for (int i = 0; i < 4; ++i)
            {
                if (sharp[0])
                    phantom(i, 0, i, 1, i, 2);
                if (sharp[2])
                    phantom(i, 3, i, 2, i, 1);
            }
            for (const auto &col : set)
                if (!std::ranges::all_of(col, std::identity{}))
                    throw std::logic_error{"subdivision: incomplete regular patch"};
            return g;
        }

        /// Uniform cubic B-spline basis on [0, 1] and its derivative of order d.
        template <std::floating_point T>
        auto basis(T t, std::size_t d) -> std::array<T, 4>
        {
            const T s = 1 - t;
            switch (d)
            {
            case 0:
                return {s * s * s / 6, (3 * t * t * t - 6 * t * t + 4) / 6, (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6, t * t * t / 6};
            case 1:
                return {-s * s / 2, (3 * t * t - 4 * t) / 2, (-3 * t * t + 2 * t + 1) / 2, t * t / 2};
            case 2:
                return {s, 3 * t - 2, 1 - 3 * t, t};
            case 3:
                return {-1, 3, -3, 1};
            default:
                return {};
            }
        }

        template <std::floating_point T>
        auto bicubic(const Grid<T> &g, T u, T v, std::size_t du, std::size_t dv) -> point<T, 3>
        {
            const auto bu = basis(u, du), bv = basis(v, dv);
            point<T, 3> p{};
            for (std::size_t i{}; i < 4; ++i)
                for (std::size_t j{}; j < 4; ++j)
                    p += (bu[i] * bv[j]) * g[i][j];
            return p;
        }

        /// p times 2^e, component by component (exact, unless the result underflows or overflows).
        template <std::floating_point T>
        auto scaled(point<T, 3> p, int e) -> point<T, 3>
        {
            for (auto &x : p)
                x = std::ldexp(x, e);
            return p;
        }

        /// End of the descent into a patch: a regular quarter, or an irregular corner hit exactly.
        template <std::floating_point T>
        struct Located
        {
            Cage<T> piece;            ///< face 0: the last quarter, with the faces around it; a point x of it is origin + 2^exponent x
            Topology topology;        ///< of the piece
            point<T, 3> origin{};
            int exponent{};
            T u{}, v{};               ///< parameters in face 0
            std::size_t depth{};      ///< quarters gone down: derivatives scale by 2^depth per order
            std::size_t corner{none}; ///< corner of face 0 hit, irregular; none: face 0 is regular

            /// A point of the piece, in the frame of the cage.
            auto point_at(const point<T, 3> &x) const -> point<T, 3> { return origin + scaled(x, exponent); }
            /// A derivative of orders summing to `order` on the piece, in the parameters of the patch.
            auto derivative(const point<T, 3> &d, std::size_t order) const -> point<T, 3> { return scaled(d, exponent + int(depth * order)); }
        };

        /**
         * @brief Goes down into face 0 of a piece, a quadrangle, to the quarter holding (u, v).
         *
         * The points of each quarter are taken relative to its first corner and scaled by a power
         * of 2 to a size of about 1, exactly (the rules are affine). Near an extraordinary vertex
         * the quarters shrink geometrically: measured from the origin, their points would keep the
         * rounding of the coordinates of the cage (at a depth of 60 from a vertex of valence 6, on a
         * cage of size 50, the normals were 3° off); at their own scale, they would leave the range
         * of the floating point numbers (at a depth of 500 from a vertex of valence 3, S_u x S_v
         * underflowed; at 2^-1000, S_u did).
         */
        template <std::floating_point T>
        auto locate(Cage<T> piece, T u, T v) -> Located<T>
        {
            constexpr std::array<std::array<T, 2>, 4> uv{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
            // every halving of (u, v) down to the smallest positive number, and some margin
            constexpr std::size_t max_depth = std::numeric_limits<T>::digits - std::numeric_limits<T>::min_exponent + 64;
            point<T, 3> origin{};
            int exponent{};
            for (std::size_t depth{};; ++depth)
            {
                const auto o = piece.vertices[piece.faces[0][0]];
                origin += scaled(o, exponent);
                T size{};
                for (auto &p : piece.vertices)
                {
                    p = p - o;
                    for (const auto x : p)
                        size = std::max(size, std::abs(x));
                }
                if (size > 0 && std::isfinite(size))
                {
                    const int e = std::ilogb(size);
                    for (auto &p : piece.vertices)
                        p = scaled(p, -e);
                    exponent += e;
                }
                Topology t{piece, false};
                std::array<bool, 4> regular;
                for (std::size_t k{}; k < 4; ++k)
                    regular[k] = regular_corner(t, k);
                if (std::ranges::all_of(regular, std::identity{}))
                    return {std::move(piece), std::move(t), origin, exponent, u, v, depth, none};
                for (std::size_t k{}; k < 4; ++k)
                    if (!regular[k] && u == uv[k][0] && v == uv[k][1])
                        return {std::move(piece), std::move(t), origin, exponent, u, v, depth, k};
                if (depth > max_depth)
                    throw std::logic_error{"subdivision: the evaluation does not reach a regular patch"};
                const auto refined = refine(piece, t);
                const int i = u >= T(0.5), j = v >= T(0.5);
                const std::size_t k = i ? (j ? 2 : 1) : (j ? 3 : 0); // the child of corner k runs along face 0
                piece = face_ring(refined, Topology{refined, false}, k, 0);
                u = 2 * u - T(i);
                v = 2 * v - T(j);
            }
        }
    }

    /**
     * @brief Limit surface of a Catmull-Clark cage, evaluated exactly (see the file comment).
     *
     * The cage is checked at construction (`CageError`). Evaluation costs one subdivision of a
     * few faces per halving of the distance to an irregular vertex; a regular patch is evaluated
     * directly.
     */
    template <std::floating_point T>
    class LimitSurface
    {
    public:
        explicit LimitSurface(Cage<T> cage) : cage_{std::move(cage)}, topo_{cage_}
        {
            patch_begin_.reserve(topo_.n_faces() + 1);
            patch_begin_.push_back(0);
            for (std::size_t f{}; f < topo_.n_faces(); ++f)
                patch_begin_.push_back(patch_begin_.back() + (topo_.face_size(f) == 4 ? 1 : topo_.face_size(f)));
        }

        auto cage() const noexcept -> const Cage<T> & { return cage_; }
        auto topology() const noexcept -> const Topology & { return topo_; }

        auto patch_count() const noexcept -> std::size_t { return patch_begin_.back(); }
        /// Patch of face f: corner of an n-gon (0 for a quadrangle).
        auto patch(std::size_t f, std::size_t corner = 0) const -> std::size_t
        {
            if (f >= topo_.n_faces() || corner >= patch_begin_[f + 1] - patch_begin_[f])
                throw std::out_of_range{"subdivision: no patch at face " + std::to_string(f) + ", corner " + std::to_string(corner)};
            return patch_begin_[f] + corner;
        }
        /// {face, corner of an n-gon (0 for a quadrangle)} of a patch.
        auto patch_face(std::size_t p) const -> std::array<std::size_t, 2>
        {
            if (p >= patch_count())
                throw std::out_of_range{"subdivision: no patch " + std::to_string(p)};
            const auto f = std::size_t(std::ranges::upper_bound(patch_begin_, p) - patch_begin_.begin()) - 1;
            return {f, p - patch_begin_[f]};
        }

        /**
         * @brief Point of patch p at (u, v), or its derivative of orders (du, dv).
         *
         * Throws OutOfBoundsEval outside [0, 1]², std::domain_error for a derivative at an
         * irregular vertex (it does not exist there).
         */
        auto value(std::size_t p, T u, T v, std::size_t du = 0, std::size_t dv = 0) const -> point<T, 3>
        {
            const auto l = locate(p, u, v);
            if (l.corner != none)
            {
                if (du + dv > 0)
                    throw std::domain_error{"subdivision: no derivative of the limit surface at an irregular vertex"};
                return l.point_at(limit_detail::vertex_limit(l.piece, l.topology, l.piece.faces[0][l.corner]));
            }
            const auto d = limit_detail::bicubic(limit_detail::regular_grid(l.piece, l.topology), l.u, l.v, du, dv);
            return du + dv == 0 ? l.point_at(d) : l.derivative(d, du + dv);
        }

        /**
         * @brief Unit normal of patch p at (u, v), S_u x S_v: outwards for faces counterclockwise
         * seen from outside. At an irregular vertex, the normal of the tangent plane on the side
         * of the patch; std::domain_error where there is none (see the file comment) or the
         * surface is degenerate.
         */
        auto normal(std::size_t p, T u, T v) const -> point<T, 3>
        {
            const auto l = locate(p, u, v);
            if (l.corner != none)
            {
                const auto n = limit_detail::vertex_normal(l.piece, l.topology, l.corner);
                if (!n)
                    throw std::domain_error{"subdivision: no tangent plane given at this irregular vertex (tied eigenvalues of its sector)"};
                // its sign: that of the normals of the patch next to the vertex, which tend to it
                const T d = std::ldexp(T(1), -16);
                return *n * normal(p, u == 0 ? d : 1 - d, v == 0 ? d : 1 - d) < 0 ? T(-1) * *n : *n;
            }
            // on the piece, of a size about 1: S_u and S_v keep their direction, without underflow
            const auto g = limit_detail::regular_grid(l.piece, l.topology);
            const auto su = limit_detail::bicubic(g, l.u, l.v, 1, 0), sv = limit_detail::bicubic(g, l.u, l.v, 0, 1);
            const auto n = su ^ sv;
            if (norm(n) <= std::numeric_limits<T>::epsilon() * norm(su) * norm(sv))
                throw std::domain_error{"subdivision: degenerate limit surface, no normal"};
            return n / norm(n);
        }

        /// Limit position of a vertex of the cage.
        auto vertex_limit(std::size_t vertex) const -> point<T, 3>
        {
            if (vertex >= cage_.vertices.size())
                throw std::out_of_range{"subdivision: no vertex " + std::to_string(vertex)};
            const auto piece = limit_detail::vertex_ring(cage_, topo_, vertex);
            return limit_detail::vertex_limit(piece, Topology{piece, false}, 0);
        }

    private:
        Cage<T> cage_;
        Topology topo_;
        std::vector<std::size_t> patch_begin_;

        auto locate(std::size_t p, T u, T v) const -> limit_detail::Located<T>
        {
            if (!(u >= 0 && u <= 1))
                throw OutOfBoundsEval<T>{u, {0, 1}, "Subdivision patch U eval "};
            if (!(v >= 0 && v <= 1))
                throw OutOfBoundsEval<T>{v, {0, 1}, "Subdivision patch V eval "};
            const auto [f, corner] = patch_face(p);
            auto piece = limit_detail::face_ring(cage_, topo_, f, 0);
            if (topo_.face_size(f) != 4)
            {
                const auto refined = refine(piece, Topology{piece, false});
                piece = limit_detail::face_ring(refined, Topology{refined, false}, corner, 0);
            }
            return limit_detail::locate(std::move(piece), u, v);
        }
    };
}
