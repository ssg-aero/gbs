#pragma once
/**
 * @file cdt.h
 * @brief Constrained Delaunay triangulation in the plane (stage 3 of the BREP,
 * `docs/sources/design/brep_palier3_mesh.md`, § 4.3).
 *
 * 1. Points are inserted one at a time (Bowyer-Watson) in a super-triangle enclosing a given
 *    box: the triangles whose circumcircle holds the point form a cavity, retriangulated as a fan
 *    around it. The cavity never crosses a constrained edge, except the one the point lies on,
 *    which is split. A point within `eps` of a vertex is that vertex.
 * 2. A constrained edge missing from the triangulation is recovered: the edges it crosses, found
 *    by walking from its first end, are flipped (Sloan), then the new edges are flipped back to
 *    Delaunay where they may. A vertex on the segment splits the constraint there.
 * 3. The triangles outside the domain are removed: from the super-triangle, the number of
 *    constrained edges crossed (its parity) tells the inside from the outside, holes included.
 *    Before that, `domain` and `place` tell the triangles and the points of the domain, so that
 *    a mesher may refine it: a point inserted inside keeps the triangles of its cavity inside
 *    (the cavity never crosses a constrained edge).
 *
 * The orientation and incircle tests are Shewchuk's exact-sign predicates
 * (`topology/robustPredicates.h`): the same decisions on every platform. Triangles are
 * counter-clockwise. Point indices are those returned by `insert`, in their order.
 *
 * The half-edge Delaunay of `inc/topology` (`tessellations.h`) is not used: its elements are
 * linked by shared pointers and its insertion visits every face (6.8 s for 20 000 random points
 * against 0.05 s here), and its constrained version handles neither holes nor inner constraints.
 */
#include <gbs/gbslib.h>
#include <topology/robustPredicates.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gbs
{
    template <std::floating_point T>
    class ConstrainedDelaunay
    {
    public:
        using index = std::uint32_t;
        static constexpr index none = std::numeric_limits<index>::max();

        /// A triangulation covering the box {x1, x2, y1, y2}: the points inserted must lie in it.
        explicit ConstrainedDelaunay(const std::array<T, 4> &box)
        {
            const T cx = (box[0] + box[1]) / T(2), cy = (box[2] + box[3]) / T(2);
            const T r = std::max({box[1] - box[0], box[3] - box[2], std::numeric_limits<T>::min()});
            // (points closer than a thousand roundings of the box are one point)
            eps_ = T(1e3) * std::numeric_limits<T>::epsilon() * (r + std::abs(cx) + std::abs(cy));
            const T big = T(20) * r;
            pts_ = {{cx - big, cy - big}, {cx + big, cy - big}, {cx, cy + big}};
            tris_.push_back({{0, 1, 2}, {none, none, none}, {false, false, false}, true});
            stamp_.push_back(0);
            inside_.push_back(0);
            vt_ = {0, 0, 0};
        }

        /// Inserts p; returns its index (an existing one when p is within eps of a vertex). Throws
        /// `std::invalid_argument` for a point that is not finite or outside the super-triangle.
        index insert(const point<T, 2> &p)
        {
            if (!std::isfinite(p[0]) || !std::isfinite(p[1]))
                throw std::invalid_argument("ConstrainedDelaunay: a point that is not finite");
            const auto [t, on_edge] = locate(p);
            for (index k : tris_[t].v)
                if (k >= 3 && dist2(pts_[k], p) <= eps_ * eps_)
                    return k - 3;
            const index v = static_cast<index>(pts_.size());
            pts_.push_back(p);
            vt_.push_back(t);
            bowyer_watson(v, t, on_edge);
            return v - 3;
        }

        /// The edge between the points a and b is kept in the triangulation (split at the points on
        /// it). Throws `std::runtime_error` when it crosses another constrained edge.
        void constrain(index a, index b)
        {
            const index ia = a + 3, ib = b + 3;
            if (ia == ib)
                return;
            classified_ = false;
            constrain_internal(ia, ib, 0);
        }

        /**
         * Removes the triangles outside the domain bounded by the constrained edges: those reached
         * from the super-triangle across an even number of them (a constrained edge that bounds
         * nothing is crossed both ways and changes nothing). After this, no point may be inserted.
         */
        void remove_outside()
        {
            classify();
            for (index t = 0; t < tris_.size(); ++t)
                if (tris_[t].alive && !inside_[t])
                    tris_[t].alive = false;
            // the boundary of the domain: no neighbour beyond
            for (auto &tr : tris_)
                if (tr.alive)
                    for (auto &n : tr.n)
                        if (n != none && !tris_[n].alive)
                            n = none;
            closed_ = true;
        }

        /// Where a point falls: see `place`.
        enum class Place
        {
            outside,     ///< outside the domain, on a constrained edge, or outside the box
            inside,      ///< inside the domain
            encroaching, ///< inside, but in the diametral circle of a constrained edge of its triangle
        };

        /**
         * Where p falls in the domain that `remove_outside` would leave: outside it (or on a
         * constrained edge), inside it, or inside it but encroaching a constrained edge of the
         * triangle holding it (in the circle of which that edge is a diameter, as Ruppert's
         * refinement tells).
         */
        Place place(const point<T, 2> &p)
        {
            if (!std::isfinite(p[0]) || !std::isfinite(p[1]))
                return Place::outside;
            index t, on;
            try
            {
                std::tie(t, on) = locate(p);
            }
            catch (const std::invalid_argument &)
            {
                return Place::outside;
            }
            if (on != none && tris_[t].c[on])
                return Place::outside;
            classify();
            if (!inside_[t])
                return Place::outside;
            for (index i = 0; i < 3; ++i)
            {
                if (!tris_[t].c[i])
                    continue;
                const auto &a = pts_[tris_[t].v[(i + 1) % 3]], &b = pts_[tris_[t].v[(i + 2) % 3]];
                if ((a[0] - p[0]) * (b[0] - p[0]) + (a[1] - p[1]) * (b[1] - p[1]) < T(0))
                    return Place::encroaching;
            }
            return Place::inside;
        }

        /// The triangles (counter-clockwise) of the domain that `remove_outside` would leave.
        [[nodiscard]] auto domain() -> std::vector<std::array<index, 3>>
        {
            classify();
            std::vector<std::array<index, 3>> r;
            for (index t = 0; t < tris_.size(); ++t)
                if (tris_[t].alive && inside_[t])
                    r.push_back({tris_[t].v[0] - 3, tris_[t].v[1] - 3, tris_[t].v[2] - 3});
            return r;
        }

        /// Is (a, b, c), counter-clockwise, a triangle of the triangulation?
        [[nodiscard]] bool has_triangle(index a, index b, index c) const
        {
            const auto e = find_edge(a + 3, b + 3);
            if (!e)
                return false;
            // (the edge as found may be seen from the other side: (b, a, x))
            index t = e->first;
            if (tris_[t].v[(e->second + 1) % 3] != a + 3)
                t = tris_[t].n[e->second];
            return t != none && tris_[t].alive && slot_of(t, a + 3) != none && slot_of(t, b + 3) != none && slot_of(t, c + 3) != none &&
                   tris_[t].v[(slot_of(t, a + 3) + 1) % 3] == b + 3;
        }

        /// The points, in the order of their insertion.
        [[nodiscard]] auto points() const -> std::vector<point<T, 2>> { return {pts_.begin() + 3, pts_.end()}; }
        [[nodiscard]] auto point_at(index i) const -> const point<T, 2> & { return pts_[i + 3]; }
        [[nodiscard]] std::size_t size() const { return pts_.size() - 3; }

        /// The triangles (counter-clockwise), those touching the super-triangle left out.
        [[nodiscard]] auto triangles() const -> std::vector<std::array<index, 3>>
        {
            std::vector<std::array<index, 3>> r;
            for (index t = 0; t < tris_.size(); ++t)
                if (tris_[t].alive && !touches_super(t))
                    r.push_back({tris_[t].v[0] - 3, tris_[t].v[1] - 3, tris_[t].v[2] - 3});
            return r;
        }

        /// The structure is consistent: neighbours reciprocal across the same edge with the same
        /// constraint, triangles counter-clockwise (for the tests).
        [[nodiscard]] bool valid() const
        {
            for (index t = 0; t < tris_.size(); ++t)
            {
                const auto &tr = tris_[t];
                if (!tr.alive)
                    continue;
                if (!(orient(tr.v[0], tr.v[1], tr.v[2]) > T(0)))
                    return false;
                for (index i = 0; i < 3; ++i)
                {
                    const index n = tr.n[i];
                    if (n == none)
                        continue;
                    if (!tris_[n].alive)
                        return false;
                    const index a = tr.v[(i + 1) % 3], b = tr.v[(i + 2) % 3];
                    bool back = false;
                    for (index j = 0; j < 3; ++j)
                        if (tris_[n].v[(j + 1) % 3] == b && tris_[n].v[(j + 2) % 3] == a)
                            back = tris_[n].n[j] == t && tris_[n].c[j] == tr.c[i];
                    if (!back)
                        return false;
                }
            }
            return true;
        }

        /// Is the edge between the points a and b constrained?
        [[nodiscard]] bool constrained(index a, index b) const
        {
            const auto e = find_edge(a + 3, b + 3);
            return e && tris_[e->first].c[e->second];
        }

    private:
        struct Triangle
        {
            std::array<index, 3> v;  ///< counter-clockwise
            std::array<index, 3> n;  ///< n[i]: across the edge opposite v[i], (v[i + 1], v[i + 2])
            std::array<bool, 3> c;   ///< c[i]: that edge is constrained
            bool alive;
        };

        std::vector<point<T, 2>> pts_;     ///< 0..2: the super-triangle
        std::vector<Triangle> tris_;
        std::vector<index> vt_;            ///< a triangle at each vertex
        std::vector<index> free_;          ///< dead triangles to reuse
        std::vector<std::uint32_t> stamp_; ///< == epoch_: in the cavity of the current insertion
        std::vector<char> inside_;         ///< in the domain (when classified_)
        std::uint32_t epoch_{};
        T eps_{};
        bool closed_{};
        bool classified_{};
        index last_{0};

        static T dist2(const point<T, 2> &a, const point<T, 2> &b)
        {
            const T dx = a[0] - b[0], dy = a[1] - b[1];
            return dx * dx + dy * dy;
        }
        /// Exact sign: > 0 when (a, b, c) is counter-clockwise.
        T orient(index a, index b, index c) const { return orient2d(pts_[a], pts_[b], pts_[c]); }
        /// Exact sign: > 0 when d is inside the circumcircle of the counter-clockwise (a, b, c).
        T in_circle(index a, index b, index c, const point<T, 2> &d) const { return incircle(pts_[a], pts_[b], pts_[c], d); }

        bool touches_super(index t) const { return tris_[t].v[0] < 3 || tris_[t].v[1] < 3 || tris_[t].v[2] < 3; }

        index slot_of(index t, index v) const
        {
            for (index i = 0; i < 3; ++i)
                if (tris_[t].v[i] == v)
                    return i;
            return none;
        }

        index new_triangle(const Triangle &t)
        {
            if (!free_.empty())
            {
                const index i = free_.back();
                free_.pop_back();
                tris_[i] = t;
                return i;
            }
            tris_.push_back(t);
            stamp_.push_back(0);
            inside_.push_back(0);
            return static_cast<index>(tris_.size() - 1);
        }

        /**
         * Tells the triangles in the domain: from the super-triangle, an odd number of constrained
         * edges crossed (0-1 breadth first: crossing a constrained edge adds one; a constrained
         * edge that bounds nothing is crossed both ways and changes nothing).
         */
        void classify()
        {
            if (classified_)
                return;
            std::vector<int> depth(tris_.size(), -1);
            std::deque<index> queue;
            for (index t = 0; t < tris_.size(); ++t)
                if (tris_[t].alive && touches_super(t))
                {
                    depth[t] = 0;
                    queue.push_back(t);
                }
            while (!queue.empty())
            {
                const index t = queue.front();
                queue.pop_front();
                for (int i = 0; i < 3; ++i)
                {
                    const index n = tris_[t].n[i];
                    if (n == none || !tris_[n].alive)
                        continue;
                    const int d = depth[t] + (tris_[t].c[i] ? 1 : 0);
                    if (depth[n] < 0 || d < depth[n])
                    {
                        depth[n] = d;
                        tris_[t].c[i] ? queue.push_back(n) : queue.push_front(n);
                    }
                }
            }
            for (index t = 0; t < tris_.size(); ++t)
                inside_[t] = tris_[t].alive && depth[t] % 2 == 1 && !touches_super(t);
            classified_ = true;
        }

        /// The triangle holding p, and the edge p lies on (none when strictly inside).
        std::pair<index, index> locate(const point<T, 2> &p)
        {
            if (closed_)
                throw std::logic_error("ConstrainedDelaunay: no insertion after remove_outside");
            index t = (last_ < tris_.size() && tris_[last_].alive) ? last_ : 0;
            while (!tris_[t].alive)
                ++t;
            for (std::size_t step = 0; step < 4 * tris_.size() + 16; ++step)
            {
                const auto &tr = tris_[t];
                index next = none, on = none;
                for (index i = 0; i < 3 && next == none; ++i)
                {
                    const T o = orient2d(pts_[tr.v[(i + 1) % 3]], pts_[tr.v[(i + 2) % 3]], p);
                    if (o < T(0))
                    {
                        if (tr.n[i] == none)
                            throw std::invalid_argument("ConstrainedDelaunay: point outside the box");
                        next = tr.n[i];
                    }
                    else if (o == T(0))
                        on = i;
                }
                if (next == none)
                {
                    last_ = t;
                    return {t, on};
                }
                t = next;
            }
            throw std::runtime_error("ConstrainedDelaunay: the walk to a point does not end");
        }

        /// Inserts the vertex v found in the triangle t (on its edge `on_edge`, if any).
        void bowyer_watson(index v, index t, index on_edge)
        {
            const auto &p = pts_[v];
            if (++epoch_ == 0) // (wrapped: the stamps start again)
            {
                std::ranges::fill(stamp_, 0);
                epoch_ = 1;
            }
            auto in = [&](index x) { return stamp_[x] == epoch_; };
            // the cavity: triangles whose circumcircle holds p, reached without crossing a
            // constrained edge, except the one p lies on
            std::vector<index> cavity{t};
            stamp_[t] = epoch_;
            const char region = inside_[t]; // (the cavity does not cross a constrained edge: one region)
            std::pair<index, index> split{none, none}; // a constrained edge that p splits
            if (on_edge != none)
            {
                const index n = tris_[t].n[on_edge];
                if (tris_[t].c[on_edge])
                    split = {tris_[t].v[(on_edge + 1) % 3], tris_[t].v[(on_edge + 2) % 3]};
                if (n != none)
                {
                    cavity.push_back(n);
                    stamp_[n] = epoch_;
                }
            }
            for (std::size_t k = 0; k < cavity.size(); ++k)
            {
                const index c = cavity[k];
                for (index i = 0; i < 3; ++i)
                {
                    const index n = tris_[c].n[i];
                    if (n == none || in(n) || tris_[c].c[i])
                        continue;
                    const auto &tn = tris_[n];
                    if (in_circle(tn.v[0], tn.v[1], tn.v[2], p) > T(0))
                    {
                        stamp_[n] = epoch_;
                        cavity.push_back(n);
                    }
                }
            }
            // its boundary: edges (a, b) counter-clockwise, with the triangle beyond and their constraint
            struct Side
            {
                index a, b, beyond;
                bool c;
            };
            std::vector<Side> sides;
            for (index c : cavity)
                for (index i = 0; i < 3; ++i)
                {
                    const index n = tris_[c].n[i];
                    if (n != none && in(n))
                        continue;
                    sides.push_back({tris_[c].v[(i + 1) % 3], tris_[c].v[(i + 2) % 3], n, tris_[c].c[i]});
                }
            for (index c : cavity)
            {
                tris_[c].alive = false;
                stamp_[c] = 0;
                free_.push_back(c);
            }
            // the fan around v: (a, b, v) for each side, linked through their shared ends
            std::unordered_map<index, index> starting_at, ending_at;
            std::vector<index> fresh;
            fresh.reserve(sides.size());
            for (const auto &s : sides)
            {
                const index nt = new_triangle({{s.a, s.b, v}, {none, none, s.beyond}, {false, false, s.c}, true});
                fresh.push_back(nt);
                starting_at[s.a] = nt;
                ending_at[s.b] = nt;
                if (s.beyond != none)
                {
                    // the triangle beyond: its edge (b, a) now faces nt
                    auto &tb = tris_[s.beyond];
                    for (index i = 0; i < 3; ++i)
                        if (tb.v[(i + 1) % 3] == s.b && tb.v[(i + 2) % 3] == s.a)
                            tb.n[i] = nt;
                }
                vt_[s.a] = nt;
                vt_[s.b] = nt;
            }
            vt_[v] = fresh.front();
            for (index nt : fresh)
            {
                auto &tr = tris_[nt];
                tr.n[0] = starting_at.at(tr.v[1]); // across (b, v): the triangle starting at b
                tr.n[1] = ending_at.at(tr.v[0]);   // across (v, a): the triangle ending at a
            }
            for (index nt : fresh)
                inside_[nt] = region;
            if (split.first != none)
            {
                classified_ = false; // (both sides of the constrained edge in the cavity)
                for (index end : {split.first, split.second})
                    if (const auto e = find_edge(v, end))
                        mark(e->first, e->second, true);
            }
            last_ = vt_[v];
        }

        /// The triangle and the slot of the edge (u, w) as seen from it (u then w counter-clockwise).
        std::optional<std::pair<index, index>> find_edge(index u, index w) const
        {
            const index start = vt_[u];
            if (start == none || !tris_[start].alive || slot_of(start, u) == none)
                return scan_edge(u, w);
            // around u both ways: a fan may stop on the boundary of the domain
            for (bool ccw : {true, false})
            {
                index t = start;
                for (std::size_t guard = 0; guard < tris_.size() + 4; ++guard)
                {
                    const index k = slot_of(t, u);
                    if (tris_[t].v[(k + 1) % 3] == w)
                        return std::pair{t, (k + 2) % 3};
                    if (tris_[t].v[(k + 2) % 3] == w)
                        return std::pair{t, (k + 1) % 3};
                    const index n = tris_[t].n[ccw ? (k + 2) % 3 : (k + 1) % 3];
                    if (n == none || n == start)
                        break;
                    t = n;
                }
            }
            return std::nullopt;
        }
        std::optional<std::pair<index, index>> scan_edge(index u, index w) const
        {
            for (index t = 0; t < tris_.size(); ++t)
            {
                if (!tris_[t].alive)
                    continue;
                const index k = slot_of(t, u);
                if (k == none)
                    continue;
                if (tris_[t].v[(k + 1) % 3] == w)
                    return std::pair{t, (k + 2) % 3};
                if (tris_[t].v[(k + 2) % 3] == w)
                    return std::pair{t, (k + 1) % 3};
            }
            return std::nullopt;
        }

        /// Marks the edge in slot i of t, and in its neighbour.
        void mark(index t, index i, bool c)
        {
            tris_[t].c[i] = c;
            const index n = tris_[t].n[i];
            if (n == none)
                return;
            const index a = tris_[t].v[(i + 1) % 3], b = tris_[t].v[(i + 2) % 3];
            for (index j = 0; j < 3; ++j)
                if (tris_[n].v[(j + 1) % 3] == b && tris_[n].v[(j + 2) % 3] == a)
                    tris_[n].c[j] = c;
        }

        /// Flips the edge in slot i of t; returns the two triangles.
        std::pair<index, index> flip(index t, index i)
        {
            const index n = tris_[t].n[i];
            auto &A = tris_[t];
            auto &B = tris_[n];
            const index x = A.v[i], u = A.v[(i + 1) % 3], w = A.v[(i + 2) % 3];
            index j = 0;
            while (B.v[j] == u || B.v[j] == w)
                ++j;
            const index y = B.v[j];
            // A = (x, u, w), B = (y, w, u): after, A = (x, u, y), B = (y, w, x)
            const index a_xu = A.n[(i + 2) % 3]; // across (x, u)
            const bool c_xu = A.c[(i + 2) % 3];
            const index a_wx = A.n[(i + 1) % 3]; // across (w, x)
            const bool c_wx = A.c[(i + 1) % 3];
            const index b_yw = B.n[(j + 2) % 3]; // across (y, w)
            const bool c_yw = B.c[(j + 2) % 3];
            const index b_uy = B.n[(j + 1) % 3]; // across (u, y)
            const bool c_uy = B.c[(j + 1) % 3];
            A = {{x, u, y}, {b_uy, n, a_xu}, {c_uy, false, c_xu}, true};
            B = {{y, w, x}, {a_wx, t, b_yw}, {c_wx, false, c_yw}, true};
            auto relink = [&](index outer, index from, index to) {
                if (outer == none)
                    return;
                for (auto &k : tris_[outer].n)
                    if (k == from)
                        k = to;
            };
            relink(b_uy, n, t);
            relink(a_wx, t, n);
            vt_[x] = t;
            vt_[u] = t;
            vt_[y] = t;
            vt_[w] = n;
            return {t, n};
        }

        /// Does the open segment (p, q) cross the open segment (a, b)?
        bool crosses(index p, index q, index a, index b) const
        {
            const T o1 = orient(p, q, a), o2 = orient(p, q, b), o3 = orient(a, b, p), o4 = orient(a, b, q);
            return ((o1 > T(0) && o2 < T(0)) || (o1 < T(0) && o2 > T(0))) && ((o3 > T(0) && o4 < T(0)) || (o3 < T(0) && o4 > T(0)));
        }

        /// Is c on the open segment (a, b)?
        bool on_segment(index a, index b, index c) const
        {
            if (orient(a, b, c) != T(0))
                return false;
            const auto &pa = pts_[a], &pb = pts_[b], &pc = pts_[c];
            const T s = (pc[0] - pa[0]) * (pb[0] - pa[0]) + (pc[1] - pa[1]) * (pb[1] - pa[1]);
            return s > T(0) && s < dist2(pa, pb);
        }

        /**
         * The edges crossed by the segment (a, b), walking from a (Sloan): first the triangle of
         * the fan of a that the segment enters, then across each crossed edge. A vertex met on
         * the segment is returned instead (the constraint is split there).
         */
        std::pair<std::vector<std::pair<index, index>>, index> walk(index a, index b) const
        {
            // the fan of a: (a, p, q) counter-clockwise with b strictly in the angle (p, q)
            index t = none, p = none, q = none;
            for (bool ccw : {true, false})
            {
                index s = vt_[a];
                for (std::size_t guard = 0; guard < tris_.size() + 4; ++guard)
                {
                    const index k = slot_of(s, a);
                    const index x = tris_[s].v[(k + 1) % 3], y = tris_[s].v[(k + 2) % 3];
                    if (on_segment(a, b, x))
                        return {{}, x};
                    if (on_segment(a, b, y))
                        return {{}, y};
                    if (orient(a, x, b) > T(0) && orient(a, y, b) < T(0))
                    {
                        t = s;
                        p = x;
                        q = y;
                        break;
                    }
                    const index n = tris_[s].n[ccw ? (k + 2) % 3 : (k + 1) % 3];
                    if (n == none || n == vt_[a])
                        break;
                    s = n;
                }
                if (t != none)
                    break;
            }
            if (t == none)
                throw std::runtime_error("ConstrainedDelaunay: no triangle at the start of a constraint");
            // p on the right of a -> b, q on its left; the edge (p, q) is crossed
            std::vector<std::pair<index, index>> crossed;
            for (std::size_t guard = 0; guard < tris_.size() + 4; ++guard)
            {
                crossed.emplace_back(p, q);
                index i = 0; // the slot of (p, q) in t: opposite its third vertex
                while (tris_[t].v[i] == p || tris_[t].v[i] == q)
                    ++i;
                if (tris_[t].c[i])
                    throw std::runtime_error("ConstrainedDelaunay: crossing constraints");
                const index n = tris_[t].n[i];
                if (n == none)
                    throw std::runtime_error("ConstrainedDelaunay: a constraint leaves the triangulation");
                index j = 0;
                while (tris_[n].v[j] == p || tris_[n].v[j] == q)
                    ++j;
                const index r = tris_[n].v[j];
                if (r == b)
                    return {crossed, none};
                if (on_segment(a, b, r))
                    return {{}, r};
                (orient(a, b, r) > T(0) ? q : p) = r;
                t = n;
            }
            throw std::runtime_error("ConstrainedDelaunay: the walk along a constraint does not end");
        }

        void constrain_internal(index a, index b, int depth)
        {
            if (depth > 64)
                throw std::runtime_error("ConstrainedDelaunay: constraint recovery does not end");
            if (const auto e = find_edge(a, b))
            {
                mark(e->first, e->second, true);
                return;
            }
            auto [edges, through] = walk(a, b);
            if (through != none)
            {
                constrain_internal(a, through, depth + 1);
                constrain_internal(through, b, depth + 1);
                return;
            }
            std::deque<std::pair<index, index>> crossing(edges.begin(), edges.end());
            std::vector<std::pair<index, index>> created;
            const std::size_t limit = 64 * (crossing.size() + 4) * (crossing.size() + 4);
            for (std::size_t guard = 0; !crossing.empty(); ++guard)
            {
                if (guard > limit)
                    throw std::runtime_error("ConstrainedDelaunay: constraint recovery does not end");
                const auto [eu, ew] = crossing.front();
                crossing.pop_front();
                const auto e = find_edge(eu, ew);
                if (!e)
                    continue;
                const auto [t, i] = *e;
                const index n = tris_[t].n[i];
                if (n == none)
                    continue;
                // (the ends as the triangle orders them: (x, u, w) counter-clockwise)
                const index x = tris_[t].v[i], u = tris_[t].v[(i + 1) % 3], w = tris_[t].v[(i + 2) % 3];
                index j = 0;
                while (tris_[n].v[j] == u || tris_[n].v[j] == w)
                    ++j;
                const index y = tris_[n].v[j];
                // the quadrilateral (x, u, y, w) must be strictly convex to flip
                if (!(orient(x, u, y) > T(0) && orient(y, w, x) > T(0)))
                {
                    crossing.emplace_back(u, w);
                    continue;
                }
                flip(t, i);
                if (crosses(a, b, x, y))
                    crossing.emplace_back(x, y);
                else
                    created.emplace_back(x, y);
            }
            const auto e = find_edge(a, b);
            if (!e)
                throw std::runtime_error("ConstrainedDelaunay: constraint not recovered");
            mark(e->first, e->second, true);
            // the new edges back to Delaunay where they may
            for (bool changed = true; changed;)
            {
                changed = false;
                for (const auto &[u, w] : created)
                {
                    if ((u == a && w == b) || (u == b && w == a))
                        continue;
                    const auto f = find_edge(u, w);
                    if (!f || tris_[f->first].c[f->second])
                        continue;
                    const auto [t, i] = *f;
                    const index n = tris_[t].n[i];
                    if (n == none)
                        continue;
                    index j = 0;
                    while (tris_[n].v[j] == u || tris_[n].v[j] == w)
                        ++j;
                    const index y = tris_[n].v[j], x = tris_[t].v[i];
                    if (in_circle(tris_[t].v[0], tris_[t].v[1], tris_[t].v[2], pts_[y]) > T(0) &&
                        orient(x, tris_[t].v[(i + 1) % 3], y) > T(0) && orient(y, tris_[t].v[(i + 2) % 3], x) > T(0))
                    {
                        flip(t, i);
                        changed = true;
                    }
                }
            }
        }
    };
}
