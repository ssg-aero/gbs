#pragma once

/**
 * @file sew.h
 * @brief Native BREP core, stage 1: simple sewing of faces into shells.
 *
 * Design: docs/sources/design/brep_core.md, section 6.3; architecture note
 * docs/sources/design/brep_pr07_sew.md.
 *
 * Two free edges are merged when their ends coincide within the sewing
 * tolerance and when each curve stays within that tolerance of the other
 * (sampled both ways). No edge is cut (T-junctions stay free), no geometry
 * is moved: the gaps are absorbed in the edge and vertex tolerances, as in
 * BRepBuilderAPI_Sewing without cutting.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gbs-brep/model.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/builders.h>
#include <gbs/bscinterp.h>

namespace gbs::brep
{
    template <std::floating_point T>
    struct SewOptions
    {
        T tol = brep_default_tolerance<T>; ///< sewing tolerance: max distance between two edges to merge
        std::size_t n_samples = 10;        ///< points per edge for the distance test (each way)
        std::size_t n_pcurve_max = 513;    ///< max samples to rebuild a pcurve on the surviving edge
    };

    struct SewReport
    {
        std::vector<ShellId> shells;                       ///< one per connected component of faces
        std::vector<std::pair<EdgeId, EdgeId>> merged;     ///< (survivor, absorbed) pairs
        std::vector<EdgeId> free_edges;                    ///< edges still used once after sewing
        std::vector<std::pair<EdgeId, EdgeId>> rejected;   ///< ends coincide but the curves part by more than tol
        std::vector<std::pair<EdgeId, EdgeId>> ambiguous;  ///< valid pair not kept because one edge was already sewn
        bool orientable{true};                             ///< false if some component cannot be consistently oriented
    };

    namespace detail
    {
        /// Gauss-Newton projection of p on crv restricted to [a, b], from the seed s. Returns {s, distance}.
        template <std::floating_point T>
        auto project_on_curve(const Curve<T, 3> &crv, T a, T b, const point<T, 3> &p, T s) -> std::pair<T, T>
        {
            s = std::clamp(s, a, b);
            for (int it = 0; it < 30; ++it)
            {
                const auto r = crv.value(s) - p;
                const auto d = crv.value(s, 1);
                const T den = d * d;
                if (!(den > T(0)))
                    break;
                const T next = std::clamp(s - (d * r) / den, a, b);
                const bool done = std::abs(next - s) <= T(1e-15) * (b - a);
                s = next;
                if (done)
                    break;
            }
            return {s, distance(crv.value(s), p)};
        }

        /// Projection with a coarse global search when the seeded one fails.
        template <std::floating_point T>
        auto project_on_curve(const Curve<T, 3> &crv, T a, T b, const point<T, 3> &p, T s, T tol) -> std::pair<T, T>
        {
            auto r = project_on_curve(crv, a, b, p, s);
            if (r.second <= tol)
                return r;
            constexpr std::size_t n = 33;
            T best = a, dbest = std::numeric_limits<T>::max();
            for (std::size_t i{}; i < n; ++i)
            {
                const T si = a + (b - a) * T(i) / T(n - 1);
                const T di = distance(crv.value(si), p);
                if (di < dbest)
                    dbest = di, best = si;
            }
            auto g = project_on_curve(crv, a, b, p, best);
            return g.second < r.second ? g : r;
        }

        /// Seed in [s1, s2] for the parameter t in [u1, u2], same or opposite sense.
        template <std::floating_point T>
        T seed(T t, T u1, T u2, T s1, T s2, bool same)
        {
            const T x = (t - u1) / (u2 - u1);
            return same ? s1 + x * (s2 - s1) : s2 - x * (s2 - s1);
        }

        /// Max distance from n samples of c1 on [u1,u2] to c2 on [s1,s2] (seeded by the parameter map).
        template <std::floating_point T>
        T one_way(const Curve<T, 3> &c1, T u1, T u2, const Curve<T, 3> &c2, T s1, T s2, bool same, std::size_t n, T tol)
        {
            T worst{};
            for (std::size_t i{}; i < n; ++i)
            {
                const T t = u1 + (u2 - u1) * T(i) / T(n - 1);
                worst = std::max(worst, project_on_curve(c2, s1, s2, c1.value(t), seed(t, u1, u2, s1, s2, same), tol).second);
            }
            return worst;
        }

        /**
         * Pcurve of the absorbed edge's co-edge, expressed in the survivor's
         * parameter: q(t) = p2(phi(t)), phi(t) = projection of C1(t) on C2.
         * Interpolated at n samples, doubled until the deviation
         * |S2(q(t)) - C1(t)| reaches `target` or `n_max` is reached.
         * Returns the pcurve and the measured deviation.
         */
        template <std::floating_point T>
        auto recompose_pcurve(const Edge<T> &e1, const Edge<T> &e2, bool same, const Curve<T, 2> &p2, const Surface<T, 3> &srf,
                              T target, std::size_t n_max, T tol) -> std::pair<std::shared_ptr<Curve<T, 2>>, T>
        {
            std::size_t n = 17;
            for (;;)
            {
                std::vector<T> t(n);
                std::vector<constrType<T, 2, 1>> q(n);
                T dev{};
                for (std::size_t i{}; i < n; ++i)
                {
                    t[i] = e1.u1 + (e1.u2 - e1.u1) * T(i) / T(n - 1);
                    const auto [s, d] = project_on_curve(*e2.curve, e2.u1, e2.u2, e1.curve->value(t[i]),
                                                         seed(t[i], e1.u1, e1.u2, e2.u1, e2.u2, same), tol);
                    q[i] = {p2.value(s)};
                }
                auto pc = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(q, t, std::min<std::size_t>(3, n - 1)));
                for (std::size_t i{}; i < n; ++i)
                {
                    const auto uv = pc->value(t[i]);
                    dev = std::max(dev, distance(srf.value(uv[0], uv[1]), e1.curve->value(t[i])));
                    if (i + 1 < n)
                    {
                        const T tm = (t[i] + t[i + 1]) / 2;
                        const auto w = pc->value(tm);
                        dev = std::max(dev, distance(srf.value(w[0], w[1]), e1.curve->value(tm)));
                    }
                }
                if (dev <= target || n >= n_max)
                    return {pc, dev};
                n = std::min(2 * n - 1, n_max);
            }
        }

        struct Use
        {
            FaceId face;
            CoEdgeRef coedge;
        };

        template <std::floating_point T>
        auto edge_uses(const Model<T> &m, std::span<const FaceId> faces) -> std::unordered_map<EdgeId, std::vector<Use>>
        {
            std::unordered_map<EdgeId, std::vector<Use>> uses;
            for (auto fid : faces)
                for (auto wid : m.face(fid).wires)
                {
                    const auto &w = m.wire(wid);
                    for (std::uint32_t i{}; i < w.coedges.size(); ++i)
                        if (!m.edge(w.coedges[i].edge).degenerate)
                            uses[w.coedges[i].edge].push_back(Use{fid, CoEdgeRef{wid, i}});
                }
            return uses;
        }
    } // namespace detail

    /**
     * @brief Sews `faces` into shells.
     *
     * 1. Free edges (used once by the faces, not degenerate) are candidates;
     *    pairs are pre-selected by overlapping bounding boxes (sweep on x).
     * 2. A pair is valid if its ends coincide within `tol` (same or opposite
     *    sense) and if each curve stays within `tol` of the other on
     *    `n_samples` points (both ways, which rejects partial overlaps).
     * 3. Valid pairs are merged from the closest to the farthest, each edge at
     *    most once (manifold result). The survivor (smallest id) keeps its
     *    curve; the absorbed edge's co-edge is redirected to it, its sense
     *    composed with the relative sense, its pcurve rebuilt in the
     *    survivor's parameter. Vertices at matched ends are merged (smallest
     *    id survives, point not moved). Tolerances are raised to the measured
     *    gaps; absorbed edges and vertices are erased.
     * 4. Faces are grouped by connected component and oriented consistently
     *    (breadth-first from the smallest face id, Forward): one shell each.
     *    A component that cannot be oriented (Möbius) is still returned and
     *    reported through `orientable = false`.
     *
     * Fails only on invalid input, then the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto sew(Model<T> &m, std::span<const FaceId> faces, std::type_identity_t<SewOptions<T>> opts = {}) -> BuildResult<SewReport>
    {
        if (faces.empty())
            return detail::fail(BuildErrc::EmptyInput, "sew needs at least one face");
        if (!detail::valid_tolerance(opts.tol))
            return detail::fail(BuildErrc::InvalidTolerance, "sewing tolerance must be finite and > 0");
        for (std::size_t i{}; i < faces.size(); ++i)
        {
            if (!m.alive(faces[i]))
                return detail::fail(BuildErrc::InvalidId, "sew: dead or invalid face", {faces[i]});
            for (std::size_t j{}; j < i; ++j)
                if (faces[j] == faces[i])
                    return detail::fail(BuildErrc::InvalidId, "sew: face given twice", {faces[i]});
            const auto &f = m.face(faces[i]);
            if (!f.surface)
                return detail::fail(BuildErrc::NullSurface, "sew: face without surface", {faces[i]});
            for (auto wid : f.wires)
            {
                if (!m.alive(wid))
                    return detail::fail(BuildErrc::InvalidId, "sew: dead wire", {wid});
                for (const auto &ce : m.wire(wid).coedges)
                    if (!m.alive(ce.edge) || !ce.pcurve)
                        return detail::fail(BuildErrc::InvalidId, "sew: face wires need live edges with pcurves", {wid});
            }
        }
        const T tol = opts.tol;
        const std::size_t ns = std::max<std::size_t>(opts.n_samples, 2);
        SewReport report;

        // ---- 1. candidates -------------------------------------------------
        auto uses = detail::edge_uses(m, faces);
        struct Cand
        {
            EdgeId e;
            BoundingBox<T> box;
            point<T, 3> a, b; // curve ends
        };
        std::vector<Cand> cands;
        for (const auto &[eid, u] : uses)
        {
            const auto &e = m.edge(eid);
            if (u.size() != 1 || !e.curve)
                continue;
            auto box = bounding_box(m, eid, ns);
            box.inflate(tol);
            cands.push_back(Cand{eid, box, e.curve->value(e.u1), e.curve->value(e.u2)});
        }
        std::ranges::sort(cands, [](const Cand &x, const Cand &y) { return x.box.min[0] < y.box.min[0] || (x.box.min[0] == y.box.min[0] && x.e < y.e); });

        // ---- 2. valid pairs ----------------------------------------------------
        struct Pair
        {
            EdgeId e1, e2; // e1 < e2
            bool same;
            T dist;
        };
        std::vector<Pair> pairs;
        for (std::size_t i{}; i < cands.size(); ++i)
            for (std::size_t j{i + 1}; j < cands.size() && cands[j].box.min[0] <= cands[i].box.max[0]; ++j)
            {
                if (!cands[i].box.intersects(cands[j].box))
                    continue;
                const auto &ci = cands[i].e < cands[j].e ? cands[i] : cands[j];
                const auto &cj = cands[i].e < cands[j].e ? cands[j] : cands[i];
                const bool same_ends = distance(ci.a, cj.a) <= tol && distance(ci.b, cj.b) <= tol;
                const bool opp_ends = distance(ci.a, cj.b) <= tol && distance(ci.b, cj.a) <= tol;
                if (!same_ends && !opp_ends)
                    continue;
                const auto &e1 = m.edge(ci.e);
                const auto &e2 = m.edge(cj.e);
                T best = std::numeric_limits<T>::max();
                bool best_same = true;
                for (bool same : {true, false})
                {
                    if ((same && !same_ends) || (!same && !opp_ends))
                        continue;
                    const T d = std::max(detail::one_way(*e1.curve, e1.u1, e1.u2, *e2.curve, e2.u1, e2.u2, same, ns, tol),
                                         detail::one_way(*e2.curve, e2.u1, e2.u2, *e1.curve, e1.u1, e1.u2, same, ns, tol));
                    if (d < best)
                        best = d, best_same = same;
                }
                if (best <= tol)
                    pairs.push_back(Pair{ci.e, cj.e, best_same, best});
                else
                    report.rejected.emplace_back(ci.e, cj.e);
            }
        std::ranges::sort(pairs, [](const Pair &x, const Pair &y) { return x.dist < y.dist || (x.dist == y.dist && std::pair{x.e1, x.e2} < std::pair{y.e1, y.e2}); });

        // ---- 3. merge ------------------------------------------------------------
        std::unordered_map<VertexId, VertexId> parent;
        auto find = [&](VertexId v) {
            while (parent.contains(v) && parent[v] != v)
                v = parent[v];
            return v;
        };
        auto unite = [&](VertexId a, VertexId b) {
            a = find(a), b = find(b);
            if (a == b)
                return;
            if (b < a)
                std::swap(a, b);
            parent[a] = a;
            parent[b] = a;
        };
        std::unordered_map<EdgeId, bool> sewn;
        std::vector<EdgeId> touched;
        for (const auto &p : pairs)
        {
            if (sewn[p.e1] || sewn[p.e2])
            {
                report.ambiguous.emplace_back(p.e1, p.e2);
                continue;
            }
            sewn[p.e1] = sewn[p.e2] = true;
            auto &e1 = m.edge(p.e1);
            const auto e2 = m.edge(p.e2); // copy: e2 is erased below
            const auto use = uses[p.e2].front();
            auto &ce = m.wire(use.coedge.wire).coedges[use.coedge.index];
            const auto &srf = *m.face(use.face).surface;
            auto [pc, dev] = detail::recompose_pcurve(e1, e2, p.same, *ce.pcurve, srf, std::max(tol, p.dist), opts.n_pcurve_max, tol);
            ce.edge = p.e1;
            ce.orient = p.same ? ce.orient : reverse(ce.orient);
            ce.pcurve = std::move(pc);
            e1.tol = std::max({e1.tol, e2.tol, p.dist, dev});
            unite(e1.v1, p.same ? e2.v1 : e2.v2);
            unite(e1.v2, p.same ? e2.v2 : e2.v1);
            m.erase(p.e2);
            uses[p.e1].push_back(use);
            uses.erase(p.e2);
            touched.push_back(p.e1);
            report.merged.emplace_back(p.e1, p.e2);
        }
        // apply the vertex merge to every live edge of the model
        if (!parent.empty())
        {
            for (auto eid : m.template ids<EdgeId>())
            {
                auto &e = m.edge(eid);
                const auto n1 = find(e.v1), n2 = find(e.v2);
                if (n1 != e.v1 || n2 != e.v2)
                {
                    e.v1 = n1;
                    e.v2 = n2;
                    touched.push_back(eid);
                }
            }
            for (const auto &[v, pv] : parent)
                if (auto r = find(v); r != v)
                {
                    auto &rv = m.vertex(r);
                    rv.tol = std::max(rv.tol, m.vertex(v).tol);
                    m.erase(v);
                }
        }
        for (auto eid : touched)
            if (m.alive(eid))
                detail::fit_vertex_tolerances(m, m.edge(eid));

        // ---- 4. components and orientation ----------------------------------------
        std::unordered_map<FaceId, std::vector<std::pair<FaceId, EdgeId>>> adj;
        for (const auto &[eid, u] : uses)
            if (u.size() == 2 && u[0].face != u[1].face)
            {
                adj[u[0].face].emplace_back(u[1].face, eid);
                adj[u[1].face].emplace_back(u[0].face, eid);
            }
        auto coedge_of = [&](FaceId f, EdgeId e) -> const CoEdge<T> & {
            for (const auto &u : uses[e])
                if (u.face == f)
                    return m.wire(u.coedge.wire).coedges[u.coedge.index];
            throw BRepError("sew: inconsistent usage table"); // unreachable
        };
        std::vector<FaceId> order(faces.begin(), faces.end());
        std::ranges::sort(order);
        std::unordered_map<FaceId, Orientation> orient;
        for (auto start : order)
        {
            if (orient.contains(start))
                continue;
            Shell sh;
            std::queue<FaceId> q;
            orient[start] = Orientation::Forward;
            q.push(start);
            while (!q.empty())
            {
                const auto f = q.front();
                q.pop();
                sh.faces.push_back(FaceUse{f, orient[f]});
                for (const auto &[g, e] : adj[f])
                {
                    const auto eff = compose(coedge_of(f, e).orient, orient[f]);
                    const auto want = compose(coedge_of(g, e).orient, reverse(eff)); // g's use must run opposite
                    if (auto it = orient.find(g); it != orient.end())
                    {
                        if (it->second != want)
                            report.orientable = false;
                        continue;
                    }
                    orient[g] = want;
                    q.push(g);
                }
            }
            const auto sid = m.add(std::move(sh));
            m.shell(sid).closed = is_closed(m, sid);
            report.shells.push_back(sid);
        }
        for (const auto &[eid, u] : uses)
            if (u.size() == 1)
                report.free_edges.push_back(eid);
        std::ranges::sort(report.free_edges);
        return report;
    }

    template <std::floating_point T>
    [[nodiscard]] auto sew(Model<T> &m, const std::vector<FaceId> &faces, std::type_identity_t<SewOptions<T>> opts = {}) -> BuildResult<SewReport>
    {
        return sew(m, std::span<const FaceId>{faces}, opts);
    }

} // namespace gbs::brep
