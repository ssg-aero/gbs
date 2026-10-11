#pragma once

/**
 * @file reader.h
 * @brief STEP reader: BREP topology of a STEP file into a gbs::brep::Model.
 *
 * Design: docs/sources/design/step_reader.md, sections 3.4, 5, 6 and 8;
 * architecture note docs/sources/design/step_pr06_topology.md.
 *
 * The topology of the file is taken as it is, without sewing: vertices and
 * edges are shared by #id, loops keep their order and senses
 * (make_wire_ordered), faces keep the senses of their loops (make_face_use).
 * Geometry goes through gbs-io/step/geometry.h: edges get the exact NURBS of
 * their curve between their vertices, faces the exact NURBS of their surface
 * over the box of their boundary, pcurves are recomputed.
 *
 * Partial import by default: a face that cannot be built is left out and
 * reported, its shell stays open and its solid is not made (the shell is
 * returned instead). In strict mode the first failure fails the reading and
 * the model is left unchanged.
 */

#include <gbs-io/step/p21.h>
#include <gbs-io/step/units.h>
#include <gbs-io/step/geometry.h>
#include <gbs-brep/brep>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <numeric>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gbs::step
{
    struct StepReadOptions
    {
        /// Length of the target unit in millimetres; by default the unit of the model (millimetre for a new model).
        std::optional<double> target_unit_mm{};
        /// Fail on the first face or entity that cannot be read, leaving the model unchanged.
        bool strict{false};
        /// Face building; `tol` is raised to the file uncertainty, `pcurve_tol` is the first accuracy tried.
        brep::MakeFaceOptions<double> face{};
        /// Largest tolerance accepted for a vertex off its curve or an edge off its surface (target units).
        double max_tolerance{1e-2};
        /// Points sampled per edge to bound the surface of a face.
        std::size_t edge_samples{16};
    };

    /// An entity of the file that was not read, with the reason.
    struct StepIssue
    {
        std::uint64_t id{};  ///< #id of the entity
        std::string type;    ///< its STEP type
        std::string message; ///< cause
    };

    struct StepReadReport
    {
        std::vector<std::string> schemas;
        std::string length_unit, angle_unit;   ///< units of the file (first representation context)
        double length_factor{1.};              ///< factor applied to the lengths of the file
        std::optional<double> file_tolerance;  ///< uncertainty of the file, in target units
        std::size_t solids{}, shells{}, faces{}, edges{}, vertices{};
        std::vector<StepIssue> unsupported;    ///< entities skipped (assemblies, faceted BREP…)
        std::vector<StepIssue> failed;         ///< faces and solids not built
        double max_tolerance{};                ///< largest edge tolerance of the result
        std::size_t raised_tolerances{};       ///< edges whose tolerance had to exceed the file uncertainty
        std::size_t sense_mismatches{};        ///< faces whose sense differs from advanced_face.same_sense
        brep::CheckReport check;               ///< check() of the result
    };

    struct StepReadResult
    {
        brep::ShapeId root; ///< the shape read, or a compound of the shapes read
        StepReadReport report;
    };

    namespace detail
    {
        using brep::BuildErrc;
        using brep::BuildError;
        using brep::EdgeId;
        using brep::FaceId;
        using brep::FaceUse;
        using brep::Orientation;
        using brep::ShapeId;
        using brep::ShellId;
        using brep::VertexId;
        using brep::WireId;

        /// Strict mode: stops the reading with this error.
        struct StrictStop
        {
            BuildError error;
        };

        class TopologyReader
        {
            const P21File &f_;
            brep::Model<double> &m_;
            const StepReadOptions &o_;
            double target_mm_;
            StepReadReport rep_;

            std::unordered_map<std::uint64_t, GeometryReader> geo_; // by representation context
            GeometryReader *g_{};                                     // current context
            double base_tol_{gbs::brep_default_tolerance<double>};

            std::unordered_map<std::uint64_t, VertexId> vertices_;
            /// edge_curve -> its edges from its v1 to its v2 (several once split at a seam)
            std::unordered_map<std::uint64_t, std::optional<std::vector<EdgeId>>> edges_;
            std::unordered_map<std::uint64_t, std::optional<FaceUse>> faces_;
            std::unordered_map<std::uint64_t, std::pair<ShellId, bool>> shells_; // shell, complete

            [[noreturn]] void stop(BuildErrc c, const std::string &msg, std::vector<ShapeId> shapes = {})
            {
                throw StrictStop{BuildError{c, msg, std::move(shapes)}};
            }

            void unsupported(std::uint64_t id, std::string_view type, std::string msg)
            {
                if (o_.strict)
                    stop(BuildErrc::UnsupportedEntity, "#" + std::to_string(id) + " " + std::string(type) + ": " + msg);
                rep_.unsupported.push_back({id, std::string(type), std::move(msg)});
            }

            void failed(std::uint64_t id, std::string_view type, std::string msg, BuildErrc code = BuildErrc::InvalidFile)
            {
                if (o_.strict)
                    stop(code, "#" + std::to_string(id) + " " + std::string(type) + ": " + msg);
                rep_.failed.push_back({id, std::string(type), std::move(msg)});
            }

            void use_context(std::optional<std::uint64_t> ctx)
            {
                const std::uint64_t key = ctx.value_or(0);
                auto it = geo_.find(key);
                if (it == geo_.end())
                {
                    Units u = ctx ? read_units(f_, *ctx, target_mm_) : Units{.length = 1. / target_mm_};
                    if (geo_.empty())
                    {
                        rep_.length_unit = u.length_name;
                        rep_.angle_unit = u.angle_name;
                        rep_.length_factor = u.length;
                        rep_.file_tolerance = u.uncertainty;
                        base_tol_ = std::max(gbs::brep_default_tolerance<double>, u.uncertainty.value_or(0.));
                    }
                    it = geo_.emplace(key, GeometryReader{f_, std::move(u)}).first;
                }
                g_ = &it->second;
            }

            void tag(const ShapeId &s, const InstanceView &in)
            {
                m_.setExternalId(s, static_cast<std::int64_t>(in.id()));
                const auto a = in.args();
                if (a.size() > 0 && a[0].is(ValueKind::String) && !a[0].as_string().empty())
                    m_.setName(s, std::string(a[0].as_string()));
            }

            // ---------------------------------------------------------------- vertices and edges

            VertexId vertex(std::uint64_t id)
            {
                if (auto it = vertices_.find(id); it != vertices_.end())
                    return it->second;
                const auto in = f_.instance(id);
                if (!in.has_type("VERTEX_POINT"))
                    throw StepUnsupported(id, in.type());
                const auto v = brep::unwrap(brep::make_vertex(m_, g_->point<3>(in.args()[1].as_ref()), base_tol_));
                tag(v, in);
                vertices_.emplace(id, v);
                return v;
            }

            static std::shared_ptr<Curve<double, 3>> reversed_nurbs(const std::shared_ptr<Curve<double, 3>> &c)
            {
                if (auto r = std::dynamic_pointer_cast<BSCurveRational<double, 3>>(c))
                    return std::make_shared<BSCurveRational<double, 3>>(r->reversed());
                if (auto p = std::dynamic_pointer_cast<BSCurve<double, 3>>(c))
                    return std::make_shared<BSCurve<double, 3>>(p->reversed());
                throw std::runtime_error("edge curve is not a NURBS");
            }

            EdgeId build_edge(std::uint64_t id)
            {
                const auto in = f_.instance(id);
                if (!in.has_type("EDGE_CURVE"))
                    throw StepUnsupported(id, in.type());
                const auto a = in.args();
                const auto v1 = vertex(a[1].as_ref()), v2 = vertex(a[2].as_ref());
                const auto c = g_->curve<3>(a[3].as_ref());
                const bool same_sense = a[4].as_bool();
                const auto p1 = m_.vertex(v1).pnt, p2 = m_.vertex(v2).pnt;

                // range [t1, t2] of the edge on the curve, increasing along the curve
                const bool forward = same_sense != c->reversed(); // v1 -> v2 along increasing parameter
                double t1, t2;
                if (c->periodic())
                {
                    double s1 = c->parameter(p1), s2 = c->parameter(p2);
                    if (!forward)
                        std::swap(s1, s2);
                    t1 = s1;
                    t2 = v1 == v2 ? s1 + 2. * std::numbers::pi : detail::wrap_to(s2, s1 + 1e-12);
                }
                else
                {
                    const auto [d1, d2] = c->domain();
                    if (v1 == v2)
                    {
                        if (!std::isfinite(d1) || !std::isfinite(d2))
                            throw StepError(id, "closed edge on an unbounded curve");
                        t1 = d1, t2 = d2;
                    }
                    else
                    {
                        t1 = c->parameter(p1), t2 = c->parameter(p2);
                        if (!forward)
                            std::swap(t1, t2);
                    }
                }
                if (!(t2 > t1))
                    throw StepError(id, "edge vertices are not in the sense of the curve");

                auto nc = c->to_nurbs(t1, t2);
                auto curve = nc.curve;
                double u1 = nc.u1, u2 = nc.u2;
                if (!forward)
                {
                    const auto [k1, k2] = curve->bounds();
                    curve = reversed_nurbs(curve);
                    std::tie(u1, u2) = std::pair{k1 + k2 - nc.u2, k1 + k2 - nc.u1};
                }
                // tolerance: the vertices must contain the curve ends
                const double d = std::max(distance(curve->value(u1), p1), distance(curve->value(u2), p2));
                if (d > o_.max_tolerance)
                    throw StepError(id, "vertex " + std::to_string(d) + " away from its curve");
                const double tol = std::max(base_tol_, 1.001 * d);
                if (d > base_tol_)
                    ++rep_.raised_tolerances;
                const auto e = brep::unwrap(brep::make_edge(m_, curve, u1, u2, v1, v2, tol));
                tag(e, in);
                return e;
            }

            /// Edges of an edge_curve, or nullopt if it failed (reported once).
            std::optional<std::vector<EdgeId>> edge(std::uint64_t id)
            {
                if (auto it = edges_.find(id); it != edges_.end())
                    return it->second;
                std::optional<std::vector<EdgeId>> r;
                try
                {
                    r = std::vector<EdgeId>{build_edge(id)};
                }
                catch (const StepUnsupported &e)
                {
                    unsupported(id, e.type(), e.what());
                }
                catch (const StrictStop &)
                {
                    throw;
                }
                catch (const std::exception &e)
                {
                    failed(id, "EDGE_CURVE", e.what());
                }
                edges_.emplace(id, r);
                return r;
            }

            // ---------------------------------------------------------------- faces

            struct Loop
            {
                std::vector<brep::OrientedEdge> coedges; // empty for a vertex loop
                std::vector<Point<3>> samples;
                bool outer{false};
            };

            void sample(const brep::OrientedEdge &ce, std::vector<Point<3>> &pts) const
            {
                const auto &e = m_.edge(ce.edge);
                const std::size_t n = std::max<std::size_t>(o_.edge_samples, 2);
                for (std::size_t i = 0; i + 1 < n; ++i) // the end is the start of the next co-edge
                {
                    double s = double(i) / double(n - 1);
                    if (ce.orient == Orientation::Reversed)
                        s = 1. - s;
                    pts.push_back(brep::edge_point(m_, ce.edge, e.u1 + s * (e.u2 - e.u1)));
                }
            }

            /// Reads the bounds of a face; nullopt (reported) if one of its edges failed.
            std::optional<std::vector<Loop>> loops(const InstanceView &face)
            {
                std::vector<Loop> r;
                for (auto bref : face.args()[1].as_list())
                {
                    const auto b = f_.instance(bref.as_ref());
                    if (!b.has_type("FACE_BOUND") && !b.has_type("FACE_OUTER_BOUND"))
                        throw StepUnsupported(b.id(), b.type());
                    const auto ba = b.args();
                    const bool orient = ba[2].as_bool();
                    const auto lp = f_.instance(ba[1].as_ref());
                    Loop l;
                    l.outer = b.has_type("FACE_OUTER_BOUND");
                    if (lp.has_type("VERTEX_LOOP"))
                    {
                        l.samples.push_back(m_.vertex(vertex(lp.args()[1].as_ref())).pnt);
                        r.push_back(std::move(l));
                        continue;
                    }
                    if (!lp.has_type("EDGE_LOOP"))
                        throw StepUnsupported(lp.id(), lp.type());
                    for (auto oref : lp.args()[1].as_list())
                    {
                        const auto oe = f_.instance(oref.as_ref());
                        if (!oe.has_type("ORIENTED_EDGE"))
                            throw StepUnsupported(oe.id(), oe.type());
                        const auto oa = oe.args();
                        const auto e = edge(oa[3].as_ref());
                        if (!e)
                            return std::nullopt;
                        if (oa[4].as_bool())
                            for (auto x : *e)
                                l.coedges.push_back({x, Orientation::Forward});
                        else
                            for (auto it = e->rbegin(); it != e->rend(); ++it)
                                l.coedges.push_back({*it, Orientation::Reversed});
                    }
                    if (!orient) // the loop is used in the opposite sense
                    {
                        std::ranges::reverse(l.coedges);
                        for (auto &ce : l.coedges)
                            ce.orient = brep::reverse(ce.orient);
                    }
                    for (const auto &ce : l.coedges)
                        sample(ce, l.samples);
                    r.push_back(std::move(l));
                }
                return r;
            }

            // ---------------------------------------------------------------- seams

            static double wrap_pi(double d) { return d - detail::two_pi * std::round(d / detail::two_pi); }

            /// u of a point of a surface closed in u, or nullopt on its axis (no angle).
            static std::optional<double> angle_of(const SurfaceDef &s, const Point<3> &p)
            {
                const auto uv = s.parameters(p);
                if (!uv)
                    return std::nullopt;
                const auto q = s.value((*uv)[0] + 1e-3, (*uv)[1]);
                if (distance(q, p) < 1e-9 * std::max(1., norm(p)))
                    return std::nullopt;
                return (*uv)[0];
            }

            /// Number of turns of a loop around the u direction (signed), from its samples.
            static int winding(const SurfaceDef &s, const std::vector<Point<3>> &pts)
            {
                std::optional<double> first, prev;
                double total = 0.;
                for (const auto &p : pts)
                    if (auto a = angle_of(s, p))
                    {
                        if (prev)
                            total += wrap_pi(*a - *prev);
                        else
                            first = a;
                        prev = a;
                    }
                if (!first)
                    return 0;
                total += wrap_pi(*first - *prev);
                return static_cast<int>(std::round(total / detail::two_pi));
            }

            /// Unwrapped angular range of a loop that does not wind.
            static std::optional<std::array<double, 2>> angular_range(const SurfaceDef &s, const std::vector<Point<3>> &pts)
            {
                std::optional<double> prev;
                double cur = 0., lo = detail::inf, hi = -detail::inf;
                for (const auto &p : pts)
                    if (auto a = angle_of(s, p))
                    {
                        cur = prev ? cur + wrap_pi(*a - *prev) : *a;
                        prev = a;
                        lo = std::min(lo, cur), hi = std::max(hi, cur);
                    }
                if (!prev)
                    return std::nullopt;
                return std::array{lo, hi};
            }

            /// Splits an edge of the model at t, updating the loops being built and the edge cache.
            void split(EdgeId e, double t, std::vector<Loop> &ls)
            {
                const auto r = brep::unwrap(brep::split_edge(m_, e, t));
                for (auto &l : ls)
                    for (std::size_t i = 0; i < l.coedges.size(); ++i)
                        if (l.coedges[i].edge == e)
                        {
                            const auto o = l.coedges[i].orient;
                            l.coedges.insert(l.coedges.begin() + static_cast<std::ptrdiff_t>(o == Orientation::Forward ? i + 1 : i), {r.second, o});
                            ++i;
                        }
                for (auto &[id, chain] : edges_)
                    if (chain)
                        if (auto it = std::ranges::find(*chain, e); it != chain->end())
                            chain->insert(std::next(it), r.second);
            }

            /// Parameter of an edge where its u crosses `us` (mod 2 pi), strictly inside the edge.
            std::optional<double> crossing(const SurfaceDef &s, EdgeId e, double us) const
            {
                const auto &ed = m_.edge(e);
                auto f = [&](double t) -> std::optional<double> {
                    auto a = angle_of(s, brep::edge_point(m_, e, t));
                    if (!a)
                        return std::nullopt;
                    return wrap_pi(*a - us);
                };
                constexpr int n = 64;
                for (int i = 0; i < n; ++i)
                {
                    double a = ed.u1 + (ed.u2 - ed.u1) * i / n, b = ed.u1 + (ed.u2 - ed.u1) * (i + 1) / n;
                    auto fa = f(a), fb = f(b);
                    if (!fa || !fb || (*fa > 0.) == (*fb > 0.) || std::abs(*fa - *fb) > std::numbers::pi)
                        continue;
                    for (int it = 0; it < 80; ++it)
                    {
                        const double c = 0.5 * (a + b);
                        auto fc = f(c);
                        if (!fc)
                            break;
                        ((*fc > 0.) == (*fa > 0.) ? a : b) = c;
                        if ((*fc > 0.) == (*fa > 0.))
                            fa = fc;
                    }
                    const double t = 0.5 * (a + b);
                    const double eps = 1e-9 * (ed.u2 - ed.u1);
                    if (t > ed.u1 + eps && t < ed.u2 - eps)
                        return t;
                }
                return std::nullopt;
            }

            /// Vertex of loop `li` at the angle us: an existing vertex, or a new one splitting an edge.
            VertexId cross(const SurfaceDef &s, std::vector<Loop> &ls, std::size_t li, double us)
            {
                for (const auto &ce : ls[li].coedges)
                {
                    const auto v = brep::coedge_start(m_, brep::CoEdge<double>{ce.edge, ce.orient, nullptr});
                    if (auto a = angle_of(s, m_.vertex(v).pnt); a && std::abs(wrap_pi(*a - us)) < 1e-9)
                        return v;
                }
                for (const auto &ce : ls[li].coedges)
                    if (auto t = crossing(s, ce.edge, us))
                    {
                        const auto e = ce.edge;
                        split(e, *t, ls);
                        return m_.edge(e).v2;
                    }
                throw StepError(0, "a loop does not reach the angle of the seam");
            }

            /// Rotates a loop so that it starts at vertex v.
            void start_at(Loop &l, VertexId v) const
            {
                auto it = std::ranges::find_if(l.coedges, [&](const brep::OrientedEdge &ce) {
                    return brep::coedge_start(m_, brep::CoEdge<double>{ce.edge, ce.orient, nullptr}) == v;
                });
                std::rotate(l.coedges.begin(), it, l.coedges.end());
            }

            /**
             * Seam of a face that turns around the u direction of its surface without a seam edge
             * (band between two circles, cap bounded by one circle, cone without seam): the seam is
             * the iso-u line at an angle us avoiding the holes, from the first winding loop to the
             * second one, or to the pole on the side of the face. Edges crossing us are split.
             * The winding loops become one loop: A, seam, B, seam reversed (or A, seam, seam reversed).
             * Returns us, or nullopt if the face needs no seam.
             */
            struct SeamPlan
            {
                double us;
                std::vector<std::size_t> winding; // indices of the winding loops
                VertexId xa;
                VertexId xb; // second loop, or the pole
                double va, vb;
                bool pole_created;
            };

            std::optional<SeamPlan> plan_seam(const SurfaceDef &s, std::vector<Loop> &ls, bool same_sense)
            {
                const auto per = s.periodic();
                if (!per[0])
                    return std::nullopt;
                std::unordered_map<EdgeId, int> uses;
                for (const auto &l : ls)
                    for (const auto &ce : l.coedges)
                        ++uses[ce.edge];
                if (std::ranges::any_of(uses, [](const auto &u) { return u.second > 1; }))
                    return std::nullopt; // the file has a seam edge
                SeamPlan p{};
                std::vector<int> turns;
                for (std::size_t i = 0; i < ls.size(); ++i)
                    if (!ls[i].coedges.empty())
                        if (int w = winding(s, ls[i].samples); w != 0)
                        {
                            p.winding.push_back(i);
                            turns.push_back(w);
                        }
                if (p.winding.empty())
                    return std::nullopt;
                if (per[1] || p.winding.size() > 2)
                    throw StepError(0, "seam insertion is supported on surfaces closed in u only, with one or two winding loops");

                // angle of the seam: a vertex of a winding loop, else a regular angle, outside the holes
                std::vector<std::array<double, 2>> holes;
                for (std::size_t i = 0; i < ls.size(); ++i)
                    if (std::ranges::find(p.winding, i) == p.winding.end())
                        if (auto r = angular_range(s, ls[i].samples))
                            holes.push_back(*r);
                std::vector<double> candidates;
                for (auto i : p.winding)
                    for (const auto &ce : ls[i].coedges)
                        if (auto a = angle_of(s, m_.vertex(brep::coedge_start(m_, brep::CoEdge<double>{ce.edge, ce.orient, nullptr})).pnt))
                            candidates.push_back(*a);
                for (int k = 0; k < 64; ++k)
                    candidates.push_back(detail::two_pi * k / 64.);
                auto free = [&](double a) {
                    return std::ranges::none_of(holes, [&](const std::array<double, 2> &h) {
                        const double m = 1e-6;
                        const double c = 0.5 * (h[0] + h[1]);
                        return std::abs(wrap_pi(a - c)) <= 0.5 * (h[1] - h[0]) + m;
                    });
                };
                auto it = std::ranges::find_if(candidates, free);
                if (it == candidates.end())
                    throw StepError(0, "no angle for a seam outside the holes");
                p.us = detail::wrap_to(*it, 0.);

                p.xa = cross(s, ls, p.winding[0], p.us);
                start_at(ls[p.winding[0]], p.xa);
                p.va = (*s.parameters(m_.vertex(p.xa).pnt))[1];
                if (p.winding.size() == 2)
                {
                    p.xb = cross(s, ls, p.winding[1], p.us);
                    start_at(ls[p.winding[1]], p.xb);
                    p.vb = (*s.parameters(m_.vertex(p.xb).pnt))[1];
                }
                else
                {
                    // cap: the face lies on the left of its loop seen from its normal, that is towards
                    // +v when the loop turns towards +u around the surface normal
                    const int side = turns[0] * (same_sense ? 1 : -1);
                    const auto dom = s.domain();
                    p.vb = side > 0 ? dom.v2 : dom.v1;
                    if (!std::isfinite(p.vb))
                        throw StepError(0, "a single winding loop on a surface without pole");
                    const auto pole = s.value(p.us, p.vb);
                    std::optional<VertexId> found;
                    for (const auto &[id, v] : vertices_)
                        if (distance(m_.vertex(v).pnt, pole) <= m_.vertex(v).tol)
                            found = v;
                    p.pole_created = !found;
                    p.xb = found ? *found : brep::unwrap(brep::make_vertex(m_, pole, base_tol_));
                }
                return p;
            }

            /**
             * A face covering a full turn in a periodic direction k gets the box [0, 2 pi] (PR 4), but
             * the seam edge of the file may lie elsewhere (OCCT torus whose axis is -z: v seam at
             * 3 pi / 2). The box then starts at the parameter of the seam edge (an edge used twice
             * in the face along an iso-k line).
             */
            void align_with_seam_edges(const SurfaceDef &s, const std::vector<Loop> &ls, UVBox &box) const
            {
                const auto per = s.periodic();
                std::unordered_map<EdgeId, int> uses;
                for (const auto &l : ls)
                    for (const auto &ce : l.coedges)
                        ++uses[ce.edge];
                std::array<double *, 2> b1{&box.u1, &box.v1}, b2{&box.u2, &box.v2};
                for (std::size_t k = 0; k < 2; ++k)
                {
                    if (!per[k] || std::abs(*b2[k] - *b1[k] - detail::two_pi) > 1e-12)
                        continue;
                    for (const auto &[e, n] : uses)
                    {
                        if (n < 2)
                            continue;
                        const auto &ed = m_.edge(e);
                        std::optional<double> val;
                        bool iso = true;
                        for (int i = 0; i <= 4 && iso; ++i)
                        {
                            const auto uv = s.parameters(brep::edge_point(m_, e, ed.u1 + (ed.u2 - ed.u1) * i / 4.));
                            if (!uv)
                            {
                                iso = false;
                                break;
                            }
                            if (!val)
                                val = (*uv)[k];
                            else
                                iso = std::abs(wrap_pi((*uv)[k] - *val)) < 1e-7;
                        }
                        if (iso && val)
                        {
                            *b1[k] = detail::wrap_to(*val, 0.);
                            if (*b1[k] > detail::two_pi - 1e-9)
                                *b1[k] = 0.;
                            *b2[k] = *b1[k] + detail::two_pi;
                            break;
                        }
                    }
                }
            }

            std::optional<FaceUse> build_face(const InstanceView &in)
            {
                const auto a = in.args();
                const bool same_sense = a[3].as_bool();
                const auto expected = same_sense ? Orientation::Forward : Orientation::Reversed;
                const auto surface = g_->surface(a[2].as_ref());
                auto ls = loops(in);
                if (!ls)
                {
                    failed(in.id(), in.type(), "an edge of the face could not be read");
                    return std::nullopt;
                }

                // seam of a face turning around its surface without seam edge
                std::vector<ShapeId> created; // erased if the face fails
                const auto seam = plan_seam(*surface, *ls, same_sense);
                if (seam && seam->pole_created)
                    created.push_back(seam->xb);
                for (auto &l : *ls) // edges may have been split: sample again
                    if (!l.coedges.empty())
                    {
                        l.samples.clear();
                        for (const auto &ce : l.coedges)
                            sample(ce, l.samples);
                    }

                std::vector<std::vector<Point<3>>> samples;
                for (const auto &l : *ls)
                    samples.push_back(l.samples);
                if (seam && seam->winding.size() == 1)
                    samples.push_back({m_.vertex(seam->xb).pnt}); // the pole bounds the box
                // a face without edge loop (none, or only vertex loops) is the whole closed surface
                const bool whole = std::ranges::none_of(*ls, [](const Loop &l) { return !l.coedges.empty(); });
                UVBox box = whole ? surface->domain() : parameter_box(*surface, samples);
                if (seam)
                    box.u1 = seam->us, box.u2 = seam->us + detail::two_pi;
                else
                    align_with_seam_edges(*surface, *ls, box);
                auto nurbs = surface->to_nurbs(box);

                if (seam)
                {
                    // seam edge: the iso-u line of the face's NURBS at its first u
                    std::shared_ptr<Curve<double, 3>> iso;
                    if (auto r = std::dynamic_pointer_cast<BSSurfaceRational<double, 3>>(nurbs))
                        iso = std::make_shared<BSCurveRational<double, 3>>(r->isoU(box.u1));
                    else if (auto q = std::dynamic_pointer_cast<BSSurface<double, 3>>(nurbs))
                        iso = std::make_shared<BSCurve<double, 3>>(q->isoU(box.u1));
                    else
                        throw StepError(in.id(), "surface is not a NURBS");
                    double a = surface->nurbs_parameters(box, seam->us, seam->va)[1];
                    double b = surface->nurbs_parameters(box, seam->us, seam->vb)[1];
                    auto va = seam->xa, vb = seam->xb;
                    const bool up = a < b; // the seam edge runs from xa to xb
                    if (!up)
                        std::swap(a, b), std::swap(va, vb);
                    const auto se = brep::make_edge(m_, iso, a, b, va, vb, base_tol_);
                    if (!se)
                    {
                        for (auto c : created)
                            m_.erase(c);
                        failed(in.id(), in.type(), "seam: " + se.error().message, se.error().code);
                        return std::nullopt;
                    }
                    created.push_back(*se);
                    const auto to_b = up ? Orientation::Forward : Orientation::Reversed;
                    Loop merged;
                    merged.outer = true;
                    merged.coedges = (*ls)[seam->winding[0]].coedges;
                    merged.coedges.push_back({*se, to_b});
                    if (seam->winding.size() == 2)
                    {
                        const auto &lb = (*ls)[seam->winding[1]].coedges;
                        merged.coedges.insert(merged.coedges.end(), lb.begin(), lb.end());
                    }
                    merged.coedges.push_back({*se, brep::reverse(to_b)});
                    std::vector<Loop> rest{std::move(merged)};
                    for (std::size_t i = 0; i < ls->size(); ++i)
                        if (std::ranges::find(seam->winding, i) == seam->winding.end())
                        {
                            (*ls)[i].outer = false;
                            rest.push_back(std::move((*ls)[i]));
                        }
                    *ls = std::move(rest);
                }
                std::erase_if(*ls, [](const Loop &l) { return l.coedges.empty(); }); // vertex loops only bound the box

                if (ls->empty()) // no edge loop: the whole surface (closed sphere or torus)
                {
                    // the vertices of vertex loops are replaced by those of the natural face
                    for (const auto &bref : in.args()[1].as_list())
                    {
                        const auto lp = f_.instance(f_.instance(bref.as_ref()).args()[1].as_ref());
                        if (!lp.has_type("VERTEX_LOOP"))
                            continue;
                        const auto vid = lp.args()[1].as_ref();
                        if (auto it = vertices_.find(vid); it != vertices_.end())
                        {
                            m_.erase(it->second);
                            vertices_.erase(it);
                        }
                    }
                    auto f = brep::make_face(m_, nurbs, base_tol_);
                    if (!f)
                    {
                        failed(in.id(), in.type(), f.error().message, f.error().code);
                        return std::nullopt;
                    }
                    tag(*f, in);
                    return FaceUse{*f, expected};
                }

                // wires, kept free until a face takes them
                std::vector<WireId> wires;
                auto drop_wires = [&] {
                    for (auto w : wires)
                        m_.erase(w);
                    for (auto c : created)
                        m_.erase(c);
                };
                for (const auto &l : *ls)
                {
                    auto w = brep::make_wire_ordered(m_, l.coedges, base_tol_);
                    if (!w)
                    {
                        drop_wires();
                        failed(in.id(), in.type(), w.error().message, w.error().code);
                        return std::nullopt;
                    }
                    wires.push_back(*w);
                }
                // outer bound first; without one, try each loop as outer
                std::vector<std::size_t> order(wires.size());
                std::iota(order.begin(), order.end(), std::size_t{0});
                std::ranges::stable_partition(order, [&](std::size_t i) { return (*ls)[i].outer; });
                const std::size_t candidates = std::ranges::any_of(*ls, [](const Loop &l) { return l.outer; }) ? 1 : wires.size();

                std::optional<BuildError> last;
                for (std::size_t c = 0; c < candidates; ++c)
                {
                    std::vector<WireId> holes;
                    for (std::size_t i = 0; i < wires.size(); ++i)
                        if (i != order[c])
                            holes.push_back(wires[i]);
                    // accuracy ladder: files whose curves are off their surfaces by more than pcurve_tol
                    auto opts = o_.face;
                    opts.tol = std::max(opts.tol, base_tol_);
                    for (;;)
                    {
                        auto r = brep::make_face_use(m_, nurbs, wires[order[c]], holes, opts);
                        if (r)
                        {
                            // the geometry gives the region of the face, same_sense its normal: the loops of
                            // some files run against it (OCCT on a partial torus), which would make their
                            // shell non-orientable; such disagreements are counted
                            if (r->orient != expected)
                                ++rep_.sense_mismatches;
                            tag(r->face, in);
                            return FaceUse{r->face, expected};
                        }
                        last = r.error();
                        const bool accuracy = r.error().code == BuildErrc::EdgeOffSurface || r.error().code == BuildErrc::PCurveApproximation;
                        if (!accuracy || opts.pcurve_tol * 10. > o_.max_tolerance)
                            break;
                        opts.pcurve_tol *= 10.;
                    }
                    if (last->code != BuildErrc::HoleOutsideOuter && last->code != BuildErrc::DegenerateWire)
                        break;
                }
                drop_wires();
                failed(in.id(), in.type(), std::string(brep::to_string(last->code)) + ": " + last->message, last->code);
                return std::nullopt;
            }

            /// Use of a face of a shell, or nullopt if it failed (reported).
            std::optional<FaceUse> face(std::uint64_t id)
            {
                if (auto it = faces_.find(id); it != faces_.end())
                    return it->second;
                std::optional<FaceUse> r;
                const auto in = f_.instance(id);
                try
                {
                    if (in.has_type("ORIENTED_FACE"))
                    {
                        const auto a = in.args();
                        r = face(a[2].as_ref());
                        if (r && !a[3].as_bool())
                            r->orient = brep::reverse(r->orient);
                    }
                    else if (in.has_type("ADVANCED_FACE") || in.has_type("FACE_SURFACE"))
                        r = build_face(in);
                    else
                        unsupported(id, in.type(), "face type not supported");
                }
                catch (const StepUnsupported &e)
                {
                    unsupported(e.id(), e.type(), e.what());
                    failed(id, in.type(), std::string("unsupported entity: ") + e.what(), BuildErrc::UnsupportedEntity);
                }
                catch (const StrictStop &)
                {
                    throw;
                }
                catch (const std::exception &e)
                {
                    failed(id, in.type(), e.what());
                }
                faces_.emplace(id, r);
                return r;
            }

            // ---------------------------------------------------------------- shells and solids

            /// Shell of a closed_shell / open_shell (possibly oriented), and whether all its faces were built.
            std::pair<ShellId, bool> shell(std::uint64_t id)
            {
                if (auto it = shells_.find(id); it != shells_.end())
                    return it->second;
                const auto in = f_.instance(id);
                std::pair<ShellId, bool> r;
                if (in.has_type("ORIENTED_CLOSED_SHELL") || in.has_type("ORIENTED_OPEN_SHELL"))
                {
                    const auto a = in.args();
                    auto [base, complete] = shell(a[2].as_ref());
                    auto s = m_.shell(base);
                    if (!a[3].as_bool())
                        for (auto &fu : s.faces)
                            fu.orient = brep::reverse(fu.orient);
                    r = {m_.add(std::move(s)), complete};
                }
                else if (in.has_type("CLOSED_SHELL") || in.has_type("OPEN_SHELL"))
                {
                    brep::Shell s;
                    bool complete = true;
                    for (auto fref : in.args()[1].as_list())
                    {
                        if (auto fu = face(fref.as_ref()))
                            s.faces.push_back(*fu);
                        else
                            complete = false;
                    }
                    r = {m_.add(std::move(s)), complete};
                }
                else
                    throw StepUnsupported(id, in.type());
                tag(r.first, in);
                shells_.emplace(id, r);
                return r;
            }

            /// Solid of a manifold_solid_brep / brep_with_voids; the outer shell alone if it is incomplete.
            ShapeId solid(const InstanceView &in)
            {
                const auto a = in.args();
                auto [outer, complete] = shell(a[1].as_ref());
                std::vector<ShellId> voids;
                if (in.has_type("BREP_WITH_VOIDS"))
                    for (auto v : a[2].as_list())
                    {
                        auto [sh, ok] = shell(v.as_ref());
                        voids.push_back(sh);
                        complete = complete && ok;
                    }
                if (!complete)
                {
                    failed(in.id(), in.type(), "faces missing: the shell is returned instead of the solid");
                    return outer;
                }
                auto so = brep::make_solid(m_, outer, voids);
                if (!so)
                {
                    failed(in.id(), in.type(), std::string(brep::to_string(so.error().code)) + ": " + so.error().message, so.error().code);
                    return outer;
                }
                tag(*so, in);
                return *so;
            }

            /// Shape of a root item, or nullopt if the item is not a BREP shape.
            std::optional<ShapeId> root(const InstanceView &in)
            {
                if (in.has_type("MANIFOLD_SOLID_BREP") || in.has_type("BREP_WITH_VOIDS"))
                    return solid(in);
                if (in.has_type("SHELL_BASED_SURFACE_MODEL"))
                {
                    std::vector<ShapeId> shells;
                    for (auto s : in.args()[1].as_list())
                        shells.push_back(shell(s.as_ref()).first);
                    if (shells.size() == 1)
                        return shells.front();
                    auto c = brep::unwrap(brep::make_compound(m_, shells));
                    tag(c, in);
                    return c;
                }
                if (in.has_type("FACETED_BREP"))
                    unsupported(in.id(), in.type(), "faceted BREP (design question 13)");
                return std::nullopt; // placements, mapped items (read by rep_geometry)…
            }

            // ---------------------------------------------------------------- products

            /// Name of the product of a product_definition ("" if none).
            std::string product_of(std::uint64_t pd) const
            {
                try
                {
                    const auto pdf = f_.instance(f_.instance(pd).args()[2].as_ref()); // product_definition_formation
                    return std::string(f_.instance(pdf.args()[2].as_ref()).args()[1].as_string());
                }
                catch (const std::exception &)
                {
                    return {};
                }
            }

            /// Representations of the shape of a product_definition (shape_definition_representation).
            std::vector<std::uint64_t> reps_of(std::uint64_t pd) const
            {
                std::vector<std::uint64_t> r;
                for (const auto &sdr : f_.instances_of("SHAPE_DEFINITION_REPRESENTATION"))
                {
                    const auto a = sdr.args();
                    const auto pds = f_.instance(a[0].as_ref());
                    if (pds.args()[2].as_ref() == pd)
                        r.push_back(a[1].as_ref());
                }
                return r;
            }

            /// Name of the product whose shape is the representation `rep` ("" if none).
            std::string product_name(std::uint64_t rep, int depth = 0) const
            {
                try
                {
                    for (const auto &sdr : f_.instances_of("SHAPE_DEFINITION_REPRESENTATION"))
                    {
                        const auto a = sdr.args();
                        if (a[1].as_ref() != rep)
                            continue;
                        const auto pds = f_.instance(a[0].as_ref()); // product_definition_shape
                        const auto def = f_.instance(pds.args()[2].as_ref());
                        if (def.has_type("PRODUCT_DEFINITION"))
                            return product_of(def.id());
                    }
                    if (depth == 0) // a BREP representation linked to the product's shape representation
                        for (const auto &[r1, r2] : plain_links())
                            if (r1 == rep || r2 == rep)
                                if (auto n = product_name(r1 == rep ? r2 : r1, 1); !n.empty())
                                    return n;
                }
                catch (const std::exception &)
                {
                }
                return {};
            }

            // ---------------------------------------------------------------- representations and instances

            using Placement = brep::Placement<double>;

            /// A shape read once from a representation, placed relatively to that representation.
            struct Instance
            {
                ShapeId shape;
                Placement place;
            };

            std::unordered_map<std::uint64_t, std::optional<ShapeId>> item_shapes_;   // root item -> shape read
            std::unordered_map<std::uint64_t, std::vector<Instance>> rep_cache_;      // representation -> instances
            std::unordered_set<std::uint64_t> reps_in_progress_;
            std::vector<ShapeId> originals_;                                            // every shape read from an item
            std::unordered_set<ShapeId> used_;                                          // originals used as an instance
            std::optional<std::vector<std::pair<std::uint64_t, std::uint64_t>>> plain_links_;

            static ListView rep_args(const InstanceView &r)
            {
                // simple instance, or complex with the REPRESENTATION partial entity
                return r.is_complex() && r.part("REPRESENTATION") ? r.args("REPRESENTATION") : r.args();
            }

            static ListView relationship_args(const InstanceView &r)
            {
                return r.is_complex() ? r.args("REPRESENTATION_RELATIONSHIP") : r.args();
            }

            /// (rep_1, rep_2) of the shape_representation_relationships without transformation.
            const std::vector<std::pair<std::uint64_t, std::uint64_t>> &plain_links() const
            {
                if (!plain_links_)
                {
                    auto &l = const_cast<TopologyReader *>(this)->plain_links_.emplace();
                    for (const auto &srr : f_.instances_of("SHAPE_REPRESENTATION_RELATIONSHIP"))
                        if (!srr.has_type("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"))
                        {
                            const auto a = relationship_args(srr);
                            l.emplace_back(a[2].as_ref(), a[3].as_ref());
                        }
                }
                return *plain_links_;
            }

            /// Placement of an axis2_placement_3d of a representation (in the units of its context).
            Placement placement(std::uint64_t axis, std::uint64_t rep)
            {
                use_context(rep_args(f_.instance(rep))[2].as_ref());
                const auto fr = g_->frame(axis);
                return Placement::frame(fr.origin, fr.x, fr.y, fr.z);
            }

            /// Shape of a root item of a representation, read once.
            std::optional<ShapeId> read_item(const InstanceView &item, std::optional<std::uint64_t> ctx, const std::string &product)
            {
                if (auto it = item_shapes_.find(item.id()); it != item_shapes_.end())
                    return it->second;
                std::optional<ShapeId> r;
                try
                {
                    use_context(ctx);
                    r = root(item);
                    if (r)
                    {
                        if (!product.empty() && m_.name(*r).empty())
                            m_.setName(*r, product);
                        originals_.push_back(*r);
                    }
                }
                catch (const StepUnsupported &e)
                {
                    unsupported(e.id(), e.type(), e.what());
                }
                catch (const StrictStop &)
                {
                    throw;
                }
                catch (const std::exception &e)
                {
                    failed(item.id(), item.type(), e.what());
                }
                item_shapes_.emplace(item.id(), r);
                return r;
            }

            /**
             * Instances of a representation in its own coordinates: its BREP items, the instances
             * of its mapped items (mapped representation placed from its mapping origin to the
             * mapping target), and those of the representations linked to it without transformation
             * (the BREP representation of a product's shape representation).
             */
            std::vector<Instance> rep_geometry(std::uint64_t rep)
            {
                if (auto it = rep_cache_.find(rep); it != rep_cache_.end())
                    return it->second;
                if (reps_in_progress_.contains(rep))
                    return {};
                reps_in_progress_.insert(rep);
                std::vector<Instance> r;
                const auto a = rep_args(f_.instance(rep));
                const auto ctx = a[2].as_ref();
                const auto product = product_name(rep);
                for (auto ref : a[1].as_list())
                {
                    const auto item = f_.instance(ref.as_ref());
                    if (item.has_type("MAPPED_ITEM"))
                    {
                        try
                        {
                            const auto ma = item.args();                      // name, mapping_source, mapping_target
                            const auto src = f_.instance(ma[1].as_ref()).args(); // representation_map: mapping_origin, mapped_representation
                            const auto mapped = src[1].as_ref();
                            const auto p = placement(ma[2].as_ref(), rep) * placement(src[0].as_ref(), mapped).inverse();
                            for (const auto &in : rep_geometry(mapped))
                                r.push_back({in.shape, p * in.place});
                        }
                        catch (const StrictStop &)
                        {
                            throw;
                        }
                        catch (const std::exception &e)
                        {
                            failed(item.id(), item.type(), e.what());
                        }
                        continue;
                    }
                    if (auto sh = read_item(item, ctx, product))
                        r.push_back({*sh, Placement{}});
                }
                for (const auto &[r1, r2] : plain_links())
                    if (r1 == rep || r2 == rep)
                        for (const auto &in : rep_geometry(r1 == rep ? r2 : r1))
                            r.push_back(in);
                reps_in_progress_.erase(rep);
                rep_cache_.emplace(rep, r);
                return r;
            }

            /// The shape itself for its first use without move, a placed copy otherwise.
            ShapeId materialize(const ShapeId &s, const Placement &p)
            {
                if (p.is_identity(1e-12) && !used_.contains(s))
                {
                    used_.insert(s);
                    return s;
                }
                return brep::copy_shape(m_, s, p);
            }

            // ---------------------------------------------------------------- assemblies

            struct Child
            {
                std::uint64_t pd, nauo;
                Placement place; // child coordinates -> parent coordinates
            };

            /// Placement of a next_assembly_usage_occurrence, from its context_dependent_shape_representation.
            Placement occurrence_placement(std::uint64_t nauo, std::uint64_t child)
            {
                for (const auto &cdsr : f_.instances_of("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"))
                {
                    const auto a = cdsr.args();
                    const auto pds = f_.instance(a[1].as_ref());
                    if (pds.args()[2].as_ref() != nauo)
                        continue;
                    const auto rel = f_.instance(a[0].as_ref());
                    const auto ra = relationship_args(rel);
                    const auto r1 = ra[2].as_ref(), r2 = ra[3].as_ref();
                    if (!rel.has_type("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"))
                        return {};
                    const auto t = f_.instance(rel.args("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")[0].as_ref());
                    if (!t.has_type("ITEM_DEFINED_TRANSFORMATION"))
                        throw StepUnsupported(t.id(), t.type());
                    const auto ta = t.args(); // name, description, transform_item_1 (in rep_1), transform_item_2 (in rep_2)
                    const auto f1 = placement(ta[2].as_ref(), r1), f2 = placement(ta[3].as_ref(), r2);
                    // the transformation maps item_1 onto item_2: from rep_1 to rep_2; rep_1 is usually the child's
                    const auto child_reps = reps_of(child);
                    const bool child_is_2 = std::ranges::find(child_reps, r2) != child_reps.end() &&
                                            std::ranges::find(child_reps, r1) == child_reps.end();
                    return child_is_2 ? f1 * f2.inverse() : f2 * f1.inverse();
                }
                return {};
            }

            std::vector<Child> children(std::uint64_t pd)
            {
                std::vector<Child> r;
                for (const auto &n : f_.instances_of("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
                {
                    const auto a = n.args(); // id, name, description, relating, related, reference_designator
                    if (a[3].as_ref() != pd)
                        continue;
                    try
                    {
                        r.push_back({a[4].as_ref(), n.id(), occurrence_placement(n.id(), a[4].as_ref())});
                    }
                    catch (const StrictStop &)
                    {
                        throw;
                    }
                    catch (const StepUnsupported &e)
                    {
                        unsupported(e.id(), e.type(), e.what());
                    }
                    catch (const std::exception &e)
                    {
                        failed(n.id(), n.type(), e.what());
                    }
                }
                return r;
            }

            /// Flattened shape of a product: its own geometry and its components, placed by `p`.
            std::optional<ShapeId> product_shape(std::uint64_t pd, const Placement &p, int depth = 0)
            {
                if (depth > 64)
                    throw StepError(pd, "assembly deeper than 64 levels (cycle?)");
                std::vector<ShapeId> parts;
                for (auto rep : reps_of(pd))
                    for (const auto &in : rep_geometry(rep))
                        parts.push_back(materialize(in.shape, p * in.place));
                const auto kids = children(pd);
                for (const auto &c : kids)
                    if (auto s = product_shape(c.pd, p * c.place, depth + 1))
                        parts.push_back(*s);
                if (parts.empty())
                    return std::nullopt;
                const auto name = product_of(pd);
                if (parts.size() == 1 && kids.empty())
                {
                    if (!name.empty() && m_.name(parts.front()).empty())
                        m_.setName(parts.front(), name);
                    return parts.front();
                }
                const auto c = brep::unwrap(brep::make_compound(m_, parts));
                if (!name.empty())
                    m_.setName(c, name);
                m_.setExternalId(c, static_cast<std::int64_t>(pd));
                return c;
            }

            /// Erases a shape never used and everything below it.
            void erase_tree(const ShapeId &s)
            {
                const auto c = brep::explore<brep::CompoundId>(m_, s);
                const auto so = brep::explore<brep::SolidId>(m_, s);
                const auto sh = brep::explore<ShellId>(m_, s);
                const auto f = brep::explore<FaceId>(m_, s);
                const auto w = brep::explore<WireId>(m_, s);
                const auto e = brep::explore<EdgeId>(m_, s);
                const auto v = brep::explore<VertexId>(m_, s);
                for (auto x : c) m_.erase(x);
                for (auto x : so) m_.erase(x);
                for (auto x : sh) m_.erase(x);
                for (auto x : f) m_.erase(x);
                for (auto x : w) m_.erase(x);
                for (auto x : e) m_.erase(x);
                for (auto x : v) m_.erase(x);
                m_.erase(s);
            }

        public:
            TopologyReader(const P21File &f, brep::Model<double> &m, const StepReadOptions &o, double target_mm)
                : f_{f}, m_{m}, o_{o}, target_mm_{target_mm}
            {
            }

            StepReadResult run()
            {
                rep_.schemas.assign(f_.schemas().begin(), f_.schemas().end());
                std::vector<ShapeId> shapes;

                // 1. assemblies: products used by an occurrence and not themselves components
                std::unordered_set<std::uint64_t> relating, related, placed_reps;
                for (const auto &n : f_.instances_of("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
                {
                    const auto a = n.args();
                    relating.insert(a[3].as_ref());
                    related.insert(a[4].as_ref());
                }
                std::vector<std::uint64_t> roots(relating.begin(), relating.end());
                std::ranges::sort(roots);
                for (auto pd : roots)
                    if (!related.contains(pd))
                    {
                        try
                        {
                            if (auto s = product_shape(pd, Placement{}))
                                shapes.push_back(*s);
                        }
                        catch (const StrictStop &)
                        {
                            throw;
                        }
                        catch (const std::exception &e)
                        {
                            failed(pd, "PRODUCT_DEFINITION", e.what());
                        }
                    }

                // 2. representations not reached: standalone parts (mapped or placed ones are instances)
                for (const auto &m : f_.instances_of("REPRESENTATION_MAP"))
                    placed_reps.insert(m.args()[1].as_ref());
                for (const auto &r : f_.instances_of("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"))
                {
                    const auto a = relationship_args(r);
                    placed_reps.insert(a[2].as_ref());
                    placed_reps.insert(a[3].as_ref());
                }
                std::vector<InstanceView> reps;
                for (const char *t : {"ADVANCED_BREP_SHAPE_REPRESENTATION", "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
                                      "SHAPE_REPRESENTATION", "FACETED_BREP_SHAPE_REPRESENTATION"})
                    for (const auto &r : f_.instances_of(t))
                        if (std::ranges::none_of(reps, [&](const InstanceView &x) { return x.id() == r.id(); }))
                            reps.push_back(r);
                std::ranges::sort(reps, {}, [](const InstanceView &x) { return x.id(); });
                for (const auto &r : reps)
                    if (!rep_cache_.contains(r.id()) && !placed_reps.contains(r.id()))
                        for (const auto &in : rep_geometry(r.id()))
                            if (!used_.contains(in.shape) || !in.place.is_identity(1e-12))
                                shapes.push_back(materialize(in.shape, in.place));

                // 3. roots outside any representation, in the first context of the file
                std::optional<std::uint64_t> first_ctx;
                if (auto c = f_.instances_of("GEOMETRIC_REPRESENTATION_CONTEXT"); !c.empty())
                    first_ctx = c.front().id();
                for (const char *t : {"MANIFOLD_SOLID_BREP", "BREP_WITH_VOIDS", "SHELL_BASED_SURFACE_MODEL"})
                    for (const auto &r : f_.instances_of(t))
                        if (!item_shapes_.contains(r.id()))
                            if (auto s = read_item(r, first_ctx, ""))
                                shapes.push_back(materialize(*s, Placement{}));

                // shapes read but only used through moved copies
                for (const auto &s : originals_)
                    if (!used_.contains(s))
                        erase_tree(s);

                if (shapes.empty())
                    throw StrictStop{BuildError{BuildErrc::InvalidFile,
                                                "no BREP shape could be read from the file (" + std::to_string(rep_.failed.size()) + " failed, " +
                                                    std::to_string(rep_.unsupported.size()) + " unsupported entities)",
                                                {}}};
                ShapeId root_shape = shapes.size() == 1 ? shapes.front() : ShapeId{brep::unwrap(brep::make_compound(m_, shapes))};

                rep_.solids = brep::explore<brep::SolidId>(m_, root_shape).size();
                rep_.shells = brep::explore<ShellId>(m_, root_shape).size();
                rep_.faces = brep::explore<FaceId>(m_, root_shape).size();
                rep_.edges = brep::explore<EdgeId>(m_, root_shape).size();
                rep_.vertices = brep::explore<VertexId>(m_, root_shape).size();
                for (auto e : brep::explore<EdgeId>(m_, root_shape))
                    rep_.max_tolerance = std::max(rep_.max_tolerance, m_.edge(e).tol);
                rep_.check = brep::check(m_, root_shape);
                return {root_shape, std::move(rep_)};
            }
        };
    } // namespace detail

    /**
     * @brief Reads the BREP shapes of a STEP file given as text into `m`.
     *
     * Solids (manifold_solid_brep, brep_with_voids), shells of surface models
     * (shell_based_surface_model) and their faces are read from every shape
     * representation of the file; one shape is returned as is, several in a
     * compound. Lengths are converted to the target unit (the model's unit by
     * default; for an empty model, `target_unit_mm` sets it). Errors that fail
     * the whole reading: syntax error, no shape, unit of a non-empty model
     * different from the target, and in strict mode any face or entity that
     * cannot be read. The model is then left unchanged.
     */
    inline auto read_step(brep::Model<double> &m, std::string_view content, const StepReadOptions &opts = {})
        -> brep::BuildResult<StepReadResult>
    {
        using brep::BuildErrc;
        auto file = parse_p21(content);
        if (!file)
            return std::unexpected(brep::BuildError{BuildErrc::InvalidFile,
                                                    "line " + std::to_string(file.error().line) + ", column " + std::to_string(file.error().column) + ": " + file.error().message,
                                                    {}});
        const bool empty_model = m.template ids<brep::VertexId>().empty() && m.template ids<brep::FaceId>().empty() &&
                                 m.template ids<brep::EdgeId>().empty();
        const double target = opts.target_unit_mm.value_or(m.unitScale());
        if (!(target > 0.) || !std::isfinite(target))
            return std::unexpected(brep::BuildError{BuildErrc::InvalidTolerance, "target unit must be positive", {}});
        if (!empty_model && target != m.unitScale())
            return std::unexpected(brep::BuildError{BuildErrc::InvalidFile, "the model unit differs from the target unit", {}});

        // snapshot for the rollback
        const std::array<std::size_t, 7> cap{m.template capacity<brep::ShapeType::Vertex>(), m.template capacity<brep::ShapeType::Edge>(),
                                             m.template capacity<brep::ShapeType::Wire>(), m.template capacity<brep::ShapeType::Face>(),
                                             m.template capacity<brep::ShapeType::Shell>(), m.template capacity<brep::ShapeType::Solid>(),
                                             m.template capacity<brep::ShapeType::Compound>()};
        const double old_scale = m.unitScale();
        auto rollback = [&] {
            auto erase_from = [&]<brep::ShapeType K>(std::size_t from) {
                for (auto i = from; i < m.template capacity<K>(); ++i)
                    m.erase(brep::Id<K>{static_cast<std::uint32_t>(i)});
            };
            erase_from.template operator()<brep::ShapeType::Compound>(cap[6]);
            erase_from.template operator()<brep::ShapeType::Solid>(cap[5]);
            erase_from.template operator()<brep::ShapeType::Shell>(cap[4]);
            erase_from.template operator()<brep::ShapeType::Face>(cap[3]);
            erase_from.template operator()<brep::ShapeType::Wire>(cap[2]);
            erase_from.template operator()<brep::ShapeType::Edge>(cap[1]);
            erase_from.template operator()<brep::ShapeType::Vertex>(cap[0]);
            m.setUnitScale(old_scale);
        };

        m.setUnitScale(target);
        try
        {
            detail::TopologyReader r{*file, m, opts, target};
            return r.run();
        }
        catch (const detail::StrictStop &s)
        {
            rollback();
            return std::unexpected(s.error);
        }
        catch (const std::exception &e) // invalid content outside any face (units, representations)
        {
            rollback();
            return std::unexpected(brep::BuildError{BuildErrc::InvalidFile, e.what(), {}});
        }
    }

    /// Reads a STEP file; see read_step().
    inline auto read_step_file(brep::Model<double> &m, const std::filesystem::path &path, const StepReadOptions &opts = {})
        -> brep::BuildResult<StepReadResult>
    {
        std::ifstream in(path, std::ios::binary);
        if (!in)
            return std::unexpected(brep::BuildError{brep::BuildErrc::InvalidFile, "cannot open " + path.string(), {}});
        std::string content{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        return read_step(m, content, opts);
    }
}
