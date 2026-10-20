#pragma once
/**
 * @file partial.h
 * @brief Faces coincident on a part of their overlap only (a B-spline face lying on another
 * surface over some of its patches, parting from it beyond), split along the border of the
 * common region before a boolean operation.
 *
 * A Bézier patch is analytic: where it lies on another analytic surface over an open set, it lies
 * on it all over. The region where a B-spline face lies on the other surface is then a union of
 * its patches, and its border runs along knot lines: exact iso lines, whatever the tolerance.
 *
 * 1. Detection: the patches of the B-spline face (the carrier) between its knot lines are sampled
 *    and projected on the other surface, those whose box reaches the other face only. A patch
 *    lies on it (every sample facing the other face within the tolerance), or not, or does not
 *    face it. The faces coincide partially when some patches lie on it and others not. Every face
 *    of the other shape is analysed before the carrier is split, once, along all the borders.
 * 2. Split: the knot lines between a patch on the other surface and one off it are the border of
 *    the common region; the carrier is split along them (iso lines of its surface, their pcurves
 *    segments). The piece on the other surface then coincides with it (the coincident faces of the
 *    boolean operation imprint each other), the pieces off it are intersected with it as usual.
 * 3. Crossings: over a patch off the other surface where the gap to it keeps one sign, the
 *    surfaces do not cross; the pieces off it are intersected with it on the patches where the
 *    gap changes sign only. Against a plane, the signs of the distances of the Bézier poles of
 *    the patch to it (the Bernstein coefficients of the gap) certify the sign; against another
 *    surface, the samples decide. Near the border, the gap stays under the tolerance over a band, where
 *    a marching would find tangencies everywhere; the surfaces touch there along the border only,
 *    an edge of both after the imprint.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/pcurve.h>
#include <gbs-brep/section.h>
#include <gbs-brep/split.h>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/intersect/surface_pieces.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    namespace partial_detail
    {
        /// A segment of a knot line on the carrier: v = `fixed` for u in [lo, hi] (`along_u`), or
        /// u = `fixed` for v in [lo, hi].
        template <std::floating_point T>
        struct Border
        {
            bool along_u{};
            T fixed{}, lo{}, hi{};
            bool operator==(const Border &) const = default;
        };

        enum class Patch
        {
            None, ///< no sample faces the other face
            On,   ///< every sample facing the other face lies on its surface
            Off
        };

        /// A partial coincidence of the faces `carrier` (a B-spline) and `other`: the patches of the
        /// carrier between the breakpoints `bu`, `bv` (u fastest), and the border between those on
        /// the other surface and those off it.
        template <std::floating_point T>
        struct Partial
        {
            FaceId carrier, other;
            std::vector<T> bu, bv;
            std::vector<Patch> patches;
            std::vector<Border<T>> borders;
            std::vector<std::array<T, 4>> crossing; ///< the patches off the other surface where the gap changes sign (carrier parameters)
            /// the surfaces do not cross beyond the border
            [[nodiscard]] bool parts() const { return crossing.empty(); }
        };

        /// The distinct knots of a B-spline surface in u and v, nullopt for other surfaces.
        template <std::floating_point T>
        auto surface_knots(const Surface<T, 3> &s) -> std::optional<std::pair<std::vector<T>, std::vector<T>>>
        {
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(&s))
                return std::pair{b->knotsU(), b->knotsV()};
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                return std::pair{b->knotsU(), b->knotsV()};
            return std::nullopt;
        }

        /// The iso line of a B-spline surface along u at v = `fixed` (`along_u`), or along v at u = `fixed`.
        template <std::floating_point T>
        auto iso_line(const Surface<T, 3> &s, bool along_u, T fixed) -> std::shared_ptr<Curve<T, 3>>
        {
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(&s))
                if (auto c = detail::row_or_iso(*b, along_u, fixed))
                    return std::make_shared<std::remove_cvref_t<decltype(*c)>>(std::move(*c));
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                if (auto c = detail::row_or_iso(*b, along_u, fixed))
                    return std::make_shared<std::remove_cvref_t<decltype(*c)>>(std::move(*c));
            return nullptr;
        }

        /// The plane of a surface (a point and the unit normal): a primitive plane, or a B-spline
        /// whose poles all lie on the plane of its corners.
        template <std::floating_point T>
        auto plane_of(const Surface<T, 3> &s) -> std::optional<std::pair<point<T, 3>, point<T, 3>>>
        {
            if (const auto *pl = dynamic_cast<const Plane<T> *>(&s))
                return std::pair{pl->frame().o, pl->frame().z};
            auto flat = [](const auto &b) -> std::optional<std::pair<point<T, 3>, point<T, 3>>> {
                const auto [u1, u2, v1, v2] = b.bounds();
                const auto o = b.value(u1, v1);
                auto n = cross(b.value(u2, v2) - o, b.value(u1, v2) - b.value(u2, v1));
                const T nn = norm(n), size = norm(b.value(u2, v2) - o) + norm(b.value(u1, v2) - b.value(u2, v1));
                if (!(nn > T(1e-12) * size * size))
                    return std::nullopt;
                n = n / nn;
                for (const auto &p : b.poles())
                {
                    point<T, 3> x;
                    for (int i = 0; i < 3; ++i)
                        x[i] = p.size() == 4 ? p[i] / p[3] : p[i];
                    if (std::abs((x - o) * n) > T(1e-12) * size)
                        return std::nullopt;
                }
                return std::pair{o, n};
            };
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(&s))
                return flat(*b);
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(&s))
                return flat(*b);
            return std::nullopt;
        }

        /// Breakpoints of [lo, hi]: its bounds and the knots strictly inside.
        template <std::floating_point T>
        auto breaks(const std::vector<T> &knots, T lo, T hi) -> std::vector<T>
        {
            const T eps = T(1e-9) * (hi - lo);
            std::vector<T> b{lo};
            for (T k : knots)
                if (k > lo + eps && k < hi - eps)
                    b.push_back(k);
            b.push_back(hi);
            return b;
        }

        /// The partial coincidence of the B-spline face `carrier` with the face `other`, if any.
        template <std::floating_point T>
        auto partial_of(const Model<T> &m, FaceId carrier, FaceId other, T tol) -> std::optional<Partial<T>>
        {
            const auto &sc = *m.face(carrier).surface, &so = *m.face(other).surface;
            const auto knots = surface_knots(sc);
            if (!knots)
                return std::nullopt;
            const auto box = classify_detail::face_uv_box(m, carrier), box_o = classify_detail::face_uv_box(m, other);
            Partial<T> res{carrier, other, breaks(knots->first, box[0], box[1]), breaks(knots->second, box[2], box[3]), {}, {}, {}};
            const auto &bu = res.bu, &bv = res.bv;
            if (bu.size() == 2 && bv.size() == 2)
                return std::nullopt; // one patch: coincident all over or not at all
            const std::size_t nu = bu.size() - 1, nv = bv.size() - 1;
            auto cell = [&](T x, const std::vector<T> &b) {
                const auto it = std::ranges::upper_bound(b, x);
                return static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(it - b.begin() - 1, 0, static_cast<std::ptrdiff_t>(b.size()) - 2));
            };
            // the Bézier patches of the carrier: their boxes skip those far from the other face;
            // their poles give the sign of the gap to a plane
            const auto reach = section_detail::face_box(m, other, tol);
            const auto plane = plane_of(so);
            std::vector<bool> near(nu * nv, true);
            std::vector<std::optional<std::array<T, 2>>> range(nu * nv); // min, max signed distance of the poles to the plane
            SurfacePieces<T> pieces{sc};
            for (const auto &pc : pieces.pieces(box))
            {
                const auto c = pc.center();
                const std::size_t k = cell(c[1], bv) * nu + cell(c[0], bu);
                Box<T, 3> b;
                b.add(reach.min);
                b.add(reach.max);
                near[k] = pc.box.overlaps(b, tol);
                if (plane)
                {
                    std::array<T, 2> r{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest()};
                    for (std::size_t j = 0; j < pc.nv; ++j)
                        for (std::size_t i = 0; i < pc.nu; ++i)
                        {
                            const T d = (pc.pole(i, j) - plane->first) * plane->second;
                            r = {std::min(r[0], d), std::max(r[1], d)};
                        }
                    range[k] = r;
                }
            }
            const T eu = T(1e-7) * (box_o[1] - box_o[0]), ev = T(1e-7) * (box_o[3] - box_o[2]);
            res.patches.assign(nu * nv, Patch::None);
            constexpr int n = 6;
            for (std::size_t j = 0; j < nv; ++j)
                for (std::size_t i = 0; i < nu; ++i)
                {
                    const std::size_t k = j * nu + i;
                    if (!near[k])
                        continue;
                    bool facing = false, off = false;
                    int sign{};
                    bool one_sign = true;
                    for (int b = 0; b < n; ++b)
                        for (int a = 0; a < n; ++a)
                        {
                            const point<T, 2> uv{bu[i] + (bu[i + 1] - bu[i]) * (T(a) + T(0.5)) / T(n),
                                                 bv[j] + (bv[j + 1] - bv[j]) * (T(b) + T(0.5)) / T(n)};
                            if (classify(m, carrier, uv, tol) != Location::In)
                                continue;
                            const auto p = sc.value(uv[0], uv[1]);
                            const auto q = intersect_detail::project_on_surface(so, box_o, p);
                            if (q[0] <= box_o[0] + eu || q[0] >= box_o[1] - eu || q[1] <= box_o[2] + ev || q[1] >= box_o[3] - ev)
                                continue;
                            if (classify(m, other, point<T, 2>{q[0], q[1]}, tol) != Location::In)
                                continue;
                            facing = true;
                            if (q[2] <= tol)
                                continue;
                            off = true;
                            const auto nq = cross(so.value(q[0], q[1], 1, 0), so.value(q[0], q[1], 0, 1));
                            const int s = (p - so.value(q[0], q[1])) * nq > T(0) ? 1 : -1;
                            if (sign == 0)
                                sign = s;
                            one_sign = one_sign && s == sign;
                        }
                    res.patches[k] = !facing ? Patch::None : (off ? Patch::Off : Patch::On);
                    // the sign over the patch: certified by the poles against a plane
                    if (off && range[k])
                        one_sign = !((*range[k])[0] < -tol && (*range[k])[1] > tol);
                    if (off && !one_sign)
                        res.crossing.push_back({bu[i], bu[i + 1], bv[j], bv[j + 1]});
                }
            // a crossing patch along the border of the common region meets the other surface on
            // that border, already an edge of both: its box is drawn back from it, by the least
            // offset where the gap exceeds 10 tol all along (a crease: little; a tangency, whose gap
            // stays under the tolerance over a band: not at all, and its intersection is singular)
            auto gap_beyond = [&](bool along_u, T fixed, T lo, T hi) {
                for (int a = 0; a <= 8; ++a)
                {
                    const T t = lo + (hi - lo) * T(a) / T(8);
                    const auto p = along_u ? sc.value(t, fixed) : sc.value(fixed, t);
                    if (intersect_detail::project_on_surface(so, box_o, p)[2] <= T(10) * tol)
                        return false;
                }
                return true;
            };
            auto patch_at = [&](T u, T v) { return res.patches[cell(v, bv) * nu + cell(u, bu)]; };
            for (auto &c : res.crossing)
            {
                const T du = c[1] - c[0], dv = c[3] - c[2];
                const T mu = (c[0] + c[1]) / T(2), mv = (c[2] + c[3]) / T(2);
                // sides: u = c[0], u = c[1] (along v), v = c[2], v = c[3] (along u), inwards
                const std::array<std::tuple<int, bool, T, T>, 4> sides{std::tuple{0, false, T(1), du}, std::tuple{1, false, T(-1), du},
                                                                       std::tuple{2, true, T(1), dv}, std::tuple{3, true, T(-1), dv}};
                const auto whole = c;
                for (const auto &[side, along_u, inward, width] : sides)
                {
                    const T at = whole[side];
                    const T eps = T(1e-9) * width;
                    const bool on_border = along_u ? (at - eps > box[2] && at + eps < box[3] && patch_at(mu, at - inward * eps) == Patch::On)
                                                   : (at - eps > box[0] && at + eps < box[1] && patch_at(at - inward * eps, mv) == Patch::On);
                    if (!on_border)
                        continue;
                    for (T d = T(1e-6) * width; d <= T(0.25) * width; d *= T(2))
                        if (along_u ? gap_beyond(true, at + inward * d, whole[0], whole[1]) : gap_beyond(false, at + inward * d, whole[2], whole[3]))
                        {
                            c[side] = at + inward * d;
                            break;
                        }
                }
            }
            const bool on = std::ranges::contains(res.patches, Patch::On);
            const bool off = std::ranges::contains(res.patches, Patch::Off);
            if (!on || !off)
                return std::nullopt;
            // the knot lines between a patch on the other surface and one off it, merged along
            auto add = [&](bool along_u, T fixed, T lo, T hi) {
                if (!res.borders.empty())
                {
                    auto &last = res.borders.back();
                    if (last.along_u == along_u && last.fixed == fixed && last.hi == lo)
                    {
                        last.hi = hi;
                        return;
                    }
                }
                res.borders.push_back({along_u, fixed, lo, hi});
            };
            auto border = [](Patch a, Patch b) { return (a == Patch::On && b == Patch::Off) || (a == Patch::Off && b == Patch::On); };
            for (std::size_t i = 0; i + 1 < nu; ++i) // u = bu[i + 1], along v
                for (std::size_t j = 0; j < nv; ++j)
                    if (border(res.patches[j * nu + i], res.patches[j * nu + i + 1]))
                        add(false, bu[i + 1], bv[j], bv[j + 1]);
            for (std::size_t j = 0; j + 1 < nv; ++j) // v = bv[j + 1], along u
                for (std::size_t i = 0; i < nu; ++i)
                    if (border(res.patches[j * nu + i], res.patches[(j + 1) * nu + i]))
                        add(true, bv[j + 1], bu[i], bu[i + 1]);
            if (res.borders.empty())
                return std::nullopt;
            return res;
        }

        /// The patches of `p` of kind `kind` inside the piece `f` of its carrier (tested at their centres).
        template <std::floating_point T>
        bool holds(const Model<T> &m, FaceId f, const Partial<T> &p, Patch kind, T tol)
        {
            const std::size_t nu = p.bu.size() - 1;
            for (std::size_t k = 0; k < p.patches.size(); ++k)
            {
                if (p.patches[k] != kind)
                    continue;
                const std::size_t i = k % nu, j = k / nu;
                const point<T, 2> c{(p.bu[i] + p.bu[i + 1]) / T(2), (p.bv[j] + p.bv[j + 1]) / T(2)};
                if (classify(m, f, c, tol) == Location::In)
                    return true;
            }
            return false;
        }

        /// Splits the carrier along the borders: edges along the iso lines, clipped to the face,
        /// with vertices from `pool`; its pieces.
        template <std::floating_point T>
        auto split_along(Model<T> &m, FaceId carrier, const std::vector<Border<T>> &borders, section_detail::VertexPool<T> &pool, T tol)
            -> BuildResult<std::vector<FaceId>>
        {
            const auto surface = m.face(carrier).surface;
            std::vector<FaceCut<T>> cuts;
            for (const auto &b : borders)
            {
                auto curve = iso_line(*surface, b.along_u, b.fixed);
                if (!curve)
                    return detail::fail(BuildErrc::NonGeneric, "partial coincidence: no iso line on the carrier");
                // the pcurve: the segment at the fixed parameter, on the parameter of the iso line
                const point<T, 2> a = b.along_u ? point<T, 2>{b.lo, b.fixed} : point<T, 2>{b.fixed, b.lo};
                const point<T, 2> c = b.along_u ? point<T, 2>{b.hi, b.fixed} : point<T, 2>{b.fixed, b.hi};
                auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{a, c}, std::vector<T>{b.lo, b.lo, b.hi, b.hi}, 1);
                std::vector<T> ts{b.lo, b.hi};
                std::vector<std::array<T, 2>> along;
                section_detail::boundary_crossings(m, carrier, *pc, b.lo, b.hi, ts, along);
                std::ranges::sort(ts);
                for (std::size_t k = 0; k + 1 < ts.size(); ++k)
                {
                    const T t0 = ts[k], t1 = ts[k + 1];
                    if (t1 - t0 <= T(1e-9) * (b.hi - b.lo))
                        continue;
                    if (classify(m, carrier, pc->value((t0 + t1) / T(2)), tol) != Location::In)
                        continue;
                    const auto v0 = pool.at(curve->value(t0)).first, v1 = pool.at(curve->value(t1)).first;
                    auto e = make_edge(m, curve, t0, t1, v0, v1, tol);
                    if (!e)
                        continue;
                    cuts.push_back({*e, pc});
                }
            }
            auto s = split_face(m, carrier, cuts, tol);
            if (!s)
                return std::unexpected(s.error());
            if (!s->dangling.empty())
                return detail::fail(BuildErrc::NonGeneric, "partial coincidence: the border does not split the face");
            return s->faces;
        }
    }

    /**
     * @brief A partial coincidence of the faces fa and fb: the first that is a B-spline face whose
     * patches lie on the surface of the other for some, off it for others; nullopt otherwise.
     */
    template <std::floating_point T>
    [[nodiscard]] auto partial_coincidence(const Model<T> &m, FaceId fa, FaceId fb, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> std::optional<partial_detail::Partial<T>>
    {
        if (auto p = partial_detail::partial_of(m, fa, fb, tol))
            return p;
        return partial_detail::partial_of(m, fb, fa, tol);
    }

    /**
     * @brief Splits the B-spline faces of the shapes `a` and `b` partially coincident with faces of
     * the other, before a section of the shapes: each carrier is analysed against every face of the
     * other shape, then split once along all the borders. Returns the restrictions of the section:
     * each piece off a face it partially coincides with is intersected with it on the crossing
     * patches only. The faces of `b` split after a restriction names them are replaced by their
     * pieces.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split_partial_coincidences(Model<T> &m, const ShapeId &a, const ShapeId &b, T tol)
        -> BuildResult<std::vector<SectionRestriction<T>>>
    {
        using namespace partial_detail;
        section_detail::VertexPool<T> pool{m, {}, tol};
        for (const auto &s : {a, b})
            for (auto v : explore<VertexId>(m, s))
                pool.ids.push_back(v);
        std::vector<SectionRestriction<T>> restricted;
        std::map<std::uint32_t, std::vector<FaceId>> split_into;
        for (const auto &[carriers, others] : {std::pair{a, b}, std::pair{b, a}})
            for (auto c : explore<FaceId>(m, carriers))
            {
                if (!m.alive(c) || !surface_knots(*m.face(c).surface))
                    continue;
                const auto box_c = section_detail::face_box(m, c, tol);
                std::vector<Partial<T>> found;
                for (auto o : explore<FaceId>(m, others))
                    if (section_detail::boxes_overlap(box_c, section_detail::face_box(m, o, tol)))
                        if (auto p = partial_of(m, c, o, tol))
                            found.push_back(std::move(*p));
                if (found.empty())
                    continue;
                std::vector<Border<T>> borders;
                for (const auto &p : found)
                    for (const auto &bd : p.borders)
                        if (std::ranges::find(borders, bd) == borders.end())
                            borders.push_back(bd);
                auto pieces = split_along(m, c, borders, pool, tol);
                if (!pieces)
                    return std::unexpected(pieces.error());
                split_into[c.index] = *pieces;
                for (const auto &p : found)
                    for (auto f : *pieces)
                        if (!holds(m, f, p, Patch::On, tol) && holds(m, f, p, Patch::Off, tol))
                            restricted.push_back({f, p.other, p.crossing});
            }
        // the faces split after a restriction named them: their pieces
        std::function<std::vector<FaceId>(FaceId)> current = [&](FaceId f) {
            const auto it = split_into.find(f.index);
            if (it == split_into.end() || m.alive(f))
                return std::vector<FaceId>{f};
            std::vector<FaceId> r;
            for (auto g : it->second)
                for (auto h : current(g))
                    r.push_back(h);
            return r;
        };
        std::vector<SectionRestriction<T>> res;
        for (const auto &r : restricted)
            for (auto f : current(r.face))
                for (auto o : current(r.other))
                    res.push_back({f, o, r.boxes});
        return res;
    }
}
