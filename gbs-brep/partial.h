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
 *    and projected on the other surface. A patch lies on it (every sample facing the other face
 *    within the tolerance), or not, or does not face it. The faces coincide partially when some
 *    patches lie on it and others not.
 * 2. Split: the knot lines between a patch on the other surface and one off it are the border of
 *    the common region; the carrier is split along them (iso lines of its surface, their pcurves
 *    segments). The piece on the other surface then coincides with it (the coincident faces of the
 *    boolean operation imprint each other), the pieces off it are intersected with it as usual.
 * 3. Crossings: over a patch off the other surface where the gap to it keeps one sign, the
 *    surfaces do not cross; the pieces off it are intersected with it on the patches where the
 *    gap changes sign only. Near the border, the gap stays under the tolerance over a band, where
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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
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
        };

        /// A partial coincidence of the faces `carrier` (a B-spline) and `other`.
        template <std::floating_point T>
        struct Partial
        {
            FaceId carrier, other;
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

        enum class Patch
        {
            None, ///< no sample faces the other face
            On,   ///< every sample facing the other face lies on its surface
            Off
        };

        /// The partial coincidence of the B-spline face `carrier` with the face `other`, if any.
        template <std::floating_point T>
        auto partial_of(const Model<T> &m, FaceId carrier, FaceId other, T tol) -> std::optional<Partial<T>>
        {
            const auto &sc = *m.face(carrier).surface, &so = *m.face(other).surface;
            const auto knots = surface_knots(sc);
            if (!knots)
                return std::nullopt;
            const auto box = classify_detail::face_uv_box(m, carrier), box_o = classify_detail::face_uv_box(m, other);
            const auto bu = breaks(knots->first, box[0], box[1]), bv = breaks(knots->second, box[2], box[3]);
            if (bu.size() == 2 && bv.size() == 2)
                return std::nullopt; // one patch: coincident all over or not at all
            const T eu = T(1e-7) * (box_o[1] - box_o[0]), ev = T(1e-7) * (box_o[3] - box_o[2]);
            const std::size_t nu = bu.size() - 1, nv = bv.size() - 1;
            std::vector<Patch> patches(nu * nv, Patch::None);
            std::vector<std::array<T, 4>> crossing;
            constexpr int n = 6;
            for (std::size_t j = 0; j < nv; ++j)
                for (std::size_t i = 0; i < nu; ++i)
                {
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
                    patches[j * nu + i] = !facing ? Patch::None : (off ? Patch::Off : Patch::On);
                    if (off && !one_sign)
                        crossing.push_back({bu[i], bu[i + 1], bv[j], bv[j + 1]});
                }
            const bool on = std::ranges::any_of(patches, [](Patch p) { return p == Patch::On; });
            const bool off = std::ranges::any_of(patches, [](Patch p) { return p == Patch::Off; });
            if (!on || !off)
                return std::nullopt;
            Partial<T> res{carrier, other, {}, std::move(crossing)};
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
                    if (border(patches[j * nu + i], patches[j * nu + i + 1]))
                        add(false, bu[i + 1], bv[j], bv[j + 1]);
            for (std::size_t j = 0; j + 1 < nv; ++j) // v = bv[j + 1], along u
                for (std::size_t i = 0; i < nu; ++i)
                    if (border(patches[j * nu + i], patches[(j + 1) * nu + i]))
                        add(true, bv[j + 1], bu[i], bu[i + 1]);
            if (res.borders.empty())
                return std::nullopt;
            return res;
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

    /// The pieces of a carrier split along the border of a partial coincidence.
    struct PartialSplit
    {
        std::vector<FaceId> on;  ///< pieces lying on the other surface
        std::vector<FaceId> off; ///< pieces off it
    };

    /**
     * @brief Splits the carrier of a partial coincidence along the border of the common region:
     * edges along the iso lines, clipped to the face, with vertices from `pool`. Returns the pieces
     * on the other surface and off it, or an error from the split.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split_partial(Model<T> &m, const partial_detail::Partial<T> &p, section_detail::VertexPool<T> &pool, T tol)
        -> BuildResult<PartialSplit>
    {
        const auto surface = m.face(p.carrier).surface;
        std::vector<FaceCut<T>> cuts;
        for (const auto &b : p.borders)
        {
            auto curve = partial_detail::iso_line(*surface, b.along_u, b.fixed);
            if (!curve)
                return detail::fail(BuildErrc::NonGeneric, "partial coincidence: no iso line on the carrier");
            // the pcurve: the segment at the fixed parameter, on the parameter of the iso line
            const point<T, 2> a = b.along_u ? point<T, 2>{b.lo, b.fixed} : point<T, 2>{b.fixed, b.lo};
            const point<T, 2> c = b.along_u ? point<T, 2>{b.hi, b.fixed} : point<T, 2>{b.fixed, b.hi};
            auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{a, c}, std::vector<T>{b.lo, b.lo, b.hi, b.hi}, 1);
            std::vector<T> ts{b.lo, b.hi};
            std::vector<std::array<T, 2>> along;
            section_detail::boundary_crossings(m, p.carrier, *pc, b.lo, b.hi, ts, along);
            std::ranges::sort(ts);
            for (std::size_t k = 0; k + 1 < ts.size(); ++k)
            {
                const T t0 = ts[k], t1 = ts[k + 1];
                if (t1 - t0 <= T(1e-9) * (b.hi - b.lo))
                    continue;
                if (classify(m, p.carrier, pc->value((t0 + t1) / T(2)), tol) != Location::In)
                    continue;
                const auto v0 = pool.at(curve->value(t0)).first, v1 = pool.at(curve->value(t1)).first;
                auto e = make_edge(m, curve, t0, t1, v0, v1, tol);
                if (!e)
                    continue;
                cuts.push_back({*e, pc});
            }
        }
        auto s = split_face(m, p.carrier, cuts, tol);
        if (!s)
            return std::unexpected(s.error());
        if (!s->dangling.empty())
            return detail::fail(BuildErrc::NonGeneric, "partial coincidence: the border does not split the face");
        PartialSplit res;
        const auto &so = *m.face(p.other).surface;
        const auto box_o = classify_detail::face_uv_box(m, p.other);
        for (auto f : s->faces)
        {
            const auto uv = interior_points(m, f, tol);
            if (uv.empty())
                return detail::fail(BuildErrc::NonGeneric, "partial coincidence: a piece without inside");
            const auto q = intersect_detail::project_on_surface(so, box_o, m.face(f).surface->value(uv.front()[0], uv.front()[1]));
            (q[2] <= tol ? res.on : res.off).push_back(f);
        }
        return res;
    }
}
