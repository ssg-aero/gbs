#pragma once
/**
 * @file fillet.h
 * @brief Geometry of a constant-radius fillet (`fillet_geometry`) and of a chamfer
 * (`chamfer_geometry`) along an edge of a solid: the surface of the blend and its two contact
 * curves on the faces of the edge. The solid is not modified; the topology (the faces cut along
 * the contacts, the ends) comes with `fillet` and `chamfer` (`blend.h`).
 *
 * **Fillet.** A ball of radius r rolling between the two faces of the edge, inside the material
 * along a convex edge, outside along a concave one (the side read on the outward normals and on
 * the directions of the faces away from the edge):
 * 1. its centre runs along the **centre line**, the intersection of the surfaces of the faces
 *    offset by r toward the ball (`offset_surface`: exact for the primitives; `intersect`: closed
 *    forms between primitives, marching otherwise), the branch nearest the edge, between the
 *    points nearest its vertices and continued past them by `extension` r (none for a closed
 *    edge: its whole loop, from its vertex, in the sense of the edge);
 * 2. it touches each face at the foot of the centre on the face, c - r n: the **contacts**, whose
 *    pcurves are those of the centre line on the offset surfaces (same parameters, but the shift
 *    of v of a cone);
 * 3. the **surface** is the envelope of the ball, between the contacts: around a straight centre
 *    line a cylinder of radius r, around a circular one a torus of minor radius r (whatever the
 *    faces); along any other, a rational B-spline: the arc of circle between the contacts, exact
 *    (quadratic, one or two segments), at stations of the centre line put back on the rolling
 *    ball by Gauss-Newton (the marched line only starts them: it loses its accuracy where the
 *    faces become nearly tangent), skinned by cubic Hermite (C1, local), and the centre line
 *    and the contacts interpolated by cubic Hermite on the same stations; each span between two
 *    stations halved until, at its middle and where the slope of a cubic Hermite strays most,
 *    the surface is within `approx_tol` / 2 of the exact arcs, its normal within
 *    `angular_tol` / 2, the centre line and the contacts within `approx_tol` / 2 of the exact
 *    ones (a span depends on its two stations only: only the spans out of tolerance are
 *    halved). A station is put where a contact crosses a knot of its face below C3 (the joins of
 *    the pieces of a section: the fillet is no smoother there). Where the faces become tangent
 *    at a vertex, the fillet vanishes there (its contacts meet) and ends there, not continued
 *    past it.
 *
 * A face closing across a seam that cannot be evaluated past it (a closed B-spline: the skinned
 * sections of a blade; an offset of one) is seen past its seam (`SeamlessSurface`: its parameters
 * brought back by the width of its bounds); a surface that repeats past its bounds already (a
 * primitive, a surface of revolution, an offset of a primitive) is left as it is. The loop of a
 * fillet round such a face starts where its contact crosses the seam, the contact then within
 * the bounds of the face (its poles within tol of them, in space, brought onto them). Round two
 * such faces, their seams must meet in the same section of the fillet (the centre line, marched
 * from seam to seam, closes only then).
 *
 * **Chamfer.** The ruled surface between two contacts, at the distances d1 on a face given and
 * d2 on the other from the edge, measured in the plane across the edge: the point of the face
 * at the distance d from the point of the edge in its normal plane, on the side of the face
 * (Gauss-Newton from the stations next to it; a pcurve interpolated through the stations,
 * doubled until within tol / 10). Between two lines in a plane it is a plane; between two
 * coaxial circles a plane, a cylinder or a cone (from its smaller circle: its semi-angle
 * positive, as STEP has it); a ruled B-spline otherwise.
 *
 * Every curve of the result (centre line or edge, contacts, their pcurves) has the same
 * parameter, on [a, b].
 *
 * Refused (`NonGeneric`): a radius beyond a radius of curvature of a face on the side of the ball
 * (its offset folds), a contact out of its face (out at the middle of the edge, or leaving it
 * across an edge other than those meeting the edge at its vertices: the radius or the distance
 * beyond the width of the face), faces tangent along the edge, an edge convex in places and
 * concave in others, no centre line near the edge, any other failure (a curve evaluated beyond
 * its bounds). `InvalidId`: an edge not of the solid, used by one face or by more than two, or
 * twice by one face (a seam); the face of d1 not a face of the edge. `InvalidBounds`: a radius or
 * a distance not beyond `tol`, tolerances not positive, an extension negative.
 *
 * Sources: Rossignac and Requicha, "Constant-radius blending in solid modelling" (Computers in
 * Mechanical Engineering, 1984); Choi and Ju, "Constant-radius blending in surface modelling"
 * (Computer-Aided Design 21(4), 1989); Vida, Martin and Várady, "A survey of blending methods
 * that use parametric surfaces" (Computer-Aided Design 26(5), 1994); Piegl and Tiller, *The NURBS
 * Book*, 2nd ed., § 7.5 (arcs of circle) and § 10.3 (skinning).
 */
#include <gbs-brep/classify.h>
#include <gbs-brep/offset.h>
#include <gbs-brep/radius_law.h>
#include <gbs-brep/thick.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscinterp.h>
#include <gbs/bsctools.h>
#include <gbs/curvecomposite.h>
#include <gbs/curveonsurface.h>
#include <gbs/curvetrimmed.h>
#include <gbs/intersect/projection.h>

#include <boost/math/tools/minima.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /**
     * A setback of a fillet at a vertex blend (stage 6, § 4.2; note `brep6_pr03_vertex_blends.md`):
     * the fillet of `edge` stops `distance` from `vertex` along the edge (arc length), farther than
     * where its contacts cross those of its neighbours (the automatic setback), for a softer corner.
     */
    template <std::floating_point T>
    struct Setback
    {
        VertexId vertex;
        EdgeId edge;
        T distance{};
    };

    /// The shape of the section of a fillet between its contacts (stage 7, `brep_palier7_blend_sections.md`, § 4).
    enum class SectionKind
    {
        Circular, ///< the arc of circle of the rolling ball (stages 4 to 6)
        Conic,    ///< the conic tangent to the faces at the contacts, of parameter rho
        G2        ///< the quartic of continuous curvature with the faces (stage 7, § 4.4)
    };

    /**
     * The section of a fillet (stage 7, § 4): the arc of circle, or the conic tangent to both
     * faces at the contacts of the ball of the radius, its shoulder point dividing the segment
     * from the middle of the contacts to the corner S of their tangents in the ratio rho (in
     * ]0, 1[: 0.5 a parabola, below an ellipse, above a hyperbola; the arc of circle where rho =
     * cos(a/2) / (1 + cos(a/2)), a the angle the arc sweeps: 0.414 at right angles). rho constant,
     * or a law along the chains of tangent edges (`rho_law`, its points on the edges filleted).
     * Or the section of continuous curvature (`g2()`, § 4.4): the quartic Bezier tangent to the
     * faces with their curvature across at its contacts, of the smoothest shape (the integral of
     * the square of the derivative of its curvature along it the least), its contacts spread so
     * that its middle is as far from S as that of the arc of radius r (1.305 r from the edge
     * between planes at right angles, where the ball touches at r): it takes as much material
     * away in the middle as the circular fillet of the same radius, its curvature without jump.
     */
    template <std::floating_point T>
    struct BlendSection
    {
        SectionKind kind{SectionKind::Circular};
        T rho{T(0.5)};
        std::optional<BlendLaw<T>> rho_law{};

        static auto circular() -> BlendSection { return {}; }
        static auto conic(T rho) -> BlendSection { return {SectionKind::Conic, rho, std::nullopt}; }
        static auto conic(const BlendLaw<T> &law) -> BlendSection { return {SectionKind::Conic, T(0.5), law}; }
        static auto g2() -> BlendSection { return {SectionKind::G2, T(0.5), std::nullopt}; }
    };

    template <std::floating_point T>
    struct FilletOptions
    {
        T tol{brep_default_tolerance<T>}; ///< tolerance of the result; the contacts of a chamfer within tol / 10 of the exact ones
        T approx_tol{T(1e-5)};            ///< approximated surfaces: deviation from the exact blend
        T angular_tol{T(1e-5)};           ///< approximated fillets: sine of the angle of their normal to that of the exact blend
        IntersectOptions<T> intersect{};  ///< intersection of the offset surfaces (its tol is raised to `tol` / 10 at least)
        T extension{T(1)};                ///< length of the blend past each end of an open edge, in radius (fillet) or largest distance (chamfer)
        std::vector<Setback<T>> setbacks{}; ///< setbacks of fillets at vertex blends (`fillet` only), each at least the automatic one
        BlendSection<T> section{};          ///< the section of a fillet (`fillet`, `fillet_geometry`; not a chamfer)
    };

    /// A contact curve of a blend on a face of the edge.
    template <std::floating_point T>
    struct BlendContact
    {
        FaceId face;
        std::shared_ptr<Surface<T, 3>> support; ///< the surface of the face, continued past the ends of the edge (same parameters)
        std::shared_ptr<Curve<T, 2>> pcurve;    ///< on `support`
        std::shared_ptr<Curve<T, 3>> curve;     ///< in space: `pcurve` on `support`
        std::shared_ptr<Curve<T, 2>> on_blend;  ///< on the surface of the blend (none for a fillet in pieces: `BlendPiece::on_blend` on each)
        T deviation{};                          ///< largest gap between `curve` and the surface of the blend (or its pieces) along `on_blend`
    };

    /// A piece of the surface of a fillet whose radius law has kinks: on [a, b], its contacts on it.
    template <std::floating_point T>
    struct BlendPiece
    {
        std::shared_ptr<Surface<T, 3>> surface; ///< u across (contact 0 at 0, contact 1 at 1), v the parameter of the blend
        T a{}, b{};
        std::array<std::shared_ptr<Curve<T, 2>>, 2> on_blend; ///< the contacts on `surface`
        T deviation{};                                         ///< from the envelope of the balls
    };

    /**
     * Where the slope of the radius of a fillet jumps (a point of a linear law), the ball there
     * between its two characteristic circles: the lune of the sphere that the rolling ball sweeps
     * from one to the other, between the pieces of the fillet on either side (stage 5, § 3.5).
     */
    template <std::floating_point T>
    struct BlendLune
    {
        T t{};                                                 ///< the parameter of the blend
        std::shared_ptr<SphericalSurface<T>> sphere;           ///< the ball
        std::array<std::shared_ptr<Curve<T, 3>>, 2> arcs;      ///< the characteristic circles before and after t, from contact 0 to contact 1, on [0, 1]: the ends of the pieces
        std::array<std::shared_ptr<Curve<T, 2>>, 2> on_sphere; ///< their pcurves on the sphere
        T deviation{};                                         ///< largest gap between an arc and its pcurve on the sphere (`sphere_pcurve`)
    };

    /**
     * Where the slope of the radius of a fillet drops (a peak of a linear law): the balls on
     * either side cover the lune of the ball there, and the envelopes on either side, continued
     * past it, meet along a crease from contact 0 to contact 1 (stage 5, § 3.5).
     */
    template <std::floating_point T>
    struct BlendCrease
    {
        T t{};                                                 ///< the parameter of the blend
        std::shared_ptr<Curve<T, 3>> curve;                    ///< on [0, 1], from contact 0 to contact 1
        std::array<std::shared_ptr<Curve<T, 2>>, 2> on_pieces; ///< its pcurves on the piece before and on the piece after (same parameter)
        T deviation{};                                         ///< largest gap between the pieces along it
    };

    /// The geometry of a fillet or of a chamfer along an edge: every curve on the parameters [a, b].
    template <std::floating_point T>
    struct EdgeBlend
    {
        EdgeId edge;
        bool convex{};                          ///< the blend removes material (adds it along a concave edge)
        bool closed{};                          ///< a closed edge: the blend runs its whole loop
        std::shared_ptr<Surface<T, 3>> surface; ///< cylinder, torus, plane, cone, or B-spline (rational for a fillet); none for a fillet in pieces
        bool exact{};                           ///< a primitive, or a B-spline built exactly (a chamfer ruled between B-spline contacts on planes), the contacts on it within tol; otherwise within `deviation` of the exact blend
        T deviation{};                          ///< largest distance of the approximated surface from the exact blend
        std::shared_ptr<Curve<T, 3>> spine;     ///< the centre line of a fillet, the edge (continued) of a chamfer
        T a{}, b{};                             ///< the blend, past the ends of the edge by the extension
        /// the parameters at the vertices v1 and v2 of the edge; t1 > t2: the spine runs the other
        /// way. A closed edge: the ends of its loop, a and b (t1 = b, t2 = a where the spine runs
        /// the other way), from its vertex or from where a contact crosses the seam of a face
        T t1{}, t2{};
        bool variable{}; ///< a fillet whose radius follows a law (where the arcs of two such fillets differ at a vertex of a chain: a lune or a crease)
        bool circular{true}; ///< a fillet whose sections are arcs of circle (a sphere at a corner possible); false: a conic or G2 section
        std::array<BlendContact<T>, 2> contacts; ///< on the faces of the two uses of the edge (a chamfer: the face of d1 first)
        /// a fillet whose radius law has kinks: its surface in pieces (`surface` null), a lune of
        /// the ball between two pieces (the contacts and the spine run on across them)
        std::vector<BlendPiece<T>> pieces;
        std::vector<BlendLune<T>> lunes;
        std::vector<BlendCrease<T>> creases; ///< between pieces where the slope of the radius drops
    };

    namespace fillet_detail
    {
        template <std::floating_point T>
        [[noreturn]] void refuse(BuildErrc c, const std::string &op, const std::string &why, std::vector<ShapeId> shapes = {})
        {
            throw BuildFailure(BuildError{c, op + ": " + why, std::move(shapes)});
        }

        template <std::floating_point T>
        point<T, 3> unit(const point<T, 3> &v)
        {
            return v / norm(v);
        }

        /// Positive and finite tolerances, a finite extension not negative.
        template <std::floating_point T>
        void check_options(const FilletOptions<T> &o, const std::string &op, EdgeId eid)
        {
            auto positive = [](T x) { return x > T(0) && std::isfinite(x); };
            if (!positive(o.tol) || !positive(o.approx_tol) || !positive(o.angular_tol))
                refuse<T>(BuildErrc::InvalidBounds, op, "the tolerances must be positive and finite", {eid});
            if (!(o.extension >= T(0)) || !std::isfinite(o.extension))
                refuse<T>(BuildErrc::InvalidBounds, op, "the extension must be finite and not negative", {eid});
            for (const auto &sb : o.setbacks)
                if (!positive(sb.distance))
                    refuse<T>(BuildErrc::InvalidBounds, op, "a setback must be positive and finite", {sb.edge});
        }

        /// A face of the edge.
        template <std::floating_point T>
        struct Side
        {
            FaceId face;
            Orientation orient{};                   ///< of the face in its shell
            std::shared_ptr<Surface<T, 3>> surface; ///< of the face
            std::shared_ptr<Curve<T, 2>> pcurve;    ///< of the edge on the face
            bool forward{};                          ///< the co-edge runs along the edge
        };

        template <std::floating_point T>
        struct EdgeSetup
        {
            std::array<Side<T>, 2> sides;
            bool convex{};
            bool closed{};
            T min_angle{}; ///< smallest angle of the wedge between the faces where the blend goes
            T max_angle{}; ///< its largest
        };

        /**
         * The faces of the edge in the solid and the side of the blend: convex when the direction
         * of the first face away from the edge, b1, goes below the second face (b1 . n2 < 0, n2
         * its outward normal), at 9 points; the angle of the wedge between b1 and b2 (its smallest
         * and its largest).
         */
        template <std::floating_point T>
        auto setup(const Model<T> &m, SolidId sid, EdgeId eid, const std::string &op) -> EdgeSetup<T>
        {
            if (!m.alive(sid) || !m.alive(eid))
                refuse<T>(BuildErrc::InvalidId, op, "dead or invalid solid or edge", {sid, eid});
            const auto &so = m.solid(sid);
            std::vector<ShellId> shells{so.outer};
            shells.insert(shells.end(), so.voids.begin(), so.voids.end());
            std::optional<std::vector<EdgeUse>> uses;
            std::map<std::uint32_t, Orientation> face_orient;
            for (auto sh : shells)
            {
                for (const auto &fu : m.shell(sh).faces)
                    face_orient[fu.face.index] = fu.orient;
                auto all = shell_edge_uses(m, sh);
                if (const auto it = all.find(eid); it != all.end())
                    uses = it->second;
            }
            const auto &ed = m.edge(eid);
            if (!uses || ed.degenerate)
                refuse<T>(BuildErrc::InvalidId, op, "not an edge of the solid", {eid});
            if (uses->size() != 2)
                refuse<T>(BuildErrc::InvalidId, op, "an edge not shared by two faces", {eid});
            if ((*uses)[0].face == (*uses)[1].face)
                refuse<T>(BuildErrc::InvalidId, op, "a seam edge (used twice by one face)", {eid});
            EdgeSetup<T> st;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &u = (*uses)[k];
                const auto &ce = m.wire(u.coedge.wire).coedges[u.coedge.index];
                if (!ce.pcurve)
                    refuse<T>(BuildErrc::MissingPCurve, op, "an edge without pcurve", {eid});
                st.sides[k] = Side<T>{u.face, face_orient.at(u.face.index), m.face(u.face).surface, ce.pcurve, ce.orient == Orientation::Forward};
            }
            st.closed = ed.v1 == ed.v2;
            // at 9 points: the outward normals and the directions of the faces away from the edge
            int convex = 0, concave = 0;
            st.min_angle = std::numbers::pi_v<T>;
            for (int i = 0; i < 9; ++i)
            {
                const T t = ed.u1 + (ed.u2 - ed.u1) * (T(i) + T(0.5)) / T(9);
                const auto te = unit(ed.curve->value(t, 1));
                std::array<point<T, 3>, 2> n, b;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto &sd = st.sides[k];
                    const auto uv = sd.pcurve->value(t);
                    const auto su = sd.surface->value(uv[0], uv[1], 1, 0), sv = sd.surface->value(uv[0], uv[1], 0, 1);
                    const auto nn = normal_direction(*sd.surface, uv[0], uv[1]);
                    n[k] = sd.orient == Orientation::Forward ? nn : T(-1) * nn;
                    // the face is on the left of its co-edge in (u, v) (outer wires counter-clockwise, holes clockwise)
                    auto d = sd.pcurve->value(t, 1);
                    if (!sd.forward)
                        d = T(-1) * d;
                    auto bk = (-d[1]) * su + d[0] * sv;
                    bk = bk - (bk * te) * te;
                    b[k] = unit(bk);
                }
                const T c = b[0] * n[1];
                if (!(std::abs(c) > T(1e-6)))
                    refuse<T>(BuildErrc::NonGeneric, op, "the faces of the edge are tangent along it", {eid});
                (c < T(0) ? convex : concave) += 1;
                const T wedge = std::acos(std::clamp(b[0] * b[1], T(-1), T(1)));
                st.min_angle = std::min(st.min_angle, wedge);
                st.max_angle = std::max(st.max_angle, wedge);
            }
            if (convex > 0 && concave > 0)
                refuse<T>(BuildErrc::NonGeneric, op, "an edge convex in places and concave in others", {eid});
            st.convex = convex > 0;
            return st;
        }

        /// The surface of a face continued around the edge, and the box where a contact may lie.
        template <std::floating_point T>
        struct Support
        {
            std::shared_ptr<Surface<T, 3>> surface;
            std::array<T, 4> box{};
        };

        /**
         * The surface of the side continued by `reach` past the face (`extended_surface`), on the
         * band of its parameters within `reach` of the edge.
         */
        template <std::floating_point T>
        auto support(const Model<T> &m, const Side<T> &sd, const Edge<T> &ed, T reach, T tol) -> Support<T>
        {
            const auto fbox = classify_detail::face_uv_box(m, sd.face);
            const auto ex = thick_detail::extended_surface(sd.surface, fbox, reach, tol);
            // the box of the pcurve, and the smallest speeds along it
            std::array<T, 4> band{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(),
                                  -std::numeric_limits<T>::infinity()};
            std::array<T, 2> speed{std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
            T fastest{};
            std::vector<std::array<T, 2>> sp;
            for (int i = 0; i <= 32; ++i)
            {
                const auto uv = sd.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(32));
                band = {std::min(band[0], uv[0]), std::max(band[1], uv[0]), std::min(band[2], uv[1]), std::max(band[3], uv[1])};
                const T a = norm(sd.surface->value(uv[0], uv[1], 1, 0)), b = norm(sd.surface->value(uv[0], uv[1], 0, 1));
                sp.push_back({a, b});
                fastest = std::max({fastest, a, b});
            }
            for (const auto &[a, b] : sp)
            {
                if (a > T(1e-9) * fastest)
                    speed[0] = std::min(speed[0], a);
                if (b > T(1e-9) * fastest)
                    speed[1] = std::min(speed[1], b);
            }
            std::array<T, 4> box = ex.box;
            for (int k = 0; k < 2; ++k)
                if (std::isfinite(speed[std::size_t(k)]))
                {
                    const T w = reach / speed[std::size_t(k)];
                    box[2 * k] = std::max(box[2 * k], band[2 * k] - w);
                    box[2 * k + 1] = std::min(box[2 * k + 1], band[2 * k + 1] + w);
                }
            if (!(box[1] > box[0]) || !(box[3] > box[2]))
                box = ex.box;
            return {ex.surface, box};
        }

        /**
         * A surface closing on itself across a seam without being periodic (a closed B-spline:
         * the skinned sections of a blade) evaluated past its seam: a parameter out of its bounds
         * brought into them by the period, the width of the bounds. The balls of a fillet along
         * an edge round such a face cross the seam: the centre line, the stations and the
         * contacts are computed on this view, their parameters running on continuously.
         */
        template <std::floating_point T>
        class SeamlessSurface : public Surface<T, 3>
        {
            std::shared_ptr<Surface<T, 3>> s_;
            std::array<T, 2> per_;

            auto in(T u, T v) const -> point<T, 2>
            {
                const auto b = s_->bounds();
                point<T, 2> w{u, v};
                for (std::size_t k = 0; k < 2; ++k)
                    if (const T p = per_[k]; p > T(0) && (w[k] < b[2 * k] || w[k] > b[2 * k + 1]))
                    {
                        w[k] = b[2 * k] + std::fmod(w[k] - b[2 * k], p);
                        if (w[k] < b[2 * k])
                            w[k] += p;
                    }
                return w;
            }

        public:
            SeamlessSurface(std::shared_ptr<Surface<T, 3>> s, const std::array<T, 2> &periods) : s_{std::move(s)}, per_{periods} {}
            [[nodiscard]] auto basis() const -> const std::shared_ptr<Surface<T, 3>> & { return s_; }
            [[nodiscard]] auto periods() const -> std::array<T, 2> { return per_; }
            auto bounds() const -> std::array<T, 4> override { return s_->bounds(); }
            auto has_second_derivatives() const -> bool override { return s_->has_second_derivatives(); }
            auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
            {
                const auto w = in(u, v);
                return s_->value(w[0], w[1], du, dv);
            }
            auto derivatives(T u, T v, size_t nu, size_t nv, std::array<T, 3> *SKL) const -> void override
            {
                const auto w = in(u, v);
                s_->derivatives(w[0], w[1], nu, nv, SKL);
            }
        };

        /// The periods of a primitive or of a seamless view (0: a parameter that does not close).
        template <std::floating_point T>
        auto view_periods(const Surface<T, 3> &s) -> std::array<T, 2>
        {
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(&s))
                return e->periods();
            if (const auto *e = dynamic_cast<const SeamlessSurface<T> *>(&s))
                return e->periods();
            return {T(0), T(0)};
        }

        /// The periods of the parameters of a surface (0: a parameter that does not close): those
        /// of a primitive, of a seamless view, else the widths of the bounds of a surface closing
        /// across a seam within tol (`surface_closure`).
        template <std::floating_point T>
        auto periods_of(const Surface<T, 3> &s, T tol) -> std::array<T, 2>
        {
            if (dynamic_cast<const ElementarySurface<T> *>(&s) || dynamic_cast<const SeamlessSurface<T> *>(&s))
                return view_periods(s);
            const auto cl = surface_closure(s, tol);
            const auto b = s.bounds();
            return {cl.closed_u ? b[1] - b[0] : T(0), cl.closed_v ? b[3] - b[2] : T(0)};
        }

        /**
         * The periods `per` of s along which it does not repeat past its bounds (a B-spline throws
         * there, an offset of one too): those a seamless view must bring back. A surface that
         * repeats already (a surface of revolution, an offset of a primitive) keeps its own
         * parameters, and what is exact about it (`offset_surface` sees an offset of a primitive).
         * Three points a quarter of a period in, against the same points a period further.
         */
        template <std::floating_point T>
        auto unrepeated(const Surface<T, 3> &s, std::array<T, 2> per, T tol) -> std::array<T, 2>
        {
            const auto b = s.bounds();
            for (std::size_t k = 0; k < 2; ++k)
                if (per[k] > T(0))
                {
                    bool repeats = true;
                    for (int i = 1; i <= 3 && repeats; ++i)
                    {
                        point<T, 2> in{};
                        in[k] = b[2 * k] + per[k] * T(i) / T(4);
                        in[1 - k] = b[2 * (1 - k)] + (b[2 * (1 - k) + 1] - b[2 * (1 - k)]) * T(i) / T(4);
                        auto out = in;
                        out[k] += per[k];
                        try
                        {
                            repeats = norm(s.value(out[0], out[1]) - s.value(in[0], in[1])) <= tol;
                        }
                        catch (const std::exception &)
                        {
                            repeats = false;
                        }
                    }
                    if (repeats)
                        per[k] = T(0);
                }
            return per;
        }

        /**
         * The parameters of the pcurve pc on [a, b] where its coordinate d crosses x + n p, n an
         * integer (the seams of a surface closing with period p; x alone where p is 0): 256
         * points, each crossing bisected, in increasing order.
         */
        template <std::floating_point T>
        auto crossings(const Curve<T, 2> &pc, T a, T b, std::size_t d, T x, T p) -> std::vector<T>
        {
            auto f = [&](T t) { return p > T(0) ? (pc.value(t)[d] - x) / p : pc.value(t)[d] - x; };
            auto level = [&](T y) { return p > T(0) ? std::floor(y) : (y < T(0) ? T(-1) : T(0)); };
            std::vector<T> r;
            T t0 = a, f0 = f(a);
            if (f0 == (p > T(0) ? std::round(f0) : T(0)))
                r.push_back(a);
            for (int i = 1; i <= 256; ++i)
            {
                const T t1 = i == 256 ? b : a + (b - a) * T(i) / T(256), f1 = f(t1);
                if (level(f0) != level(f1) && !(r.size() && r.back() == t0))
                {
                    const T n = p > T(0) ? std::max(level(f0), level(f1)) : T(0);
                    const bool up = f1 > f0;
                    T lo = t0, hi = t1;
                    for (int j = 0; j < 80 && hi - lo > T(1e-15) * (T(1) + std::abs(a) + std::abs(b)); ++j)
                    {
                        const T mid = (lo + hi) / T(2);
                        ((f(mid) < n) == up ? lo : hi) = mid;
                    }
                    r.push_back((lo + hi) / T(2));
                }
                t0 = t1, f0 = f1;
            }
            return r;
        }

        /// The interior knots of the flat knots k of degree p where the B-spline is less than C^c
        /// (a multiplicity m with p - m < c).
        template <std::floating_point T>
        auto weak_knots(const std::vector<T> &k, std::size_t p, std::size_t c) -> std::vector<T>
        {
            std::vector<T> r;
            for (std::size_t i = p + 1; i + p + 1 < k.size();)
            {
                std::size_t j = i;
                while (j + p + 1 < k.size() && k[j] == k[i])
                    ++j;
                if (p < j - i + c)
                    r.push_back(k[i]);
                i = j;
            }
            return r;
        }

        /// The interior knots of a B-spline surface (or of its seamless view) in u (d = 0) or v
        /// (d = 1) where it is less than C3: a fillet along it is no longer smooth across them.
        template <std::floating_point T>
        auto weak_knots(const Surface<T, 3> &s, std::size_t d) -> std::vector<T>
        {
            const Surface<T, 3> *base = &s;
            if (const auto *w = dynamic_cast<const SeamlessSurface<T> *>(base))
                base = w->basis().get();
            auto go = [&](const auto &bs) {
                return d == 0 ? weak_knots(bs.knotsFlatsU(), bs.degreeU(), 3) : weak_knots(bs.knotsFlatsV(), bs.degreeV(), 3);
            };
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(base))
                return go(*b);
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(base))
                return go(*b);
            return {};
        }

        /// The interior knots of the curve of an edge, a B-spline (none for another curve), where
        /// it is less than C^c.
        template <std::floating_point T>
        auto edge_weak_knots(const Curve<T, 3> &curve, std::size_t c) -> std::vector<T>
        {
            if (const auto *b = dynamic_cast<const BSCurve<T, 3> *>(&curve))
                return weak_knots(b->knotsFlats(), b->degree(), c);
            if (const auto *b = dynamic_cast<const BSCurveRational<T, 3> *>(&curve))
                return weak_knots(b->knotsFlats(), b->degree(), c);
            return {};
        }

        /**
         * The pcurve pc of a contact running once round a face that closes across a seam
         * without being periodic (from the seam to the seam), on its seamless view `sv`, moved by
         * whole periods into the bounds of the face, its poles within tol of them in space (a
         * parameter step times the fastest speed of the surface along that parameter on pc, and
         * 1e-9 of a period at least) brought onto them (the B-spline throws past its bounds; the
         * contact moves with its pcurve, the gap to the blend measured afterwards by `on_blend`):
         * nullptr where it does not fit in.
         */
        template <std::floating_point T>
        auto into_bounds(const std::shared_ptr<Curve<T, 2>> &pc, const SeamlessSurface<T> &sv, T tol) -> std::shared_ptr<Curve<T, 2>>
        {
            const auto b = sv.bounds();
            const auto per = sv.periods();
            const auto [a, z] = pc->bounds();
            point<T, 2> lo{std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()}, hi{-lo[0], -lo[1]};
            std::array<T, 2> speed{};
            for (int i = 0; i <= 64; ++i)
            {
                const auto q = pc->value(a + (z - a) * T(i) / T(64));
                for (std::size_t d = 0; d < 2; ++d)
                    lo[d] = std::min(lo[d], q[d]), hi[d] = std::max(hi[d], q[d]);
                speed[0] = std::max(speed[0], norm(sv.value(q[0], q[1], 1, 0)));
                speed[1] = std::max(speed[1], norm(sv.value(q[0], q[1], 0, 1)));
            }
            point<T, 2> shift{};
            for (std::size_t d = 0; d < 2; ++d)
                if (per[d] > T(0))
                    shift[d] = per[d] * std::round(((b[2 * d] + b[2 * d + 1]) - (lo[d] + hi[d])) / (T(2) * per[d]));
            auto moved = shifted_pcurve(pc, shift);
            auto *bs = dynamic_cast<BSCurve<T, 2> *>(moved.get());
            if (!bs)
                return nullptr;
            auto poles = bs->poles();
            for (auto &q : poles)
                for (std::size_t d = 0; d < 2; ++d)
                    if (per[d] > T(0))
                    {
                        const T eps = std::max(T(1e-9) * per[d], speed[d] > T(0) ? tol / speed[d] : T(0));
                        if (q[d] < b[2 * d] - eps || q[d] > b[2 * d + 1] + eps)
                            return nullptr;
                        q[d] = std::clamp(q[d], b[2 * d], b[2 * d + 1]);
                    }
            return std::make_shared<BSCurve<T, 2>>(poles, bs->knotsFlats(), bs->degree());
        }

        /// u and v of a pcurve point brought into the box of the face by its periods `per` (`periods_of`).
        template <std::floating_point T>
        auto into_face(const std::array<T, 2> &per, const std::array<T, 4> &fbox, point<T, 2> uv) -> point<T, 2>
        {
            for (std::size_t k = 0; k < 2; ++k)
                if (per[k] > T(0))
                {
                    const T lo = fbox[2 * k], eps = T(1e-9) * per[k];
                    while (uv[k] < lo - eps)
                        uv[k] += per[k];
                    while (uv[k] > fbox[2 * k + 1] + eps)
                        uv[k] -= per[k];
                }
            return uv;
        }

        /// The point of the edge e nearest p, {its parameter, its distance}: the nearest of 65 points, refined by Newton on [u1, u2].
        template <std::floating_point T>
        auto edge_point(const Edge<T> &e, const point<T, 3> &p) -> std::array<T, 2>
        {
            const std::array<T, 2> range{std::min(e.u1, e.u2), std::max(e.u1, e.u2)};
            T u0 = range[0], d0 = std::numeric_limits<T>::infinity();
            for (int i = 0; i <= 64; ++i)
            {
                const T u = range[0] + (range[1] - range[0]) * T(i) / T(64);
                if (const T d = norm(e.curve->value(u) - p); d < d0)
                    d0 = d, u0 = u;
            }
            const auto x = newton_point_curve(*e.curve, p, u0, range, T(1e-14) * (T(1) + range[1] - range[0]));
            return x[1] < d0 ? x : std::array<T, 2>{u0, d0};
        }

        /// The edge of the face f nearest the point p, but `eid` (`edge_point` on each).
        template <std::floating_point T>
        auto nearest_edge(const Model<T> &m, FaceId f, EdgeId eid, const point<T, 3> &p) -> std::optional<EdgeId>
        {
            std::optional<EdgeId> best;
            T gap = std::numeric_limits<T>::infinity();
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    const auto &e = m.edge(ce.edge);
                    if (ce.edge == eid || e.degenerate || !e.curve)
                        continue;
                    if (const T d = edge_point(e, p)[1]; d < gap)
                        gap = d, best = ce.edge;
                }
            return best;
        }

        /**
         * The contact stays on its face along the edge: in at the middle of the edge; from there
         * toward each vertex (t1 at v1, t2 at v2), it may leave the face only across an edge
         * meeting the edge at that vertex (past a corner more acute than the blend: the face at
         * that end cuts it, PR 6), and a closed edge nowhere. 16 points each way, the way out
         * found by bisection.
         */
        template <std::floating_point T>
        void check_on_face(const Model<T> &m, const Side<T> &sd, const Curve<T, 2> &pc, const Edge<T> &ed, T t1, T t2, T tol, const std::string &op,
                           EdgeId eid)
        {
            const auto fbox = classify_detail::face_uv_box(m, sd.face);
            const auto per = periods_of(*sd.surface, tol);
            auto uv_at = [&](T t) { return into_face(per, fbox, pc.value(t)); };
            auto inside = [&](T t) { return classify(m, sd.face, uv_at(t), tol) != Location::Out; };
            auto off = [&] { refuse<T>(BuildErrc::NonGeneric, op, "a contact runs off its face (beyond the width of the face)", {eid, sd.face}); };
            const T tm = (t1 + t2) / T(2);
            if (!inside(tm))
                off();
            const bool closed = ed.v1 == ed.v2;
            for (int k = 0; k < 2; ++k)
            {
                const T te = k == 0 ? t1 : t2;
                const VertexId v = k == 0 ? ed.v1 : ed.v2;
                T in = tm;
                for (int i = 1; i < 16; ++i)
                {
                    const T t = tm + (te - tm) * T(i) / T(16);
                    if (inside(t))
                    {
                        in = t;
                        continue;
                    }
                    if (closed)
                        off();
                    T out = t;
                    for (int j = 0; j < 60 && std::abs(out - in) > T(1e-12) * (T(1) + std::abs(t1) + std::abs(t2)); ++j)
                    {
                        const T mid = (in + out) / T(2);
                        (inside(mid) ? in : out) = mid;
                    }
                    const auto uv = uv_at((in + out) / T(2));
                    const auto across = nearest_edge(m, sd.face, eid, sd.surface->value(uv[0], uv[1]));
                    if (!across || !(m.edge(*across).v1 == v || m.edge(*across).v2 == v))
                        off();
                    break;
                }
            }
        }

        /// The quadratic rational arc of circle from p1 to p2 around c (one or two segments, on [0, 1]).
        template <std::floating_point T>
        auto arc(const point<T, 3> &c, const point<T, 3> &p1, const point<T, 3> &p2, int segments) -> BSCurveRational<T, 3>
        {
            const auto e1 = p1 - c, e2 = p2 - c;
            const T rho = (norm(e1) + norm(e2)) / T(2);
            const auto u1 = e1 / norm(e1), u2 = e2 / norm(e2);
            const T phi = std::atan2(norm(cross(u1, u2)), u1 * u2);
            const auto s = u1 + u2;
            if (!(norm(s) > T(1e-6)))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: contacts opposite on the ball", {}});
            const auto bis = s / norm(s);
            if (segments == 1)
            {
                const T w = std::cos(phi / T(2));
                return BSCurveRational<T, 3>(points_vector<T, 3>{p1, c + (rho / w) * bis, p2}, std::vector<T>{0, 0, 0, 1, 1, 1}, std::vector<T>{1, w, 1}, 2);
            }
            const T w = std::cos(phi / T(4));
            const auto q1 = u1 + bis, q2 = bis + u2;
            return BSCurveRational<T, 3>(points_vector<T, 3>{p1, c + (rho / w) * (q1 / norm(q1)), c + rho * bis, c + (rho / w) * (q2 / norm(q2)), p2},
                                         std::vector<T>{0, 0, 0, T(0.5), T(0.5), 1, 1, 1}, std::vector<T>{1, w, 1, w, 1}, 2);
        }

        /**
         * The quadratic rational arc of circle from p1 to p2 around c that passes by `via`, a
         * point of the circle (two segments, on [0, 1], joined at its middle): whatever it sweeps,
         * a half turn and more (the section of a full round fillet between flanks that are
         * parallel or come nearer away from the middle face, `full_round.h`), where `arc` takes
         * the shorter arc and has none between opposite points. Refused past 300 degrees.
         */
        template <std::floating_point T>
        auto arc_through(const point<T, 3> &c, const point<T, 3> &p1, const point<T, 3> &p2, const point<T, 3> &via) -> BSCurveRational<T, 3>
        {
            const auto e1 = p1 - c, e2 = p2 - c;
            const T rho = (norm(e1) + norm(e2)) / T(2);
            const auto u1 = e1 / norm(e1), u2 = e2 / norm(e2);
            // the middle of the arc: across its chord, on the side of via
            auto d = u2 - u1;
            auto w = via - c;
            if (norm(d) > T(1e-9))
                w = w - ((w * d) / (d * d)) * d;
            if (!(norm(w) > T(1e-9) * rho))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: an arc through a point at its centre", {}});
            const auto bis = w / norm(w);
            const T half = std::atan2(norm(cross(u1, bis)), u1 * bis);
            if (!(half < T(150) * std::numbers::pi_v<T> / T(180)))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: an arc of more than 300 degrees between its contacts", {}});
            const T wt = std::cos(half / T(2));
            const auto q1 = u1 + bis, q2 = bis + u2;
            return BSCurveRational<T, 3>(points_vector<T, 3>{p1, c + (rho / wt) * (q1 / norm(q1)), c + rho * bis, c + (rho / wt) * (q2 / norm(q2)), p2},
                                         std::vector<T>{0, 0, 0, T(0.5), T(0.5), 1, 1, 1}, std::vector<T>{1, wt, 1, wt, 1}, 2);
        }

        /// The section of a fillet for `skinned`: the arc of circle, a conic of rho(t) (and where its law kinks), or the G2 quartic.
        template <std::floating_point T>
        struct SectionOf
        {
            bool conic{};
            bool g2{};
            bool constant{true}; ///< rho the same everywhere
            std::function<T(T)> rho;
            std::vector<T> kinks;
            auto weight(T t) const -> T
            {
                const T r = rho(t);
                return r / (T(1) - r);
            }
        };

        /**
         * The conic tangent to the faces at the contacts p1 and p2 of the station (stage 7, § 4.2):
         * the rational quadratic of poles p1, S, p2 and weights 1, w, 1 on [0, 1], S the corner
         * of the tangents of the arc of circle of the station (its middle pole), w = rho / (1 -
         * rho). Refused where the arc sweeps more than 170 degrees (S too far: faces nearly
         * tangent).
         */
        template <std::floating_point T>
        auto conic_section(const point<T, 3> &c, const point<T, 3> &p1, const point<T, 3> &p2, T w) -> BSCurveRational<T, 3>
        {
            const auto e1 = p1 - c, e2 = p2 - c;
            const T rr = (norm(e1) + norm(e2)) / T(2);
            const auto u1 = e1 / norm(e1), u2 = e2 / norm(e2);
            const T phi = std::atan2(norm(cross(u1, u2)), u1 * u2);
            if (!(phi < T(170) * std::numbers::pi_v<T> / T(180)))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: a conic section between faces nearly tangent (the arc sweeps more than 170 degrees)", {}});
            const auto s = u1 + u2;
            const auto S = c + (rr / std::cos(phi / T(2))) * (s / norm(s));
            return BSCurveRational<T, 3>(points_vector<T, 3>{p1, S, p2}, std::vector<T>{0, 0, 0, 1, 1, 1}, std::vector<T>{1, w, 1}, 2);
        }

        /// A segment of (u, v) from (u1, v1) at t = a to (u2, v2) at t = b.
        template <std::floating_point T>
        auto segment(const point<T, 2> &p, const point<T, 2> &q, T a, T b) -> std::shared_ptr<Curve<T, 2>>
        {
            return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{p, q}, std::vector<T>{a, a, b, b}, 1);
        }

        /// A station of the ball: its centre, its contacts and their parameters on the surfaces of the faces.
        template <std::floating_point T>
        struct Station
        {
            point<T, 3> c, p1, p2;
            point<T, 2> uv1, uv2;
            /// the centre of the arc between the contacts where it is not c: the characteristic
            /// circle of a ball whose radius varies, in the plane normal to the centre line
            std::optional<point<T, 3>> m;
            /// a G2 section: its quartic Bezier between the contacts p1 and p2, which are then
            /// its own (`g2_station`), no longer those of the ball
            std::optional<std::array<point<T, 3>, 5>> g2;
            /// a point the arc between the contacts passes by (the contact of the ball with the
            /// middle face of a full round fillet): the arc is `arc_through`, a half turn and more
            std::optional<point<T, 3>> via;
            auto centre() const -> point<T, 3> { return m ? *m : c; }
        };

        /// The rational B-spline of a fillet and the curves fitted on its stations.
        template <std::floating_point T>
        struct Skinned
        {
            std::shared_ptr<Surface<T, 3>> surface;
            std::shared_ptr<Curve<T, 3>> spine;
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pcurves; ///< of the contacts, on the surfaces of the faces
            T deviation{};
        };

        /// The quartic Bezier of poles P on [0, 1] (Bernstein), or its derivative of order k.
        template <std::floating_point T>
        auto quartic(const std::array<point<T, 3>, 5> &P, T t, int k = 0) -> point<T, 3>
        {
            std::array<point<T, 3>, 5> Q = P;
            std::size_t n = 4;
            for (int j = 0; j < k; ++j, --n)
                for (std::size_t i = 0; i < n; ++i)
                    Q[i] = T(n) * (Q[i + 1] - Q[i]);
            // de Casteljau on the n + 1 poles left
            for (std::size_t j = n; j > 0; --j)
                for (std::size_t i = 0; i < j; ++i)
                    Q[i] = (T(1) - t) * Q[i] + t * Q[i + 1];
            return Q[0];
        }

        /**
         * The fairness of a planar quartic (normal n of its plane): the integral along it of the
         * square of the derivative of its curvature in arc length, \int (dk/ds)^2 ds = \int
         * (dk/dt)^2 / |B'| dt, by Gauss-Legendre (8 points) on 16 equal intervals of [0, 1] (the
         * lambda of the least within 1e-15 of it at right angles; 1e-7 with 4 intervals).
         */
        template <std::floating_point T>
        auto fairness(const std::array<point<T, 3>, 5> &P, const point<T, 3> &n) -> T
        {
            static constexpr std::array<double, 4> x{0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
            static constexpr std::array<double, 4> w{0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
            T e{};
            for (int q = 0; q < 16; ++q)
                for (std::size_t i = 0; i < 4; ++i)
                    for (const T sgn : {T(-1), T(1)})
                    {
                        const T t = (T(q) + T(0.5) + sgn * T(x[i]) / T(2)) / T(16);
                        const auto d1 = quartic(P, t, 1), d2 = quartic(P, t, 2), d3 = quartic(P, t, 3);
                        const T s = norm(d1);
                        const T dk = (cross(d1, d3) * n) / (s * s * s) - T(3) * (cross(d1, d2) * n) * (d1 * d2) / (s * s * s * s * s);
                        e += T(w[i]) / T(32) * dk * dk / s;
                    }
            return e;
        }

        /// The curvature vector across, in the plane of normal n, of the section of the surface s
        /// by that plane at (u, v): its normal curvature along the tangent of the section, over
        /// the cosine of the angle of the normal of s to the plane; with that tangent, toward +n x
        /// the normal of s.
        template <std::floating_point T>
        auto section_curvature(const Surface<T, 3> &s, const point<T, 2> &uv, const point<T, 3> &n) -> std::pair<point<T, 3>, point<T, 3>>
        {
            const auto su = s.value(uv[0], uv[1], 1, 0), sv = s.value(uv[0], uv[1], 0, 1);
            const auto suu = s.value(uv[0], uv[1], 2, 0), suv = s.value(uv[0], uv[1], 1, 1), svv = s.value(uv[0], uv[1], 0, 2);
            const auto N = unit(cross(su, sv));
            const auto t = unit(cross(n, N));
            // t = a su + b sv (the first fundamental form)
            const T E = su * su, F = su * sv, G = sv * sv, det = E * G - F * F;
            const T a = (G * (t * su) - F * (t * sv)) / det, b = (E * (t * sv) - F * (t * su)) / det;
            const T kn = (a * a * suu + T(2) * a * b * suv + b * b * svv) * N;
            const auto nu = cross(n, t); // in the plane, across t
            const T c = nu * N;
            return {t, (kn / c) * nu};
        }

        /**
         * The G2 section of a station of the ball (stage 7, § 4.4, rule C), in the plane of the arc
         * of the ball (normal n, through its centre o and its contacts), S_b the corner of the
         * tangents of that arc at its contacts and L their distance to it:
         *
         * - its contacts P_k on the faces in that plane, at m L from S_b along each tangent of the
         *   arc (projected on it): on two planes, m L along the faces;
         * - its poles P_1, Q_1, Q_2, Q_3, P_2: Q_1 and Q_3 on the tangents of the faces in the plane
         *   at P_1 and P_2, at lambda of the way to their corner S; Q_2 at the heights over those
         *   tangents that give the quartic the curvature across of the faces (h = 4/3 k a^2, a =
         *   |P Q_1|; Q_2 = S between planes);
         * - lambda that of the smoothest quartic (`fairness` the least: Brent, then Newton on its
         *   derivative by 5 points, for a lambda accurate and smooth from station to station;
         *   along the secant, Newton from the lambda of the shape before);
         * - m by the secant so that its middle B(1/2) is as far from S_b as that of the arc
         *   (r / cos(psi / 2) - r, psi the angle the arc sweeps): exact at the first step between
         *   planes (the section scales with m), m = 1.305 at right angles.
         *
         * The station then has the contacts of the section (p1, p2, uv1, uv2) and its poles;
         * its centres stay those of the ball. An arc reduced to a point (faces tangent): the
         * contacts of the ball, the poles on the segment. Refused (`NonGeneric`, without shapes:
         * `with_edge` gives the edge): an arc of more than 170 degrees, contacts that cannot be
         * found on the faces, a section whose tangents do not meet, a secant that does not reach
         * the sag of the arc (within 1e-10).
         */
        template <std::floating_point T>
        auto g2_station(Station<T> st, const std::array<const Surface<T, 3> *, 2> &supports) -> Station<T>
        {
            const auto o = st.centre();
            const auto e1 = st.p1 - o, e2 = st.p2 - o;
            const T r = (norm(e1) + norm(e2)) / T(2);
            const T psi = std::atan2(norm(cross(e1, e2)), e1 * e2);
            if (!(psi < T(170) * std::numbers::pi_v<T> / T(180)))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: a G2 section between faces nearly tangent (the arc sweeps more than 170 degrees)", {}});
            if (!(psi > T(1e-7)))
            {
                st.g2 = std::array<point<T, 3>, 5>{};
                for (std::size_t i = 0; i < 5; ++i)
                    (*st.g2)[i] = st.p1 + (T(i) / T(4)) * (st.p2 - st.p1);
                return st;
            }
            const auto n = unit(cross(e1, e2));
            const auto bis = unit(unit(e1) + unit(e2));
            const auto Sb = o + (r / std::cos(psi / T(2))) * bis;
            const T sag_arc = r / std::cos(psi / T(2)) - r;
            const std::array<point<T, 3>, 2> pb{st.p1, st.p2};
            std::array<point<T, 3>, 2> d;
            std::array<T, 2> L;
            for (std::size_t k = 0; k < 2; ++k)
                L[k] = norm(Sb - pb[k]), d[k] = (Sb - pb[k]) / L[k];
            std::array<point<T, 2>, 2> uv{st.uv1, st.uv2};
            const T scale = T(1) + norm(o) + r;
            struct Shape
            {
                std::array<point<T, 3>, 5> P;
                std::array<point<T, 3>, 2> p;
                std::array<point<T, 2>, 2> uv;
                T sag{};
                T lambda{};
            };
            // (from: the lambda of the shape before in the secant, where Newton starts)
            auto shape = [&](T m, std::optional<T> from = std::nullopt) -> Shape {
                Shape x;
                // the contacts: on face k, in the plane, at m L from S_b along d_k (Newton on (u, v))
                for (std::size_t k = 0; k < 2; ++k)
                {
                    auto w = uv[k];
                    bool ok = false;
                    for (int it = 0; it < 50; ++it)
                    {
                        const auto P = supports[k]->value(w[0], w[1]);
                        const auto su = supports[k]->value(w[0], w[1], 1, 0), sv = supports[k]->value(w[0], w[1], 0, 1);
                        const T f0 = (P - o) * n, f1 = (P - Sb) * d[k] + m * L[k];
                        if (std::abs(f0) + std::abs(f1) <= T(1e-14) * scale)
                        {
                            ok = true;
                            break;
                        }
                        const T a = su * n, b = sv * n, c = su * d[k], e = sv * d[k];
                        const T det = a * e - b * c;
                        if (!(std::abs(det) > T(0)))
                            break;
                        w[0] -= (e * f0 - b * f1) / det;
                        w[1] -= (a * f1 - c * f0) / det;
                    }
                    if (!ok)
                    {
                        const auto P = supports[k]->value(w[0], w[1]);
                        if (!(std::abs((P - o) * n) + std::abs((P - Sb) * d[k] + m * L[k]) <= T(1e-10) * scale))
                            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: the contacts of a G2 section cannot be found on the faces", {}});
                    }
                    x.uv[k] = w;
                    x.p[k] = supports[k]->value(w[0], w[1]);
                }
                // the tangents in the plane (toward S_b) and the curvature vectors across
                std::array<point<T, 3>, 2> t, kv, nu;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    auto [tk, kk] = section_curvature(*supports[k], x.uv[k], n);
                    if (tk * d[k] < T(0))
                        tk = T(-1) * tk;
                    t[k] = tk, kv[k] = kk, nu[k] = cross(n, tk);
                }
                // the corner S of the tangents: P_1 + s t_1 = P_2 + u t_2 in the plane
                const T c12 = t[0] * t[1], den = T(1) - c12 * c12;
                if (!(den > T(1e-12)))
                    throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: the tangents of a G2 section do not meet", {}});
                const auto dp = x.p[1] - x.p[0];
                const T sa = (dp * t[0] - c12 * (dp * t[1])) / den;
                const auto S = x.p[0] + sa * t[0];
                auto poles = [&](T lam) {
                    std::array<point<T, 3>, 5> P;
                    P[0] = x.p[0], P[4] = x.p[1];
                    P[1] = x.p[0] + lam * (S - x.p[0]);
                    P[3] = x.p[1] + lam * (S - x.p[1]);
                    // Q_2 = S + y: y . nu_k = h_k = 4/3 a_k^2 (k_k . nu_k), y = alpha nu_1 + beta nu_2
                    std::array<T, 2> h;
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        const T a = norm(P[k == 0 ? 1 : 3] - x.p[k]);
                        h[k] = T(4) / T(3) * a * a * (kv[k] * nu[k]);
                    }
                    const T g = nu[0] * nu[1], dd = T(1) - g * g;
                    const T alpha = (h[0] - g * h[1]) / dd, beta = (h[1] - g * h[0]) / dd;
                    P[2] = S + alpha * nu[0] + beta * nu[1];
                    return P;
                };
                auto E = [&](T lam) { return fairness(poles(lam), n); };
                // Newton on the derivative of E by 5 points (smooth in the geometry and accurate,
                // unlike the argmin of a flat minimum, within sqrt(epsilon) only); none where E is
                // not convex where it starts
                auto newton = [&](T lam) -> std::optional<T> {
                    const T dl = T(1e-3);
                    for (int it = 0; it < 8; ++it)
                    {
                        const T em2 = E(lam - T(2) * dl), em = E(lam - dl), ec = E(lam), ep = E(lam + dl), ep2 = E(lam + T(2) * dl);
                        const T g = (em2 - T(8) * em + T(8) * ep - ep2) / (T(12) * dl);
                        const T gg = (-em2 + T(16) * em - T(30) * ec + T(16) * ep - ep2) / (T(12) * dl * dl);
                        if (!(gg > T(0)))
                            return it == 0 ? std::nullopt : std::optional<T>{lam};
                        const T step = g / gg;
                        lam = std::clamp(lam - step, T(0.02), T(0.98));
                        if (std::abs(step) < T(1e-13))
                            break;
                    }
                    return lam;
                };
                // (Brent first, unless the lambda of the shape before starts Newton)
                std::optional<T> lam = from ? newton(*from) : std::nullopt;
                if (!lam)
                {
                    const T l0 = boost::math::tools::brent_find_minima(E, T(0.02), T(0.98), std::numeric_limits<T>::digits / 2).first;
                    lam = newton(l0).value_or(l0);
                }
                x.lambda = *lam;
                x.P = poles(*lam);
                x.sag = norm(quartic(x.P, T(0.5)) - Sb);
                return x;
            };
            // the secant on m for the sag of the arc
            T m0 = T(1.3), m1 = T(1.31);
            Shape s0 = shape(m0), s1 = shape(m1, s0.lambda);
            for (int it = 0; it < 30; ++it)
            {
                const T f0 = s0.sag - sag_arc, f1 = s1.sag - sag_arc;
                if (std::abs(f1) <= T(1e-14) * scale || f1 == f0)
                    break;
                const T m2 = m1 - f1 * (m1 - m0) / (f1 - f0);
                if (!(m2 > T(0)) || !std::isfinite(m2))
                    break;
                const T before = s1.lambda;
                m0 = m1, s0 = std::move(s1);
                m1 = m2, s1 = shape(m1, before);
                uv = s1.uv;
            }
            // (not reached: refused, not a section of another sag)
            if (!(std::abs(s1.sag - sag_arc) <= T(1e-10) * scale))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet_geometry: no spread of the contacts of a G2 section gives it the sag of the arc", {}});
            st.p1 = s1.p[0], st.p2 = s1.p[1];
            st.uv1 = s1.uv[0], st.uv2 = s1.uv[1];
            st.g2 = s1.P;
            return st;
        }

        /// The quartic Bezier of a G2 section, as a rational B-spline (weights 1) on [0, 1].
        template <std::floating_point T>
        auto g2_curve(const std::array<point<T, 3>, 5> &P) -> BSCurveRational<T, 3>
        {
            return BSCurveRational<T, 3>(points_vector<T, 3>{P[0], P[1], P[2], P[3], P[4]}, std::vector<T>{0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, std::vector<T>(5, T(1)), 4);
        }

        template <std::floating_point T, std::size_t d>
        auto joined(const std::vector<BSCurve<T, d>> &cs) -> BSCurve<T, d>;

        /**
         * The rational B-spline of the arcs of circle around the centres between the contacts of
         * the stations of `station` on [a, b] (u across, from p1 at 0 to p2 at 1; v along): the
         * homogeneous poles of the arcs at the stations and their derivatives along the centre
         * line (central differences), interpolated by cubic Hermite (C1, the tangent planes of the
         * arcs at the stations); the centres and the parameters of the contacts (on `supports`,
         * continuous across a seam) interpolated by cubic Hermite at the same stations, from their
         * values and derivatives. The stations start at 8 spans, cut at `breaks` (where a contact
         * crosses a knot of its face below C3; one within 1e-3 (b - a) of a station is dropped);
         * each span depends on its two stations only, and those out of tolerance are halved (up to
         * 2048 spans) until, at the middle of each span and at 1/2 -+ 1/(2 sqrt 3) (where the
         * slope of a cubic Hermite strays most), the surface is within approx_tol / 2 of the exact
         * arcs, its normal within angular_tol / 2 of theirs (where the arc is not reduced to a
         * point), and the centre line and the contacts within approx_tol / 2 of the stations (half
         * the tolerances: a margin for the points between those checked). `deviation`: the largest
         * gap found at those points. The arcs are around `Station::centre()` (the characteristic
         * circle of a ball of varying radius). A span whose points cannot be evaluated (the pcurve
         * of a coarse span out of the support) is halved too, down to 1e-6 (b - a); a refusal of a
         * station (`BuildFailure`) goes on up. A G2 section (`g2_station`, each kept once
         * computed) creases where its contacts cross a knot of their face below C3 or the seam of
         * a face closing without being periodic (its curvature kinks there): stations there,
         * found on its own contacts, where the skin is C0 (the derivatives on either side).
         */
        template <std::floating_point T, typename F>
        auto skinned(F &&ball, const std::array<const Surface<T, 3> *, 2> &supports, T a, T b, T approx_tol, T angular_tol,
                     const std::vector<T> &breaks = {}, const SectionOf<T> &sec = {}) -> Skinned<T>
        {
            // a G2 section: the station of the ball turned into that of the section (its own
            // contacts), each kept (the nodes, the checks and the creases ask for some again)
            std::map<T, Station<T>> g2_done;
            auto station = [&](T s) -> Station<T> {
                if (!sec.g2)
                    return ball(s);
                if (const auto it = g2_done.find(s); it != g2_done.end())
                    return it->second;
                return g2_done.emplace(s, g2_station(ball(s), supports)).first->second;
            };
            const bool circular = !sec.conic && !sec.g2;
            // (the arcs in two segments past 120 degrees)
            T phi{};
            if (circular)
                for (int i = 0; i <= 16; ++i)
                {
                    const auto st = station(a + (b - a) * T(i) / T(16));
                    const auto e1 = st.p1 - st.centre(), e2 = st.p2 - st.centre();
                    phi = std::max(phi, std::atan2(norm(cross(e1, e2)), e1 * e2));
                }
            const int segments = circular && phi > T(2) * std::numbers::pi_v<T> / T(3) ? 2 : 1;
            // (a conic: one segment, its weight that of rho at the station; G2: the quartic of the station)
            auto section = [&](T s, const Station<T> &st) {
                return sec.g2      ? g2_curve(*st.g2)
                       : sec.conic ? conic_section(st.centre(), st.p1, st.p2, sec.weight(s))
                       : st.via    ? arc_through(st.centre(), st.p1, st.p2, *st.via)
                                   : arc(st.centre(), st.p1, st.p2, segments);
            };
            const auto first = section(a, station(a));
            const auto ku = first.knotsFlats();
            const std::size_t pu = first.degree();
            const std::size_t nu = ku.size() - pu - 1;
            std::array<std::array<T, 2>, 2> per{};
            for (std::size_t k = 0; k < 2; ++k)
                per[k] = view_periods(*supports[k]);
            // a G2 section takes the curvature of the faces at its contacts: where a contact crosses
            // a knot of its face below C3, or the seam of a face closing without being periodic,
            // the sections kink along the fillet (its surface creases). Those stations found on the
            // contacts of the section itself (the contacts of the ball, `breaks`, cross them
            // elsewhere): at 64 stations, then by bisection
            std::vector<T> creases;
            if (sec.g2)
            {
                constexpr int ns = 64;
                auto at = [&](int i) { return i == ns ? b : a + (b - a) * T(i) / T(ns); };
                std::vector<Station<T>> sts;
                for (int i = 0; i <= ns; ++i)
                    sts.push_back(station(at(i)));
                for (std::size_t k = 0; k < 2; ++k)
                    for (std::size_t d = 0; d < 2; ++d)
                    {
                        auto xs = weak_knots(*supports[k], d);
                        if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(supports[k]); sv && sv->periods()[d] > T(0))
                            xs.push_back(sv->bounds()[2 * d]);
                        const T p = per[k][d];
                        for (const T x : xs)
                        {
                            auto g = [&](const Station<T> &q) {
                                const T y = (k == 0 ? q.uv1 : q.uv2)[d] - x;
                                return p > T(0) ? y - p * std::round(y / p) : y;
                            };
                            for (int i = 0; i < ns; ++i)
                            {
                                T glo = g(sts[std::size_t(i)]);
                                const T ghi = g(sts[std::size_t(i + 1)]);
                                // (a jump of half a period: the parameter brought back, no crossing)
                                if (!(glo * ghi < T(0)) || (p > T(0) && std::abs(ghi - glo) > p / T(2)))
                                    continue;
                                T lo = at(i), hi = at(i + 1);
                                while (hi - lo > T(1e-12) * (b - a))
                                {
                                    const T mid = (lo + hi) / T(2), gm = g(station(mid));
                                    if (gm * glo > T(0))
                                        lo = mid, glo = gm;
                                    else
                                        hi = mid;
                                }
                                // (not within 4 hs of an end: the stencils of the derivatives there)
                                if (const T x0 = (lo + hi) / T(2); x0 > a + T(4e-6) * (b - a) && x0 < b - T(4e-6) * (b - a))
                                    creases.push_back(x0);
                            }
                        }
                    }
            }
            // a node: the station, the homogeneous poles of its arc, its centre and the parameters
            // of its contacts, with their derivatives along the centre line (at a crease, those
            // before it, and those after it apart)
            struct Node
            {
                Station<T> st;
                std::vector<constrType<T, 4, 2>> h;
                point<T, 3> dc{};
                std::array<point<T, 2>, 2> duv{};
                bool crease{};
                std::vector<point<T, 4>> hr;
                point<T, 3> dcr{};
                std::array<point<T, 2>, 2> duvr{};
            };
            const T hs = T(1e-6) * (b - a);
            auto node = [&](T s) {
                Node x{station(s), std::vector<constrType<T, 4, 2>>(nu)};
                x.crease = std::ranges::contains(creases, s);
                const auto ax = section(s, x.st);
                // the derivatives from s0 to s1 (central differences, one-sided at an end)
                auto slopes = [&](T s0, T s1, std::vector<point<T, 4>> &dh, point<T, 3> &dc, std::array<point<T, 2>, 2> &duv) {
                    const auto st0 = s0 == s ? x.st : station(s0), st1 = s1 == s ? x.st : station(s1);
                    const auto h0 = section(s0, st0), h1 = section(s1, st1);
                    dh.resize(nu);
                    for (std::size_t j = 0; j < nu; ++j)
                        dh[j] = (h1.poles()[j] - h0.poles()[j]) / (s1 - s0);
                    dc = (st1.c - st0.c) / (s1 - s0);
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        auto d = k == 0 ? st1.uv1 - st0.uv1 : st1.uv2 - st0.uv2;
                        for (std::size_t i = 0; i < 2; ++i)
                            if (per[k][i] > T(0))
                                d[i] -= per[k][i] * std::round(d[i] / per[k][i]);
                        duv[k] = d / (s1 - s0);
                    }
                };
                // (a G2 section at an end or a crease: one-sided of the second order, from s and the
                // stations hs and 2 hs away on the side dir)
                auto one_sided = [&](T dir, std::vector<point<T, 4>> &dh, point<T, 3> &dc, std::array<point<T, 2>, 2> &duv) {
                    const T t1 = s + dir * hs, t2 = s + T(2) * dir * hs;
                    const auto st1 = station(t1), st2 = station(t2);
                    const auto h1 = section(t1, st1), h2 = section(t2, st2);
                    const T c = dir / (T(2) * hs);
                    dh.resize(nu);
                    for (std::size_t j = 0; j < nu; ++j)
                        dh[j] = c * (T(4) * (h1.poles()[j] - ax.poles()[j]) - (h2.poles()[j] - ax.poles()[j]));
                    dc = c * (T(4) * (st1.c - x.st.c) - (st2.c - x.st.c));
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        const auto u0 = k == 0 ? x.st.uv1 : x.st.uv2;
                        auto d1 = (k == 0 ? st1.uv1 : st1.uv2) - u0, d2 = (k == 0 ? st2.uv1 : st2.uv2) - u0;
                        for (std::size_t i = 0; i < 2; ++i)
                            if (per[k][i] > T(0))
                            {
                                d1[i] -= per[k][i] * std::round(d1[i] / per[k][i]);
                                d2[i] -= per[k][i] * std::round(d2[i] / per[k][i]);
                            }
                        duv[k] = c * (T(4) * d1 - d2);
                    }
                };
                const T s0 = std::max(a, s - hs), s1 = std::min(b, s + hs);
                std::vector<point<T, 4>> dh;
                if (sec.g2 && (x.crease || s == a || s == b))
                {
                    one_sided(s == a ? T(1) : T(-1), dh, x.dc, x.duv);
                    if (x.crease)
                        one_sided(T(1), x.hr, x.dcr, x.duvr);
                }
                else
                    slopes(s0, s1, dh, x.dc, x.duv);
                for (std::size_t j = 0; j < nu; ++j)
                    x.h[j] = {ax.poles()[j], dh[j]};
                return x;
            };
            // 8 spans, cut at the creases, then at the breaks (where a contact crosses a knot of
            // its face below C3; a station of the 8 spans or a break within 1e-3 (b - a) of one
            // before dropped)
            std::vector<T> v{a, b};
            auto add = [&](T x, T near) {
                if (x > a && x < b && std::ranges::none_of(v, [&](T y) { return std::abs(y - x) < near; }))
                    v.push_back(x);
            };
            const T gap = T(1e-3) * (b - a);
            for (const T x : creases)
                add(x, T(1e-9) * (b - a));
            for (int i = 1; i < 8; ++i)
                add(a + (b - a) * T(i) / T(8), gap);
            for (const auto *xs : {&breaks, &sec.kinks})
                for (const T x : *xs)
                    add(x, gap);
            std::ranges::sort(v);
            std::vector<Node> N(v.size());
            for (std::size_t i = 0; i < v.size(); ++i)
                N[i] = node(v[i]);
            T worst_d{}, worst_ang{}; // the worst span out of the tolerances, relative to each
            for (;;)
            {
                const std::size_t n = v.size() - 1;
                std::vector<constrType<T, 3, 2>> centres(n + 1);
                std::array<std::vector<constrType<T, 2, 2>>, 2> uv{std::vector<constrType<T, 2, 2>>(n + 1), std::vector<constrType<T, 2, 2>>(n + 1)};
                for (std::size_t i = 0; i <= n; ++i)
                {
                    centres[i] = {N[i].st.c, N[i].dc};
                    uv[0][i] = {N[i].st.uv1, N[i].duv[0]};
                    uv[1][i] = {N[i].st.uv2, N[i].duv[1]};
                }
                // the parameters of the contacts continuous from station to station (a periodic
                // parameter brought back by its period across the seam of its surface)
                for (std::size_t k = 0; k < 2; ++k)
                    for (std::size_t d = 0; d < 2; ++d)
                        if (const T p = per[k][d]; p > T(0))
                            for (std::size_t i = 1; i <= n; ++i)
                                uv[k][i][0][d] += p * std::round((uv[k][i - 1][0][d] - uv[k][i][0][d]) / p);
                // the cubic Hermite of the nodes: C1, but C0 at a crease (the pieces between the
                // creases on their own, from the derivatives after the crease at their start, joined)
                const bool creased = std::ranges::any_of(N, [](const Node &x) { return x.crease; });
                auto hermite = [&]<std::size_t dim>(const std::vector<constrType<T, dim, 2>> &row, auto &&after) -> BSCurve<T, dim> {
                    if (!creased)
                        return interpolate<T, dim, 2>(row, v);
                    std::vector<BSCurve<T, dim>> parts;
                    std::size_t i0 = 0;
                    for (std::size_t i = 1; i <= n; ++i)
                        if (i == n || N[i].crease)
                        {
                            std::vector<constrType<T, dim, 2>> piece(row.begin() + std::ptrdiff_t(i0), row.begin() + std::ptrdiff_t(i + 1));
                            if (N[i0].crease)
                                piece.front()[1] = after(i0);
                            parts.push_back(interpolate<T, dim, 2>(piece, std::vector<T>(v.begin() + std::ptrdiff_t(i0), v.begin() + std::ptrdiff_t(i + 1))));
                            i0 = i;
                        }
                    return joined(parts);
                };
                std::vector<BSCurve<T, 4>> along;
                for (std::size_t j = 0; j < nu; ++j)
                {
                    std::vector<constrType<T, 4, 2>> row(n + 1);
                    for (std::size_t i = 0; i <= n; ++i)
                        row[i] = N[i].h[j];
                    along.push_back(hermite(row, [&](std::size_t i) { return N[i].hr[j]; }));
                }
                const auto kv = along.front().knotsFlats();
                const std::size_t nv = along.front().poles().size();
                std::vector<std::array<T, 4>> poles(nu * nv);
                for (std::size_t i = 0; i < nv; ++i)
                    for (std::size_t j = 0; j < nu; ++j)
                        poles[j + i * nu] = along[j].poles()[i];
                Skinned<T> r;
                r.surface = std::make_shared<BSSurfaceRational<T, 3>>(poles, ku, kv, pu, std::size_t{3});
                r.spine = std::make_shared<BSCurve<T, 3>>(hermite(centres, [&](std::size_t i) { return N[i].dcr; }));
                for (std::size_t k = 0; k < 2; ++k)
                    r.pcurves[k] = std::make_shared<BSCurve<T, 2>>(hermite(uv[k], [&](std::size_t i) { return N[i].duvr[k]; }));
                // each span checked on its own (the interpolation is local: a span depends on its
                // two stations only), the spans out of the tolerances halved
                T dev{};
                std::vector<std::size_t> out;
                worst_d = worst_ang = T(0);
                // the piece between the creases (or the ends) round each span
                std::vector<T> piece_lo(n), piece_hi(n);
                for (std::size_t i = 0; i < n; ++i)
                    piece_lo[i] = i == 0 || N[i].crease ? v[i] : piece_lo[i - 1];
                for (std::size_t i = n; i-- > 0;)
                    piece_hi[i] = i + 1 == n || N[i + 1].crease ? v[i + 1] : piece_hi[i + 1];
                for (std::size_t i = 0; i < n; ++i)
                {
                    T d{}, ang{};
                    // (a pcurve of a coarse span strayed out of the support continued: that span
                    // halved, down to 1e-6 of [a, b]; a refusal of a station goes on up)
                    try
                    {
                        // the middle (where a cubic Hermite strays most) and the points where its slope
                        // does, 1/2 -+ 1/(2 sqrt 3) (the normals)
                        for (const T w : {T(0.5) - T(0.5) / std::numbers::sqrt3_v<T>, T(0.5), T(0.5) + T(0.5) / std::numbers::sqrt3_v<T>})
                        {
                            const T s = v[i] + (v[i + 1] - v[i]) * w;
                            const auto st = station(s);
                            const auto exact = section(s, st);
                            // a conic or G2: the normal of the exact surface from the sections on
                            // either side (not across a crease)
                            std::optional<BSCurveRational<T, 3>> em, ep;
                            if (!circular)
                            {
                                const T sm = std::max(piece_lo[i], s - hs), sp = std::min(piece_hi[i], s + hs);
                                em = section(sm, station(sm)), ep = section(sp, station(sp));
                            }
                            d = std::max(d, norm(r.spine->value(s) - st.c));
                            for (std::size_t k = 0; k < 2; ++k)
                            {
                                const auto q = r.pcurves[k]->value(s);
                                d = std::max(d, norm(supports[k]->value(q[0], q[1]) - (k == 0 ? st.p1 : st.p2)));
                            }
                            const auto o = st.centre();
                            const T arc_length = norm(st.p1 - o) * std::atan2(norm(cross(st.p1 - o, st.p2 - o)), (st.p1 - o) * (st.p2 - o));
                            for (int k = 0; k <= 8; ++k)
                            {
                                const T u = T(k) / T(8);
                                d = std::max(d, norm(r.surface->value(u, s) - exact.value(u)));
                                // the normal of the envelope, along the radius (undefined where the arc
                                // is reduced to a point: a fillet vanishing where its faces become tangent)
                                if (arc_length > T(1e3) * approx_tol)
                                {
                                    point<T, 3> nx;
                                    if (!circular)
                                        nx = unit(cross(exact.value(u, 1), ep->value(u) - em->value(u)));
                                    else
                                    {
                                        const auto e = exact.value(u) - st.c;
                                        nx = e / norm(e);
                                    }
                                    ang = std::max(ang, norm(cross(normal_direction(*r.surface, u, s), nx)));
                                }
                            }
                        }
                    }
                    catch (const BuildFailure &)
                    {
                        throw; // a refusal of a station (a ball of a varying radius): not a coarse span
                    }
                    catch (const std::exception &)
                    {
                        // a span already narrow that still cannot be evaluated: the contact does
                        // leave the support, not a coarse pcurve
                        if (!(v[i + 1] - v[i] > T(1e-6) * (b - a)))
                            throw;
                        d = std::numeric_limits<T>::infinity();
                    }
                    if (std::isfinite(d))
                        dev = std::max(dev, d);
                    // half the tolerances at three points: a margin for the points between them
                    if (!(d <= approx_tol / T(2)) || !(ang <= angular_tol / T(2)))
                    {
                        out.push_back(i);
                        worst_d = std::max(worst_d, d / approx_tol), worst_ang = std::max(worst_ang, ang / angular_tol);
                    }
                }
                if (out.empty())
                {
                    r.deviation = dev;
                    return r;
                }
                if (n + out.size() > 2048)
                    break;
                std::vector<T> v2;
                std::vector<Node> N2;
                for (std::size_t i = 0, o = 0; i <= n; ++i)
                {
                    v2.push_back(v[i]);
                    N2.push_back(std::move(N[i]));
                    if (o < out.size() && out[o] == i)
                    {
                        v2.push_back((v[i] + v[i + 1]) / T(2));
                        N2.push_back(node(v2.back()));
                        ++o;
                    }
                }
                v = std::move(v2), N = std::move(N2);
            }
            // the tolerance actually missed (relative to it, the worse one if both)
            throw BuildFailure(BuildError{BuildErrc::NonGeneric,
                                          worst_ang > worst_d ? "fillet_geometry: the normal of the surface cannot be approximated within angular_tol"
                                                              : "fillet_geometry: the surface cannot be approximated within approx_tol",
                                          {}});
        }

        /**
         * The blend surface as a primitive around a straight (cylinder) or circular (torus) spine,
         * on the box of the contacts; nullptr for another spine.
         */
        template <std::floating_point T>
        auto canal(const Curve<T, 3> &spine, const Curve<T, 3> &p1, const Curve<T, 3> &p2, T a, T b, T r) -> std::shared_ptr<ElementarySurface<T>>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            const T sm = (a + b) / T(2);
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(&spine))
            {
                const auto o = spine.value(sm);
                auto z = unit(l->vector());
                const auto e1 = p1.value(sm) - o;
                const auto x = unit(e1 - (e1 * z) * z);
                if ((p2.value(sm) - o) * cross(z, x) < T(0))
                    z = T(-1) * z;
                const Frame3<T> f{ax2<T, 3>{o, z, x}};
                T u1 = std::numeric_limits<T>::infinity(), u2 = -u1;
                for (int i = 0; i <= 32; ++i)
                {
                    const T s = a + (b - a) * T(i) / T(32);
                    for (const auto *p : {&p1, &p2})
                    {
                        const auto q = p->value(s) - o;
                        const T u = std::atan2(q * f.y, q * f.x);
                        u1 = std::min(u1, u), u2 = std::max(u2, u);
                    }
                }
                const T v1 = (spine.value(a) - o) * z, v2 = (spine.value(b) - o) * z;
                return std::make_shared<CylindricalSurface<T>>(r, f, std::array<T, 4>{u1, u2, std::min(v1, v2), std::max(v1, v2)});
            }
            if (const auto *e = dynamic_cast<const EllipseCurve<T> *>(&spine))
            {
                if (!(std::abs(e->radius1() - e->radius2()) <= T(1e-12) * e->radius1()))
                    return nullptr;
                const auto &f = e->frame();
                auto minor_angle = [&](const Curve<T, 3> &p, T s) {
                    const auto q = p.value(s) - spine.value(s);
                    return std::atan2(q * f.z, q * f.radial(s));
                };
                // the arc between the contacts shorter than a half turn, each contact continuous along
                const T w1 = minor_angle(p1, sm);
                T w2 = minor_angle(p2, sm);
                while (w2 > w1 + pi)
                    w2 -= T(2) * pi;
                while (w2 < w1 - pi)
                    w2 += T(2) * pi;
                T v1 = std::numeric_limits<T>::infinity(), v2 = -v1;
                for (int i = 0; i <= 32; ++i)
                {
                    const T s = a + (b - a) * T(i) / T(32);
                    for (const auto &[p, w] : {std::pair{&p1, w1}, std::pair{&p2, w2}})
                    {
                        T v = minor_angle(*p, s);
                        v += T(2) * pi * std::round((w - v) / (T(2) * pi));
                        v1 = std::min(v1, v), v2 = std::max(v2, v);
                    }
                }
                return std::make_shared<ToroidalSurface<T>>(e->radius1(), r, f, std::array<T, 4>{a, b, v1, v2});
            }
            return nullptr;
        }

        /**
         * The contacts' pcurves on the blend and the gap of each contact to its surface: on a
         * primitive fitted (`primitive_pcurve`); on an approximated surface its sides, u = 0 and
         * u = 1 across a fillet (`across_u`), v = 0 and v = 1 across a chamfer.
         */
        template <std::floating_point T>
        void on_blend(EdgeBlend<T> &r, T tol, bool across_u)
        {
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(r.surface.get()))
            {
                for (auto &c : r.contacts)
                {
                    const auto fit = intersect_detail::primitive_pcurve(*c.curve, r.a, r.b, *e, e->bounds(), tol);
                    c.on_blend = fit.pcurve;
                    c.deviation = fit.deviation;
                }
                return;
            }
            for (std::size_t k = 0; k < 2; ++k)
            {
                auto &c = r.contacts[k];
                const T w = T(k);
                c.on_blend = across_u ? segment<T>({w, r.a}, {w, r.b}, r.a, r.b) : segment<T>({r.a, w}, {r.b, w}, r.a, r.b);
                T dev{};
                for (int i = 0; i <= 256; ++i)
                {
                    const T s = r.a + (r.b - r.a) * T(i) / T(256);
                    const auto q = c.on_blend->value(s);
                    dev = std::max(dev, norm(r.surface->value(q[0], q[1]) - c.curve->value(s)));
                }
                c.deviation = dev;
            }
        }

        /**
         * The parameters on [a, b] where a contact crosses a knot of its face below C3 (the joins
         * of the pieces of a section: the fillet is no smoother there), found on the pcurves pc
         * of the contacts on the supports, then bisected on the stations themselves.
         */
        template <std::floating_point T, typename F>
        auto knot_breaks(F &&station, const std::array<std::shared_ptr<Curve<T, 2>>, 2> &pc, const std::array<Support<T>, 2> &sup, T a, T b) -> std::vector<T>
        {
            std::vector<T> breaks;
            for (std::size_t k = 0; k < 2; ++k)
                for (std::size_t d = 0; d < 2; ++d)
                    for (const T x : weak_knots(*sup[k].surface, d))
                    {
                        const T p = view_periods(*sup[k].surface)[d];
                        for (const T t : crossings(*pc[k], a, b, d, x, p))
                        {
                            auto g = [&](T u) {
                                const auto st = station(u);
                                const T y = (k == 0 ? st.uv1 : st.uv2)[d] - x;
                                return p > T(0) ? y - p * std::round(y / p) : y;
                            };
                            const T w = T(1e-3) * (b - a);
                            T lo = std::max(a, t - w), hi = std::min(b, t + w);
                            T glo = g(lo);
                            if (!(glo * g(hi) < T(0)))
                            {
                                breaks.push_back(t);
                                continue;
                            }
                            for (int j = 0; j < 60 && hi - lo > T(1e-14) * (T(1) + std::abs(a) + std::abs(b)); ++j)
                            {
                                const T mid = (lo + hi) / T(2), gm = g(mid);
                                if (gm * glo > T(0))
                                    lo = mid, glo = gm;
                                else
                                    hi = mid;
                            }
                            breaks.push_back((lo + hi) / T(2));
                        }
                    }
            return breaks;
        }

        template <std::floating_point T, std::size_t d>
        auto joined(const std::vector<BSCurve<T, d>> &cs) -> BSCurve<T, d>;

        /// f(), a refusal without shapes (a station of a G2 section, a skin: they know no edge) given the edge e.
        template <typename F>
        auto with_edge(EdgeId e, F &&f) -> decltype(f())
        {
            try
            {
                return f();
            }
            catch (const BuildFailure &x)
            {
                if (!x.error().shapes.empty())
                    throw;
                auto err = x.error();
                err.shapes = {ShapeId{e}};
                throw BuildFailure(std::move(err));
            }
        }

        template <std::floating_point T>
        auto fillet(const Model<T> &m, SolidId sid, EdgeId eid, T r, const FilletOptions<T> &o, const SectionOf<T> &sec = {}) -> EdgeBlend<T>
        {
            const std::string op = "fillet_geometry";
            check_options(o, op, eid);
            if (!(r > o.tol) || !std::isfinite(r))
                refuse<T>(BuildErrc::InvalidBounds, op, "the radius must be larger than the tolerance", {eid});
            const auto st = setup(m, sid, eid, op);
            const auto &ed = m.edge(eid);
            const T tol = o.tol;
            const T ext = st.closed ? T(0) : o.extension * r;
            // the contacts within r / tan(a / 2) of the edge, a the angle of the wedge
            const T reach = r / std::tan(std::max(st.min_angle, T(1e-3)) / T(2)) + r + ext;
            const T side = st.convex ? T(-1) : T(1); // the ball along the outward normals
            std::array<Support<T>, 2> sup;
            std::array<ExactOffset<T>, 2> off;
            std::array<std::array<T, 4>, 2> obox;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &sd = st.sides[k];
                sup[k] = support(m, sd, ed, reach, tol);
                // a face closing across a seam that cannot be evaluated past it: seen past its seam
                if (!dynamic_cast<const ElementarySurface<T> *>(sup[k].surface.get()))
                    if (const auto per = unrepeated(*sup[k].surface, periods_of(*sup[k].surface, tol), tol); per[0] > T(0) || per[1] > T(0))
                        sup[k].surface = std::make_shared<SeamlessSurface<T>>(sup[k].surface, per);
                const T d = side * r * (sd.orient == Orientation::Forward ? T(1) : T(-1));
                const auto x = offset_surface<T>(sup[k].surface, d, tol, sup[k].box);
                if (!x)
                    refuse<T>(BuildErrc::NonGeneric, op, "the radius exceeds a radius of curvature of a face", {eid, sd.face});
                off[k] = *x;
                obox[k] = {sup[k].box[0], sup[k].box[1], sup[k].box[2] + x->dv, sup[k].box[3] + x->dv};
            }
            IntersectOptions<T> io = o.intersect;
            io.tol = std::max(io.tol, tol / T(10));
            const auto res = intersect(*off[0].surface, obox[0], *off[1].surface, obox[1], io);

            // ---- the branch of the centre line nearest the edge
            const T um = (ed.u1 + ed.u2) / T(2);
            const auto pm = ed.curve->value(um);
            const auto q1 = m.vertex(ed.v1).pnt, q2 = m.vertex(ed.v2).pnt;
            const IntersectionCurve<T> *best = nullptr;
            T best_gap = std::numeric_limits<T>::infinity();
            for (const auto &c : res.curves)
            {
                const auto [s, d] = extrema_curve_point(*c.curve, pm, T(1e-12) * (T(1) + std::abs(c.b - c.a)));
                if (d < best_gap)
                    best_gap = d, best = &c;
            }
            if (!best)
                refuse<T>(BuildErrc::NonGeneric, op, "no centre line near the edge (the offsets of its faces do not meet)", {eid});
            const auto &c = *best;
            const auto *elementary = dynamic_cast<const ElementaryCurve<T> *>(c.curve.get());
            const T span = c.b - c.a, eps = T(1e-12) * (T(1) + std::abs(c.a) + std::abs(c.b));
            const bool loop = c.closed || norm(c.curve->value(c.a) - c.curve->value(c.b)) <= T(10) * std::max(tol, io.tol);
            auto proj = [&](const point<T, 3> &p) { return extrema_curve_point(*c.curve, p, T(1e-12) * (T(1) + std::abs(span)))[0]; };
            auto wrap = [&](T t) { return loop && t >= c.b - eps ? t - span : t; };
            T start{}, end{}, t1{}, t2{};
            if (st.closed)
            {
                if (!loop)
                {
                    // seen past their seams, two faces stop the marching at theirs: it closes only
                    // where both seams meet in the same section of the fillet
                    if (dynamic_cast<const SeamlessSurface<T> *>(sup[0].surface.get()) && dynamic_cast<const SeamlessSurface<T> *>(sup[1].surface.get()))
                        refuse<T>(BuildErrc::NonGeneric, op,
                                  "a closed edge between two closed B-spline faces whose seams do not meet in the same section of the fillet (not handled)",
                                  {eid, st.sides[0].face, st.sides[1].face});
                    refuse<T>(BuildErrc::NonGeneric, op, "the centre line of a closed edge is not closed", {eid});
                }
                start = wrap(proj(q1));
                // round a face closing across a seam without being periodic: from where the
                // contact on it crosses the seam, the contact then within the bounds of the face
                for (std::size_t k = 0; k < 2; ++k)
                    if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(sup[k].surface.get()))
                    {
                        const auto per = sv->periods();
                        const auto b = sv->bounds();
                        const std::size_t d = per[0] > T(0) ? 0 : 1;
                        if (const auto x = crossings(*(k == 0 ? c.pc1 : c.pc2), c.a, c.b, d, b[2 * d] + (d == 1 ? off[k].dv : T(0)), per[d]); !x.empty())
                        {
                            start = wrap(x.front());
                            // a G2 section: from where its own contact crosses the seam (wider than
                            // that of the ball: elsewhere), found on the stations of the marched
                            // centre line, outward from there by steps of span / 1024 (up to an
                            // eighth of it), then by bisection
                            if (sec.g2)
                            {
                                const std::array<const Surface<T, 3> *, 2> on{sup[0].surface.get(), sup[1].surface.get()};
                                auto g = [&](T t) {
                                    t = c.a + std::fmod(t - c.a, span);
                                    if (t < c.a)
                                        t += span;
                                    Station<T> x0;
                                    x0.c = c.curve->value(t);
                                    const auto w1 = c.pc1->value(t), w2 = c.pc2->value(t);
                                    x0.uv1 = {w1[0], w1[1] - off[0].dv}, x0.uv2 = {w2[0], w2[1] - off[1].dv};
                                    x0.p1 = on[0]->value(x0.uv1[0], x0.uv1[1]), x0.p2 = on[1]->value(x0.uv2[0], x0.uv2[1]);
                                    const auto q = g2_station(x0, on);
                                    const T y = (k == 0 ? q.uv1 : q.uv2)[d] - b[2 * d];
                                    return y - per[d] * std::round(y / per[d]);
                                };
                                const T g0 = g(start), step = span / T(1024);
                                std::optional<std::array<T, 2>> bracket;
                                for (int j = 1; j <= 128 && !bracket; ++j)
                                    for (const T dir : {T(-1), T(1)})
                                    {
                                        const T t0 = start + dir * step * T(j - 1), t1 = start + dir * step * T(j);
                                        const T g1 = g(t1);
                                        if (g1 * (j == 1 ? g0 : g(t0)) <= T(0) && std::abs(g1 - (j == 1 ? g0 : g(t0))) < per[d] / T(2))
                                        {
                                            bracket = std::array<T, 2>{t0, t1};
                                            break;
                                        }
                                    }
                                if (!bracket)
                                    refuse<T>(BuildErrc::NonGeneric, op, "the contact of a G2 section does not cross the seam of a closed face near that of the ball", {eid, st.sides[k].face});
                                T lo = (*bracket)[0], hi = (*bracket)[1], glo = g(lo);
                                while (std::abs(hi - lo) > T(1e-13) * (T(1) + std::abs(span)))
                                {
                                    const T mid = (lo + hi) / T(2), gm = g(mid);
                                    if (gm * glo > T(0))
                                        lo = mid, glo = gm;
                                    else
                                        hi = mid;
                                }
                                start = wrap(c.a + std::fmod((lo + hi) / T(2) - c.a + span, span));
                            }
                            break;
                        }
                    }
                // within eps of the start of a B-spline: from its start (it throws beyond knot_eps of its bounds)
                if (!elementary && start < c.a + eps)
                    start = c.a;
                end = start + span;
                // the sense of the centre line along the edge, at the middle of the edge
                const bool along = c.curve->value(proj(pm), 1) * ed.curve->value(um, 1) >= T(0);
                t1 = along ? start : end;
                t2 = along ? end : start;
            }
            else
            {
                const T s1 = proj(q1), s2 = proj(q2);
                bool along = s1 <= s2; // the centre line runs from v1 to v2
                if (!loop)
                    start = std::min(s1, s2), end = std::max(s1, s2);
                else
                {
                    // the arc between the vertices whose middle is nearest that of the edge
                    const T w1 = wrap(s1), w2 = wrap(s2);
                    auto at = [&](T t) {
                        t = c.a + std::fmod(t - c.a, span);
                        return c.curve->value(t < c.a ? t + span : t);
                    };
                    const std::array<std::array<T, 2>, 2> arcs{{{w1, w2 > w1 ? w2 : w2 + span}, {w2, w1 > w2 ? w1 : w1 + span}}};
                    const std::size_t i = norm(at((arcs[0][0] + arcs[0][1]) / T(2)) - pm) <= norm(at((arcs[1][0] + arcs[1][1]) / T(2)) - pm) ? 0 : 1;
                    start = arcs[i][0], end = arcs[i][1];
                    along = i == 0;
                }
                t1 = along ? start : end;
                t2 = along ? end : start;
                // continued past the vertices
                const T room = loop ? std::max(T(0), span - (end - start)) / T(2) : std::numeric_limits<T>::infinity();
                const T da = std::min(room, ext / norm(c.curve->value(std::clamp(start, c.a, c.b), 1)));
                const T db = std::min(room, ext / norm(c.curve->value(std::clamp(end <= c.b ? end : end - span, c.a, c.b), 1)));
                start -= da, end += db;
                if (!loop && !elementary)
                    start = std::max(start, c.a), end = std::min(end, c.b);
                if (loop && start < c.a - eps)
                    start += span, end += span, t1 += span, t2 += span;
                if (loop && !elementary && start < c.a + eps)
                {
                    // from the start of the B-spline, within its bounds (it throws beyond knot_eps of them)
                    start = c.a, end = std::min(end, c.b);
                    t1 = std::clamp(t1, start, end), t2 = std::clamp(t2, start, end);
                }
            }
            if (!(end - start > T(0)))
                refuse<T>(BuildErrc::NonGeneric, op, "the centre line reduces to a point", {eid});

            // ---- the centre line on [start, end] and its pcurves on the supports of the faces
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.a = start, rb.b = end, rb.t1 = t1, rb.t2 = t2;
            // a pcurve on an offset brought onto the support: (u, v - dv) on a cone, before it is joined from the vertex
            auto on_support = [&](const std::shared_ptr<Curve<T, 2>> &p, std::size_t k) {
                auto s = offset_detail::shifted(p, -off[k].dv);
                if (!s)
                    refuse<T>(BuildErrc::NonGeneric, op, "a pcurve on a cone that cannot be moved", {eid});
                return s;
            };
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pc;
            if (elementary)
            {
                rb.spine = intersect_detail::restricted(*elementary, start, end);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto *e = dynamic_cast<const ElementarySurface<T> *>(off[k].surface.get());
                    if (!e)
                        refuse<T>(BuildErrc::NonGeneric, op, "a closed-form centre line on a surface that is not a primitive", {eid});
                    pc[k] = on_support(intersect_detail::primitive_pcurve(*rb.spine, start, end, *e, obox[k], io.tol).pcurve, k);
                }
            }
            else
            {
                const std::array<std::shared_ptr<Curve<T, 2>>, 2> pcs{on_support(c.pc1, 0), on_support(c.pc2, 1)};
                if (end > c.b)
                {
                    // past the end of a loop: run from start (start > c.a + eps here)
                    rb.spine = thick_detail::from_parameter(c.curve, c.a, c.b, start);
                    for (std::size_t k = 0; k < 2; ++k)
                        pc[k] = thick_detail::from_parameter(pcs[k], c.a, c.b, start);
                }
                else
                    rb.spine = c.curve, pc = pcs;
            }

            // ---- the contacts
            for (std::size_t k = 0; k < 2; ++k)
            {
                auto &ct = rb.contacts[k];
                ct.face = st.sides[k].face;
                ct.support = sup[k].surface;
                ct.pcurve = pc[k];
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, st.sides[k], *ct.pcurve, ed, t1, t2, tol, op, eid);
            }

            // a fillet vanishing at a vertex (its faces tangent there: its contacts meet) ends
            // there, not continued past it: no face is to cut it at that end
            if (!st.closed)
                for (const T t : {t1, t2})
                    if (norm(rb.contacts[0].curve->value(t) - rb.contacts[1].curve->value(t)) <= T(10) * tol)
                        (t == std::min(t1, t2) ? rb.a : rb.b) = t;

            // ---- the surface
            rb.circular = !sec.conic && !sec.g2;
            std::shared_ptr<ElementarySurface<T>> primitive;
            if (rb.circular)
                primitive = canal(*rb.spine, *rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, r);
            const std::array<const Surface<T, 3> *, 2> supports{sup[0].surface.get(), sup[1].surface.get()};
            // a G2 section around a straight centre line: the same quartic all along, translated,
            // where its contacts run straight on (uv) between those at the ends (planes, cylinders
            // along the edge): its extrusion, exact (degree (4, 1)), the contacts those of the quartic
            auto g2_extruded = [&]() -> bool {
                if (!sec.g2 || !dynamic_cast<const LineCurve<T> *>(rb.spine.get()))
                    return false;
                auto at = [&](T t) {
                    Station<T> x;
                    x.c = rb.spine->value(t);
                    x.uv1 = pc[0]->value(t), x.uv2 = pc[1]->value(t);
                    x.p1 = rb.contacts[0].curve->value(t), x.p2 = rb.contacts[1].curve->value(t);
                    return g2_station(x, supports);
                };
                const auto sa = at(rb.a), sb = at(rb.b), sm = at((rb.a + rb.b) / T(2));
                for (std::size_t i = 0; i < 5; ++i)
                    if (!(norm((*sm.g2)[i] - ((*sa.g2)[i] + (*sb.g2)[i]) / T(2)) <= tol / T(10)))
                        return false;
                std::array<std::shared_ptr<Curve<T, 2>>, 2> g2pc{segment<T>(sa.uv1, sb.uv1, rb.a, rb.b), segment<T>(sa.uv2, sb.uv2, rb.a, rb.b)};
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto q = g2pc[k]->value((rb.a + rb.b) / T(2));
                    if (!(norm(sup[k].surface->value(q[0], q[1]) - (k == 0 ? sm.p1 : sm.p2)) <= tol / T(10)))
                        return false;
                }
                std::vector<std::array<T, 4>> poles;
                for (const auto *x : {&sa, &sb})
                    for (const auto &q : *x->g2)
                        poles.push_back({q[0], q[1], q[2], T(1)});
                rb.surface = std::make_shared<BSSurfaceRational<T, 3>>(poles, std::vector<T>{0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, std::vector<T>{rb.a, rb.a, rb.b, rb.b},
                                                                       std::size_t{4}, std::size_t{1});
                rb.exact = true;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    rb.contacts[k].pcurve = g2pc[k];
                    rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(g2pc[k], rb.contacts[k].support);
                }
                return true;
            };
            if (primitive)
            {
                rb.surface = primitive;
                rb.exact = true;
            }
            else if (with_edge(eid, g2_extruded))
            {
            }
            else if (sec.conic && sec.constant && dynamic_cast<const LineCurve<T> *>(rb.spine.get()))
            {
                // a conic of constant rho around a straight centre line: the same section all
                // along, translated: its extrusion, exact (rational, degree (2, 1); u across, v along)
                std::vector<std::array<T, 4>> poles;
                for (const T t : {rb.a, rb.b})
                {
                    const auto sc = conic_section(rb.spine->value(t), rb.contacts[0].curve->value(t), rb.contacts[1].curve->value(t), sec.weight(t));
                    for (const auto &q : sc.poles())
                        poles.push_back(q);
                }
                rb.surface = std::make_shared<BSSurfaceRational<T, 3>>(poles, std::vector<T>{0, 0, 0, 1, 1, 1}, std::vector<T>{rb.a, rb.a, rb.b, rb.b}, std::size_t{2},
                                                                       std::size_t{1});
                rb.exact = true;
            }
            else
            {
                // the stations of the ball, exact: from the marched centre line and its pcurves on
                // the offsets, Gauss-Newton on O1(U1, V1) = O2(U2, V2) and (O1 - c) . c' = 0 (the
                // centre in the plane across the marched line at s; Choi and Ju). The marched
                // centre line is accurate to its tolerance divided by the sine of the angle of the
                // offsets, poor where the faces become tangent; its contacts would no longer be at
                // the same distance from it and the arcs no circles. The marched station is kept
                // where Gauss-Newton does not converge or goes further than r / 10 (another root).
                // The pcurves pc are on the supports: (u, v + dv) on the offsets
                const T scale = T(1) + norm(rb.spine->value(rb.a)) + r;
                auto station = [&](T s) {
                    const auto c0 = rb.spine->value(s);
                    const auto t0 = unit(rb.spine->value(s, 1));
                    const auto w1 = pc[0]->value(s), w2 = pc[1]->value(s);
                    thick_detail::Vec<T> x(4), x0(4);
                    x0 << w1[0], w1[1] + off[0].dv, w2[0], w2[1] + off[1].dv;
                    x = x0;
                    const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                        std::array<point<T, 3>, 4> D1, D2;
                        off[0].surface->derivatives(y[0], y[1], 1, 1, D1.data());
                        off[1].surface->derivatives(y[2], y[3], 1, 1, D2.data());
                        // (derivatives: [0] value, [1] d/dv, [2] d/du, [3] d2/dudv)
                        const auto d = D1[0] - D2[0];
                        f.resize(4);
                        J.resize(4, 4);
                        f << d[0], d[1], d[2], (D1[0] - c0) * t0;
                        for (int i = 0; i < 3; ++i)
                            J.row(i) << D1[2][i], D1[1][i], -D2[2][i], -D2[1][i];
                        J.row(3) << D1[2] * t0, D1[1] * t0, T(0), T(0);
                    });
                    if (const T jump = norm(off[0].surface->value(x[0], x[1]) - c0); !(res <= T(1e-12) * scale) || jump > T(0.1) * r)
                        x = x0;
                    Station<T> st;
                    st.c = off[0].surface->value(x[0], x[1]);
                    st.uv1 = {x[0], x[1] - off[0].dv};
                    st.uv2 = {x[2], x[3] - off[1].dv};
                    st.p1 = sup[0].surface->value(st.uv1[0], st.uv1[1]);
                    st.p2 = sup[1].surface->value(st.uv2[0], st.uv2[1]);
                    return st;
                };
                // where a contact crosses a knot of its face below C3, a station
                const auto breaks = knot_breaks<T>(station, pc, sup, rb.a, rb.b);
                // (a law of rho on the parameter of the edge: brought onto that of the centre line,
                // linearly between the vertices, where they match)
                auto on_spine = sec;
                if (sec.conic && !sec.constant && rb.t2 != rb.t1)
                {
                    const T u1 = ed.u1, u2 = ed.u2, s1 = rb.t1, s2 = rb.t2;
                    on_spine.rho = [f = sec.rho, u1, u2, s1, s2](T s) { return f(u1 + (s - s1) * (u2 - u1) / (s2 - s1)); };
                    for (auto &k : on_spine.kinks)
                        k = s1 + (k - u1) * (s2 - s1) / (u2 - u1);
                }
                // where the law of rho kinks the surface has a crease (the sections, the same on
                // either side, sweep on otherwise), that one smooth skin cannot follow: a skin on
                // each interval between the kinks, the pieces sharing the section at each (an arc
                // common to both, `blend`), the centre line and the contacts joined across them
                std::vector<T> cuts{rb.a};
                for (const T k : on_spine.kinks)
                    if (k > rb.a + T(1e-9) * (rb.b - rb.a) && k < rb.b - T(1e-9) * (rb.b - rb.a))
                        cuts.push_back(k);
                std::ranges::sort(cuts);
                cuts.push_back(rb.b);
                if (cuts.size() > 2)
                {
                    std::vector<BSCurve<T, 3>> spines;
                    std::array<std::vector<BSCurve<T, 2>>, 2> pcs;
                    for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
                    {
                        const auto sk = with_edge(eid, [&] {
                            return skinned(station, std::array<const Surface<T, 3> *, 2>{sup[0].surface.get(), sup[1].surface.get()}, cuts[i], cuts[i + 1],
                                           o.approx_tol, o.angular_tol, breaks, on_spine);
                        });
                        spines.push_back(*std::dynamic_pointer_cast<BSCurve<T, 3>>(sk.spine));
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            auto c = *std::dynamic_pointer_cast<BSCurve<T, 2>>(sk.pcurves[k]);
                            if (!pcs[k].empty())
                            {
                                // on the periods of the face, on from the piece before
                                const auto per = view_periods(*sup[k].surface);
                                const auto prev = pcs[k].back().value(pcs[k].back().bounds()[1]), first = c.value(c.bounds()[0]);
                                point<T, 2> shift{};
                                for (std::size_t d = 0; d < 2; ++d)
                                    if (per[d] > T(0))
                                        shift[d] = per[d] * std::round((prev[d] - first[d]) / per[d]);
                                translate(c, shift);
                            }
                            pcs[k].push_back(c);
                        }
                        BlendPiece<T> pi;
                        pi.surface = sk.surface;
                        pi.a = cuts[i], pi.b = cuts[i + 1];
                        pi.deviation = sk.deviation;
                        rb.deviation = std::max(rb.deviation, pi.deviation);
                        rb.pieces.push_back(std::move(pi));
                    }
                    rb.spine = std::make_shared<BSCurve<T, 3>>(joined(spines));
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        std::shared_ptr<Curve<T, 2>> pk = std::make_shared<BSCurve<T, 2>>(joined(pcs[k]));
                        if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(sup[k].surface.get()))
                        {
                            pk = into_bounds(pk, *sv, tol);
                            if (!pk)
                                refuse<T>(BuildErrc::NonGeneric, op, "a contact crossing the seam of a face that is not periodic (a closed B-spline)", {eid, st.sides[k].face});
                        }
                        rb.contacts[k].pcurve = pk;
                        rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(pk, rb.contacts[k].support);
                    }
                    for (auto &pi : rb.pieces)
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            pi.on_blend[k] = segment<T>({T(k), pi.a}, {T(k), pi.b}, pi.a, pi.b);
                            for (int j = 0; j <= 64; ++j)
                            {
                                const T t = pi.a + (pi.b - pi.a) * T(j) / T(64);
                                const auto q = pi.on_blend[k]->value(t);
                                rb.contacts[k].deviation = std::max(rb.contacts[k].deviation, norm(pi.surface->value(q[0], q[1]) - rb.contacts[k].curve->value(t)));
                            }
                        }
                    rb.surface = nullptr;
                    return rb;
                }
                auto sk = with_edge(eid, [&] {
                    return skinned(station, std::array<const Surface<T, 3> *, 2>{sup[0].surface.get(), sup[1].surface.get()}, rb.a, rb.b, o.approx_tol,
                                   o.angular_tol, breaks, on_spine);
                });
                rb.surface = sk.surface;
                rb.deviation = sk.deviation;
                rb.spine = sk.spine;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(sup[k].surface.get()))
                    {
                        sk.pcurves[k] = into_bounds(sk.pcurves[k], *sv, tol);
                        if (!sk.pcurves[k])
                            refuse<T>(BuildErrc::NonGeneric, op, "a contact crossing the seam of a face that is not periodic (a closed B-spline)", {eid, st.sides[k].face});
                    }
                    rb.contacts[k].pcurve = sk.pcurves[k];
                    rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(sk.pcurves[k], rb.contacts[k].support);
                }
            }
            // the contacts of a G2 section, wider than those of the ball: on the faces too
            if (sec.g2)
                for (std::size_t k = 0; k < 2; ++k)
                    check_on_face(m, st.sides[k], *rb.contacts[k].pcurve, ed, t1, t2, tol, op, eid);
            on_blend(rb, io.tol, true);
            return rb;
        }

        /// B-splines end to end, each starting where the one before ends: one B-spline, C0 at the
        /// joins (`gbs::join`: the pole there the middle of both).
        template <std::floating_point T, std::size_t d>
        auto joined(const std::vector<BSCurve<T, d>> &cs) -> BSCurve<T, d>
        {
            return std::accumulate(std::next(cs.begin()), cs.end(), cs.front(), [](const BSCurve<T, d> &x, const BSCurve<T, d> &y) { return BSCurve<T, d>{join(x, y)}; });
        }

        /**
         * The crease where the surfaces s1 and s2 of two pieces of fillets meet (u across, 0 at
         * contact 0 and 1 at contact 1 on s1, the same on s2 where `same`, swapped otherwise; v
         * along), from their contacts at v = t1 on s1 and v = t2 on s2, each continued past them:
         * for w across s1, Gauss-Newton on s1(w, v1) = s2(u2, v2), marched from the middle out by
         * continuation, the ends the contacts themselves; v1, u2, v2 interpolated (cubic) on w,
         * the points doubled from 32 until within tol / 10 of the exact points between them.
         * nullopt where a point cannot be found.
         */
        template <std::floating_point T>
        auto crease_curve(const std::shared_ptr<Surface<T, 3>> &s1, T t1, const std::shared_ptr<Surface<T, 3>> &s2, T t2, bool same, T tol)
            -> std::optional<BlendCrease<T>>
        {
            const T scale = T(1) + norm(s1->value(T(0.5), t1));
            auto solve = [&](T w, std::array<T, 3> x0) -> std::optional<std::array<T, 3>> {
                thick_detail::Vec<T> x(3);
                x << x0[0], x0[1], x0[2];
                const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                    const auto d = s1->value(w, y[0]) - s2->value(y[1], y[2]);
                    const auto a = s1->value(w, y[0], 0, 1), bu = s2->value(y[1], y[2], 1, 0), bv = s2->value(y[1], y[2], 0, 1);
                    f.resize(3);
                    J.resize(3, 3);
                    for (int i = 0; i < 3; ++i)
                    {
                        f[i] = d[std::size_t(i)];
                        J.row(i) << a[std::size_t(i)], -bu[std::size_t(i)], -bv[std::size_t(i)];
                    }
                });
                // (near the contacts the surfaces become tangent: a point within tol / 100 of both)
                if (!(res <= std::max(T(1e-11) * scale, tol / T(100))))
                    return std::nullopt;
                return std::array<T, 3>{x[0], x[1], x[2]};
            };
            auto across = [&](T w) { return same ? w : T(1) - w; };
            // the points w clustered toward the contacts (Chebyshev), where the surfaces become
            // tangent and the crease turns fastest
            auto wi = [](std::size_t i, std::size_t n) { return (T(1) - std::cos(std::numbers::pi_v<T> * T(i) / T(n))) / T(2); };
            for (std::size_t n = 32; n <= 1024; n *= 2)
            {
                std::vector<std::array<T, 3>> x(n + 1);
                x[0] = {t1, across(0), t2};
                x[n] = {t1, across(1), t2};
                const std::size_t mid = n / 2;
                auto xm = solve(wi(mid, n), {t1, across(wi(mid, n)), t2});
                if (!xm)
                    return std::nullopt;
                x[mid] = *xm;
                // from the middle out; near the contacts, a point not found is left to the
                // interpolation between its neighbours
                std::vector<bool> found(n + 1, true);
                for (int dir : {-1, 1})
                {
                    auto prev = x[mid];
                    for (std::size_t i = mid + dir; i > 0 && i < n; i += dir)
                    {
                        if (const auto y = solve(wi(i, n), prev))
                            x[i] = *y, prev = *y;
                        else
                            found[i] = false;
                    }
                }
                std::vector<T> ws;
                std::array<points_vector<T, 2>, 2> uv;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    if (!found[i])
                        continue;
                    ws.push_back(wi(i, n));
                    uv[0].push_back({ws.back(), x[i][0]});
                    uv[1].push_back({x[i][1], x[i][2]});
                }
                if (ws.size() < 5)
                    return std::nullopt;
                BlendCrease<T> c;
                for (std::size_t k = 0; k < 2; ++k)
                    c.on_pieces[k] = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(uv[k], ws, std::size_t{3}));
                // between the points: the curves on both pieces within tol / 10 of each other
                T gap{};
                for (std::size_t i = 0; i < n; ++i)
                    for (int j = 1; j < 4; ++j)
                    {
                        const T w = wi(i, n) + (wi(i + 1, n) - wi(i, n)) * T(j) / T(4);
                        const auto p = c.on_pieces[0]->value(w), q = c.on_pieces[1]->value(w);
                        gap = std::max(gap, norm(s1->value(p[0], p[1]) - s2->value(q[0], q[1])));
                    }
                c.deviation = gap;
                c.curve = std::make_shared<CurveOnSurface<T, 3>>(c.on_pieces[0], s1);
                if (gap <= tol / T(10) || n == 1024)
                    return gap <= tol ? std::optional<BlendCrease<T>>{c} : std::nullopt;
            }
            return std::nullopt;
        }

        /// A closed curve on [u1, u2] read on [lo, hi], its parameter brought back into [u1, u2] by
        /// the period u2 - u1: a loop run from another of its points.
        template <std::floating_point T>
        class LoopCurve : public Curve<T, 3>
        {
            std::shared_ptr<Curve<T, 3>> c_;
            T u1_, u2_, lo_, hi_;

        public:
            LoopCurve(std::shared_ptr<Curve<T, 3>> c, T u1, T u2, T lo, T hi) : c_{std::move(c)}, u1_{u1}, u2_{u2}, lo_{lo}, hi_{hi} {}
            auto value(T u, size_t d = 0) const -> std::array<T, 3> override
            {
                const T p = u2_ - u1_;
                if (u < u1_ || u > u2_)
                    u = u1_ + std::fmod(std::fmod(u - u1_, p) + p, p);
                return c_->value(u, d);
            }
            auto bounds() const -> std::array<T, 2> override { return {lo_, hi_}; }
        };

        /// A radius along an edge: r(t) on the parameter t of the edge (past its first and last
        /// points the law goes on with its slope there: constant for a smooth law, straight on for
        /// a linear one; past the vertices of an open edge, the law of the edge with its slope at
        /// the vertex), the parameters where the law has a knot (its slope may jump there), its
        /// values at the vertices and its largest value on the edge (the vertices included; past
        /// them, `fillet_law` adds the ends of the blend and the pieces continued past a peak).
        template <std::floating_point T>
        struct EdgeRadius
        {
            std::function<T(T)> r;
            std::vector<T> knots; ///< where the law is C1 but no smoother (a smooth law; a closed one: where it closes too)
            std::vector<T> kinks; ///< where its slope jumps (a linear law)
            T r1{}, r2{}, rmax{};
            /// whether a chain goes on at the vertex v1, v2: the blend continued there by 1.5 radius
            /// at least, room for the crease with the blend of the next edge (`joint`)
            std::array<bool, 2> joined{};
        };

        /**
         * The pcurve on the sphere s of the curve c on [0, 1] lying on it (an arc of a corner or of
         * a lune): the parameters of n + 1 of its points interpolated (cubic, u continuous), n
         * doubled from 32 while the interpolation strays more than tol / 10 from them between the
         * points (up to 1024: 33 points left a pcurve 1.4e-4 off its arc near a pole), and the
         * largest gap in space between the pcurve on s and c, at the points and between them (what
         * the tolerance of an edge on c must cover: `settle` and `check` sample at the points
         * interpolated).
         */
        template <std::floating_point T>
        auto sphere_pcurve(const SphericalSurface<T> &s, const Curve<T, 3> &c, T tol) -> std::pair<std::shared_ptr<BSCurve<T, 2>>, T>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            std::shared_ptr<BSCurve<T, 2>> on;
            std::size_t n = 32;
            for (;; n *= 2)
            {
                std::vector<point<T, 2>> uv;
                std::vector<T> w;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    w.push_back(T(i) / T(n));
                    auto q = s.parameters(c.value(w.back()));
                    if (!uv.empty())
                        q[0] += T(2) * pi * std::round((uv.back()[0] - q[0]) / (T(2) * pi));
                    uv.push_back(q);
                }
                on = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(uv, w, std::size_t{3}));
                T stray{};
                for (std::size_t i = 0; i < n; ++i)
                {
                    const T x = (T(i) + T(0.5)) / T(n);
                    const auto a = on->value(x), b = s.parameters(c.value(x));
                    stray = std::max(stray, norm(s.value(a[0], a[1]) - s.value(b[0], b[1])));
                }
                if (stray <= tol / T(10) || n >= 1024)
                    break;
            }
            T gap{};
            for (std::size_t i = 0; i <= 2 * n; ++i)
            {
                const T x = T(i) / T(2 * n);
                const auto a = on->value(x);
                gap = std::max(gap, norm(s.value(a[0], a[1]) - c.value(x)));
            }
            return {on, gap};
        }

        /**
         * The fillet of a radius varying along the edge (stage 5, § 3): at the parameter t of the
         * edge, the ball of radius r(t) in the plane normal to the edge at e(t), touching both
         * faces (Gauss-Newton on the parameters of its contacts, from the nearest station solved,
         * or from the edge with the radius grown in eight steps); its centres c(t) the centre line.
         * The surface is the envelope of the balls: at each station the arc between the contacts
         * of the characteristic circle, in the plane normal to c', its centre c + l c'/|c'| with
         * l the mean of (p - c) . c'/|c'| over the contacts (both are -r r'/|c'|), c' the
         * derivative of the station's equations (implicit: one linear solve, no other station);
         * skinned as the fillets of constant radius, a station at each knot of the law. Every
         * curve on the parameter of the edge, continued at each vertex by `extension` times the
         * larger of the radii at the vertices (1.5 times at least where a chain goes on: room for
         * the crease with the next blend). The supports reach, and the offsets are checked
         * against folding, for the largest radius of the blend: its ends and the pieces continued
         * past a peak included (the fillet done again at their radius where it is larger,
         * `rmax_floor`).
         */
        template <std::floating_point T>
        auto fillet_law(const Model<T> &m, SolidId sid, EdgeId eid, const EdgeRadius<T> &R, const FilletOptions<T> &o, T rmax_floor = T(0),
                        const SectionOf<T> &sec = {}) -> EdgeBlend<T>
        {
            const std::string op = "fillet_geometry";
            check_options(o, op, eid);
            const auto st = setup(m, sid, eid, op);
            const auto &ed = m.edge(eid);
            const T tol = o.tol;
            const T rv = std::max(R.r1, R.r2);
            const T ext = st.closed ? T(0) : std::max(o.extension * rv, R.joined[0] || R.joined[1] ? T(1.5) * rv : T(0));
            const T side = st.convex ? T(-1) : T(1);
            // the knots of the edge below C3 (the planes of the stations turn with the edge, the
            // characteristic circles with its second derivative) and below C2 (its curvature jumps)
            const auto weak3 = edge_weak_knots(*ed.curve, 3), weak2 = edge_weak_knots(*ed.curve, 2);
            // a piece may be continued past a peak: at a kink of the law, or where the curvature of
            // the edge jumps within it
            const bool may_peak = !R.kinks.empty() || std::ranges::any_of(weak2, [&](T k) { return k > ed.u1 && k < ed.u2; });
            // the largest radius of the blend: on the edge, and at its ends (a linear law goes on
            // with its slope past the vertices): first where the edge continued at its speed at
            // each vertex ends, then, if larger, at the ends of the continued edge
            const T v1 = norm(ed.curve->value(ed.u1, 1)), v2 = norm(ed.curve->value(ed.u2, 1));
            T rmax = st.closed ? std::max(R.rmax, rmax_floor) : std::max({R.rmax, rmax_floor, R.r(ed.u1 - ext / v1), R.r(ed.u2 + ext / v2)});
            std::array<Support<T>, 2> sup;
            std::array<T, 2> sign{};
            T a = ed.u1, b = ed.u2, reach{}, edge_lo = ed.u1, edge_hi = ed.u2;
            std::shared_ptr<Curve<T, 3>> edge;
            for (int pass = 0; pass < 2; ++pass)
            {
                reach = rmax / std::tan(std::max(st.min_angle, T(1e-3)) / T(2)) + rmax + ext;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto &sd = st.sides[k];
                    sup[k] = support(m, sd, ed, reach, tol);
                    if (!dynamic_cast<const ElementarySurface<T> *>(sup[k].surface.get()))
                        if (const auto per = unrepeated(*sup[k].surface, periods_of(*sup[k].surface, tol), tol); per[0] > T(0) || per[1] > T(0))
                            sup[k].surface = std::make_shared<SeamlessSurface<T>>(sup[k].surface, per);
                    sign[k] = side * (sd.orient == Orientation::Forward ? T(1) : T(-1));
                    if (!offset_surface<T>(sup[k].surface, sign[k] * rmax, tol, sup[k].box))
                        refuse<T>(BuildErrc::NonGeneric, op, "the radius exceeds a radius of curvature of a face", {eid, sd.face});
                }
                // the edge, continued past its vertices: the planes of the stations; further than
                // the blend (by 1.5 radius) where a piece may be continued past a peak near a vertex
                a = ed.u1, b = ed.u2;
                if (st.closed)
                    edge = std::make_shared<LoopCurve<T>>(ed.curve, ed.u1, ed.u2, ed.u1, ed.u2);
                else
                {
                    T lo{}, hi{};
                    edge = continued_edge(ed, may_peak ? ext + T(1.5) * rmax : ext, st.sides, {sup[0].surface, sup[1].surface}, tol, lo, hi);
                    a = std::max(lo, ed.u1 - ext / v1), b = std::min(hi, ed.u2 + ext / v2);
                    edge_lo = lo, edge_hi = hi;
                }
                if (st.closed || !(std::max(R.r(a), R.r(b)) > rmax * (T(1) + T(1e-9))))
                    break;
                rmax = std::max(R.r(a), R.r(b));
            }
            for (const T t : {a, b})
                if (!(R.r(t) > tol))
                    refuse<T>(BuildErrc::NonGeneric, op, "the radius law, going on with its slope past its first or last point, is not positive at an end of the blend",
                              {eid});

            // ---- the stations
            const T scale = T(1) + norm(ed.curve->value(ed.u1)) + rmax;
            const T period = ed.u2 - ed.u1;
            auto on_edge = [&](T t) { return st.closed ? ed.u1 + std::fmod(std::fmod(t - ed.u1, period) + period, period) : std::clamp(t, ed.u1, ed.u2); };
            // the stations of a law (R.r, or that of a piece continued past a peak of the law),
            // each law with the stations it solved
            using Law = std::function<T(T)>;
            using Solved = std::map<T, std::array<T, 4>>;
            Solved solved;
            // (the plane at t normal to the edge just inside [lo, hi]: at a corner of the edge, that
            // of the side of the piece)
            auto solve = [&](T t, const std::array<T, 4> &x0, T r, T lo, T hi) -> std::optional<std::array<T, 4>> {
                const auto e = edge->value(t);
                const auto te = unit(edge->value(std::clamp(t, lo + T(1e-9) * (hi - lo), hi - T(1e-9) * (hi - lo)), 1));
                const OffsetSurface<T> o1{sup[0].surface, sign[0] * r}, o2{sup[1].surface, sign[1] * r};
                thick_detail::Vec<T> x(4);
                x << x0[0], x0[1], x0[2], x0[3];
                const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                    std::array<point<T, 3>, 4> D1, D2;
                    o1.derivatives(y[0], y[1], 1, 1, D1.data());
                    o2.derivatives(y[2], y[3], 1, 1, D2.data());
                    // (derivatives: [0] value, [1] d/dv, [2] d/du, [3] d2/dudv)
                    const auto d = D1[0] - D2[0];
                    f.resize(4);
                    J.resize(4, 4);
                    f << d[0], d[1], d[2], (D1[0] - e) * te;
                    for (int i = 0; i < 3; ++i)
                        J.row(i) << D1[2][i], D1[1][i], -D2[2][i], -D2[1][i];
                    J.row(3) << D1[2] * te, D1[1] * te, T(0), T(0);
                });
                // (another root: a ball touching the faces far from the edge)
                if (!(res <= T(1e-12) * scale) || !(norm(o1.value(x[0], x[1]) - e) <= reach))
                    return std::nullopt;
                return std::array<T, 4>{x[0], x[1], x[2], x[3]};
            };
            auto raw_of = [&](T t, const Law &rf, Solved &solved, T lo, T hi) {
                const T r = rf(t);
                std::optional<std::array<T, 4>> x;
                if (const auto it = solved.find(t); it != solved.end())
                    x = it->second;
                if (!x && !solved.empty())
                {
                    auto it = solved.lower_bound(t);
                    if (it == solved.end() || (it != solved.begin() && t - std::prev(it)->first < it->first - t))
                        --it;
                    x = solve(t, it->second, r, lo, hi);
                    // a ball that jumped to another root (past a sharp turn of the edge, the far side
                    // of a thin face in the plane): its contacts moved further than the edge
                    // could carry them; from the edge instead
                    if (x)
                    {
                        const auto &y = it->second;
                        const T dt = std::abs(t - it->first), speed = norm(edge->value(t, 1));
                        const T moved = std::max(norm(sup[0].surface->value((*x)[0], (*x)[1]) - sup[0].surface->value(y[0], y[1])),
                                                 norm(sup[1].surface->value((*x)[2], (*x)[3]) - sup[1].surface->value(y[2], y[3])));
                        if (moved > T(20) * speed * dt + T(0.1) * r)
                            x.reset();
                    }
                }
                if (!x)
                {
                    // from the edge (the ball of radius 0), the radius grown
                    const T tc = on_edge(t);
                    const auto w1 = st.sides[0].pcurve->value(tc), w2 = st.sides[1].pcurve->value(tc);
                    std::optional<std::array<T, 4>> y = std::array<T, 4>{w1[0], w1[1], w2[0], w2[1]};
                    for (int k = 1; k <= 8 && y; ++k)
                        y = solve(t, *y, r * T(k) / T(8), lo, hi);
                    x = y;
                }
                if (!x)
                    refuse<T>(BuildErrc::NonGeneric, op, "no ball of the radius touches both faces in the plane normal to the edge", {eid});
                solved[t] = *x;
                Station<T> s;
                s.uv1 = {(*x)[0], (*x)[1]};
                s.uv2 = {(*x)[2], (*x)[3]};
                s.p1 = sup[0].surface->value(s.uv1[0], s.uv1[1]);
                s.p2 = sup[1].surface->value(s.uv2[0], s.uv2[1]);
                s.c = OffsetSurface<T>{sup[0].surface, sign[0] * r}.value(s.uv1[0], s.uv1[1]);
                return s;
            };
            auto raw = [&](T t) { return raw_of(t, R.r, solved, a, b); };
            // the stations of a law marched from lo to hi (512 steps, halved where Gauss-Newton
            // fails or the contacts jump further than the edge carries them: another root, the far
            // side of a thin face in the plane), so that any station starts from one next to it
            // on the same branch
            // (from the station at `from`, already in the cache, to `to`, either way; [lo, hi] the
            // piece, where the planes are taken)
            auto march = [&](const Law &rf, Solved &cache, T from, T to, T lo, T hi) {
                auto prev = cache.at(from);
                const T full = std::abs(to - from), dir = to < from ? T(-1) : T(1);
                T t = from, step = full / T(512);
                while (dir * (to - t) > T(0))
                {
                    const T next = dir > T(0) ? std::min(to, t + step) : std::max(to, t - step);
                    auto x = solve(next, prev, rf(next), lo, hi);
                    bool ok = x.has_value();
                    if (ok)
                    {
                        const T speed = norm(edge->value(next, 1));
                        const T moved = std::max(norm(sup[0].surface->value((*x)[0], (*x)[1]) - sup[0].surface->value(prev[0], prev[1])),
                                                 norm(sup[1].surface->value((*x)[2], (*x)[3]) - sup[1].surface->value(prev[2], prev[3])));
                        ok = moved <= T(20) * speed * std::abs(next - t) + T(1e-9) * scale;
                    }
                    if (!ok)
                    {
                        if (!(step > T(1e-9) * (hi - lo)))
                            refuse<T>(BuildErrc::NonGeneric, op, "the balls of the fillet cannot be followed along the edge", {eid});
                        step /= T(2);
                        continue;
                    }
                    cache[next] = *x;
                    prev = *x;
                    t = next;
                    step = std::min(full / T(512), step * T(2));
                }
            };
            auto prime = [&](const Law &rf, Solved &cache, T lo, T hi) {
                (void)raw_of(lo, rf, cache, lo, hi);
                march(rf, cache, lo, hi, lo, hi);
            };
            // the characteristic circle: in the plane normal to the centre line, c' the derivative
            // of the equations of the station F(x, t) = 0 (x the parameters of the contacts): J dx
            // = -dF/dt, J their Jacobian, dF/dt from r' (one-sided within [lo, hi]: at a kink of
            // the law, the circle on that side), the direction of the edge and its derivative;
            // c' = dO1/dt, O1 the offset of the first face (dO/dr = (O - S) / r)
            const T hs = T(1e-6) * (b - a);
            auto tilted = [&](T t, T lo, T hi, const Law &rf, Solved &solved) {
                auto s = raw_of(t, rf, solved, lo, hi);
                const T r = rf(t);
                const T t0 = std::max(lo, t - hs), t1 = std::min(hi, t + hs);
                const T dr = (rf(t1) - rf(t0)) / (t1 - t0);
                const OffsetSurface<T> o1{sup[0].surface, sign[0] * r}, o2{sup[1].surface, sign[1] * r};
                std::array<point<T, 3>, 4> D1, D2;
                o1.derivatives(s.uv1[0], s.uv1[1], 1, 1, D1.data());
                o2.derivatives(s.uv2[0], s.uv2[1], 1, 1, D2.data());
                // (derivatives: [0] value, [1] d/dv, [2] d/du)
                const auto n1 = (D1[0] - s.p1) / r, n2 = (D2[0] - s.p2) / r;
                // (the derivatives of the edge just inside [lo, hi]: at a vertex where its curvature
                // jumps, those of the side of the piece)
                const T ti = std::clamp(t, lo + T(1e-9) * (hi - lo), hi - T(1e-9) * (hi - lo));
                const auto e = edge->value(t), de = edge->value(ti, 1), dde = edge->value(ti, 2);
                const T speed = norm(de);
                const auto te = de / speed, dte = (dde - (dde * te) * te) / speed;
                thick_detail::Mat<T> J(4, 4);
                thick_detail::Vec<T> g(4);
                for (int i = 0; i < 3; ++i)
                {
                    J.row(i) << D1[2][i], D1[1][i], -D2[2][i], -D2[1][i];
                    g[i] = dr * (n1[i] - n2[i]);
                }
                J.row(3) << D1[2] * te, D1[1] * te, T(0), T(0);
                g[3] = dr * (n1 * te) - speed + (D1[0] - e) * dte;
                const thick_detail::Vec<T> dx = J.colPivHouseholderQr().solve(-g);
                const auto dc = dx[0] * D1[2] + dx[1] * D1[1] + dr * n1;
                if (!dx.allFinite() || !(norm(dc) > T(0)))
                    refuse<T>(BuildErrc::NonGeneric, op, "the centre line of the balls stops (a station of singular equations)", {eid});
                const auto tc = dc / norm(dc);
                const T l = ((s.p1 - s.c) * tc + (s.p2 - s.c) * tc) / T(2);
                if (!(std::abs(l) < (T(1) - T(1e-6)) * norm(s.p1 - s.c)))
                    refuse<T>(BuildErrc::NonGeneric, op, "the radius varies faster than the centre of the ball moves (no envelope)", {eid});
                s.m = s.c + l * tc;
                return s;
            };

            prime(R.r, solved, a, b);
            // round a face that closes across a seam without being periodic: the loop from where the
            // contact on it crosses the seam (its contact then within the bounds of the face), the
            // edge read past its vertex
            std::vector<T> kinks = R.kinks, knots = R.knots;
            if (st.closed)
                for (std::size_t k = 0; k < 2; ++k)
                    if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(sup[k].surface.get()))
                    {
                        std::vector<T> ts;
                        points_vector<T, 2> uv;
                        const auto per = sv->periods();
                        for (int i = 0; i <= 64; ++i)
                        {
                            ts.push_back(a + (b - a) * T(i) / T(64));
                            const auto s0 = raw(ts.back());
                            auto q = k == 0 ? s0.uv1 : s0.uv2;
                            if (!uv.empty())
                                for (std::size_t d = 0; d < 2; ++d)
                                    if (per[d] > T(0))
                                        q[d] += per[d] * std::round((uv.back()[d] - q[d]) / per[d]);
                            uv.push_back(q);
                        }
                        const BSCurve<T, 2> c{interpolate<T, 2>(uv, ts, std::size_t{3})};
                        const std::size_t d = per[0] > T(0) ? 0 : 1;
                        const T seam = sv->bounds()[2 * d];
                        const auto x = crossings(c, a, b, d, seam, per[d]);
                        if (x.empty() || !(x.front() > a + T(1e-9) * period))
                            break;
                        // the crossing put back on the stations themselves (bisection)
                        auto g = [&](T u) {
                            const auto st0 = raw(u);
                            const T y = (k == 0 ? st0.uv1 : st0.uv2)[d] - seam;
                            return y - per[d] * std::round(y / per[d]);
                        };
                        T lo = std::max(a, x.front() - T(1e-3) * period), hi = std::min(b, x.front() + T(1e-3) * period);
                        T glo = g(lo);
                        if (glo * g(hi) < T(0))
                            for (int j = 0; j < 80 && hi - lo > T(1e-15) * (T(1) + std::abs(lo)); ++j)
                            {
                                const T mid = (lo + hi) / T(2), gm = g(mid);
                                if (gm * glo > T(0))
                                    lo = mid, glo = gm;
                                else
                                    hi = mid;
                            }
                        const T t0 = glo * g(hi) <= T(0) ? (lo + hi) / T(2) : x.front();
                        a = t0, b = t0 + period;
                        edge = std::make_shared<LoopCurve<T>>(ed.curve, ed.u1, ed.u2, a, b);
                        for (auto *v : {&kinks, &knots})
                            for (auto &t : *v)
                                if (t < a)
                                    t += period;
                        // the stations marched, round the loop from there (the same a period on)
                        Solved shifted;
                        for (const auto &[t, x] : solved)
                            shifted[t < a ? t + period : t] = x;
                        solved = std::move(shifted);
                        break;
                    }
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.variable = true;
            rb.a = a, rb.b = b;
            rb.t1 = st.closed ? a : ed.u1, rb.t2 = st.closed ? b : ed.u2;
            // the contacts at 65 stations, for the knots of the faces they cross
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pc;
            {
                std::vector<T> ts;
                std::array<points_vector<T, 2>, 2> uv;
                for (int i = 0; i <= 64; ++i)
                {
                    ts.push_back(a + (b - a) * T(i) / T(64));
                    const auto s = raw(ts.back());
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        auto q = k == 0 ? s.uv1 : s.uv2;
                        const auto per = view_periods(*sup[k].surface);
                        if (!uv[k].empty())
                            for (std::size_t d = 0; d < 2; ++d)
                                if (per[d] > T(0))
                                    q[d] += per[d] * std::round((uv[k].back()[d] - q[d]) / per[d]);
                        uv[k].push_back(q);
                    }
                }
                for (std::size_t k = 0; k < 2; ++k)
                    pc[k] = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(uv[k], ts, std::size_t{3}));
            }
            auto breaks = knot_breaks<T>(raw, pc, sup, a, b);
            breaks.insert(breaks.end(), knots.begin(), knots.end());
            // the knots of the edge below C3: the planes of the stations turn with the edge, the
            // characteristic circles with the second derivative of the edge (a radius varying)
            for (T t : weak3)
            {
                if (st.closed)
                    while (t < a)
                        t += period;
                if (t > a && t < b)
                    breaks.push_back(t);
            }

            // ---- the pieces, each skinned on its own, between the points where the balls jump
            // from one characteristic circle to another: kinks of the law, and points where the
            // curvature of the edge jumps (the vertex of a closed edge, knots of the edge below
            // C2) while the radius varies; the centre line has a corner there. At such a point
            // the balls next to it cover the middle of the lune of its ball between the two
            // circles (the envelopes on either side meet along a crease: a peak of the law), or
            // not (the lune is part of the fillet: a hollow of the law); the same circle on either
            // side (within tol): no cut. The stations of either side and of the balls next to the
            // point from those of the law marched (`solved`), on the branch of the fillet.
            enum class Cut
            {
                None,
                Lune,
                Crease
            };
            const T eps = T(1e-9) * (b - a);
            auto kind = [&](T t) {
                const T w = T(1e-3) * (b - a), e1 = T(1e-9) * (b - a);
                const auto L = tilted(t - e1, t - w, t - e1, R.r, solved), Rt = tilted(t + e1, t + e1, t + w, R.r, solved);
                const auto al = arc(*L.m, L.p1, L.p2, 1), ar = arc(*Rt.m, Rt.p1, Rt.p2, 1);
                const auto ml = al.value(T(0.5)), mr = ar.value(T(0.5));
                if (!(norm(ml - mr) > tol))
                    return Cut::None;
                const auto c = (L.c + Rt.c) / T(2);
                const T r = norm(L.p1 - L.c);
                const auto q = c + r * unit((ml + mr) / T(2) - c);
                auto covers = [&](T u) {
                    const auto st0 = raw_of(u, R.r, solved, a, b);
                    return norm(q - st0.c) < norm(st0.p1 - st0.c) - T(1e-12) * (T(1) + r);
                };
                const T d = T(1e-4) * (b - a);
                return covers(t - d) || covers(t + d) ? Cut::Crease : Cut::Lune;
            };
            // the candidates: the kinks of the law, the vertex of a closed edge, the knots of the
            // edge below C2; within an open edge (a kink at its vertex: the joint of two blends of
            // a chain, `joint`)
            std::vector<T> candidates = kinks;
            if (st.closed)
                candidates.push_back(ed.u1);
            candidates.insert(candidates.end(), weak2.begin(), weak2.end());
            const T lo_c = st.closed ? a : ed.u1, hi_c = st.closed ? b : ed.u2;
            std::vector<T> cuts{a};
            std::map<T, Cut> kinds;
            Cut at_ends = Cut::None; // (a cut at the start of a loop: between its last piece and its first)
            for (T k : candidates)
            {
                if (st.closed)
                    while (k < a - eps)
                        k += period;
                if (st.closed && (std::abs(k - a) <= eps || std::abs(k - b) <= eps))
                {
                    if (const auto c = kind(a); c != Cut::None)
                        at_ends = c;
                    continue;
                }
                if (!(k > lo_c + eps && k < hi_c - eps) || std::ranges::any_of(cuts, [&](T c) { return std::abs(c - k) <= eps; }))
                    continue;
                if (const auto c = kind(k); c != Cut::None)
                {
                    cuts.push_back(k);
                    kinds[k] = c;
                }
            }
            const bool lune_at_ends = at_ends != Cut::None;
            // (a law of rho that kinks: the pieces of a fillet of constant radius only, `fillet`)
            if (sec.conic && std::ranges::any_of(sec.kinks, [&](T k) { return k > a + eps && k < b - eps; }))
                refuse<T>(BuildErrc::NonGeneric, op,
                          "a conic fillet whose law of rho kinks along a fillet whose radius varies (not handled: a smooth law of rho, or a constant radius)",
                          {eid});
            std::ranges::sort(cuts);
            cuts.push_back(b);
            const T hk = T(1e-7) * (b - a);
            auto slope_before = [&](T t) { return (R.r(t) - R.r(t - hk)) / hk; };
            auto slope_after = [&](T t) { return (R.r(t + hk) - R.r(t)) / hk; };
            auto peak = [&](T t) {
                if (const auto it = kinds.find(t); it != kinds.end())
                    return it->second == Cut::Crease;
                return (std::abs(t - a) <= eps || std::abs(t - b) <= eps) && at_ends == Cut::Crease;
            };
            // each piece continued past a peak at its ends (by 1.5 radius) with its slope there
            const std::size_t np = cuts.size() - 1;
            // pieces with lunes of the ball between them meet arcs of circle only: refused for a
            // conic or G2 section before the pieces are skinned
            if (np > 1 || lune_at_ends)
            {
                if (sec.conic)
                    refuse<T>(BuildErrc::NonGeneric, op, "a conic fillet whose radius law kinks (the lunes of the balls meet arcs of circle only: not handled)", {eid});
                if (sec.g2)
                    refuse<T>(BuildErrc::NonGeneric, op, "a G2 fillet whose radius law kinks (the lunes of the balls meet arcs of circle only: not handled)", {eid});
            }
            std::vector<T> lo_x(np), hi_x(np);
            std::vector<Law> laws(np);
            std::vector<Solved> caches(np);
            for (std::size_t i = 0; i < np; ++i)
            {
                const T c0 = cuts[i], c1 = cuts[i + 1];
                const bool p0 = (i > 0 || (st.closed && lune_at_ends)) && peak(c0), p1 = (i + 1 < np || (st.closed && lune_at_ends)) && peak(c1);
                lo_x[i] = p0 ? std::max(st.closed ? c0 - (b - a) / T(4) : edge_lo, c0 - T(1.5) * R.r(c0) / norm(edge->value(c0, 1))) : c0;
                hi_x[i] = p1 ? std::min(st.closed ? c1 + (b - a) / T(4) : edge_hi, c1 + T(1.5) * R.r(c1) / norm(edge->value(c1, 1))) : c1;
                const T r0 = R.r(c0), r1 = R.r(c1), d0 = slope_after(c0), d1 = slope_before(c1);
                laws[i] = [=, &R](T t) { return t < c0 ? r0 + d0 * (t - c0) : t > c1 ? r1 + d1 * (t - c1) : R.r(t); };
            }
            // past a peak the radius of a piece continued with its slope there grows beyond the
            // peak (0.41 for a peak of 0.3 on the edge of the box): the supports, their reach and
            // the offsets for that radius, the fillet done again (5 % more, the edge then
            // continued a little further)
            {
                T rx{};
                for (std::size_t i = 0; i < np; ++i)
                {
                    if (lo_x[i] < cuts[i])
                        rx = std::max(rx, laws[i](lo_x[i]));
                    if (hi_x[i] > cuts[i + 1])
                        rx = std::max(rx, laws[i](hi_x[i]));
                }
                if (rx > rmax * (T(1) + T(1e-9)))
                    return fillet_law(m, sid, eid, R, o, T(1.05) * rx, sec);
            }
            std::vector<Skinned<T>> sks;
            for (std::size_t i = 0; i < np; ++i)
            {
                const T lo = lo_x[i], hi = hi_x[i];
                std::vector<T> within;
                for (const T t : breaks)
                    if (t > lo + eps && t < hi - eps && t != cuts[i] && t != cuts[i + 1])
                        within.push_back(t);
                if (lo < cuts[i])
                    within.push_back(cuts[i]);
                if (hi > cuts[i + 1])
                    within.push_back(cuts[i + 1]);
                // its stations: those of R.r strictly within it (its law there), its ends solved
                // from them, the piece continued past a peak marched from its ends
                for (auto it = solved.upper_bound(cuts[i] + eps); it != solved.end() && it->first < cuts[i + 1] - eps; ++it)
                    caches[i].insert(*it);
                if (caches[i].empty())
                    prime(laws[i], caches[i], lo, hi);
                else
                {
                    (void)raw_of(cuts[i], laws[i], caches[i], lo, hi);
                    (void)raw_of(cuts[i + 1], laws[i], caches[i], lo, hi);
                    if (lo < cuts[i])
                        march(laws[i], caches[i], cuts[i], lo, lo, hi);
                    if (hi > cuts[i + 1])
                        march(laws[i], caches[i], cuts[i + 1], hi, lo, hi);
                }
                sks.push_back(with_edge(eid, [&] {
                    return skinned([&, lo, hi, i](T t) { return tilted(t, lo, hi, laws[i], caches[i]); },
                                   std::array<const Surface<T, 3> *, 2>{sup[0].surface.get(), sup[1].surface.get()}, lo, hi, o.approx_tol, o.angular_tol, within, sec);
                }));
            }
            // the centre line and the contacts joined across the pieces (C0 at the kinks)
            auto as_bs3 = [](const std::shared_ptr<Curve<T, 3>> &c) { return *std::dynamic_pointer_cast<BSCurve<T, 3>>(c); };
            auto as_bs2 = [](const std::shared_ptr<Curve<T, 2>> &c) { return *std::dynamic_pointer_cast<BSCurve<T, 2>>(c); };
            std::vector<BSCurve<T, 3>> spines;
            std::array<std::vector<BSCurve<T, 2>>, 2> pcs;
            for (std::size_t i = 0; i < np; ++i)
            {
                const auto &sk = sks[i];
                auto sp = as_bs3(sk.spine);
                if (lo_x[i] < cuts[i] || hi_x[i] > cuts[i + 1])
                    sp.trim(cuts[i], cuts[i + 1]);
                spines.push_back(sp);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    auto c = as_bs2(sk.pcurves[k]);
                    if (lo_x[i] < cuts[i] || hi_x[i] > cuts[i + 1])
                        c.trim(cuts[i], cuts[i + 1]);
                    if (!pcs[k].empty())
                    {
                        // on the periods of the face, on from the piece before
                        const auto per = view_periods(*sup[k].surface);
                        const auto prev = pcs[k].back().value(pcs[k].back().bounds()[1]), first = c.value(c.bounds()[0]);
                        point<T, 2> shift{};
                        for (std::size_t d = 0; d < 2; ++d)
                            if (per[d] > T(0))
                                shift[d] = per[d] * std::round((prev[d] - first[d]) / per[d]);
                        translate(c, shift);
                    }
                    pcs[k].push_back(c);
                }
            }
            rb.spine = std::make_shared<BSCurve<T, 3>>(joined(spines));
            for (std::size_t k = 0; k < 2; ++k)
            {
                auto &ct = rb.contacts[k];
                ct.face = st.sides[k].face;
                ct.support = sup[k].surface;
                std::shared_ptr<Curve<T, 2>> p = std::make_shared<BSCurve<T, 2>>(joined(pcs[k]));
                if (const auto *sv = dynamic_cast<const SeamlessSurface<T> *>(sup[k].surface.get()))
                {
                    p = into_bounds(p, *sv, tol);
                    if (!p)
                        refuse<T>(BuildErrc::NonGeneric, op, "a contact crossing the seam of a face that is not periodic (a closed B-spline)", {eid, ct.face});
                }
                ct.pcurve = p;
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, st.sides[k], *ct.pcurve, ed, rb.t1, rb.t2, tol, op, eid);
            }
            rb.circular = !sec.conic && !sec.g2;
            if (sks.size() == 1 && !lune_at_ends)
            {
                rb.surface = sks.front().surface;
                rb.deviation = sks.front().deviation;
                on_blend(rb, std::max(o.intersect.tol, tol / T(10)), true);
                return rb;
            }

            // ---- the pieces and the lunes between them (arcs of circle: conic and G2 sections
            // refused above)
            for (std::size_t i = 0; i < sks.size(); ++i)
            {
                BlendPiece<T> pi;
                pi.surface = sks[i].surface;
                pi.a = cuts[i], pi.b = cuts[i + 1];
                pi.deviation = sks[i].deviation;
                rb.deviation = std::max(rb.deviation, pi.deviation);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    pi.on_blend[k] = segment<T>({T(k), pi.a}, {T(k), pi.b}, pi.a, pi.b);
                    T dev{};
                    for (int j = 0; j <= 64; ++j)
                    {
                        const T t = pi.a + (pi.b - pi.a) * T(j) / T(64);
                        const auto q = pi.on_blend[k]->value(t);
                        dev = std::max(dev, norm(pi.surface->value(q[0], q[1]) - rb.contacts[k].curve->value(t)));
                    }
                    rb.contacts[k].deviation = std::max(rb.contacts[k].deviation, dev);
                }
                rb.pieces.push_back(std::move(pi));
            }
            auto lune = [&](T t, std::size_t before, std::size_t after) {
                constexpr T pi = std::numbers::pi_v<T>;
                BlendLune<T> l;
                l.t = t;
                const auto s0 = raw(t);
                const T r = norm(s0.p1 - s0.c);
                // the end of the piece before and the start of the piece after (at b and a for the
                // lune at the vertex of a closed edge)
                const auto &pb = rb.pieces[before], &pa = rb.pieces[after];
                l.arcs = {std::make_shared<CurveOnSurface<T, 3>>(segment<T>({T(0), pb.b}, {T(1), pb.b}, T(0), T(1)), pb.surface),
                          std::make_shared<CurveOnSurface<T, 3>>(segment<T>({T(0), pa.a}, {T(1), pa.a}, T(0), T(1)), pa.surface)};
                // the sphere, its seam and poles away from the lune: x toward its middle, z across the contacts
                auto x = l.arcs[0]->value(T(0.5)) + l.arcs[1]->value(T(0.5)) - T(2) * s0.c;
                x = x / norm(x);
                auto z = cross(x, s0.p2 - s0.p1);
                z = z / norm(z);
                l.sphere = std::make_shared<SphericalSurface<T>>(r, ax2<T, 3>{s0.c, z, x}, std::array<T, 4>{-pi, pi, -pi / T(2), pi / T(2)});
                for (std::size_t j = 0; j < 2; ++j)
                {
                    const auto [on, gap] = sphere_pcurve(*l.sphere, *l.arcs[j], tol);
                    l.on_sphere[j] = on;
                    l.deviation = std::max(l.deviation, gap);
                }
                rb.lunes.push_back(std::move(l));
            };
            auto crease = [&](T t, std::size_t before, std::size_t after) {
                const auto &pb = rb.pieces[before], &pa = rb.pieces[after];
                auto c = crease_curve(pb.surface, pb.b, pa.surface, pa.a, true, tol);
                if (!c)
                    refuse<T>(BuildErrc::NonGeneric, op, "the fillets on either side of a peak of the radius law do not meet along a crease", {eid});
                c->t = t;
                rb.creases.push_back(std::move(*c));
            };
            auto between = [&](T t, std::size_t before, std::size_t after) {
                if (peak(t))
                    crease(t, before, after);
                else
                    lune(t, before, after);
            };
            for (std::size_t i = 0; i + 1 < rb.pieces.size(); ++i)
                between(rb.pieces[i].b, i, i + 1);
            if (lune_at_ends)
                between(b, rb.pieces.size() - 1, 0);
            return rb;
        }

        /// The arc length along an edge, s(t), and its inverse t(s): interpolated through the
        /// lengths of 64 spans (`abs_curv` gives t(s) the same way), s continued past the vertices
        /// at the speed of the edge there.
        template <std::floating_point T>
        struct ArcLength
        {
            T u1{}, u2{}, L{}, v1{}, v2{};
            BSCfunction<T> s_of, t_of;

            auto s(T t) const -> T { return t <= u1 ? (t - u1) * v1 : t >= u2 ? L + (t - u2) * v2 : s_of(t); }
            /// t(s) on [0, L], put back on s(t) by Newton
            auto t(T s) const -> T
            {
                T x = t_of(std::clamp(s, T(0), L));
                for (int j = 0; j < 8; ++j)
                    x = std::clamp(x - (s_of(x) - s) / s_of(x, 1), u1, u2);
                return x;
            }
        };

        template <std::floating_point T>
        auto arc_length(const Edge<T> &ed) -> ArcLength<T>
        {
            ArcLength<T> a;
            a.u1 = ed.u1, a.u2 = ed.u2;
            a.L = length(*ed.curve, ed.u1, ed.u2);
            a.v1 = norm(ed.curve->value(ed.u1, 1)), a.v2 = norm(ed.curve->value(ed.u2, 1));
            points_vector<T, 1> ts, ss;
            T acc{};
            for (int i = 0; i <= 64; ++i)
            {
                const T t = ed.u1 + (ed.u2 - ed.u1) * T(i) / T(64);
                if (i > 0)
                    acc += length(*ed.curve, ts.back()[0], t);
                ts.push_back({t});
                ss.push_back({acc});
            }
            std::vector<T> tv, sv;
            for (std::size_t i = 0; i < ts.size(); ++i)
            {
                ss[i][0] *= a.L / acc;
                tv.push_back(ts[i][0]), sv.push_back(ss[i][0]);
            }
            a.s_of = BSCfunction<T>{interpolate<T, 1>(ss, tv, std::size_t{3})};
            a.t_of = BSCfunction<T>{interpolate<T, 1>(ts, sv, std::size_t{3})};
            return a;
        }

        /// A law of the arc length S along a chain of edges: f through its points, going on with
        /// its slope past its first and last points (constant for a smooth law, straight on for a
        /// linear one), periodic along a closed chain of length L.
        template <std::floating_point T>
        struct ChainRadius
        {
            BSCfunction<T> f;
            T f1{}, f2{}, d1{}, d2{}, L{};
            bool closed{};
            RadiusVariation variation{};

            auto operator()(T S) const -> T
            {
                if (closed)
                {
                    while (S < f1)
                        S += L;
                    while (S > f1 + L)
                        S -= L;
                }
                return S < f1 ? f(f1) + d1 * (S - f1) : S > f2 ? f(f2) + d2 * (S - f2) : f(S);
            }
            /// its slope at S
            auto slope(T S) const -> T
            {
                if (closed)
                {
                    while (S < f1)
                        S += L;
                    while (S > f1 + L)
                        S -= L;
                }
                return S < f1 ? d1 : S > f2 ? d2 : f(S, 1);
            }
            /// the abscissae where its slope jumps (linear) or where it is no smoother than C1 (smooth)
            auto knots() const -> std::vector<T>
            {
                std::vector<T> r;
                const auto &kf = f.basisCurve().knotsFlats();
                for (std::size_t i = 1; i + 1 < kf.size(); ++i)
                {
                    const T s = kf[i];
                    if (s == kf[i - 1] || s <= f1 || s >= f2)
                        continue;
                    if (variation == RadiusVariation::Linear && !(std::abs(f(s + T(1e-9) * L, 1) - f(s - T(1e-9) * L, 1)) > T(1e-12)))
                        continue;
                    r.push_back(s);
                }
                // where a closed law closes: a kink of a linear law whose slopes differ there, a knot of a smooth one
                if (closed && (variation == RadiusVariation::Smooth || std::abs(d1 - d2) > T(1e-12)))
                    r.push_back(f1);
                return r;
            }
        };

        /**
         * The law of the radius of the edge `ed` of a chain, as a function of the parameter of the
         * edge: S(t) = offset + s(t) along the chain (offset + L - s(t) where the edge runs against
         * it). Past the vertices of an open edge, the law of the edge goes on with its slope there,
         * taken on the edge (a kink of the law at a vertex of a chain: the lune or the crease
         * between the blends of its two edges, decided at their joint); the knots and kinks of the
         * law strictly inside an open edge, its vertex too along a closed edge. The largest radius
         * on the edge: at its vertices or at a knot of the law within it (a linear law is straight
         * between its kinks, a smooth one monotone between its points).
         */
        template <std::floating_point T>
        auto edge_radius_in(const Edge<T> &ed, const ArcLength<T> &al, const ChainRadius<T> &F, T offset, bool forward) -> EdgeRadius<T>
        {
            EdgeRadius<T> R;
            const bool closed = ed.v1 == ed.v2;
            const T u1 = ed.u1, u2 = ed.u2, L = al.L;
            auto at = [=](T s) { return offset + (forward ? s : L - s); };
            auto on_chain = [=](T s) { return F(at(s)); };
            // the slopes dr/ds at the vertices, taken just inside the edge
            const T h = T(1e-9) * L, f1 = on_chain(T(0)), f2 = on_chain(L);
            const T g1 = (forward ? T(1) : T(-1)) * F.slope(at(h)), g2 = (forward ? T(1) : T(-1)) * F.slope(at(L - h));
            R.r = [=](T t) {
                if (closed)
                    return on_chain(al.s(u1 + std::fmod(std::fmod(t - u1, u2 - u1) + (u2 - u1), u2 - u1)));
                const T s = al.s(t);
                return s < T(0) ? f1 + g1 * s : s > L ? f2 + g2 * (s - L) : on_chain(s);
            };
            R.r1 = R.r(u1), R.r2 = R.r(u2);
            R.rmax = std::max(R.r1, R.r2);
            const T margin = T(1e-12) * F.L;
            for (T S : F.knots())
            {
                // onto the edge: S within [offset, offset + L], a period away along a closed chain
                T x = S - offset;
                if (F.closed)
                    x -= F.L * std::floor((x + margin) / F.L);
                if (closed ? !(x >= -margin && x <= L + margin) : !(x > margin && x < L - margin))
                    continue;
                const T t = al.t(forward ? x : L - x);
                (F.variation == RadiusVariation::Linear ? R.kinks : R.knots).push_back(t);
                R.rmax = std::max(R.rmax, F(S));
            }
            return R;
        }

        /// An edge of a chain: the edge, and whether it runs along the chain (from its v1 to its v2).
        struct ChainEdge
        {
            EdgeId edge;
            bool forward{true};
        };

        /**
         * The law with its points given by position (`BlendPoint::where`, stage 8, PR 2,
         * `brep8_pr02_partial_strips_laws.md`, § 4.5) put on the edges `edges`: each at its foot
         * on the nearest of them, at the fraction of the length of that edge there. Refused
         * (`InvalidBounds`): a point farther from every edge than ten times its value and the
         * length of the edge of its foot (a point put elsewhere), two points at the same foot
         * with different values.
         */
        template <std::floating_point T>
        auto law_by_position(const Model<T> &m, std::span<const EdgeId> edges, BlendLaw<T> law, const std::string &op) -> BlendLaw<T>
        {
            if (std::ranges::none_of(law.points, [](const BlendPoint<T> &p) { return p.where.has_value(); }))
                return law;
            std::map<EdgeId, ArcLength<T>> als;
            std::vector<std::pair<point<T, 3>, const BlendPoint<T> *>> feet;
            for (auto &p : law.points)
            {
                if (!p.where)
                    continue;
                const auto &q = *p.where;
                if (!std::ranges::all_of(q, [](T x) { return std::isfinite(x); }))
                    refuse<T>(BuildErrc::InvalidPoint, op, "a point of the law with a coordinate not finite");
                EdgeId best{};
                T bu{}, bd = std::numeric_limits<T>::infinity();
                for (auto e : edges)
                {
                    if (!m.alive(e) || m.edge(e).degenerate || !m.edge(e).curve)
                        continue;
                    const auto &ed = m.edge(e);
                    // the nearest of 64 samples, then Newton within the edge
                    T u = ed.u1, d = std::numeric_limits<T>::infinity();
                    for (int i = 0; i <= 64; ++i)
                    {
                        const T t = ed.u1 + (ed.u2 - ed.u1) * T(i) / T(64);
                        if (const T x = norm(ed.curve->value(t) - q); x < d)
                            d = x, u = t;
                    }
                    try
                    {
                        const auto x = extrema_curve_point(*ed.curve, q, u, T(1e-14) * (T(1) + std::abs(ed.u2 - ed.u1)));
                        const T t = std::clamp(x[0], ed.u1, ed.u2);
                        if (const T y = norm(ed.curve->value(t) - q); y < d)
                            d = y, u = t;
                    }
                    catch (const std::exception &)
                    {
                    }
                    if (d < bd)
                        bd = d, bu = u, best = e;
                }
                if (!best.valid())
                    refuse<T>(BuildErrc::InvalidBounds, op, "a point of the law given by position without an edge to put it on");
                auto it = als.find(best);
                if (it == als.end())
                    it = als.emplace(best, arc_length(m.edge(best))).first;
                const auto &al = it->second;
                if (bd > std::max(T(10) * std::abs(p.value), al.L))
                    refuse<T>(BuildErrc::InvalidBounds, op, "a point of the law given by position away from the edges (farther than ten times its value and the length of the edge)",
                              {best});
                p.edge = best;
                p.at = al.L > T(0) ? std::clamp(al.s(bu) / al.L, T(0), T(1)) : T(0);
                feet.emplace_back(m.edge(best).curve->value(bu), &p);
                p.where.reset();
            }
            for (std::size_t i = 0; i < feet.size(); ++i)
                for (std::size_t j = i + 1; j < feet.size(); ++j)
                    if (norm(feet[i].first - feet[j].first) <= T(1e-9) * (T(1) + norm(feet[i].first)) && feet[i].second->value != feet[j].second->value)
                        refuse<T>(BuildErrc::InvalidBounds, op, "two points of the law given by position at the same foot on the edges, with different values",
                                  {feet[i].second->edge});
            return law;
        }

        /**
         * The radius law along the chain of edges `chain` (in order, from one end to the other, or
         * round a loop: `closed`), from the points of `law` on its edges: the abscissa of a point
         * at the fraction `at` of the length of its edge from the first vertex of the edge; two
         * points at the same place (a vertex of two edges) with the same radius. One `EdgeRadius`
         * per edge; nullopt where all the radii of the chain are the same (a fillet of constant
         * radius). Refused (`InvalidBounds`): no point on the chain, a point out of its edge, a
         * radius not larger than tol, two radii at a point.
         */
        template <std::floating_point T>
        auto chain_radius(const Model<T> &m, const std::vector<ChainEdge> &chain, bool closed, const RadiusLaw<T> &law, T tol, const std::string &op,
                          const std::string &what = "radius") -> std::optional<std::vector<EdgeRadius<T>>>
        {
            std::vector<ArcLength<T>> als;
            std::vector<T> offsets;
            T total{};
            for (const auto &ce : chain)
            {
                if (!m.alive(ce.edge))
                    refuse<T>(BuildErrc::InvalidId, op, "dead or invalid edge", {ce.edge});
                const auto &ed = m.edge(ce.edge);
                if (ed.degenerate || !ed.curve)
                    refuse<T>(BuildErrc::InvalidId, op, "a degenerate edge (no curve)", {ce.edge});
                als.push_back(arc_length(ed));
                offsets.push_back(total);
                total += als.back().L;
            }
            std::vector<std::pair<T, T>> pts;
            for (const auto &p : law.points)
                for (std::size_t i = 0; i < chain.size(); ++i)
                {
                    if (p.edge != chain[i].edge)
                        continue;
                    if (!(p.at >= T(0) && p.at <= T(1)))
                        refuse<T>(BuildErrc::InvalidBounds, op, "a point of the law out of its edge (at not in [0, 1])", {p.edge});
                    if (!(p.value > tol) || !std::isfinite(p.value))
                        refuse<T>(BuildErrc::InvalidBounds, op, "a " + what + " of the law not larger than the tolerance", {p.edge});
                    T S = offsets[i] + (chain[i].forward ? p.at : T(1) - p.at) * als[i].L;
                    if (closed && S >= total * (T(1) - T(1e-12)))
                        S = T(0);
                    pts.emplace_back(S, p.value);
                }
            std::vector<ShapeId> ids;
            for (const auto &ce : chain)
                ids.push_back(ce.edge);
            if (pts.empty())
                refuse<T>(BuildErrc::InvalidBounds, op, chain.size() > 1 ? "no " + what + " on a chain of edges" : "no " + what + " on the edge", ids);
            std::ranges::sort(pts);
            std::vector<std::pair<T, T>> u;
            for (const auto &p : pts)
            {
                if (!u.empty() && std::abs(u.back().first - p.first) <= T(1e-12) * total)
                {
                    if (u.back().second != p.second)
                        refuse<T>(BuildErrc::InvalidBounds, op,
                                  "two values of the " + what + " at " + (closed && chain.size() == 1 && p.first == T(0) ? "the vertex of a closed edge" : "the same point of the edge"),
                                  ids);
                    continue;
                }
                u.push_back(p);
            }
            if (std::ranges::all_of(u, [&](const auto &p) { return p.second == u.front().second; }))
                return std::nullopt;
            ChainRadius<T> F;
            F.f = radius_detail::radius_function(u, law.variation, closed, total);
            std::tie(F.f1, F.f2) = std::pair{F.f.bounds()[0], F.f.bounds()[1]};
            F.d1 = F.f(F.f1, 1), F.d2 = F.f(F.f2, 1);
            F.L = total;
            F.closed = closed;
            F.variation = law.variation;
            std::vector<EdgeRadius<T>> r;
            for (std::size_t i = 0; i < chain.size(); ++i)
            {
                r.push_back(edge_radius_in(m.edge(chain[i].edge), als[i], F, offsets[i], chain[i].forward));
                // the chain going on before and after the edge, at its v1 and v2 as it runs
                const bool before = closed || i > 0, after = closed || i + 1 < chain.size();
                r.back().joined = chain[i].forward ? std::array<bool, 2>{before, after} : std::array<bool, 2>{after, before};
            }
            return r;
        }

        /// The radius law on the edge e alone (its points on e): `chain_radius` of the chain of e.
        template <std::floating_point T>
        auto edge_radius(const Model<T> &m, EdgeId eid, const RadiusLaw<T> &law, T tol, const std::string &op, const std::string &what = "radius")
            -> std::optional<EdgeRadius<T>>
        {
            if (!m.alive(eid))
                refuse<T>(BuildErrc::InvalidId, op, "dead or invalid edge", {eid});
            const auto &ed = m.edge(eid);
            const auto r = chain_radius(m, std::vector<ChainEdge>{{eid, true}}, ed.v1 == ed.v2, law_by_position(m, std::span<const EdgeId>{&eid, 1}, law, op), tol, op, what);
            if (!r)
                return std::nullopt;
            return r->front();
        }

        /**
         * The value of `law` on the edge e alone where it is constant (`edge_radius` gives none):
         * that of its points on e, those given by position put there first (`law_by_position`;
         * the law as given may have no point on e).
         */
        template <std::floating_point T>
        auto edge_value(const Model<T> &m, EdgeId eid, const BlendLaw<T> &law, const std::string &op) -> T
        {
            const auto placed = law_by_position(m, std::span<const EdgeId>{&eid, 1}, law, op);
            const auto it = std::ranges::find(placed.points, eid, &BlendPoint<T>::edge);
            if (it == placed.points.end())
                refuse<T>(BuildErrc::InvalidBounds, op, "no value of the law on the edge", {eid});
            return it->value;
        }

        // ---- chamfers

        /**
         * The edge continued by ext past each end, on its parameter [lo, hi]: a line or a circle
         * on a wider range (an arc short of a turn), a B-spline by the polynomials of its end
         * spans (`extrapolated`); a curve on the surface of one of its faces (the edges of a
         * boolean) by its pcurve continued the same way, on the support of that face (the surface
         * continued); any other curve is the edge on [u1, u2] joined at its ends to the pcurve of
         * a face continued (the fewest poles, the lowest degree) on its support. The edge as it
         * is, on [u1, u2], when it cannot be continued.
         */
        template <std::floating_point T>
        auto continued_edge(const Edge<T> &ed, T ext, const std::array<Side<T>, 2> &sides, const std::array<std::shared_ptr<Surface<T, 3>>, 2> &sup, T tol,
                            T &lo, T &hi) -> std::shared_ptr<Curve<T, 3>>
        {
            lo = ed.u1, hi = ed.u2;
            if (!(ext > T(0)))
                return ed.curve;
            T d1 = ext / norm(ed.curve->value(ed.u1, 1)), d2 = ext / norm(ed.curve->value(ed.u2, 1));
            if (const auto *ec = dynamic_cast<const ElementaryCurve<T> *>(ed.curve.get()))
            {
                if (dynamic_cast<const EllipseCurve<T> *>(ec))
                {
                    // at most half of the rest of the turn
                    const T room = T(2) * std::numbers::pi_v<T> - (ed.u2 - ed.u1);
                    if (!(room > T(0)))
                        return ed.curve;
                    const T k = std::min(T(1), room / (T(2) * (d1 + d2)));
                    d1 *= k, d2 *= k;
                }
                lo -= d1, hi += d2;
                return intersect_detail::restricted(*ec, lo, hi);
            }
            // a B-spline (in 3D, or in 2D for a pcurve) gone on by the polynomials of its end spans over [u1 - d1, u2 + d2]
            auto go_on = [&]<std::size_t dim, bool rational>(const BSCurveGeneral<T, dim, rational> &b) -> std::shared_ptr<Curve<T, dim>> {
                const auto bb = b.bounds();
                const std::array<T, 2> dt{std::max(T(0), bb[0] - (ed.u1 - d1)), std::max(T(0), (ed.u2 + d2) - bb[1])};
                try
                {
                    return std::make_shared<std::conditional_t<rational, BSCurveRational<T, dim>, BSCurve<T, dim>>>(extrapolated(b, dt));
                }
                catch (const std::invalid_argument &)
                {
                    return nullptr;
                }
            };
            std::shared_ptr<Curve<T, 3>> x;
            if (const auto *b = dynamic_cast<const BSCurve<T, 3> *>(ed.curve.get()))
                x = go_on(*b);
            else if (const auto *b = dynamic_cast<const BSCurveRational<T, 3> *>(ed.curve.get()))
                x = go_on(*b);
            if (!x)
            {
                // the pcurve of the edge itself when it lies on the surface of a face (its pcurve on the support of
                // that face runs through the edge within tol at its ends and middle); otherwise the pcurve of a face
                // with the fewest poles, of the lowest degree (the polynomial of an end span gone on the furthest)
                auto poles = [](const Curve<T, 2> &p) {
                    if (const auto *b = dynamic_cast<const BSCurve<T, 2> *>(&p))
                        return std::pair{b->degree(), b->poles().size()};
                    if (const auto *b = dynamic_cast<const BSCurveRational<T, 2> *>(&p))
                        return std::pair{b->degree(), b->poles().size()};
                    return std::pair{std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()};
                };
                std::size_t k = poles(*sides[1].pcurve) < poles(*sides[0].pcurve) ? 1 : 0;
                std::shared_ptr<Curve<T, 2>> base = sides[k].pcurve;
                bool own = false;
                if (const auto *cs = dynamic_cast<const CurveOnSurface<T, 3> *>(ed.curve.get()))
                    for (std::size_t j = 0; j < 2 && !own; ++j)
                        if (std::ranges::all_of(std::array{ed.u1, (ed.u1 + ed.u2) / T(2), ed.u2}, [&](T u) {
                                const auto uv = cs->p_basisCurve()->value(u);
                                return norm(sup[j]->value(uv[0], uv[1]) - ed.curve->value(u)) <= tol;
                            }))
                            base = cs->p_basisCurve(), k = j, own = true;
                std::shared_ptr<Curve<T, 2>> pc;
                if (const auto *p = dynamic_cast<const BSCurve<T, 2> *>(base.get()))
                    pc = go_on(*p);
                else if (const auto *p = dynamic_cast<const BSCurveRational<T, 2> *>(base.get()))
                    pc = go_on(*p);
                if (!pc)
                    return ed.curve;
                const auto on = std::make_shared<CurveOnSurface<T, 3>>(pc, sup[k]);
                try
                {
                    // within the support continued
                    (void)on->value(ed.u1 - d1), (void)on->value(ed.u2 + d2);
                }
                catch (const std::exception &)
                {
                    return ed.curve;
                }
                if (own)
                    x = on; // the edge itself on [u1, u2]
                else
                    x = std::make_shared<CurveComposite<T, 3>>(std::vector<std::shared_ptr<Curve<T, 3>>>{
                        std::make_shared<CurveTrimmed<T, 3>>(on, ed.u1 - d1, ed.u1), std::make_shared<CurveTrimmed<T, 3>>(ed.curve, ed.u1, ed.u2),
                        std::make_shared<CurveTrimmed<T, 3>>(on, ed.u2, ed.u2 + d2)});
            }
            lo -= d1, hi += d2;
            return x;
        }

        /**
         * The point of `s` at the distance d from e(t) in the plane across the curve e at t,
         * from uv (Gauss-Newton on (|S - e|^2 - d^2) / 2d and (S - e) . e' / |e'|); nullopt
         * when it does not converge within tol / 100 or lands on the far side of the edge
         * ((S - e) . side <= 0, side: toward the face from the edge).
         */
        template <std::floating_point T>
        auto foot(const Surface<T, 3> &s, const Curve<T, 3> &e, T t, T d, point<T, 2> uv, const point<T, 3> &side, T tol) -> std::optional<point<T, 2>>
        {
            using thick_detail::Vec;
            using thick_detail::Mat;
            const auto p = e.value(t);
            const auto te = unit(e.value(t, 1));
            Vec<T> x(2);
            x << uv[0], uv[1];
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &r, Mat<T> &J) {
                const auto q = s.value(y[0], y[1]) - p;
                const auto su = s.value(y[0], y[1], 1, 0), sv = s.value(y[0], y[1], 0, 1);
                r.resize(2);
                J.resize(2, 2);
                r << (q * q - d * d) / (T(2) * d), q * te;
                J << (q * su) / d, (q * sv) / d, su * te, sv * te;
            });
            if (!(res <= tol / T(100)))
                return std::nullopt;
            if (!((s.value(x[0], x[1]) - p) * side > T(0)))
                return std::nullopt;
            return point<T, 2>{x[0], x[1]};
        }

        /**
         * The contact at the distance d(t) on a side, on [lo, hi]: the feet at n + 1 stations, each
         * from the one before (the first, at the middle of the edge, from its pcurve moved by d
         * toward the face), the move of the edge brought onto the tangent plane, the step halved
         * where it does not converge; interpolated (cubic); n doubled from 8 until the pcurve is
         * within `acc` of the feet at the middles, which are then the new stations.
         */
        template <std::floating_point T>
        auto contact_pcurve(const Side<T> &sd, const Surface<T, 3> &s, const Curve<T, 3> &e, const Edge<T> &ed, T lo, T hi,
                            const std::function<T(T)> &dist, T acc, T tol, const std::string &op, EdgeId eid) -> std::shared_ptr<Curve<T, 2>>
        {
            // the seed: the pcurve of the edge at its middle moved by d along the face, away from it
            const T tm = (ed.u1 + ed.u2) / T(2), d = dist(tm);
            const auto uv0 = sd.pcurve->value(tm);
            const auto su = s.value(uv0[0], uv0[1], 1, 0), sv = s.value(uv0[0], uv0[1], 0, 1);
            auto dd = sd.pcurve->value(tm, 1);
            if (!sd.forward)
                dd = T(-1) * dd;
            const point<T, 2> left{-dd[1], dd[0]};
            auto b = left[0] * su + left[1] * sv;
            const auto te = unit(e.value(tm, 1));
            b = unit(b - (b * te) * te);
            // (u, v) of d b on the tangent plane, by the metric
            const T E = su * su, F = su * sv, G = sv * sv, det = E * G - F * F;
            const T bu = d * (G * (b * su) - F * (b * sv)) / det, bv = d * (E * (b * sv) - F * (b * su)) / det;
            const auto first = foot(s, e, tm, d, point<T, 2>{uv0[0] + bu, uv0[1] + bv}, b, tol);
            if (!first || (s.value((*first)[0], (*first)[1]) - e.value(tm)) * b < T(0.5) * d)
                refuse<T>(BuildErrc::NonGeneric, op, "no point of a face at the distance asked from the edge", {eid, sd.face});
            // the foot at t from the one at t0: the move of the edge brought onto the tangent plane at
            // the foot (first-order predictor), the step halved where Gauss-Newton does not converge,
            // each foot on the side of the face of the one before
            auto predicted = [&](T t, T t0, const point<T, 2> &uv) {
                const auto a = s.value(uv[0], uv[1], 1, 0), c = s.value(uv[0], uv[1], 0, 1);
                const auto dp = e.value(t) - e.value(t0);
                const T E = a * a, F = a * c, G = c * c, det = E * G - F * F;
                if (!(std::abs(det) > T(0)))
                    return uv;
                return point<T, 2>{uv[0] + (G * (dp * a) - F * (dp * c)) / det, uv[1] + (E * (dp * c) - F * (dp * a)) / det};
            };
            auto march = [&](auto &&self, T t, T t0, const point<T, 2> &uv0, int depth) -> std::optional<point<T, 2>> {
                const auto side = s.value(uv0[0], uv0[1]) - e.value(t0);
                if (auto x = foot(s, e, t, dist(t), predicted(t, t0, uv0), side, tol))
                    return x;
                if (depth == 0)
                    return std::nullopt;
                const T tm = (t + t0) / T(2);
                const auto m = self(self, tm, t0, uv0, depth - 1);
                return m ? self(self, t, tm, *m, depth - 1) : std::nullopt;
            };
            auto solve = [&](T t, T t0, const point<T, 2> &seed) {
                const auto x = march(march, t, t0, seed, 10);
                if (!x)
                    refuse<T>(BuildErrc::NonGeneric, op, "no point of a face at the distance asked from the edge", {eid, sd.face});
                return *x;
            };
            std::size_t n = 8;
            std::vector<T> t(n + 1);
            std::vector<point<T, 2>> uv(n + 1);
            for (std::size_t i = 0; i <= n; ++i)
                t[i] = i == n ? hi : lo + (hi - lo) * T(i) / T(n);
            // from the station nearest the middle of the edge, both ways
            const std::size_t i0 = std::size_t(std::clamp<T>(std::round((tm - lo) / (hi - lo) * T(n)), T(0), T(n)));
            uv[i0] = solve(t[i0], tm, *first);
            for (std::size_t i = i0 + 1; i <= n; ++i)
                uv[i] = solve(t[i], t[i - 1], uv[i - 1]);
            for (std::size_t i = i0; i-- > 0;)
                uv[i] = solve(t[i], t[i + 1], uv[i + 1]);
            for (;;)
            {
                std::vector<constrType<T, 2, 1>> q(n + 1);
                for (std::size_t i = 0; i <= n; ++i)
                    q[i] = {uv[i]};
                auto pc = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(q, t, std::size_t{3}));
                // the feet at the middles: the check, and the new stations
                std::vector<T> tt(n);
                std::vector<point<T, 2>> mid(n);
                T dev{};
                for (std::size_t i = 0; i < n; ++i)
                {
                    tt[i] = (t[i] + t[i + 1]) / T(2);
                    mid[i] = solve(tt[i], t[i], uv[i]);
                    const auto y = pc->value(tt[i]);
                    dev = std::max(dev, norm(s.value(y[0], y[1]) - s.value(mid[i][0], mid[i][1])));
                }
                if (dev <= acc)
                    return pc;
                if (2 * n > 4096)
                    break;
                std::vector<T> t2(2 * n + 1);
                std::vector<point<T, 2>> uv2(2 * n + 1);
                for (std::size_t i = 0; i <= n; ++i)
                {
                    t2[2 * i] = t[i], uv2[2 * i] = uv[i];
                    if (i < n)
                        t2[2 * i + 1] = tt[i], uv2[2 * i + 1] = mid[i];
                }
                t = std::move(t2), uv = std::move(uv2), n *= 2;
            }
            refuse<T>(BuildErrc::NonGeneric, op, "a contact cannot be approximated within the tolerance", {eid, sd.face});
        }

        /// The contact at the constant distance d (`contact_pcurve`).
        template <std::floating_point T>
        auto contact_pcurve(const Side<T> &sd, const Surface<T, 3> &s, const Curve<T, 3> &e, const Edge<T> &ed, T lo, T hi, T d, T acc, T tol,
                            const std::string &op, EdgeId eid) -> std::shared_ptr<Curve<T, 2>>
        {
            return contact_pcurve(sd, s, e, ed, lo, hi, std::function<T(T)>{[d](T) { return d; }}, acc, tol, op, eid);
        }

        /// The circle through p, q, w: centre, unit normal, radius; nullopt for aligned points.
        template <std::floating_point T>
        auto circle3(const point<T, 3> &p, const point<T, 3> &q, const point<T, 3> &w) -> std::optional<std::tuple<point<T, 3>, point<T, 3>, T>>
        {
            const auto a = q - p, b = w - p;
            const auto n = cross(a, b);
            const T nn = n * n;
            if (!(nn > T(1e-24) * (a * a) * (b * b)))
                return std::nullopt;
            const auto c = p + ((b * b) * cross(n, a) + (a * a) * cross(b, n)) / (T(2) * nn);
            return std::tuple{c, n / std::sqrt(nn), norm(p - c)};
        }

        /**
         * The ruled surface between the contacts c1 and c2 on [a, b] as a primitive: a plane
         * between two lines in a plane, a plane, a cylinder or a cone between two coaxial circles;
         * nullptr otherwise. Tested at 34 points within tol, away from the stations of the
         * contacts (a power of 2) but at the ends; a circle is the one through the stations at
         * a, a + (b - a) / 4 and (a + b) / 2, where the contacts are the exact feet.
         */
        template <std::floating_point T>
        auto ruled_primitive(const Curve<T, 3> &c1, const Curve<T, 3> &c2, T a, T b, T tol) -> std::shared_ptr<ElementarySurface<T>>
        {
            constexpr int n = 33;
            std::array<std::vector<point<T, 3>>, 2> P;
            for (int i = 0; i <= n; ++i)
            {
                const T s = i == n ? b : a + (b - a) * T(i) / T(n);
                P[0].push_back(c1.value(s));
                P[1].push_back(c2.value(s));
            }
            auto plane_box = [&](const Frame3<T> &f) {
                std::array<T, 4> box{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(),
                                     -std::numeric_limits<T>::infinity()};
                for (const auto &Q : P)
                    for (const auto &p : Q)
                    {
                        const T u = (p - f.o) * f.x, v = (p - f.o) * f.y;
                        box = {std::min(box[0], u), std::max(box[1], u), std::min(box[2], v), std::max(box[3], v)};
                    }
                return box;
            };
            // two lines in a plane
            auto straight = [&](const std::vector<point<T, 3>> &Q) {
                const auto d = Q.back() - Q.front();
                if (!(norm(d) > tol))
                    return false;
                const auto u = unit(d);
                return std::ranges::all_of(Q, [&](const point<T, 3> &p) { return norm(cross(p - Q.front(), u)) <= tol; });
            };
            if (straight(P[0]) && straight(P[1]))
            {
                const auto x = unit(P[0].back() - P[0].front());
                const auto z = cross(x, P[1].front() - P[0].front()), z2 = cross(x, P[1].back() - P[0].front());
                const auto nz = norm(z) >= norm(z2) ? z : z2;
                if (!(norm(nz) > tol))
                    return nullptr;
                const Frame3<T> f{ax2<T, 3>{P[0].front(), unit(nz), x}};
                if (std::ranges::all_of(P[1], [&](const point<T, 3> &p) { return std::abs((p - f.o) * f.z) <= tol; }))
                    return std::make_shared<Plane<T>>(f, plane_box(f));
                return nullptr;
            }
            // two coaxial circles
            std::array<std::tuple<point<T, 3>, point<T, 3>, T>, 2> circ;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &ck_curve = k == 0 ? c1 : c2;
                const auto ck = circle3(ck_curve.value(a), ck_curve.value(a + (b - a) / T(4)), ck_curve.value((a + b) / T(2)));
                if (!ck)
                    return nullptr;
                const auto &[cc, cn, cr] = *ck;
                for (const auto &p : P[k])
                    if (std::abs(norm(p - cc) - cr) > tol || std::abs((p - cc) * cn) > tol)
                        return nullptr;
                circ[k] = *ck;
            }
            const auto &[o1, n1, r1] = circ[0];
            const auto &[o2, n2, r2] = circ[1];
            if (norm(cross(n1, n2)) > T(1e-9) || norm(cross(o2 - o1, n1)) > tol)
                return nullptr;
            const T h = (o2 - o1) * n1;
            const auto z = h >= T(0) ? n1 : T(-1) * n1;
            if (std::abs(h) <= tol)
            {
                const Frame3<T> f{ax2<T, 3>{o1, z, P[0].front() - o1}};
                return std::make_shared<Plane<T>>(f, plane_box(f));
            }
            // a cone from its smaller circle, its semi-angle positive (as STEP has it)
            const bool shrinks = std::abs(r1 - r2) > tol && r2 < r1;
            const auto &ob = shrinks ? o2 : o1;
            const Frame3<T> f{ax2<T, 3>{ob, shrinks ? T(-1) * z : z, P[shrinks ? 1 : 0].front() - ob}};
            // the angles of the contacts, continuous along them
            T u1 = std::numeric_limits<T>::infinity(), u2 = -u1;
            for (const auto &Q : P)
            {
                T prev{};
                for (std::size_t i = 0; i < Q.size(); ++i)
                {
                    const auto q = Q[i] - f.o;
                    T u = std::atan2(q * f.y, q * f.x);
                    if (i > 0)
                        u += T(2) * std::numbers::pi_v<T> * std::round((prev - u) / (T(2) * std::numbers::pi_v<T>));
                    prev = u;
                    u1 = std::min(u1, u), u2 = std::max(u2, u);
                }
            }
            const std::array<T, 4> box{u1, u2, T(0), std::abs(h)};
            if (std::abs(r1 - r2) <= tol)
                return std::make_shared<CylindricalSurface<T>>(r1, f, box);
            return std::make_shared<ConicalSurface<T>>(std::min(r1, r2), std::atan(std::abs(r2 - r1) / std::abs(h)), f, box);
        }

        /**
         * The pcurve pc on the surface s as a B-spline of space, exactly, where it is one: a
         * B-spline pcurve on a plane (the plane affine in its parameters, the poles mapped by it);
         * nullopt otherwise.
         */
        template <std::floating_point T>
        auto on_plane(const Curve<T, 2> &pc, const Surface<T, 3> &s) -> std::optional<BSCurve<T, 3>>
        {
            const auto *pl = dynamic_cast<const Plane<T> *>(&s);
            const auto *b = dynamic_cast<const BSCurve<T, 2> *>(&pc);
            if (!pl || !b)
                return std::nullopt;
            points_vector<T, 3> poles;
            for (const auto &q : b->poles())
                poles.push_back(pl->value(q[0], q[1]));
            return BSCurve<T, 3>(poles, b->knotsFlats(), b->degree());
        }

        /// The curve c on [a, b] as a B-spline, exactly, where it is one: a B-spline pcurve on [a, b] on a plane (`on_plane`).
        template <std::floating_point T>
        auto exact_spline(const Curve<T, 3> &c, T a, T b) -> std::optional<BSCurve<T, 3>>
        {
            const auto *cs = dynamic_cast<const CurveOnSurface<T, 3> *>(&c);
            if (!cs)
                return std::nullopt;
            const auto bd = cs->basisCurve().bounds();
            const T eps = T(1e-12) * (T(1) + std::abs(b - a));
            if (std::abs(bd[0] - a) > eps || std::abs(bd[1] - b) > eps)
                return std::nullopt;
            return on_plane(cs->basisCurve(), cs->basisSurface());
        }

        /// The ruled B-spline between the B-splines c1 (v = 0) and c2 (v = 1) of space, exactly (their degrees and knots unified).
        template <std::floating_point T>
        auto ruled_exact(const BSCurve<T, 3> &c1, const BSCurve<T, 3> &c2) -> std::shared_ptr<BSSurface<T, 3>>
        {
            std::vector<BSCurve<T, 3>> both{c1, c2};
            unify_curves(both);
            auto poles = both[0].poles();
            const auto p2 = both[1].poles();
            poles.insert(poles.end(), p2.begin(), p2.end());
            return std::make_shared<BSSurface<T, 3>>(poles, both[0].knotsFlats(), std::vector<T>{0, 0, 1, 1}, both[0].degree(), std::size_t{1});
        }

        /**
         * The ruled B-spline between c1 (v = 0) and c2 (v = 1), u along on [a, b]: both
         * interpolated at the same stations, doubled from 8 (the stations kept, their middles
         * added) until within approx_tol at three points between each pair.
         */
        template <std::floating_point T>
        auto ruled_spline(const Curve<T, 3> &c1, const Curve<T, 3> &c2, T a, T b, T approx_tol, T &deviation) -> std::shared_ptr<Surface<T, 3>>
        {
            // two contacts that are B-splines exactly (on planes): the ruled surface between them, exact
            if (auto e1 = exact_spline(c1, a, b), e2 = exact_spline(c2, a, b); e1 && e2)
            {
                deviation = T(0);
                return ruled_exact(*e1, *e2);
            }
            std::size_t n = 8;
            std::vector<T> t(n + 1);
            std::array<points_vector<T, 3>, 2> q;
            for (std::size_t i = 0; i <= n; ++i)
            {
                t[i] = i == n ? b : a + (b - a) * T(i) / T(n);
                q[0].push_back(c1.value(t[i]));
                q[1].push_back(c2.value(t[i]));
            }
            for (;;)
            {
                const auto b1 = interpolate<T, 3>(q[0], t, std::size_t{3}), b2 = interpolate<T, 3>(q[1], t, std::size_t{3});
                T dev{};
                for (std::size_t i = 0; i < n && dev <= approx_tol; ++i)
                    for (int j = 1; j < 4; ++j)
                    {
                        const T s = t[i] + (t[i + 1] - t[i]) * T(j) / T(4);
                        dev = std::max({dev, norm(b1.value(s) - c1.value(s)), norm(b2.value(s) - c2.value(s))});
                    }
                if (dev <= approx_tol)
                {
                    deviation = dev;
                    auto poles = b1.poles();
                    const auto p2 = b2.poles();
                    poles.insert(poles.end(), p2.begin(), p2.end());
                    return std::make_shared<BSSurface<T, 3>>(poles, b1.knotsFlats(), std::vector<T>{0, 0, 1, 1}, b1.degree(), std::size_t{1});
                }
                if (2 * n > 2048)
                    break;
                std::vector<T> t2(2 * n + 1);
                std::array<points_vector<T, 3>, 2> q2{points_vector<T, 3>(2 * n + 1), points_vector<T, 3>(2 * n + 1)};
                for (std::size_t i = 0; i <= n; ++i)
                {
                    t2[2 * i] = t[i], q2[0][2 * i] = q[0][i], q2[1][2 * i] = q[1][i];
                    if (i < n)
                    {
                        const T s = (t[i] + t[i + 1]) / T(2);
                        t2[2 * i + 1] = s, q2[0][2 * i + 1] = c1.value(s), q2[1][2 * i + 1] = c2.value(s);
                    }
                }
                t = std::move(t2), q = std::move(q2), n *= 2;
            }
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "chamfer_geometry: the surface cannot be approximated within approx_tol", {}});
        }

        template <std::floating_point T>
        auto chamfer(const Model<T> &m, SolidId sid, EdgeId eid, FaceId fid, T d1, T d2, const FilletOptions<T> &o) -> EdgeBlend<T>
        {
            const std::string op = "chamfer_geometry";
            check_options(o, op, eid);
            if (!(d1 > o.tol) || !(d2 > o.tol) || !std::isfinite(d1) || !std::isfinite(d2))
                refuse<T>(BuildErrc::InvalidBounds, op, "the distances must be larger than the tolerance", {eid});
            auto st = setup(m, sid, eid, op);
            // d1 on the face fid: its side first
            if (st.sides[1].face == fid)
                std::swap(st.sides[0], st.sides[1]);
            else if (st.sides[0].face != fid)
                refuse<T>(BuildErrc::InvalidId, op, "the face of d1 is not a face of the edge", {eid, fid});
            const auto &ed = m.edge(eid);
            const T tol = o.tol, io_tol = std::max(o.intersect.tol, tol / T(10));
            const T ext = st.closed ? T(0) : o.extension * std::max(d1, d2);
            const std::array<T, 2> d{d1, d2};
            std::array<std::shared_ptr<Surface<T, 3>>, 2> sup;
            for (std::size_t k = 0; k < 2; ++k)
                sup[k] = support(m, st.sides[k], ed, T(2) * d[k] + ext, tol).surface;
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.spine = continued_edge(ed, ext, st.sides, sup, tol, rb.a, rb.b);
            rb.t1 = ed.u1, rb.t2 = ed.u2;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &sd = st.sides[k];
                auto &ct = rb.contacts[k];
                ct.face = sd.face;
                ct.support = sup[k];
                // the contacts as accurate as those of a fillet on primitives
                ct.pcurve = contact_pcurve(sd, *ct.support, *rb.spine, ed, rb.a, rb.b, d[k], std::min(o.approx_tol, io_tol), tol, op, eid);
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, sd, *ct.pcurve, ed, ed.u1, ed.u2, tol, op, eid);
            }
            // a primitive when the contacts lie on it within tol, a ruled B-spline otherwise
            if (auto e = ruled_primitive(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, tol))
            {
                rb.surface = e;
                on_blend(rb, io_tol, false);
                rb.exact = std::ranges::all_of(rb.contacts, [&](const BlendContact<T> &c) { return c.deviation <= tol; });
            }
            if (!rb.exact)
            {
                rb.surface = ruled_spline(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, o.approx_tol, rb.deviation);
                on_blend(rb, io_tol, false);
                // (built exactly between B-spline contacts on planes)
                rb.exact = rb.deviation == T(0) && std::ranges::all_of(rb.contacts, [&](const BlendContact<T> &c) { return c.deviation <= tol; });
            }
            return rb;
        }

        /// A distance of a chamfer along its edge (stage 7, § 3): d(t) at the parameter t of the edge (continued past its vertices), where its slope jumps, its largest value on the edge.
        template <std::floating_point T>
        struct EdgeDistance
        {
            std::function<T(T)> d;
            std::vector<T> kinks;
            T dmax{};
        };

        /// The distance of a law on an edge (`edge_radius`, its knots smooth): the same function.
        template <std::floating_point T>
        auto edge_distance(const EdgeRadius<T> &R) -> EdgeDistance<T>
        {
            return {R.r, R.kinks, R.rmax};
        }

        /// The same distance d along the edge.
        template <std::floating_point T>
        auto edge_distance(T d) -> EdgeDistance<T>
        {
            return {[d](T) { return d; }, {}, d};
        }

        /**
         * The chamfer whose distances vary along the edge (stage 7, § 3; note
         * `brep7_pr01_variable_chamfers.md`): its contacts at the distances d1(t) on the face fid
         * and d2(t) on the other in the plane normal to the edge (`contact_pcurve`), or, with no
         * d2, the other contact where the line from the first at the angle `angle` from the
         * tangent plane of its face there, in that plane, into the solid along a convex edge and
         * out of it along a concave one, meets the other face (Gauss-Newton at each station, from
         * the station before). The surface ruled between the contacts: a primitive where they
         * allow it, a ruled B-spline otherwise, exact between B-spline contacts on planes (a
         * linear law along a straight edge: a hyperbolic paraboloid; `exact` then). Where the
         * slope of a law jumps (a point of a linear law), the chamfer comes in pieces, ruled
         * each, a crease between them along the ruling there (`EdgeBlend::pieces`, `creases`),
         * the contacts joined across it. Refused: a distance not beyond tol on the chamfer
         * continued, an angle out of ]0, 90[ degrees, a line at the angle that does not meet the
         * other face; the refusals of `chamfer`.
         */
        template <std::floating_point T>
        auto chamfer_law(const Model<T> &m, SolidId sid, EdgeId eid, FaceId fid, const EdgeDistance<T> &D1,
                         const std::type_identity_t<std::optional<EdgeDistance<T>>> &D2, std::type_identity_t<T> angle, const FilletOptions<T> &o) -> EdgeBlend<T>
        {
            const std::string op = "chamfer_geometry";
            check_options(o, op, eid);
            constexpr T pi = std::numbers::pi_v<T>;
            if (!D2 && !(angle > T(0) && angle < pi / T(2)))
                refuse<T>(BuildErrc::InvalidBounds, op, "the angle of a chamfer must be within ]0, 90[ degrees", {eid});
            auto st = setup(m, sid, eid, op);
            if (st.sides[1].face == fid)
                std::swap(st.sides[0], st.sides[1]);
            else if (st.sides[0].face != fid)
                refuse<T>(BuildErrc::InvalidId, op, "the face of d1 is not a face of the edge", {eid, fid});
            const auto &ed = m.edge(eid);
            const T tol = o.tol, io_tol = std::max(o.intersect.tol, tol / T(10));
            // the other distance at an angle, for the reach of the supports and the extension: by
            // the law of sines in the triangle of the edge and the contacts, the faces their
            // tangent planes, the wedge at its widest along the edge (a sine under 0.05: twenty
            // times d1), half as much again for their curvature; twice d1 max(1, tan) at least
            const T dm1 = D1.dmax;
            T dm2 = D2 ? D2->dmax : dm1 * T(2) * std::max(T(1), std::tan(angle));
            if (!D2)
            {
                const T sn = std::sin(st.max_angle + angle);
                dm2 = std::max(dm2, sn > T(0.05) ? T(1.5) * dm1 * std::sin(angle) / sn : T(20) * dm1);
            }
            if (!(dm1 > tol) || !(dm2 > tol))
                refuse<T>(BuildErrc::InvalidBounds, op, "the distances must be larger than the tolerance", {eid});
            const T ext = st.closed ? T(0) : o.extension * std::max(dm1, dm2);
            std::array<std::shared_ptr<Surface<T, 3>>, 2> sup;
            sup[0] = support(m, st.sides[0], ed, T(2) * (dm1 + ext), tol).surface;
            sup[1] = support(m, st.sides[1], ed, T(2) * (dm2 + ext), tol).surface;
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.spine = continued_edge(ed, ext, st.sides, sup, tol, rb.a, rb.b);
            rb.t1 = ed.u1, rb.t2 = ed.u2;
            for (int i = 0; i <= 64; ++i)
            {
                const T t = rb.a + (rb.b - rb.a) * T(i) / T(64);
                if (!(D1.d(t) > tol) || (D2 && !(D2->d(t) > tol)))
                    refuse<T>(BuildErrc::InvalidBounds, op,
                              "a distance of the chamfer not larger than the tolerance on the chamfer continued past the vertices with the slope of "
                              "its law (FilletOptions::extension shorter, or the law less steep there)",
                              {eid});
            }
            // the pieces, between the kinks of the laws within the chamfer
            std::vector<T> cuts;
            for (const auto *D : {&D1, D2 ? &*D2 : nullptr})
                if (D)
                    for (T k : D->kinks)
                        if (k > rb.a && k < rb.b)
                            cuts.push_back(k);
            std::ranges::sort(cuts);
            {
                std::vector<T> u{rb.a};
                for (T k : cuts)
                    if (k - u.back() > T(1e-9) * (rb.b - rb.a))
                        u.push_back(k);
                if (rb.b - u.back() <= T(1e-9) * (rb.b - rb.a))
                    u.pop_back();
                u.push_back(rb.b);
                cuts = std::move(u);
            }
            const T acc = std::min(o.approx_tol, io_tol);
            std::array<std::vector<BSCurve<T, 2>>, 2> parts; // the pcurve of each contact piece by piece
            auto contact = [&](std::size_t k, const std::function<T(T)> &dist) -> std::shared_ptr<Curve<T, 2>> {
                for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
                {
                    const auto pc = contact_pcurve(st.sides[k], *sup[k], *rb.spine, ed, cuts[i], cuts[i + 1], dist, acc, tol, op, eid);
                    parts[k].push_back(*std::dynamic_pointer_cast<BSCurve<T, 2>>(pc));
                }
                return std::make_shared<BSCurve<T, 2>>(parts[k].size() == 1 ? parts[k].front() : joined(parts[k]));
            };
            auto set_contact = [&](std::size_t k, std::shared_ptr<Curve<T, 2>> pc) {
                auto &ct = rb.contacts[k];
                ct.face = st.sides[k].face;
                ct.support = sup[k];
                ct.pcurve = std::move(pc);
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, st.sides[k], *ct.pcurve, ed, ed.u1, ed.u2, tol, op, eid);
            };
            set_contact(0, contact(0, D1.d));
            if (D2)
                set_contact(1, contact(1, D2->d));
            else
            {
                // the other distance: where the line at the angle from the tangent plane of the
                // face at the first contact meets the other face, in the plane normal to the edge
                // (from the nearest point found); the tangent, not the chord to the edge, which
                // turned the line on a face curved across the edge (36 degrees for 30, d1 = 0.3 on
                // a cylinder of radius 1)
                const auto c0 = rb.contacts[0].curve;
                const auto s0 = sup[0], s1 = sup[1];
                const bool fwd0 = st.sides[0].orient == Orientation::Forward;
                const T into = st.convex ? T(-1) : T(1);
                auto found = std::make_shared<std::map<T, std::array<T, 3>>>(); // t -> (lambda, u, v)
                const auto pc1 = st.sides[1].pcurve;
                const T u1 = ed.u1, u2 = ed.u2;
                const auto spine = rb.spine;
                std::function<T(T)> d2 = [=](T t) {
                    const auto e = spine->value(t), te = unit(spine->value(t, 1));
                    const auto p1 = c0->value(t);
                    const auto uv0 = rb.contacts[0].pcurve->value(t);
                    auto n = normal_direction(*s0, uv0[0], uv0[1]);
                    if (!fwd0)
                        n = T(-1) * n;
                    n = unit(n - (n * te) * te);
                    // the tangent of the face in the plane normal to the edge, toward the edge
                    auto u = unit(cross(n, te));
                    if (u * (e - p1) < T(0))
                        u = T(-1) * u;
                    const auto dir = std::cos(angle) * u + into * std::sin(angle) * n;
                    thick_detail::Vec<T> x(3);
                    if (found->empty())
                    {
                        const auto uv = pc1->value(std::clamp(t, u1, u2));
                        x << norm(e - p1), uv[0], uv[1];
                    }
                    else
                    {
                        auto it = found->lower_bound(t);
                        if (it == found->end() || (it != found->begin() && t - std::prev(it)->first < it->first - t))
                            it = std::prev(it);
                        x << it->second[0], it->second[1], it->second[2];
                    }
                    const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                        const auto q = p1 + y[0] * dir - s1->value(y[1], y[2]);
                        const auto su = s1->value(y[1], y[2], 1, 0), sv = s1->value(y[1], y[2], 0, 1);
                        f.resize(3);
                        J.resize(3, 3);
                        for (int i = 0; i < 3; ++i)
                        {
                            f[i] = q[std::size_t(i)];
                            J.row(i) << dir[std::size_t(i)], -su[std::size_t(i)], -sv[std::size_t(i)];
                        }
                    });
                    if (!(res <= tol / T(100)) || !(x[0] > T(0)))
                        refuse<T>(BuildErrc::NonGeneric, op, "the line at the angle of the chamfer does not meet the other face", {eid, st.sides[1].face});
                    (*found)[t] = {x[0], x[1], x[2]};
                    return norm(p1 + x[0] * dir - e);
                };
                set_contact(1, contact(1, d2));
            }
            if (cuts.size() == 2)
            {
                // one piece: a primitive when the contacts lie on it within tol, a ruled B-spline otherwise
                if (auto e = ruled_primitive(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, tol))
                {
                    rb.surface = e;
                    on_blend(rb, io_tol, false);
                    rb.exact = std::ranges::all_of(rb.contacts, [&](const BlendContact<T> &c) { return c.deviation <= tol; });
                }
                if (!rb.exact)
                {
                    rb.surface = ruled_spline(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, o.approx_tol, rb.deviation);
                    on_blend(rb, io_tol, false);
                    // (built exactly between B-spline contacts on planes)
                    rb.exact = rb.deviation == T(0) && std::ranges::all_of(rb.contacts, [&](const BlendContact<T> &c) { return c.deviation <= tol; });
                }
                return rb;
            }
            // in pieces: each ruled, (u, v) = (t, across), exactly between the pieces of B-spline
            // contacts on planes; a crease along the ruling at each cut
            for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
            {
                BlendPiece<T> pc;
                pc.a = cuts[i], pc.b = cuts[i + 1];
                if (const auto e1 = on_plane(parts[0][i], *sup[0]), e2 = on_plane(parts[1][i], *sup[1]); e1 && e2)
                {
                    pc.surface = ruled_exact(*e1, *e2);
                    pc.deviation = T(0);
                }
                else
                    pc.surface = ruled_spline(*rb.contacts[0].curve, *rb.contacts[1].curve, pc.a, pc.b, o.approx_tol, pc.deviation);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    pc.on_blend[k] = segment<T>({pc.a, T(k)}, {pc.b, T(k)}, pc.a, pc.b);
                    for (int j = 0; j <= 64; ++j)
                    {
                        const T s = pc.a + (pc.b - pc.a) * T(j) / T(64);
                        rb.contacts[k].deviation =
                            std::max(rb.contacts[k].deviation, norm(pc.surface->value(s, T(k)) - rb.contacts[k].curve->value(s)));
                    }
                }
                rb.deviation = std::max(rb.deviation, pc.deviation);
                rb.pieces.push_back(std::move(pc));
            }
            for (std::size_t i = 1; i + 1 < cuts.size(); ++i)
            {
                BlendCrease<T> cr;
                cr.t = cuts[i];
                const auto p = rb.contacts[0].curve->value(cr.t), q = rb.contacts[1].curve->value(cr.t);
                cr.curve = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{p, q}, std::vector<T>{T(0), T(0), T(1), T(1)}, std::size_t{1});
                const auto w = segment<T>({cr.t, T(0)}, {cr.t, T(1)}, T(0), T(1));
                cr.on_pieces = {w, w};
                for (int j = 0; j <= 16; ++j)
                {
                    const T s = T(j) / T(16);
                    const auto a = rb.pieces[i - 1].surface->value(cr.t, s), b = rb.pieces[i].surface->value(cr.t, s);
                    cr.deviation = std::max({cr.deviation, norm(a - b), norm(a - cr.curve->value(s))});
                }
                rb.creases.push_back(std::move(cr));
            }
            rb.exact = rb.deviation == T(0) && std::ranges::all_of(rb.contacts, [&](const BlendContact<T> &c) { return c.deviation <= tol; });
            return rb;
        }

        /// The checks of a fillet of chord c, before its edge is set up.
        template <std::floating_point T>
        void check_chord(const FilletOptions<T> &o, T c, EdgeId eid, const std::string &op)
        {
            check_options(o, op, eid);
            if (!(c > o.tol) || !std::isfinite(c))
                refuse<T>(BuildErrc::InvalidBounds, op, "the chord must be larger than the tolerance", {eid});
            // (the contacts of a G2 section are not those of its ball: its chord would not be c)
            if (o.section.kind == SectionKind::G2)
                refuse<T>(BuildErrc::NonGeneric, op, "a fillet of constant chord with a G2 section (not handled: a circular or conic section)", {eid});
        }

        /**
         * The radius of the fillet of chord c along an edge between two planes (stage 7, § 4.3):
         * their angle constant, the radius c / (2 sin(phi / 2)) (phi the angle of their outward
         * normals), the fillet of that constant radius exact (a cylinder for a circular section);
         * nullopt for other faces. `st` the setup of the edge `ed` (faces tangent along it refused there).
         */
        template <std::floating_point T>
        auto planar_chord_radius(const EdgeSetup<T> &st, const Edge<T> &ed, T c) -> std::optional<T>
        {
            std::array<point<T, 3>, 2> n;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto *pl = dynamic_cast<const Plane<T> *>(st.sides[k].surface.get());
                if (!pl)
                    return std::nullopt;
                const auto uv = st.sides[k].pcurve->value((ed.u1 + ed.u2) / T(2));
                n[k] = normal_direction(*pl, uv[0], uv[1]);
                if (st.sides[k].orient != Orientation::Forward)
                    n[k] = T(-1) * n[k];
            }
            const T phi = std::atan2(norm(cross(n[0], n[1])), n[0] * n[1]);
            return c / (T(2) * std::sin(phi / T(2)));
        }

        /**
         * The radius of the fillet of chord c along the edge (stage 7, § 4.3; note
         * `brep7_pr02_conic_chordal.md`): at the parameter t of the edge (continued past its
         * vertices as far as `fillet_law` continues the fillet: max(`extension`, 1.5) times the
         * larger of c and the radius, a quarter more for the radius to grow; held beyond), the
         * radius r of the ball that, in the plane normal to the edge at t, touches both faces at
         * points c apart. The secant on r from the radius between
         * two planes, c / (2 sin(phi / 2)) (phi the angle of the outward normals of the faces at
         * the edge); for each r, the ball by Gauss-Newton on its contacts (u1, v1, u2, v2):
         * S1 + s r N1 = S2 + s r N2, the centre in the plane normal to the edge (s: into the solid
         * along a convex edge, out of it along a concave one), from the nearest station solved, the
         * first ones from the contacts at the distance of two planes (`contact_pcurve`). Between
         * two planes at a constant angle, the radius is constant: the fillet a cylinder. Refused:
         * a ball that cannot be put on both faces at that chord, a chord not beyond tol.
         */
        template <std::floating_point T>
        auto chordal_radius(const Model<T> &m, const EdgeSetup<T> &st, EdgeId eid, T c, const FilletOptions<T> &o) -> EdgeRadius<T>
        {
            const std::string op = "fillet_geometry";
            const auto &ed = m.edge(eid);
            const T tol = o.tol;
            // the planar radius and contact distance at t (the normals of the faces at the edge, t within it)
            // (copies: the radius is evaluated after this returns)
            auto outward = [sides = st.sides, u1 = ed.u1, u2 = ed.u2](std::size_t k, T t) {
                const auto &sd = sides[k];
                const auto uv = sd.pcurve->value(std::clamp(t, u1, u2));
                const auto n = normal_direction(*sd.surface, uv[0], uv[1]);
                return sd.orient == Orientation::Forward ? n : T(-1) * n;
            };
            auto planar = [outward, c](T t) {
                const auto n1 = outward(0, t), n2 = outward(1, t);
                const T phi = std::atan2(norm(cross(n1, n2)), n1 * n2);
                const T r = c / (T(2) * std::sin(std::max(phi, T(1e-3)) / T(2)));
                return std::pair{r, r * std::tan(std::max(phi, T(1e-3)) / T(2))}; // (radius, distance of the contacts from the edge)
            };
            const T rmax0 = std::max(planar(ed.u1).first, planar(ed.u2).first);
            // (as far as `fillet_law` continues the fillet: `extension` radii, 1.5 at a vertex
            // where a chain goes on; `extension` c only stopped short of that along a chain)
            const T ext = st.closed ? T(0) : T(1.25) * std::max(o.extension, T(1.5)) * std::max(c, rmax0);
            const T reach = T(4) * (rmax0 + c) + ext;
            std::array<std::shared_ptr<Surface<T, 3>>, 2> sup;
            for (std::size_t k = 0; k < 2; ++k)
                sup[k] = support(m, st.sides[k], ed, reach, tol).surface;
            T a{}, b{};
            const auto spine = continued_edge(ed, ext, st.sides, sup, tol, a, b);
            // the seeds: the contacts at the distance of two planes
            std::array<std::shared_ptr<Curve<T, 2>>, 2> seed;
            for (std::size_t k = 0; k < 2; ++k)
                seed[k] = contact_pcurve(st.sides[k], *sup[k], *spine, ed, a, b, std::function<T(T)>{[&](T t) { return planar(t).second; }},
                                         T(1e-3) * c, tol, op, eid);
            const T side = st.convex ? T(-1) : T(1);
            const std::array<bool, 2> fwd{st.sides[0].orient == Orientation::Forward, st.sides[1].orient == Orientation::Forward};
            auto cache = std::make_shared<std::map<T, std::array<T, 5>>>(); // t -> (u1, v1, u2, v2, r)
            auto ball = [=](T t, T r, std::array<T, 4> x0) -> std::optional<std::array<T, 4>> {
                const auto e = spine->value(t), te = unit(spine->value(t, 1));
                auto centre = [&](std::size_t k, T u, T v) {
                    auto n = normal_direction(*sup[k], u, v);
                    if (!fwd[k])
                        n = T(-1) * n;
                    return sup[k]->value(u, v) + side * r * n;
                };
                thick_detail::Vec<T> x(4);
                x << x0[0], x0[1], x0[2], x0[3];
                auto eval = [&](const thick_detail::Vec<T> &y) {
                    const auto c1 = centre(0, y[0], y[1]), c2 = centre(1, y[2], y[3]);
                    const auto d = c1 - c2;
                    thick_detail::Vec<T> f(4);
                    f << d[0], d[1], d[2], (c1 - e) * te;
                    return f;
                };
                const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                    f = eval(y);
                    J.resize(4, 4);
                    for (int j = 0; j < 4; ++j)
                    {
                        auto yp = y, ym = y;
                        const T h = T(1e-7) * (T(1) + std::abs(y[j]));
                        yp[j] += h, ym[j] -= h;
                        J.col(j) = (eval(yp) - eval(ym)) / (T(2) * h);
                    }
                });
                if (!(res <= T(1e-12) * (T(1) + norm(e) + r)))
                    return std::nullopt;
                return std::array<T, 4>{x[0], x[1], x[2], x[3]};
            };
            EdgeRadius<T> R;
            const bool closed = st.closed;
            const T u1 = ed.u1, u2 = ed.u2, step = (u2 - u1) / T(64);
            // the ball at t whose contacts are c apart, from (x0, r0): the secant on r
            auto solve = [=](T t, std::array<T, 4> x, T r) -> std::pair<std::array<T, 4>, T> {
                auto chord = [&](T rr, std::array<T, 4> &y) -> T {
                    const auto z = ball(t, rr, y);
                    if (!z)
                        refuse<T>(BuildErrc::NonGeneric, op, "a ball of the chord cannot be put on both faces", {eid});
                    y = *z;
                    return norm(sup[0]->value(y[0], y[1]) - sup[1]->value(y[2], y[3])) - c;
                };
                T r0 = r, r1 = r * (T(1) + T(1e-3));
                T g0 = chord(r0, x);
                auto x1 = x;
                T g1 = chord(r1, x1);
                for (int i = 0; i < 50 && std::abs(g1) > T(1e-13) * c; ++i)
                {
                    if (!(std::abs(g1 - g0) > T(0)))
                        break;
                    const T r2 = r1 - g1 * (r1 - r0) / (g1 - g0);
                    if (!(r2 > T(0)))
                        refuse<T>(BuildErrc::NonGeneric, op, "a ball of the chord cannot be put on both faces", {eid});
                    r0 = r1, g0 = g1, r1 = r2;
                    g1 = chord(r1, x1);
                }
                if (!(std::abs(g1) <= T(1e-9) * c))
                    refuse<T>(BuildErrc::NonGeneric, op, "a ball of the chord cannot be put on both faces", {eid});
                return {x1, r1};
            };
            R.r = [=](T t) {
                // (a closed edge: t brought back onto its loop; an open one: within the edge
                // continued, the radius held past its ends)
                if (closed)
                    t = u1 + std::fmod(std::fmod(t - u1, u2 - u1) + (u2 - u1), u2 - u1);
                else
                    t = std::clamp(t, a, b);
                if (const auto it = cache->find(t); it != cache->end())
                    return it->second[4];
                // the first ball at the middle of the edge, from the seeds; then marched from the
                // nearest ball found, by steps of a 64th of the edge (the same branch all along)
                if (cache->empty())
                {
                    const T tm = (u1 + u2) / T(2);
                    const auto w1 = seed[0]->value(tm), w2 = seed[1]->value(tm);
                    const auto [x, r] = solve(tm, {w1[0], w1[1], w2[0], w2[1]}, planar(tm).first);
                    (*cache)[tm] = {x[0], x[1], x[2], x[3], r};
                }
                auto it = cache->lower_bound(t);
                if (it == cache->end() || (it != cache->begin() && t - std::prev(it)->first < it->first - t))
                    it = std::prev(it);
                T at = it->first;
                std::array<T, 4> x{it->second[0], it->second[1], it->second[2], it->second[3]};
                T r = it->second[4];
                for (;;)
                {
                    const T next = std::abs(t - at) <= step ? t : at + (t > at ? step : -step);
                    std::tie(x, r) = solve(next, x, r);
                    (*cache)[next] = {x[0], x[1], x[2], x[3], r};
                    if (next == t)
                        return r;
                    at = next;
                }
            };
            R.r1 = R.r(ed.u1), R.r2 = R.r(ed.u2);
            R.rmax = std::max(R.r1, R.r2);
            for (int i = 1; i < 32; ++i)
                R.rmax = std::max(R.rmax, R.r(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(32)));
            return R;
        }

        /// The same, the chord checked and the edge set up (`check_chord`, `setup`).
        template <std::floating_point T>
        auto chordal_radius(const Model<T> &m, SolidId sid, EdgeId eid, T c, const FilletOptions<T> &o) -> EdgeRadius<T>
        {
            const std::string op = "fillet_geometry";
            check_chord(o, c, eid, op);
            return chordal_radius(m, setup(m, sid, eid, op), eid, c, o);
        }

        /**
         * The fillet of chord c along the edge: between two planes, the radius of the fillet of
         * constant radius, exact (`planar_chord_radius`); otherwise the law of its radius
         * (`chordal_radius`). The chord checked and the edge set up once.
         */
        template <std::floating_point T>
        auto chord_radius(const Model<T> &m, SolidId sid, EdgeId eid, T c, const FilletOptions<T> &o) -> std::variant<T, EdgeRadius<T>>
        {
            const std::string op = "fillet_geometry";
            check_chord(o, c, eid, op);
            const auto st = setup(m, sid, eid, op);
            if (const auto r = planar_chord_radius(st, m.edge(eid), c))
                return *r;
            return chordal_radius(m, st, eid, c, o);
        }

        /// The section of the fillet of the edge e alone (`BlendSection`; a law of rho by its points on e; the G2 quartic). Refused (`InvalidBounds`): rho out of ]0, 1[.
        template <std::floating_point T>
        auto section_of_edge(const Model<T> &m, EdgeId e, const BlendSection<T> &sec, T tol, const std::string &op) -> SectionOf<T>
        {
            SectionOf<T> s;
            s.g2 = sec.kind == SectionKind::G2;
            if (sec.kind != SectionKind::Conic)
                return s;
            s.conic = true;
            auto inside = [](T x) { return x > T(0) && x < T(1); };
            if (!sec.rho_law)
            {
                if (!inside(sec.rho))
                    refuse<T>(BuildErrc::InvalidBounds, op, "the rho of a conic section must be within ]0, 1[", {e});
                s.rho = [rho = sec.rho](T) { return rho; };
                return s;
            }
            for (const auto &p : sec.rho_law->points)
                if (!inside(p.value))
                    refuse<T>(BuildErrc::InvalidBounds, op, "the rho of a conic section must be within ]0, 1[", {p.edge});
            const auto R = edge_radius(m, e, *sec.rho_law, std::min(tol, T(1e-9)), op, "rho");
            if (!R)
            {
                const T rho = edge_value(m, e, *sec.rho_law, op);
                s.rho = [rho](T) { return rho; };
                return s;
            }
            s.constant = false;
            s.kinks = R->kinks;
            s.rho = [f = R->r](T t) { return std::clamp(f(t), T(1e-3), T(1) - T(1e-3)); };
            return s;
        }

        template <typename F>
        auto guarded(const std::string &op, F &&f) -> BuildResult<std::invoke_result_t<F>>
        {
            try
            {
                return f();
            }
            catch (const BuildFailure &e)
            {
                return std::unexpected(e.error());
            }
            catch (const std::exception &e)
            {
                return detail::fail(BuildErrc::NonGeneric, op + ": " + e.what());
            }
        }
    }

    /**
     * @brief The geometry of the fillet of radius r along the edge `e` of the solid `s` (file
     * header): its centre line, its contacts on the two faces of the edge and its surface,
     * continued by `extension` r past the vertices of an open edge. The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_geometry(const Model<T> &m, SolidId s, EdgeId e, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("fillet_geometry",
                                      [&] { return fillet_detail::fillet(m, s, e, r, o, fillet_detail::section_of_edge(m, e, o.section, o.tol, "fillet_geometry")); });
    }

    /**
     * @brief The geometry of the fillet along the edge `e` of the solid `s` whose radius follows
     * `law` (its points on `e`; stage 5, § 3): the balls of radius r(s) in the planes normal to
     * the edge, their envelope approximated within `approx_tol` and `angular_tol`; every curve on
     * the parameter of the edge, continued past each vertex of an open edge by `extension` times
     * the larger of the radii at the vertices. Past its first and last points the law goes on
     * with its slope there: constant for a smooth law, straight on for a linear one (refused
     * where it is no longer positive at an end of the blend). The same radius everywhere on the
     * edge: the fillet of constant radius. The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_geometry(const Model<T> &m, SolidId s, EdgeId e, const RadiusLaw<T> &law, const FilletOptions<T> &o = {})
        -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("fillet_geometry", [&] {
            const auto R = fillet_detail::edge_radius(m, e, law, o.tol, "fillet_geometry");
            if (!R)
                return fillet_detail::fillet(m, s, e, fillet_detail::edge_value(m, e, law, "fillet_geometry"), o,
                                             fillet_detail::section_of_edge(m, e, o.section, o.tol, "fillet_geometry"));
            return fillet_detail::fillet_law(m, s, e, *R, o, T(0), fillet_detail::section_of_edge(m, e, o.section, o.tol, "fillet_geometry"));
        });
    }

    /**
     * @brief The geometry of the fillet of chord c along the edge `e` of the solid `s` (stage 7,
     * § 4.3): its contacts c apart at every station, the radius that follows, the envelope of
     * its balls (as a fillet whose radius follows a law), the section of `o`. Between two planes,
     * the fillet of constant radius c / (2 sin(phi / 2)), exact (a cylinder for a circular
     * section). The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_chordal_geometry(const Model<T> &m, SolidId s, EdgeId e, std::type_identity_t<T> c, const FilletOptions<T> &o = {})
        -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("fillet_geometry", [&] {
            const auto R = fillet_detail::chord_radius(m, s, e, c, o);
            const auto sec = fillet_detail::section_of_edge(m, e, o.section, o.tol, "fillet_geometry");
            if (const auto *r = std::get_if<T>(&R))
                return fillet_detail::fillet(m, s, e, *r, o, sec); // exact, as a fillet of that radius
            return fillet_detail::fillet_law(m, s, e, std::get<fillet_detail::EdgeRadius<T>>(R), o, T(0), sec);
        });
    }

    /**
     * @brief The geometry of the chamfer along the edge `e` of the solid `s` at distances following
     * laws (stage 7, § 3): `d1` on the face `f`, `d2` on the other (their points on `e`); in
     * pieces where the slope of a law jumps (`EdgeBlend::pieces`, `creases`). The model is not
     * modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_geometry(const Model<T> &m, SolidId s, EdgeId e, FaceId f, const BlendLaw<T> &d1, const BlendLaw<T> &d2,
                                        const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("chamfer_geometry", [&] {
            auto dist = [&](const BlendLaw<T> &law) {
                const auto R = fillet_detail::edge_radius(m, e, law, o.tol, "chamfer_geometry", "distance");
                if (R)
                    return fillet_detail::edge_distance(*R);
                return fillet_detail::edge_distance<T>(fillet_detail::edge_value(m, e, law, "chamfer_geometry"));
            };
            return fillet_detail::chamfer_law(m, s, e, f, dist(d1), std::optional{dist(d2)}, T(0), o);
        });
    }

    /**
     * @brief The geometry of the chamfer along the edge `e` of the solid `s` at the distance `d1`
     * on the face `f` and the angle `angle` (radians) from it (stage 7, § 3.2). The model is not
     * modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_angle_geometry(const Model<T> &m, SolidId s, EdgeId e, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> angle,
                                              const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("chamfer_geometry",
                                      [&] { return fillet_detail::chamfer_law(m, s, e, f, fillet_detail::edge_distance<T>(d1), std::nullopt, angle, o); });
    }

    /**
     * @brief The geometry of the chamfer along the edge `e` of the solid `s`, at the distances d1
     * on the face `f` of the edge and d2 on the other (file header): its contacts (on `f` first)
     * and its ruled surface, continued by `extension` max(d1, d2) past the vertices of an open
     * edge. The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_geometry(const Model<T> &m, SolidId s, EdgeId e, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> d2,
                                        const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded("chamfer_geometry", [&] { return fillet_detail::chamfer(m, s, e, f, d1, d2, o); });
    }
}
