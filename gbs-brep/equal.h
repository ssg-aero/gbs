#pragma once

/**
 * @file equal.h
 * @brief Native BREP core: exact equality of two models (`equal`, `first_difference`).
 *
 * Architecture note docs/sources/design/io_brep_json.md, § 7. Two models are equal when one
 * could stand for the other wherever an id or a geometry is used:
 * - the same tables: the same size, the same entities alive, every field the same, every real
 *   bit for bit (-0 and 0 differ), the ids stored the same;
 * - the same attributes: names, external ids, instances (part and placement bit for bit), the
 *   next part key, the unit;
 * - the same geometry: of the same exact class, with the same defining data bit for bit (poles,
 *   weights, knots, degrees, bounds, frames, radii, placements...) and the same values and
 *   derivatives at sample parameters, bit for bit;
 * - the same sharing: the pairs of geometry objects met form a bijection, so that an object
 *   shared by several entities (or inside several geometries) on one side is one object shared
 *   the same way on the other, and two distinct objects stay two.
 *
 * A model and its copy (`Model` copies share the geometry) are equal; a model written and read
 * back by gbs-io/brep_json.h is equal to it. A geometry of a class gbs does not know is compared
 * by its bounds and its values at the samples only.
 */

#include <gbs-brep/edit.h>
#include <gbs-brep/model.h>
#include <gbs-brep/seamless.h>
#include <gbs/curves>
#include <gbs/offsets.h>
#include <gbs/primitives.h>
#include <gbs/surfaces>

#include <array>
#include <cmath>
#include <cstring>
#include <format>
#include <optional>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace gbs::brep
{
    namespace equal_detail
    {
        template <std::floating_point T>
        bool same(T a, T b) noexcept { return std::memcmp(&a, &b, sizeof(T)) == 0; }
        template <std::floating_point T, std::size_t N>
        bool same(const std::array<T, N> &a, const std::array<T, N> &b) noexcept { return std::memcmp(a.data(), b.data(), sizeof(T) * N) == 0; }
        template <typename X>
        bool same(const std::vector<X> &a, const std::vector<X> &b) noexcept
        {
            if (a.size() != b.size())
                return false;
            for (std::size_t i = 0; i < a.size(); ++i)
                if (!same(a[i], b[i]))
                    return false;
            return true;
        }
        template <std::floating_point T>
        bool same(const Placement<T> &a, const Placement<T> &b) noexcept
        {
            return same(a.r[0], b.r[0]) && same(a.r[1], b.r[1]) && same(a.r[2], b.r[2]) && same(a.t, b.t);
        }
        template <std::floating_point T>
        bool same(const Frame3<T> &a, const Frame3<T> &b) noexcept { return same(a.o, b.o) && same(a.x, b.x) && same(a.y, b.y) && same(a.z, b.z); }

        /// The first difference found, thrown to the top.
        struct Difference
        {
            std::string what;
        };

        template <std::floating_point T>
        class Comparer
        {
            std::unordered_map<const void *, const void *> a2b_, b2a_;
            std::string ctx_; ///< the entity being compared

            [[noreturn]] void differ(std::string_view what) const { throw Difference{ctx_.empty() ? std::string{what} : std::format("{}: {}", ctx_, what)}; }
            void expect(bool ok, std::string_view what) const
            {
                if (!ok)
                    differ(what);
            }

            /// true the first time the pair is met; a pair of objects met before must be the same pair
            bool pair(const void *a, const void *b)
            {
                const auto ia = a2b_.find(a);
                const auto ib = b2a_.find(b);
                if (ia == a2b_.end() && ib == b2a_.end())
                {
                    a2b_.emplace(a, b);
                    b2a_.emplace(b, a);
                    return true;
                }
                if (ia == a2b_.end() || ib == b2a_.end() || ia->second != b || ib->second != a)
                    differ("a geometry shared on one side and not the same way on the other");
                return false;
            }

            template <std::size_t D>
            void sample(const Curve<T, D> &a, const Curve<T, D> &b) const
            {
                const auto [u1, u2] = a.bounds();
                if (!std::isfinite(u2 - u1))
                    return;
                for (int i = 0; i <= 4; ++i)
                    for (std::size_t d = 0; d <= 2; ++d)
                    {
                        const T u = u1 + (u2 - u1) * T(i) / T(4);
                        std::optional<std::array<T, D>> pa, pb;
                        try { pa = a.value(u, d); } catch (const std::exception &) {}
                        try { pb = b.value(u, d); } catch (const std::exception &) {}
                        expect(pa.has_value() == pb.has_value() && (!pa || same(*pa, *pb)), std::format("curve values differ at u = {} (derivative {})", u, d));
                    }
            }
            template <std::size_t D>
            void sample(const Surface<T, D> &a, const Surface<T, D> &b) const
            {
                const auto bx = a.bounds();
                if (!std::isfinite(bx[1] - bx[0]) || !std::isfinite(bx[3] - bx[2]))
                    return;
                for (int i = 0; i <= 2; ++i)
                    for (int j = 0; j <= 2; ++j)
                        for (std::size_t d = 0; d <= 2; ++d)
                        {
                            const T u = bx[0] + (bx[1] - bx[0]) * T(i) / T(2), v = bx[2] + (bx[3] - bx[2]) * T(j) / T(2);
                            std::optional<std::array<T, D>> pa, pb;
                            try { pa = a.value(u, v, d == 1, d == 2); } catch (const std::exception &) {}
                            try { pb = b.value(u, v, d == 1, d == 2); } catch (const std::exception &) {}
                            expect(pa.has_value() == pb.has_value() && (!pa || same(*pa, *pb)),
                                   std::format("surface values differ at (u, v) = ({}, {}) (derivative {})", u, v, d));
                        }
            }

            template <typename X, typename B>
            static const X *as(const B &b) { return typeid(b) == typeid(X) ? static_cast<const X *>(&b) : nullptr; }

            template <std::size_t N, bool R>
            void bspline(const BSCurveGeneral<T, N, R> &a, const BSCurveGeneral<T, N, R> &b) const
            {
                expect(a.degree() == b.degree(), "B-spline degrees differ");
                expect(same(a.knotsFlats(), b.knotsFlats()), "B-spline knots differ");
                expect(same(a.poles(), b.poles()), "B-spline poles differ");
                expect(same(a.bounds(), b.bounds()), "B-spline bounds differ");
            }
            template <std::size_t N, bool R>
            void bspline(const BSSurfaceGeneral<T, N, R> &a, const BSSurfaceGeneral<T, N, R> &b) const
            {
                expect(a.degreeU() == b.degreeU() && a.degreeV() == b.degreeV(), "B-spline surface degrees differ");
                expect(same(a.knotsFlatsU(), b.knotsFlatsU()) && same(a.knotsFlatsV(), b.knotsFlatsV()), "B-spline surface knots differ");
                expect(same(a.poles(), b.poles()), "B-spline surface poles differ");
            }

        public:
            void context(std::string c) { ctx_ = std::move(c); }

            template <std::size_t D>
            void curve(const Curve<T, D> *a, const Curve<T, D> *b)
            {
                if (!a || !b)
                {
                    expect(a == b, "a curve on one side only");
                    return;
                }
                if (!pair(dynamic_cast<const void *>(a), dynamic_cast<const void *>(b)) || a == b)
                    return;
                expect(typeid(*a) == typeid(*b), std::format("curves of classes {} and {}", typeid(*a).name(), typeid(*b).name()));
                expect(same(a->bounds(), b->bounds()), "curve bounds differ");
                if (auto x = as<BSCurve<T, D>>(*a))
                    bspline(*x, *as<BSCurve<T, D>>(*b));
                else if (auto x = as<BSCurveRational<T, D>>(*a))
                    bspline(*x, *as<BSCurveRational<T, D>>(*b));
                else if (auto x = as<CurveTrimmed<T, D>>(*a))
                    curve<D>(x->basisCurve().get(), as<CurveTrimmed<T, D>>(*b)->basisCurve().get());
                else if (auto x = as<CurveComposite<T, D>>(*a))
                {
                    const auto &ca = x->curves();
                    const auto &cb = as<CurveComposite<T, D>>(*b)->curves();
                    expect(ca.size() == cb.size(), "composite curves of different numbers of curves");
                    for (std::size_t i = 0; i < ca.size(); ++i)
                        curve<D>(ca[i].get(), cb[i].get());
                }
                else if (auto x = as<CurveReversed<T, D>>(*a))
                    curve<D>(x->basisCurve().get(), as<CurveReversed<T, D>>(*b)->basisCurve().get());
                else if (auto x = as<CurveReparametrized<T, D>>(*a))
                {
                    const auto y = as<CurveReparametrized<T, D>>(*b);
                    bspline(x->paramFunction().basisCurve(), y->paramFunction().basisCurve());
                    curve<D>(x->basisCurve().get(), y->basisCurve().get());
                }
                else if (auto x = as<CurveExtended<T, D>>(*a))
                {
                    const auto y = as<CurveExtended<T, D>>(*b);
                    expect(x->isClamped() == y->isClamped(), "extended curves, clamped on one side only");
                    curve<D>(x->basisCurve().get(), y->basisCurve().get());
                }
                else if (auto x = as<CurveOnSurface<T, D>>(*a))
                {
                    const auto y = as<CurveOnSurface<T, D>>(*b);
                    curve<2>(x->p_basisCurve().get(), y->p_basisCurve().get());
                    surface<D>(x->p_basisSurface().get(), y->p_basisSurface().get());
                }
                if constexpr (D == 3)
                {
                    if (auto x = as<LineCurve<T>>(*a))
                        expect(same(x->origin(), as<LineCurve<T>>(*b)->origin()) && same(x->vector(), as<LineCurve<T>>(*b)->vector()), "lines differ");
                    else if (auto x = as<EllipseCurve<T>>(*a))
                    {
                        const auto y = as<EllipseCurve<T>>(*b);
                        expect(same(x->radius1(), y->radius1()) && same(x->radius2(), y->radius2()) && same(x->frame(), y->frame()), "ellipses differ");
                    }
                    else if (auto x = as<MovedCurve<T>>(*a))
                    {
                        expect(same(x->placement(), as<MovedCurve<T>>(*b)->placement()), "placements of moved curves differ");
                        curve<3>(x->basis().get(), as<MovedCurve<T>>(*b)->basis().get());
                    }
                    else if (auto x = as<CurveOffset3D<T, BSCfunction<T>>>(*a))
                    {
                        const auto y = as<CurveOffset3D<T, BSCfunction<T>>>(*b);
                        expect(same(x->direction(), y->direction()), "directions of offset curves differ");
                        bspline(x->offset().basisCurve(), y->offset().basisCurve());
                        curve<3>(&x->basisCurve(), &y->basisCurve());
                    }
                }
                else if (auto x = as<CurveOffset2D<T, BSCfunction<T>>>(*a))
                {
                    const auto y = as<CurveOffset2D<T, BSCfunction<T>>>(*b);
                    bspline(x->offset().basisCurve(), y->offset().basisCurve());
                    curve<2>(&x->basisCurve(), &y->basisCurve());
                }
                sample(*a, *b);
            }

            template <std::size_t D>
            void surface(const Surface<T, D> *a, const Surface<T, D> *b)
            {
                if (!a || !b)
                {
                    expect(a == b, "a surface on one side only");
                    return;
                }
                if (!pair(dynamic_cast<const void *>(a), dynamic_cast<const void *>(b)) || a == b)
                    return;
                expect(typeid(*a) == typeid(*b), std::format("surfaces of classes {} and {}", typeid(*a).name(), typeid(*b).name()));
                expect(same(a->bounds(), b->bounds()), "surface bounds differ");
                if (auto x = as<BSSurface<T, D>>(*a))
                    bspline(*x, *as<BSSurface<T, D>>(*b));
                else if (auto x = as<BSSurfaceRational<T, D>>(*a))
                    bspline(*x, *as<BSSurfaceRational<T, D>>(*b));
                if constexpr (D == 3)
                {
                    if (auto x = dynamic_cast<const ElementarySurface<T> *>(a))
                    {
                        const auto y = dynamic_cast<const ElementarySurface<T> *>(b);
                        expect(same(x->frame(), y->frame()), "frames of primitives differ");
                        if (auto c = as<CylindricalSurface<T>>(*a))
                            expect(same(c->radius(), as<CylindricalSurface<T>>(*b)->radius()), "radii of cylinders differ");
                        else if (auto c = as<ConicalSurface<T>>(*a))
                            expect(same(c->radius(), as<ConicalSurface<T>>(*b)->radius()) && same(c->semi_angle(), as<ConicalSurface<T>>(*b)->semi_angle()),
                                   "cones differ");
                        else if (auto c = as<SphericalSurface<T>>(*a))
                            expect(same(c->radius(), as<SphericalSurface<T>>(*b)->radius()), "radii of spheres differ");
                        else if (auto c = as<ToroidalSurface<T>>(*a))
                            expect(same(c->major_radius(), as<ToroidalSurface<T>>(*b)->major_radius()) &&
                                       same(c->minor_radius(), as<ToroidalSurface<T>>(*b)->minor_radius()),
                                   "radii of tori differ");
                    }
                    else if (auto x = as<SurfaceOfRevolution<T>>(*a))
                    {
                        const auto y = as<SurfaceOfRevolution<T>>(*b);
                        expect(same(x->axis2()[0], y->axis2()[0]) && same(x->axis2()[1], y->axis2()[1]) && same(x->axis2()[2], y->axis2()[2]),
                               "axes of surfaces of revolution differ");
                        curve<2>(x->basisCurve().get(), y->basisCurve().get());
                    }
                    else if (auto x = as<OffsetSurface<T>>(*a))
                    {
                        expect(same(x->distance(), as<OffsetSurface<T>>(*b)->distance()), "distances of offset surfaces differ");
                        surface<3>(x->basis().get(), as<OffsetSurface<T>>(*b)->basis().get());
                    }
                    else if (auto x = as<MovedSurface<T>>(*a))
                    {
                        expect(same(x->placement(), as<MovedSurface<T>>(*b)->placement()), "placements of moved surfaces differ");
                        surface<3>(x->basis().get(), as<MovedSurface<T>>(*b)->basis().get());
                    }
                    else if (auto x = as<SeamlessSurface<T>>(*a))
                    {
                        expect(same(x->periods(), as<SeamlessSurface<T>>(*b)->periods()), "periods of seamless surfaces differ");
                        surface<3>(x->basis().get(), as<SeamlessSurface<T>>(*b)->basis().get());
                    }
                }
                sample(*a, *b);
            }

            template <ShapeType K>
            void attributes(const Model<T> &a, const Model<T> &b, Id<K> h) const
            {
                expect(a.name(h) == b.name(h), "names differ");
                expect(a.externalId(h) == b.externalId(h), "external ids differ");
                const auto ia = a.instance(h), ib = b.instance(h);
                expect(ia.has_value() == ib.has_value() && (!ia || (ia->part == ib->part && same(ia->place, ib->place))), "instances differ");
            }

            template <ShapeType K>
            void table(const Model<T> &a, const Model<T> &b)
            {
                static constexpr std::array<const char *, 7> names{"vertex", "edge", "wire", "face", "shell", "solid", "compound"};
                const char *name = names[static_cast<std::size_t>(K)];
                ctx_.clear();
                expect(a.template capacity<K>() == b.template capacity<K>(),
                       std::format("{} tables of {} and {} entities", name, a.template capacity<K>(), b.template capacity<K>()));
                for (std::uint32_t i = 0; i < a.template capacity<K>(); ++i)
                {
                    const Id<K> h{i};
                    ctx_ = std::format("{} {}", name, i);
                    expect(a.alive(h) == b.alive(h), "alive on one side only");
                    if (!a.alive(h))
                        continue;
                    if constexpr (K == ShapeType::Vertex)
                        expect(same(a.vertex(h).pnt, b.vertex(h).pnt) && same(a.vertex(h).tol, b.vertex(h).tol), "points or tolerances differ");
                    else if constexpr (K == ShapeType::Edge)
                    {
                        const auto &x = a.edge(h), &y = b.edge(h);
                        expect(same(x.u1, y.u1) && same(x.u2, y.u2), "bounds differ");
                        expect(x.v1 == y.v1 && x.v2 == y.v2, "vertices differ");
                        expect(same(x.tol, y.tol), "tolerances differ");
                        expect(x.degenerate == y.degenerate && x.same_parameter == y.same_parameter, "flags differ");
                        curve<3>(x.curve.get(), y.curve.get());
                    }
                    else if constexpr (K == ShapeType::Wire)
                    {
                        const auto &x = a.wire(h), &y = b.wire(h);
                        expect(x.closed == y.closed, "closed on one side only");
                        expect(x.coedges.size() == y.coedges.size(), "numbers of co-edges differ");
                        for (std::size_t k = 0; k < x.coedges.size(); ++k)
                        {
                            ctx_ = std::format("{} {}, co-edge {}", name, i, k);
                            expect(x.coedges[k].edge == y.coedges[k].edge && x.coedges[k].orient == y.coedges[k].orient, "edges or orientations differ");
                            curve<2>(x.coedges[k].pcurve.get(), y.coedges[k].pcurve.get());
                        }
                        ctx_ = std::format("{} {}", name, i);
                    }
                    else if constexpr (K == ShapeType::Face)
                    {
                        const auto &x = a.face(h), &y = b.face(h);
                        expect(x.wires == y.wires, "wires differ");
                        expect(same(x.tol, y.tol), "tolerances differ");
                        expect(x.natural_bounds == y.natural_bounds, "natural bounds on one side only");
                        surface<3>(x.surface.get(), y.surface.get());
                    }
                    else if constexpr (K == ShapeType::Shell)
                    {
                        const auto &x = a.shell(h), &y = b.shell(h);
                        expect(x.closed == y.closed, "closed on one side only");
                        expect(x.faces.size() == y.faces.size(), "numbers of faces differ");
                        for (std::size_t k = 0; k < x.faces.size(); ++k)
                            expect(x.faces[k].face == y.faces[k].face && x.faces[k].orient == y.faces[k].orient, "face uses differ");
                    }
                    else if constexpr (K == ShapeType::Solid)
                        expect(a.solid(h).outer == b.solid(h).outer && a.solid(h).voids == b.solid(h).voids, "shells differ");
                    else
                        expect(a.compound(h).shapes == b.compound(h).shapes, "shapes differ");
                    attributes(a, b, h);
                }
                ctx_.clear();
            }

            void run(const Model<T> &a, const Model<T> &b)
            {
                expect(same(a.unitScale(), b.unitScale()), "units differ");
                expect(a.nextPart() == b.nextPart(), "next part keys differ");
                table<ShapeType::Vertex>(a, b);
                table<ShapeType::Edge>(a, b);
                table<ShapeType::Wire>(a, b);
                table<ShapeType::Face>(a, b);
                table<ShapeType::Shell>(a, b);
                table<ShapeType::Solid>(a, b);
                table<ShapeType::Compound>(a, b);
            }
        };
    }

    /**
     * @brief The first difference between two models (as "edge 12: B-spline poles differ"), none
     * if they are equal (see the file comment). Never modifies them.
     */
    template <std::floating_point T>
    [[nodiscard]] auto first_difference(const Model<T> &a, const Model<T> &b) -> std::optional<std::string>
    {
        try
        {
            equal_detail::Comparer<T>{}.run(a, b);
        }
        catch (const equal_detail::Difference &d)
        {
            return d.what;
        }
        return std::nullopt;
    }

    /// Whether two models are equal: the same tables, ids, attributes, geometry bit for bit and sharing.
    template <std::floating_point T>
    [[nodiscard]] bool equal(const Model<T> &a, const Model<T> &b)
    {
        return !first_difference(a, b);
    }
}
