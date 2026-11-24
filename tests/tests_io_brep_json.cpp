// Exact writing and reading of a BREP model in JSON (gbs-io/brep_json.h).
// Design: docs/sources/design/io_brep_json.md.
#include <doctest_gtest.hpp>
#include <gbs-io/brep_json.h>
#include <gbs-io/step/reader.h>
#include <gbs-brep/brep>
#include <blade_fixture.h>
#include "tests_brep_geom.h"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <numbers>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;
using namespace gbs::io;
using rapidjson::SizeType;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    // ---- bit for bit ---------------------------------------------------------

    template <std::floating_point R>
    bool same(R a, R b) { return std::memcmp(&a, &b, sizeof(R)) == 0; }
    template <std::floating_point R, std::size_t N>
    bool same(const std::array<R, N> &a, const std::array<R, N> &b) { return std::memcmp(a.data(), b.data(), sizeof(R) * N) == 0; }
    template <typename X>
    bool same(const std::vector<X> &a, const std::vector<X> &b)
    {
        if (a.size() != b.size())
            return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (!same(a[i], b[i]))
                return false;
        return true;
    }
    template <std::floating_point R>
    bool same(const Placement<R> &a, const Placement<R> &b)
    {
        return same(a.r[0], b.r[0]) && same(a.r[1], b.r[1]) && same(a.r[2], b.r[2]) && same(a.t, b.t);
    }
    template <std::floating_point R>
    bool same(const Frame3<R> &a, const Frame3<R> &b) { return same(a.o, b.o) && same(a.x, b.x) && same(a.y, b.y) && same(a.z, b.z); }

    /**
     * Two geometries compared bit for bit, without the writer: the same exact class, the same
     * data (poles, knots, frames, radii...), the same values at samples; and the same sharing:
     * the pairs of objects met form a bijection, so an object shared by several entities (or
     * inside several geometries) on one side is one object shared the same way on the other.
     */
    template <std::floating_point R>
    class SameGeometry
    {
        std::unordered_map<const void *, const void *> a2b_, b2a_;

        /// true the first time the pair is met
        bool pair(const void *a, const void *b)
        {
            const auto ia = a2b_.find(a);
            const auto ib = b2a_.find(b);
            if (ia != a2b_.end() || ib != b2a_.end())
            {
                REQUIRE(ia != a2b_.end());
                REQUIRE(ib != b2a_.end());
                REQUIRE(ia->second == b);
                REQUIRE(ib->second == a);
                ++shared;
                return false;
            }
            REQUIRE(a != b); // a copy, not the object itself
            a2b_.emplace(a, b);
            b2a_.emplace(b, a);
            return true;
        }

        template <std::size_t D>
        static void sample(const Curve<R, D> &a, const Curve<R, D> &b)
        {
            const auto [u1, u2] = a.bounds();
            if (!std::isfinite(u2 - u1))
                return;
            for (int i = 0; i <= 8; ++i)
                for (std::size_t d = 0; d <= 1; ++d)
                {
                    const R u = u1 + (u2 - u1) * R(i) / R(8);
                    std::optional<std::array<R, D>> pa, pb;
                    try { pa = a.value(u, d); } catch (const std::exception &) {}
                    try { pb = b.value(u, d); } catch (const std::exception &) {}
                    REQUIRE(pa.has_value() == pb.has_value());
                    if (pa)
                        REQUIRE(same(*pa, *pb));
                }
        }
        template <std::size_t D>
        static void sample(const Surface<R, D> &a, const Surface<R, D> &b)
        {
            const auto bx = a.bounds();
            if (!std::isfinite(bx[1] - bx[0]) || !std::isfinite(bx[3] - bx[2]))
                return;
            for (int i = 0; i <= 4; ++i)
                for (int j = 0; j <= 4; ++j)
                    for (std::size_t d = 0; d <= 2; ++d)
                    {
                        const R u = bx[0] + (bx[1] - bx[0]) * R(i) / R(4), v = bx[2] + (bx[3] - bx[2]) * R(j) / R(4);
                        std::optional<std::array<R, D>> pa, pb;
                        try { pa = a.value(u, v, d == 1, d == 2); } catch (const std::exception &) {}
                        try { pb = b.value(u, v, d == 1, d == 2); } catch (const std::exception &) {}
                        REQUIRE(pa.has_value() == pb.has_value());
                        if (pa)
                            REQUIRE(same(*pa, *pb));
                    }
        }
        template <typename X, typename B>
        static auto as(const B &b) { return dynamic_cast<const X *>(&b); }

        template <std::size_t N, bool Rat>
        static void bspline(const BSCurveGeneral<R, N, Rat> &a, const BSCurveGeneral<R, N, Rat> &b)
        {
            REQUIRE(a.degree() == b.degree());
            REQUIRE(same(a.knotsFlats(), b.knotsFlats()));
            REQUIRE(same(a.poles(), b.poles()));
            REQUIRE(same(a.bounds(), b.bounds()));
        }

    public:
        std::size_t shared{}; ///< pairs met again: the sharing found on both sides
        std::size_t objects() const { return a2b_.size(); }

        template <std::size_t D>
        void curve(const Curve<R, D> *a, const Curve<R, D> *b)
        {
            REQUIRE((a == nullptr) == (b == nullptr));
            if (!a)
                return;
            REQUIRE(typeid(*a) == typeid(*b));
            if (!pair(dynamic_cast<const void *>(a), dynamic_cast<const void *>(b)))
                return;
            REQUIRE(same(a->bounds(), b->bounds()));
            if (auto x = as<BSCurve<R, D>>(*a))
                bspline(*x, *as<BSCurve<R, D>>(*b));
            else if (auto x = as<BSCurveRational<R, D>>(*a))
                bspline(*x, *as<BSCurveRational<R, D>>(*b));
            else if (auto x = as<CurveTrimmed<R, D>>(*a))
                curve<D>(x->basisCurve().get(), as<CurveTrimmed<R, D>>(*b)->basisCurve().get());
            else if (auto x = as<CurveComposite<R, D>>(*a))
            {
                const auto &ca = x->curves();
                const auto &cb = as<CurveComposite<R, D>>(*b)->curves();
                REQUIRE(ca.size() == cb.size());
                for (std::size_t i = 0; i < ca.size(); ++i)
                    curve<D>(ca[i].get(), cb[i].get());
            }
            else if (auto x = as<CurveReversed<R, D>>(*a))
                curve<D>(x->basisCurve().get(), as<CurveReversed<R, D>>(*b)->basisCurve().get());
            else if (auto x = as<CurveReparametrized<R, D>>(*a))
            {
                const auto y = as<CurveReparametrized<R, D>>(*b);
                curve<D>(x->basisCurve().get(), y->basisCurve().get());
                bspline(x->paramFunction().basisCurve(), y->paramFunction().basisCurve());
            }
            else if (auto x = as<CurveExtended<R, D>>(*a))
            {
                const auto y = as<CurveExtended<R, D>>(*b);
                REQUIRE(x->isClamped() == y->isClamped());
                curve<D>(x->basisCurve().get(), y->basisCurve().get());
            }
            else if (auto x = as<CurveOnSurface<R, D>>(*a))
            {
                const auto y = as<CurveOnSurface<R, D>>(*b);
                curve<2>(x->p_basisCurve().get(), y->p_basisCurve().get());
                surface<D>(x->p_basisSurface().get(), y->p_basisSurface().get());
            }
            if constexpr (D == 3)
            {
                if (auto x = as<LineCurve<R>>(*a))
                {
                    REQUIRE(same(x->origin(), as<LineCurve<R>>(*b)->origin()));
                    REQUIRE(same(x->vector(), as<LineCurve<R>>(*b)->vector()));
                }
                else if (auto x = as<EllipseCurve<R>>(*a))
                {
                    const auto y = as<EllipseCurve<R>>(*b);
                    REQUIRE(same(x->radius1(), y->radius1()));
                    REQUIRE(same(x->radius2(), y->radius2()));
                    REQUIRE(same(x->frame(), y->frame()));
                }
                else if (auto x = as<MovedCurve<R>>(*a))
                {
                    REQUIRE(same(x->placement(), as<MovedCurve<R>>(*b)->placement()));
                    curve<3>(x->basis().get(), as<MovedCurve<R>>(*b)->basis().get());
                }
                else if (auto x = as<CurveOffset3D<R, BSCfunction<R>>>(*a))
                {
                    const auto y = as<CurveOffset3D<R, BSCfunction<R>>>(*b);
                    REQUIRE(same(x->direction(), y->direction()));
                    bspline(x->offset().basisCurve(), y->offset().basisCurve());
                    curve<3>(&x->basisCurve(), &y->basisCurve());
                }
            }
            else if (auto x = as<CurveOffset2D<R, BSCfunction<R>>>(*a))
            {
                const auto y = as<CurveOffset2D<R, BSCfunction<R>>>(*b);
                bspline(x->offset().basisCurve(), y->offset().basisCurve());
                curve<2>(&x->basisCurve(), &y->basisCurve());
            }
            sample(*a, *b);
        }

        template <std::size_t D>
        void surface(const Surface<R, D> *a, const Surface<R, D> *b)
        {
            REQUIRE((a == nullptr) == (b == nullptr));
            if (!a)
                return;
            REQUIRE(typeid(*a) == typeid(*b));
            if (!pair(dynamic_cast<const void *>(a), dynamic_cast<const void *>(b)))
                return;
            REQUIRE(same(a->bounds(), b->bounds()));
            auto bspline_surface = [](const auto &x, const auto &y) {
                REQUIRE(x.degreeU() == y.degreeU());
                REQUIRE(x.degreeV() == y.degreeV());
                REQUIRE(same(x.knotsFlatsU(), y.knotsFlatsU()));
                REQUIRE(same(x.knotsFlatsV(), y.knotsFlatsV()));
                REQUIRE(same(x.poles(), y.poles()));
            };
            if (auto x = as<BSSurface<R, D>>(*a))
                bspline_surface(*x, *as<BSSurface<R, D>>(*b));
            else if (auto x = as<BSSurfaceRational<R, D>>(*a))
                bspline_surface(*x, *as<BSSurfaceRational<R, D>>(*b));
            if constexpr (D == 3)
            {
                if (auto x = as<ElementarySurface<R>>(*a))
                {
                    REQUIRE(same(x->frame(), as<ElementarySurface<R>>(*b)->frame()));
                    if (auto c = as<CylindricalSurface<R>>(*a))
                        REQUIRE(same(c->radius(), as<CylindricalSurface<R>>(*b)->radius()));
                    if (auto c = as<ConicalSurface<R>>(*a))
                    {
                        REQUIRE(same(c->radius(), as<ConicalSurface<R>>(*b)->radius()));
                        REQUIRE(same(c->semi_angle(), as<ConicalSurface<R>>(*b)->semi_angle()));
                    }
                    if (auto c = as<SphericalSurface<R>>(*a))
                        REQUIRE(same(c->radius(), as<SphericalSurface<R>>(*b)->radius()));
                    if (auto c = as<ToroidalSurface<R>>(*a))
                    {
                        REQUIRE(same(c->major_radius(), as<ToroidalSurface<R>>(*b)->major_radius()));
                        REQUIRE(same(c->minor_radius(), as<ToroidalSurface<R>>(*b)->minor_radius()));
                    }
                }
                else if (auto x = as<SurfaceOfRevolution<R>>(*a))
                {
                    const auto y = as<SurfaceOfRevolution<R>>(*b);
                    REQUIRE(same(x->axis2()[0], y->axis2()[0]));
                    REQUIRE(same(x->axis2()[1], y->axis2()[1]));
                    REQUIRE(same(x->axis2()[2], y->axis2()[2]));
                    curve<2>(x->basisCurve().get(), y->basisCurve().get());
                }
                else if (auto x = as<OffsetSurface<R>>(*a))
                {
                    REQUIRE(same(x->distance(), as<OffsetSurface<R>>(*b)->distance()));
                    surface<3>(x->basis().get(), as<OffsetSurface<R>>(*b)->basis().get());
                }
                else if (auto x = as<MovedSurface<R>>(*a))
                {
                    REQUIRE(same(x->placement(), as<MovedSurface<R>>(*b)->placement()));
                    surface<3>(x->basis().get(), as<MovedSurface<R>>(*b)->basis().get());
                }
                else if (auto x = as<SeamlessSurface<R>>(*a))
                {
                    REQUIRE(same(x->periods(), as<SeamlessSurface<R>>(*b)->periods()));
                    surface<3>(x->basis().get(), as<SeamlessSurface<R>>(*b)->basis().get());
                }
            }
            sample(*a, *b);
        }
    };

    template <ShapeType K, std::floating_point R>
    void same_attributes(const Model<R> &a, const Model<R> &b, Id<K> h)
    {
        REQUIRE(a.name(h) == b.name(h));
        REQUIRE(a.externalId(h) == b.externalId(h));
        const auto ia = a.instance(h), ib = b.instance(h);
        REQUIRE(ia.has_value() == ib.has_value());
        if (ia)
        {
            REQUIRE(ia->part == ib->part);
            REQUIRE(same(ia->place, ib->place));
        }
    }

    template <ShapeType K, std::floating_point R>
    void same_table(const Model<R> &a, const Model<R> &b, SameGeometry<R> &g)
    {
        REQUIRE(a.template capacity<K>() == b.template capacity<K>());
        REQUIRE(a.template count<K>() == b.template count<K>());
        for (std::uint32_t i = 0; i < a.template capacity<K>(); ++i)
        {
            const Id<K> h{i};
            REQUIRE(a.alive(h) == b.alive(h));
            if (!a.alive(h))
                continue;
            if constexpr (K == ShapeType::Vertex)
            {
                REQUIRE(same(a.vertex(h).pnt, b.vertex(h).pnt));
                REQUIRE(same(a.vertex(h).tol, b.vertex(h).tol));
            }
            else if constexpr (K == ShapeType::Edge)
            {
                const auto &x = a.edge(h), &y = b.edge(h);
                g.template curve<3>(x.curve.get(), y.curve.get());
                REQUIRE(same(x.u1, y.u1));
                REQUIRE(same(x.u2, y.u2));
                REQUIRE(x.v1 == y.v1);
                REQUIRE(x.v2 == y.v2);
                REQUIRE(same(x.tol, y.tol));
                REQUIRE(x.degenerate == y.degenerate);
                REQUIRE(x.same_parameter == y.same_parameter);
            }
            else if constexpr (K == ShapeType::Wire)
            {
                const auto &x = a.wire(h), &y = b.wire(h);
                REQUIRE(x.closed == y.closed);
                REQUIRE(x.coedges.size() == y.coedges.size());
                for (std::size_t k = 0; k < x.coedges.size(); ++k)
                {
                    REQUIRE(x.coedges[k].edge == y.coedges[k].edge);
                    REQUIRE(x.coedges[k].orient == y.coedges[k].orient);
                    g.template curve<2>(x.coedges[k].pcurve.get(), y.coedges[k].pcurve.get());
                }
            }
            else if constexpr (K == ShapeType::Face)
            {
                const auto &x = a.face(h), &y = b.face(h);
                g.template surface<3>(x.surface.get(), y.surface.get());
                REQUIRE(x.wires == y.wires);
                REQUIRE(same(x.tol, y.tol));
                REQUIRE(x.natural_bounds == y.natural_bounds);
            }
            else if constexpr (K == ShapeType::Shell)
            {
                const auto &x = a.shell(h), &y = b.shell(h);
                REQUIRE(x.closed == y.closed);
                REQUIRE(x.faces.size() == y.faces.size());
                for (std::size_t k = 0; k < x.faces.size(); ++k)
                {
                    REQUIRE(x.faces[k].face == y.faces[k].face);
                    REQUIRE(x.faces[k].orient == y.faces[k].orient);
                }
            }
            else if constexpr (K == ShapeType::Solid)
            {
                REQUIRE(a.solid(h).outer == b.solid(h).outer);
                REQUIRE(a.solid(h).voids == b.solid(h).voids);
            }
            else
                REQUIRE(a.compound(h).shapes == b.compound(h).shapes);
            same_attributes(a, b, h);
        }
    }

    /// The model b is the model a, entity by entity and bit for bit, sharing included.
    template <std::floating_point R>
    auto expect_same(const Model<R> &a, const Model<R> &b) -> SameGeometry<R>
    {
        SameGeometry<R> g;
        REQUIRE(same(a.unitScale(), b.unitScale()));
        REQUIRE(a.nextPart() == b.nextPart());
        same_table<ShapeType::Vertex>(a, b, g);
        same_table<ShapeType::Edge>(a, b, g);
        same_table<ShapeType::Wire>(a, b, g);
        same_table<ShapeType::Face>(a, b, g);
        same_table<ShapeType::Shell>(a, b, g);
        same_table<ShapeType::Solid>(a, b, g);
        same_table<ShapeType::Compound>(a, b, g);
        return g;
    }

    /// Written, read back, compared; written again: the same text.
    template <std::floating_point R>
    auto round_trip(const Model<R> &m) -> SameGeometry<R>
    {
        const auto text = write_brep(m);
        const auto back = read_brep<R>(text);
        auto g = expect_same(m, back);
        REQUIRE(write_brep(back) == text);
        // the indented text is the same model
        const auto pretty = write_brep(m, {.pretty = true});
        REQUIRE(pretty != text);
        REQUIRE(write_brep(read_brep<R>(pretty)) == text);
        return g;
    }

    /// The number of geometries the file holds.
    std::size_t pool_size(const std::string &text)
    {
        rapidjson::Document d;
        d.Parse(text.data(), text.size());
        return d["geometry"].Size();
    }

    std::string step_dir() { return std::string(GBS_TESTS_IN_DIR) + "/step"; }

    // ---- models ----------------------------------------------------------------

    std::shared_ptr<BSCurve<T, 2>> bs2(T y = 0.)
    {
        return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0., y}, {0.3, y + 0.2}, {0.7, y - 0.1}, {1., y}}, std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, 2);
    }
    std::shared_ptr<BSCurve<T, 3>> bs3()
    {
        return std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {1., 0.3, 0.1}, {2., -0.2, 0.3}, {3., 0., 0.}},
                                               std::vector<T>{0., 0., 0., 0., 1., 1., 1., 1.}, 3);
    }
    BSCfunction<T> law() { return BSCfunction<T>{std::vector<T>{0.05, 0.1, 0.07}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2}; }
    std::shared_ptr<BSSurface<T, 3>> bss3()
    {
        points_vector<T, 3> p;
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 4; ++i)
                p.push_back({T(i), T(j), 0.1 * T(i * j) - 0.05 * T(i)});
        return std::make_shared<BSSurface<T, 3>>(p, std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, std::vector<T>{0., 0., 1., 2., 2.}, 2, 1);
    }

    /// One object of every class the format holds, in edges, co-edges and faces, some shared,
    /// some nested: the entities hold them, the topology does not need to make sense.
    Model<T> every_geometry()
    {
        Model<T> m;
        const auto v0 = m.add(Vertex<T>{{0., 0., 0.}, 1e-7});
        const auto v1 = m.add(Vertex<T>{{1., 0., 0.}, 2e-7});
        auto edge = [&](std::shared_ptr<Curve<T, 3>> c) {
            const auto [u1, u2] = c->bounds();
            const bool finite = std::isfinite(u2 - u1);
            return m.add(Edge<T>{c, finite ? u1 : 0., finite ? u2 : 1., v0, v1, 1e-7});
        };
        // 2D curves
        auto b2 = bs2();
        std::vector<std::shared_ptr<Curve<T, 2>>> c2{
            b2,
            std::make_shared<BSCurveRational<T, 2>>(points_vector<T, 2>{{0., 0.}, {1., 1.}, {2., 0.}}, std::vector<T>{0., 0., 0., 1., 1., 1.},
                                                    std::vector<T>{1., std::sqrt(0.5), 1.}, 2),
            std::make_shared<CurveTrimmed<T, 2>>(b2, 0.2, 0.7),
            std::make_shared<CurveReversed<T, 2>>(b2),
            std::make_shared<CurveReparametrized<T, 2>>(b2, BSCfunction<T>{std::vector<T>{0., 0.4, 1.}, std::vector<T>{0., 0., 0., 2., 2., 2.}, 2}),
            std::make_shared<CurveExtended<T, 2>>(b2, true),
            std::make_shared<CurveOffset2D<T, BSCfunction<T>>>(b2, law()),
        };
        c2.push_back(std::make_shared<CurveComposite<T, 2>>(std::vector<std::shared_ptr<Curve<T, 2>>>{c2[2], c2[3]}));
        {
            points_vector<T, 2> p{{0., 0.}, {1., 0.}, {0., 1.}, {1.2, 1.1}};
            auto s2 = std::make_shared<BSSurface<T, 2>>(p, std::vector<T>{0., 0., 1., 1.}, std::vector<T>{0., 0., 1., 1.}, 1, 1);
            c2.push_back(std::make_shared<CurveOnSurface<T, 2>>(bs2(0.5), s2));
        }
        // 3D surfaces
        Frame3<T> fr{ax2<T, 3>{{{1., 2., 3.}, {0.1, 0.2, 1.}, {1., 0., 0.}}}};
        auto plane = std::make_shared<Plane<T>>(fr, std::array<T, 4>{-1., 2., -3., 4.});
        auto b3s = bss3();
        auto trimmed_u = bss3(); // U bounds an ulp inside the end knots, kept by trimU
        trimmed_u->trimU(std::nextafter(0., 1.), std::nextafter(1., 0.));
        REQUIRE(trimmed_u->bounds()[0] > 0.);
        auto nurbs_s = std::make_shared<BSSurfaceRational<T, 3>>(
            points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 1.}}, std::vector<T>{1., 0.8, 0.9, 1.2}, std::vector<T>{0., 0., 1., 1.},
            std::vector<T>{0., 0., 1., 1.}, 1, 1);
        auto revolution = std::make_shared<SurfaceOfRevolution<T>>(bs2(1.), ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, 0.25, 2.5);
        Placement<T> turn = Placement<T>::frame({1., -2., 0.5}, {std::cos(0.3), std::sin(0.3), 0.}, {-std::sin(0.3), std::cos(0.3), 0.}, {0., 0., 1.});
        auto moved_revolution = std::make_shared<MovedSurface<T>>(revolution, turn);
        std::vector<std::shared_ptr<Surface<T, 3>>> s3{
            plane,
            std::make_shared<CylindricalSurface<T>>(2., fr, std::array<T, 4>{-pi, 3. * pi, 0., 5.})->with_box({0.5, 7., -1., 6.}),
            std::make_shared<ConicalSurface<T>>(1.5, 0.3, fr, std::array<T, 4>{0., 2. * pi, -1., 2.}),
            // a sphere a little past its poles, as with_box leaves it (its constructor refuses it)
            SphericalSurface<T>(3., fr, std::array<T, 4>{0., pi, -1., 1.}).with_box({0., 2. * pi, -0.5 * pi * (1. + 1e-11), 0.5 * pi}),
            std::make_shared<ToroidalSurface<T>>(4., 1., fr, std::array<T, 4>{0., pi, -pi, pi}),
            b3s,
            trimmed_u,
            nurbs_s,
            revolution,
            std::make_shared<OffsetSurface<T>>(b3s, 0.2),
            std::make_shared<OffsetSurface<T>>(moved_revolution, -0.1),
            moved_revolution,
            std::make_shared<SeamlessSurface<T>>(b3s, std::array<T, 2>{1., 0.}),
        };
        // 3D curves
        auto b3 = bs3();
        auto unclamped_bounds = bs3();
        unclamped_bounds->trim(0.1, 0.9, false); // the bounds only
        Frame3<T> moved_frame = fr.moved([&](const point<T, 3> &p) { return turn.apply(p); }); // not quite orthonormal
        auto ellipse = std::make_shared<EllipseCurve<T>>(2., 1., moved_frame, 5.9, 6.6);
        auto composite3 = std::make_shared<CurveComposite<T, 3>>(
            std::vector<std::shared_ptr<Curve<T, 3>>>{std::make_shared<CurveTrimmed<T, 3>>(b3, 0., 0.5), std::make_shared<CurveTrimmed<T, 3>>(b3, 0.5, 1.)});
        std::vector<std::shared_ptr<Curve<T, 3>>> c3{
            b3,
            unclamped_bounds,
            std::make_shared<BSCurveRational<T, 3>>(build_ellipse_arc<T, 3>(1., 2., 0.2, 2.9, {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.})),
            std::make_shared<LineCurve<T>>(point<T, 3>{1., 2., 3.}, point<T, 3>{0.1, -0.7, 0.3}, -2., 5.),
            ellipse,
            composite3,
            std::make_shared<MovedCurve<T>>(composite3, turn),
            std::make_shared<CurveReversed<T, 3>>(ellipse),
            std::make_shared<CurveExtended<T, 3>>(b3, false),
            std::make_shared<CurveReparametrized<T, 3>>(b3, law()),
            std::make_shared<CurveOffset3D<T, BSCfunction<T>>>(b3, law(), std::array<T, 3>{0., 0., 1.}),
            std::make_shared<CurveOnSurface<T, 3>>(c2[0], plane),     // the surface of a face, the pcurve of a co-edge
            std::make_shared<CurveOnSurface<T, 3>>(bs2(0.3), s3.back()), // a fillet's contact on a seamless view
        };
        for (const auto &c : c3)
            edge(c);
        edge(b3); // shared by two edges
        // a wire holding every 2D curve, one twice
        Wire<T> w;
        for (std::size_t k = 0; k < c2.size(); ++k)
            w.coedges.push_back({EdgeId{std::uint32_t(k % 4)}, k % 2 ? Orientation::Reversed : Orientation::Forward, c2[k]});
        w.coedges.push_back({EdgeId{0}, Orientation::Forward, c2[0]});
        w.coedges.push_back({EdgeId{1}, Orientation::Forward, nullptr}); // a free co-edge
        const auto wid = m.add(w);
        for (const auto &s : s3)
            m.add(Face<T>{s, {wid}, 1e-6, false});
        m.add(Face<T>{plane, {}, 1e-6, true}); // a surface shared by two faces
        return m;
    }

    /// The edge of the solid whose middle is nearest to p.
    EdgeId edge_near(const Model<T> &m, SolidId s, point<T, 3> p)
    {
        EdgeId best;
        T d = std::numeric_limits<T>::infinity();
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (!ed.curve)
                continue;
            if (const T x = norm(ed.curve->value(0.5 * (ed.u1 + ed.u2)) - p); x < d)
                d = x, best = e;
        }
        return best;
    }

    SolidId block(Model<T> &m, const point<T, 3> &a, const point<T, 3> &b)
    {
        const auto f = unwrap(make_face(m, brep_tests::quad({a[0], a[1], a[2]}, {b[0], a[1], a[2]}, {a[0], b[1], a[2]}, {b[0], b[1], a[2]})));
        return explore<SolidId>(m, unwrap(make_prism(m, ShapeId{f}, point<T, 3>{0., 0., b[2] - a[2]})).shape).front();
    }

    /// A text of a small model, with geometry of several kinds, for the corrupt files.
    std::string small_text()
    {
        Model<T> m;
        const auto s = block(m, {0., 0., 0.}, {2., 1., 1.});
        m.setName(ShapeId{s}, "block");
        return write_brep(m);
    }

    /// The text with one change, made on its JSON document.
    template <typename F>
    std::string edited(const std::string &text, F &&f)
    {
        rapidjson::Document d;
        d.Parse<rapidjson::kParseFullPrecisionFlag>(text.data(), text.size());
        REQUIRE(!d.HasParseError());
        f(d, d.GetAllocator());
        rapidjson::StringBuffer b;
        rapidjson::Writer<rapidjson::StringBuffer> w(b);
        d.Accept(w);
        return {b.GetString(), b.GetSize()};
    }

    /// read_brep refuses the text with this code, at this place.
    void expect_refused(const std::string &text, BrepReadErrc code, std::string_view where, std::string_view reason_part = {})
    {
        try
        {
            (void)read_brep<T>(text, "file.brep");
        }
        catch (const BrepReadError &e)
        {
            INFO(e.what());
            REQUIRE(e.code() == code);
            REQUIRE(e.where() == where);
            REQUIRE(std::string_view{e.what()}.starts_with("file.brep: "));
            if (!reason_part.empty())
                REQUIRE(e.reason().find(reason_part) != std::string_view::npos);
            return;
        }
        FAIL("the text was read");
    }
}

TEST(tests_io_brep_json, empty_model)
{
    Model<T> m;
    const auto text = write_brep(m);
    ASSERT_EQ(text, R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":1.0,"next_part":1,"geometry":[],"vertices":[],"edges":[],)"
                    R"("wires":[],"faces":[],"shells":[],"solids":[],"compounds":[]})");
    round_trip(m);
}

TEST(tests_io_brep_json, every_geometry_class)
{
    const auto m = every_geometry();
    const auto text = write_brep(m);
    const auto g = round_trip(m);
    // each object written once: the plane of two faces and of a curve on surface, the B-spline of
    // two edges, of the trimmed and composite curves, the pcurve of two co-edges...
    ASSERT_GT(g.shared, 10u);
    ASSERT_EQ(pool_size(text), g.objects());
    // read back, the first and last edges share their curve, as written
    const auto back = read_brep<T>(text);
    const auto last = EdgeId{static_cast<std::uint32_t>(back.capacity<ShapeType::Edge>() - 1)};
    ASSERT_EQ(back.edge(EdgeId{0}).curve.get(), back.edge(last).curve.get());
}

TEST(tests_io_brep_json, dead_entities_ids_and_attributes)
{
    Model<T> m;
    const auto s = block(m, {0., 0., 0.}, {2., 1., 1.});
    const auto t = block(m, {3., 0., 0.}, {4., 1., 2.});
    const auto c = m.add(Compound{{ShapeId{s}, ShapeId{t}, ShapeId{VertexId{2}}}});
    // attributes of every kind, bytes of names as they are
    m.setName(ShapeId{s}, "pièce \"A\"\\ \n\t\x01 é");
    m.setName(ShapeId{c}, std::string("nul\0inside", 10));
    m.setName(ShapeId{FaceId{3}}, "a face");
    m.setExternalId(ShapeId{s}, std::numeric_limits<std::int64_t>::min());
    m.setExternalId(ShapeId{EdgeId{5}}, std::numeric_limits<std::int64_t>::max());
    const auto part = m.newPart();
    (void)m.newPart(); // a key given, not used
    m.setInstance(ShapeId{s}, InstanceOf<T>{part, Placement<T>::frame({1., 2., 3.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.})});
    m.setUnitScale(25.4);
    // dead entities everywhere, some still referred to (erase does not cascade), invalid ids
    m.erase(ShapeId{t});
    m.erase(ShapeId{VertexId{2}});
    m.erase(ShapeId{EdgeId{7}});
    m.erase(ShapeId{FaceId{1}});
    m.add(Edge<T>{nullptr, 0., 1., VertexId{}, VertexId{0}, 1e-6, true}); // degenerate, an invalid id
    m.add(Shell{});
    m.erase(ShapeId{ShellId{static_cast<std::uint32_t>(m.capacity<ShapeType::Shell>() - 1)}});
    m.add(Solid{ShellId{}, {}});
    const auto text = write_brep(m);
    round_trip(m);
    const auto back = read_brep<T>(text);
    ASSERT_FALSE(back.alive(ShapeId{t}));
    ASSERT_TRUE(back.alive(ShapeId{s}));
    ASSERT_EQ(back.compound(c).shapes[1], ShapeId{t}); // still referred to, dead
    ASSERT_EQ(back.name(ShapeId{c}), std::string("nul\0inside", 10));
    ASSERT_EQ(back.unitScale(), 25.4);
    // the part keys go on from where they were
    auto copy = back;
    ASSERT_EQ(copy.newPart(), 3u);
    // compacted before writing: no dead entity, the ids renumbered by the remap
    auto compacted = m;
    const auto remap = compacted.compact();
    const auto back2 = read_brep<T>(write_brep(compacted));
    expect_same(compacted, back2);
    ASSERT_EQ(back2.capacity<ShapeType::Solid>(), back2.count<ShapeType::Solid>());
    ASSERT_EQ(back2.name(remap.map(ShapeId{s})), m.name(ShapeId{s}));
}

TEST(tests_io_brep_json, numbers_bit_for_bit)
{
    // reals of every kind: random bits (subnormals included), -0, infinities, the extremes
    Model<T> m;
    std::mt19937_64 g{12345};
    std::vector<T> x{-0., 0., std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::denorm_min(),
                     std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::min(), 0.1, 1. / 3.};
    while (x.size() < 30000)
    {
        const auto d = std::bit_cast<T>(g());
        if (!std::isnan(d))
            x.push_back(d);
        x.push_back(std::bit_cast<T>(g() & 0x000FFFFFFFFFFFFFull)); // subnormal
    }
    for (std::size_t i = 0; i + 3 < x.size(); i += 4)
        m.add(Vertex<T>{{x[i], x[i + 1], x[i + 2]}, x[i + 3]});
    m.add(Vertex<T>{{std::numeric_limits<T>::quiet_NaN(), 1., 2.}, 1e-7});
    const auto back = read_brep<T>(write_brep(m));
    for (std::uint32_t i = 0; i + 1 < m.capacity<ShapeType::Vertex>(); ++i)
    {
        ASSERT_TRUE(same(m.vertex(VertexId{i}).pnt, back.vertex(VertexId{i}).pnt));
        ASSERT_TRUE(same(m.vertex(VertexId{i}).tol, back.vertex(VertexId{i}).tol));
    }
    // a NaN comes back a NaN (its sign and payload are not kept)
    ASSERT_TRUE(std::isnan(back.vertex(VertexId{static_cast<std::uint32_t>(m.capacity<ShapeType::Vertex>() - 1)}).pnt[0]));
}

TEST(tests_io_brep_json, float_model)
{
    Model<float> m;
    const auto v = m.add(Vertex<float>{{0.1f, std::nextafter(1.f, 2.f), -0.f}, 1e-5f});
    auto c = std::make_shared<BSCurve<float, 3>>(points_vector<float, 3>{{0.f, 0.f, 0.f}, {1.f / 3.f, 0.2f, 0.f}, {1.f, 0.f, 0.f}},
                                                 std::vector<float>{0.f, 0.f, 0.f, 1.f, 1.f, 1.f}, 2);
    m.add(Edge<float>{c, 0.f, 1.f, v, v, 1e-5f});
    const auto text = write_brep(m);
    ASSERT_NE(text.find(R"("real":"float")"), std::string::npos);
    round_trip(m);
    // a float model is not read as a double one, nor the reverse
    expect_refused(text, BrepReadErrc::Invalid, "/real");
    try
    {
        (void)read_brep<float>(write_brep(Model<T>{}));
        FAIL("read");
    }
    catch (const BrepReadError &e)
    {
        ASSERT_EQ(e.where(), "/real");
    }
}

TEST(tests_io_brep_json, booleans_and_fillets)
{
    Model<T> m;
    // a block minus a cylinder (STEP) and a sphere: primitives, B-splines, section curves on surfaces
    const auto bx = ShapeId{block(m, {-6., -6., 1.}, {6., 6., 4.})};
    const auto cyl = unwrap(step::read_step_file(m, step_dir() + "/cylinder.stp")).root;
    const auto cut1 = unwrap(cut(m, bx, cyl));
    // a box 0..4 minus the sphere R 3 at its corner, the axis of the sphere along (1, 1, 1) so that
    // its poles and seam stay off the edges of the box (as in tests_brep_boolean)
    const auto bx2 = ShapeId{unwrap(make_solid(m, unwrap(sew(m, brep_tests::box_faces(m, {20., 0., 0.}, 4.))).shells.front()))};
    const auto sph = unwrap(step::read_step_file(m, step_dir() + "/sphere.stp")).root;
    const point<T, 3> z{1. / std::sqrt(3.), 1. / std::sqrt(3.), 1. / std::sqrt(3.)}, x{1. / std::sqrt(2.), -1. / std::sqrt(2.), 0.};
    const auto cut2 = unwrap(cut(m, bx2, copy_shape(m, sph, Placement<T>::frame({24., 4., 4.}, x, cross(z, x), z))));
    // a fillet along an edge of a block
    const auto s = block(m, {10., 0., 0.}, {12., 1.5, 1.});
    const auto f = unwrap(fillet(m, s, std::vector<EdgeId>{edge_near(m, s, {12., 0.75, 1.})}, 0.3));
    m.setName(ShapeId{f}, "filleted");
    // the root of a blade on its hub: the fillet runs across the seam of a closed B-spline, its
    // contacts on a seamless view of the blade (SeamlessSurface)
    const auto [blade, root] = blade_fixture::blade_on_hub(m);
    const auto bf = unwrap(fillet(m, blade, root, 0.03));
    std::size_t seamless = 0;
    for (auto e : explore<EdgeId>(m, ShapeId{bf}))
        if (auto c = dynamic_cast<const CurveOnSurface<T, 3> *>(m.edge(e).curve.get()))
            seamless += dynamic_cast<const SeamlessSurface<T> *>(c->p_basisSurface().get()) != nullptr;
    ASSERT_GT(seamless, 0u);
    const auto g = round_trip(m);
    ASSERT_GT(g.shared, 0u);
    // the results read back are valid solids of the same volume
    const auto back = read_brep<T>(write_brep(m));
    for (const ShapeId r : {cut1.shape, cut2.shape, ShapeId{f}, ShapeId{bf}})
    {
        ASSERT_EQ(check(back, r).entries.size(), check(m, r).entries.size());
        for (auto so : explore<SolidId>(back, r))
            ASSERT_TRUE(same(signed_volume(back, back.solid(so).outer, 64), signed_volume(m, m.solid(so).outer, 64)));
    }
}

TEST(tests_io_brep_json, step_files)
{
    // every STEP file of the tests, each in a model of its own, then all in one model with
    // some of them erased
    Model<T> all;
    std::vector<ShapeId> roots;
    std::size_t files{};
    for (const auto &entry : std::filesystem::directory_iterator(step_dir()))
    {
        if (entry.path().extension() != ".stp")
            continue;
        INFO(entry.path().filename().string());
        Model<T> m;
        auto r = step::read_step_file(m, entry.path(), step::StepReadOptions{.target_unit_mm = 1.});
        if (!r)
            continue; // a file the reader refuses holds no model
        ++files;
        round_trip(m);
        roots.push_back(all.append(m).map(r->root));
    }
    ASSERT_GE(files, 20u);
    for (std::size_t k = 0; k < roots.size(); k += 3)
        all.erase(roots[k]);
    round_trip(all);
}

TEST(tests_io_brep_json, files_and_streams)
{
    const auto m = every_geometry();
    const auto dir = std::filesystem::temp_directory_path() / ("gbs_brep_json_" + std::to_string(std::random_device{}()));
    std::filesystem::create_directories(dir);
    const auto path = dir / "modèle.brep";
    write_brep_file(m, path);
    expect_same(m, read_brep_file<T>(path));
    write_brep_file(m, path, {.pretty = true}); // replaced whole
    expect_same(m, read_brep_file<T>(path));
    std::stringstream ss;
    write_brep(m, ss);
    expect_same(m, read_brep<T>(ss));
    // a file missing: filesystem_error; a file cut short: BrepReadError, with its name
    ASSERT_THROW((void)read_brep_file<T>(dir / "missing.brep"), std::filesystem::filesystem_error);
    write_text_atomic(dir / "cut.brep", write_brep(m).substr(0, 100));
    try
    {
        (void)read_brep_file<T>(dir / "cut.brep");
        FAIL("read");
    }
    catch (const BrepReadError &e)
    {
        ASSERT_EQ(e.code(), BrepReadErrc::Syntax);
        ASSERT_NE(std::string_view{e.what()}.find("cut.brep"), std::string_view::npos);
    }
    std::filesystem::remove_all(dir);
}

TEST(tests_io_brep_json, version_and_format)
{
    const auto text = small_text();
    const auto newer = edited(text, [](auto &d, auto &) { d["version"].SetUint(2); });
    expect_refused(newer, BrepReadErrc::NewerVersion, "/version", "version 2");
    // a newer version is reported as such even if what follows is unknown
    expect_refused(edited(text, [](auto &d, auto &a) {
                       d["version"].SetUint(7);
                       d.AddMember("new_table", rapidjson::Value{rapidjson::kArrayType}, a);
                   }),
                   BrepReadErrc::NewerVersion, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d["version"].SetUint(0); }), BrepReadErrc::Invalid, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d.RemoveMember("version"); }), BrepReadErrc::Invalid, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d["format"].SetString("step"); }), BrepReadErrc::NotBrep, "/format");
    expect_refused(edited(text, [](auto &d, auto &) { d.RemoveMember("format"); }), BrepReadErrc::NotBrep, "/format");
    expect_refused("[1, 2]", BrepReadErrc::NotBrep, "");
    expect_refused(R"({"a": 1})", BrepReadErrc::NotBrep, "/format");
}

TEST(tests_io_brep_json, truncated_and_malformed)
{
    const auto text = small_text();
    // every text cut short is refused as such, never read
    for (std::size_t n = 0; n < text.size(); ++n)
        expect_refused(text.substr(0, n), BrepReadErrc::Syntax, "");
    expect_refused(text + "x", BrepReadErrc::Syntax, "");
    expect_refused(text.substr(0, 40) + '\0' + text.substr(40), BrepReadErrc::Syntax, "");
    // a deep nesting: refused without exhausting the stack (the parser is not recursive)
    expect_refused(std::string(1'000'000, '['), BrepReadErrc::Syntax, "");
    auto deep = text;
    const std::string real = R"("real":"double")";
    deep.replace(deep.find(real), real.size(), R"("real":)" + std::string(1'000'000, '[') + std::string(1'000'000, ']'));
    expect_refused(deep, BrepReadErrc::Invalid, "/real");
}

TEST(tests_io_brep_json, corrupt_values)
{
    const auto text = small_text();
    using V = rapidjson::Value;
    rapidjson::Document doc;
    doc.Parse(text.data(), text.size());
    const auto n_geometry = doc["geometry"].Size(), n_compounds = doc["compounds"].Size();
    const auto edge_curve = doc["edges"][0]["curve"].GetUint(), face_surface = doc["faces"][0]["surface"].GetUint();
    SizeType bspline = n_geometry; // the first B-spline curve
    for (SizeType i = 0; i < n_geometry && bspline == n_geometry; ++i)
        if (std::string_view{doc["geometry"][i]["type"].GetString()} == "bspline_curve")
            bspline = i;
    REQUIRE(bspline < n_geometry);
    const auto at_bspline = "/geometry/" + std::to_string(bspline);
    auto refused = [&](auto &&edit, std::string_view where, std::string_view reason = {}) {
        expect_refused(edited(text, edit), BrepReadErrc::Invalid, where, reason);
    };
    refused([](auto &d, auto &) { d["edges"][0]["v"][0].SetUint(999); }, "/edges/0/v/0", "vertex index below");
    refused([](auto &d, auto &) { d["edges"][0]["v"][1].SetInt(-1); }, "/edges/0/v/1");
    refused([](auto &d, auto &) { d["wires"][0]["coedges"][0][0].SetUint(10'000); }, "/wires/0/coedges/0/0");
    refused([](auto &d, auto &) { d["wires"][0]["coedges"][0][1].SetString("X"); }, "/wires/0/coedges/0/1");
    refused([](auto &d, auto &a) { d["compounds"].PushBack(V{rapidjson::kObjectType}, a); }, "/compounds/" + std::to_string(n_compounds),
            "member \"shapes\" missing");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("colour", V{1}, a); }, "/vertices/0", "unknown member \"colour\"");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("tol", V{1.}, a); }, "/vertices/0", "given twice");
    refused([](auto &d, auto &) { d["vertices"][0]["tol"].SetString("small"); }, "/vertices/0/tol", "a number");
    refused([](auto &d, auto &) { d["vertices"][0]["p"].PopBack(); }, "/vertices/0/p", "an array of 3");
    refused([](auto &d, auto &) { d["faces"][0]["natural_bounds"].SetInt(1); }, "/faces/0/natural_bounds");
    refused([](auto &d, auto &) { d["unit_scale"].SetDouble(-1.); }, "/unit_scale");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("name", V{"", 0}, a); }, "/vertices/0/name");
    refused([](auto &d, auto &a) {
        V inst{rapidjson::kObjectType};
        V place{rapidjson::kObjectType};
        V r{rapidjson::kArrayType}, t{rapidjson::kArrayType};
        for (int i = 0; i < 3; ++i)
        {
            V row{rapidjson::kArrayType};
            for (int j = 0; j < 3; ++j)
                row.PushBack(V{i == j ? 1. : 0.}, a);
            r.PushBack(row, a);
            t.PushBack(V{0.}, a);
        }
        place.AddMember("r", r, a).AddMember("t", t, a);
        inst.AddMember("part", V{1u}, a).AddMember("placement", place, a);
        d["solids"][0].AddMember("instance", inst, a);
    }, "/solids/0/instance/part", "below next_part");
    // geometry: of the wrong kind, referring to itself or to one after it, out of its range
    refused([&](auto &d, auto &) { d["edges"][0]["curve"].SetUint(face_surface); }, "/edges/0/curve", "a 3D curve expected");
    refused([&](auto &d, auto &) { d["faces"][0]["surface"].SetUint(edge_curve); }, "/faces/0/surface", "a 3D surface expected");
    refused([&](auto &d, auto &) { d["edges"][0]["curve"].SetUint(d["geometry"].Size()); }, "/edges/0/curve", "below");
    refused(
        [&](auto &d, auto &a) {
            V c{rapidjson::kObjectType};
            V list{rapidjson::kArrayType};
            list.PushBack(V{d["geometry"].Size()}, a); // itself
            c.AddMember("type", "composite_curve", a).AddMember("dim", V{3}, a).AddMember("curves", list, a);
            d["geometry"].PushBack(c, a);
        },
        "/geometry/" + std::to_string(n_geometry) + "/curves/0", "below " + std::to_string(n_geometry));
    refused([&](auto &d, auto &) { d["geometry"][bspline]["type"].SetString("spiral"); }, at_bspline + "/type", "curve of type \"spiral\"");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["dim"].SetUint(4); }, at_bspline + "/dim");
    refused([&](auto &d, auto &) {
        auto &k = d["geometry"][bspline]["knots"];
        std::swap(k[0], k[k.Size() - 1]); // out of order
    }, at_bspline + "/knots", "in order");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["degree"].SetUint(40); }, at_bspline + "/knots");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["poles"].PopBack(); }, at_bspline + "/knots");
}

TEST(tests_io_brep_json, forward_reference_refused)
{
    // a geometry refers only to those before it: no cycle can be read
    const auto text = small_text();
    rapidjson::Document d;
    d.Parse(text.data(), text.size());
    const auto n = d["geometry"].Size();
    const auto t = edited(text, [n](auto &doc, auto &a) {
        rapidjson::Value c{rapidjson::kObjectType};
        c.AddMember("type", "reversed_curve", a).AddMember("dim", rapidjson::Value{3}, a).AddMember("basis", rapidjson::Value{n}, a);
        doc["geometry"].PushBack(c, a);
    });
    expect_refused(t, BrepReadErrc::Invalid, "/geometry/" + std::to_string(n) + "/basis", "below " + std::to_string(n));
}

TEST(tests_io_brep_json, constructors_refusals_reported)
{
    // values a constructor of gbs refuses: reported at the geometry, as invalid
    Model<T> m;
    const auto v = m.add(Vertex<T>{});
    m.add(Edge<T>{std::make_shared<EllipseCurve<T>>(1., 2., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}), 0., 1., v, v});
    const auto text = edited(write_brep(m), [](auto &d, auto &) { d["geometry"][0]["radii"][0].SetDouble(-1.); });
    expect_refused(text, BrepReadErrc::Invalid, "/geometry/0", "radii must be positive");
}

TEST(tests_io_brep_json, random_damage)
{
    // bytes changed at random: the text is read or refused with BrepReadError, nothing else
    const auto text = small_text();
    std::mt19937 g{7};
    const std::string alphabet = "0123456789-+.eE,:[]{}\"nultrfaseNIy x";
    std::size_t read{}, refused{};
    for (int k = 0; k < 3000; ++k)
    {
        auto t = text;
        for (int j = 0; j < 1 + k % 3; ++j)
            t[g() % t.size()] = alphabet[g() % alphabet.size()];
        try
        {
            (void)read_brep<T>(t);
            ++read;
        }
        catch (const BrepReadError &)
        {
            ++refused;
        }
        catch (const std::exception &e)
        {
            INFO(e.what());
            FAIL("an exception other than BrepReadError");
        }
    }
    ASSERT_GT(read, 0u);
    ASSERT_GT(refused, 0u);
}

TEST(tests_io_brep_json, refused_writes)
{
    // a geometry that cannot be read back exactly: an old Line (its constructor normalizes), a
    // class of an application (even derived from a B-spline): refused, with the entity
    struct MyCurve : BSCurve<T, 3>
    {
        MyCurve(const points_vector<T, 3> &p, const std::vector<T> &k, std::size_t d) : BSCurve<T, 3>(p, k, d) {}
    };
    {
        Model<T> m;
        const auto v = m.add(Vertex<T>{});
        m.add(Edge<T>{std::make_shared<Line<T, 3>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 1., 0.}), 0., 1., v, v});
        try
        {
            (void)write_brep(m);
            FAIL("written");
        }
        catch (const BrepWriteError &e)
        {
            ASSERT_NE(std::string_view{e.what()}.find("edge 0"), std::string_view::npos);
        }
        m.edge(EdgeId{0}).curve = std::make_shared<MyCurve>(points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}}, std::vector<T>{0., 0., 1., 1.}, 1);
        ASSERT_THROW((void)write_brep(m), BrepWriteError);
        // inside another geometry too
        m.edge(EdgeId{0}).curve = std::make_shared<CurveReversed<T, 3>>(std::make_shared<Circle3d<T>>(1., ax1<T, 3>{{{0., 0., 0.}, {0., 0., 1.}}}));
        ASSERT_THROW((void)write_brep(m), BrepWriteError);
    }
    {
        // an id beyond its table: nothing written
        Model<T> m;
        m.add(Edge<T>{nullptr, 0., 1., VertexId{3}, VertexId{}, 1e-7, true});
        try
        {
            (void)write_brep(m);
            FAIL("written");
        }
        catch (const BrepWriteError &e)
        {
            ASSERT_NE(std::string_view{e.what()}.find("vertex id 3 beyond its table"), std::string_view::npos);
        }
        const auto dir = std::filesystem::temp_directory_path() / ("gbs_brep_json_w_" + std::to_string(std::random_device{}()));
        std::filesystem::create_directories(dir);
        write_text_atomic(dir / "kept.brep", "old");
        ASSERT_THROW(write_brep_file(m, dir / "kept.brep"), BrepWriteError);
        ASSERT_EQ(read_text(dir / "kept.brep").value(), "old");
        std::filesystem::remove_all(dir);
    }
}
