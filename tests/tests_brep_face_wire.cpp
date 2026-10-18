#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include "tests_brep_geom.h"

#include <cmath>
#include <numbers>
#include <random>
#include <vector>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using namespace brep_tests;
    const T tol = 1e-6;

    // Plane z = 0 whose parameters are x and y on [-2, 2]^2.
    std::shared_ptr<Surface<T, 3>> xy_plane()
    {
        points_vector<T, 3> poles{{-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}};
        std::vector<T> k{-2., -2., 2., 2.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);
    }

    // Bicubic bump on [-2, 2]^2.
    std::shared_ptr<BSSurface<T, 3>> bump()
    {
        points_vector<T, 3> poles;
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 4; ++i)
            {
                const T x = -2. + 4. * i / 3., y = -2. + 4. * j / 3.;
                const T z = (i == 1 || i == 2) && (j == 1 || j == 2) ? 1.5 : 0.;
                poles.push_back({x, y, z});
            }
        std::vector<T> k{-2., -2., -2., -2., 2., 2., 2., 2.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 3, 3);
    }

    std::shared_ptr<Curve<T, 2>> line2d(point<T, 2> a, point<T, 2> b)
    {
        return std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>(a, b));
    }

    std::shared_ptr<Curve<T, 3>> circle3d(T r, point<T, 3> c)
    {
        return std::make_shared<AnalyticCurve>([r, c](T t) { return point<T, 3>{c[0] + r * std::cos(t), c[1] + r * std::sin(t), c[2]}; },
                                               std::array<T, 2>{0., 2. * pi});
    }

    WireId square_wire(Model<T> &m, T h, bool clockwise = false)
    {
        point<T, 3> a{-h, -h, 0.}, b{h, -h, 0.}, c{h, h, 0.}, d{-h, h, 0.};
        std::vector<EdgeId> e = clockwise
                                    ? std::vector{unwrap(make_edge(m, a, d)), unwrap(make_edge(m, d, c)), unwrap(make_edge(m, c, b)), unwrap(make_edge(m, b, a))}
                                    : std::vector{unwrap(make_edge(m, a, b)), unwrap(make_edge(m, b, c)), unwrap(make_edge(m, c, d)), unwrap(make_edge(m, d, a))};
        return unwrap(make_wire(m, e));
    }

    // S(pcurve(t)) agrees with the 3D co-edge point within the edge tolerance, for every co-edge of the face.
    void check_pcurves(const Model<T> &m, FaceId fid, size_t n = 41)
    {
        const auto &f = m.face(fid);
        for (auto wid : f.wires)
        {
            REQUIRE(is_closed(m, wid));
            const auto &w = m.wire(wid);
            for (size_t i = 0; i < w.coedges.size(); ++i)
            {
                const auto &ce = w.coedges[i];
                const auto &e = m.edge(ce.edge);
                REQUIRE(ce.pcurve);
                CHECK(e.tol >= f.tol);
                CHECK(m.vertex(e.v1).tol >= e.tol);
                CHECK(m.vertex(e.v2).tol >= e.tol);
                for (size_t j = 0; j < n; ++j)
                {
                    const T t = e.u1 + (e.u2 - e.u1) * T(j) / T(n - 1);
                    const auto uv = coedge_uv(m, ce, t);
                    CHECK(distance(f.surface->value(uv[0], uv[1]), coedge_point(m, wid, i, t)) <= e.tol * (1. + 1e-9));
                }
            }
        }
    }

    size_t alive_wires_with_pcurves(const Model<T> &m)
    {
        size_t n = 0;
        for (auto wid : m.ids<WireId>())
            for (const auto &ce : m.wire(wid).coedges)
                n += ce.pcurve != nullptr;
        return n;
    }
}

TEST(tests_brep_face_wire, disk_on_plane)
{
    Model<T> m;
    auto srf = xy_plane();
    auto e = unwrap(make_edge(m, std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(1.))));
    auto w = unwrap(make_wire(m, std::vector{e}));
    auto f = unwrap(make_face(m, srf, w));

    const auto &face = m.face(f);
    ASSERT_FALSE(face.natural_bounds);
    ASSERT_EQ(face.wires.size(), 1);
    ASSERT_TRUE(face.wires.front() == w); // the free wire becomes the face's wire
    ASSERT_TRUE(face.surface == srf);
    ASSERT_NEAR(uv_signed_area(m, w, 2048), pi, 1e-5); // (u,v) = (x,y): the unit disk
    ASSERT_LE(m.edge(e).tol, MakeFaceOptions<T>{}.pcurve_tol);
    check_pcurves(m, f);
}

TEST(tests_brep_face_wire, outer_reoriented_counter_clockwise)
{
    Model<T> m;
    auto w = square_wire(m, 1., true); // clockwise in (x, y)
    const auto first_edge = m.wire(w).coedges.front().edge;
    auto f = unwrap(make_face(m, xy_plane(), w));
    ASSERT_NEAR(uv_signed_area(m, w), 4., 1e-12);
    ASSERT_TRUE(m.wire(w).coedges.back().edge == first_edge); // order reversed…
    check_pcurves(m, f);                                       // …and senses flipped consistently
    // straight edges on a plane: exact up to round-off
    for (const auto &ce : m.wire(w).coedges)
        ASSERT_LE(m.edge(ce.edge).tol, 1e-6);
}

TEST(tests_brep_face_wire, face_with_holes)
{
    Model<T> m;
    auto outer = square_wire(m, 1.5);
    auto h1 = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.4, {-0.7, 0., 0.})))}));        // counter-clockwise: flipped
    auto h2 = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.3, {0.7, 0.5, 0.})))}));
    auto f = unwrap(make_face(m, xy_plane(), outer, std::vector{h1, h2}));
    const auto &face = m.face(f);
    ASSERT_EQ(face.wires.size(), 3);
    ASSERT_TRUE(face.wires[0] == outer);
    ASSERT_NEAR(uv_signed_area(m, outer), 9., 1e-12);
    ASSERT_NEAR(uv_signed_area(m, h1, 2048), -pi * 0.16, 1e-5);
    ASSERT_NEAR(uv_signed_area(m, h2, 2048), -pi * 0.09, 1e-5);
    check_pcurves(m, f);
}

TEST(tests_brep_face_wire, hole_outside_outer_leaves_model_unchanged)
{
    Model<T> m;
    auto outer = square_wire(m, 0.5);
    auto hole = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.3, {1.2, 0., 0.})))}));
    const auto tol_before = m.edge(m.wire(hole).coedges.front().edge).tol;
    auto r = make_face(m, xy_plane(), outer, std::vector{hole});
    ASSERT_FALSE(r.has_value());
    ASSERT_TRUE(r.error().code == BuildErrc::HoleOutsideOuter);
    ASSERT_EQ(m.count<ShapeType::Face>(), 0);
    ASSERT_EQ(alive_wires_with_pcurves(m), 0);
    ASSERT_DOUBLE_EQ(m.edge(m.wire(hole).coedges.front().edge).tol, tol_before);
    // the same wires are still usable
    ASSERT_TRUE(make_face(m, xy_plane(), outer).has_value());
}

TEST(tests_brep_face_wire, projection_on_curved_surface)
{
    Model<T> m;
    auto srf = bump();
    auto other = std::make_shared<BSSurface<T, 3>>(*srf); // same geometry, another object: no exact extraction
    std::vector<EdgeId> e;
    const std::array<point<T, 2>, 4> c{{{-1., -1.}, {1., -1.}, {1., 1.}, {-1., 1.}}};
    for (size_t i = 0; i < 4; ++i)
        e.push_back(unwrap(make_edge(m, std::make_shared<CurveOnSurface<T, 3>>(line2d(c[i], c[(i + 1) % 4]), other))));
    auto w = unwrap(make_wire(m, e));
    auto f = unwrap(make_face(m, srf, w));
    ASSERT_NEAR(uv_signed_area(m, w), 4., 1e-6);
    for (const auto &ce : m.wire(w).coedges)
    {
        ASSERT_GT(m.edge(ce.edge).tol, 0.);
        ASSERT_LE(m.edge(ce.edge).tol, MakeFaceOptions<T>{}.pcurve_tol);
        // the projected pcurve is the original 2D segment
        const auto &ed = m.edge(ce.edge);
        auto cos = std::dynamic_pointer_cast<CurveOnSurface<T, 3>>(ed.curve);
        for (T s : {0., 0.3, 0.7, 1.})
        {
            const T t = ed.u1 + s * (ed.u2 - ed.u1);
            ASSERT_LE(distance(ce.pcurve->value(t), cos->basisCurve().value(t)), 1e-4);
        }
    }
    check_pcurves(m, f);
}

TEST(tests_brep_face_wire, exact_pcurves_from_2d_loops)
{
    Model<T> m;
    auto srf = bump();
    std::vector<std::shared_ptr<Curve<T, 2>>> outer, hole;
    const std::array<point<T, 2>, 4> c{{{-1.5, -1.}, {1.5, -1.}, {1.5, 1.}, {-1.5, 1.}}};
    for (size_t i = 0; i < 4; ++i)
        outer.push_back(line2d(c[i], c[(i + 1) % 4]));
    const std::array<point<T, 2>, 3> h{{{0., 0.5}, {0.5, -0.5}, {-0.5, -0.5}}}; // clockwise triangle
    for (size_t i = 0; i < 3; ++i)
        hole.push_back(line2d(h[i], h[(i + 1) % 3]));

    auto f = unwrap(make_face(m, srf, outer, {hole}));
    const auto &face = m.face(f);
    ASSERT_EQ(face.wires.size(), 2);
    ASSERT_NEAR(uv_signed_area(m, face.wires[0]), 6., 1e-12);
    ASSERT_NEAR(uv_signed_area(m, face.wires[1]), -0.5, 1e-12);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 7); // shared corners merged
    for (auto wid : face.wires)
        for (const auto &ce : m.wire(wid).coedges)
        {
            // exact: the pcurve is the given 2D curve itself, no deviation
            ASSERT_TRUE(std::ranges::find(outer, ce.pcurve) != outer.end() || std::ranges::find(hole, ce.pcurve) != hole.end());
            ASSERT_DOUBLE_EQ(m.edge(ce.edge).tol, brep_default_tolerance<T>);
        }
    check_pcurves(m, f);

    // failure: everything created is erased
    const auto nv = m.count<ShapeType::Vertex>(), ne = m.count<ShapeType::Edge>(), nw = m.count<ShapeType::Wire>();
    std::vector<std::shared_ptr<Curve<T, 2>>> far;
    const std::array<point<T, 2>, 3> g{{{1.6, 1.6}, {1.9, 1.6}, {1.9, 1.9}}};
    for (size_t i = 0; i < 3; ++i)
        far.push_back(line2d(g[i], g[(i + 1) % 3]));
    auto r = make_face(m, srf, outer, {far});
    ASSERT_TRUE(r.error().code == BuildErrc::HoleOutsideOuter);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), nv);
    ASSERT_EQ(m.count<ShapeType::Edge>(), ne);
    ASSERT_EQ(m.count<ShapeType::Wire>(), nw);
    ASSERT_EQ(m.count<ShapeType::Face>(), 1);
}

TEST(tests_brep_face_wire, seam_patches_on_cylinder)
{
    auto srf = cylinder(1., 1.);
    // patch [0, pi]: its left side lies on the seam u = 0 == 2 pi; patch [pi, 2 pi]: right side on the seam
    for (auto [a, b, seam_u] : {std::array<T, 3>{0., pi, 0.}, std::array<T, 3>{pi, 2. * pi, 2. * pi}})
    {
        Model<T> m;
        auto arc = [&](T z) {
            return unwrap(make_edge(m, std::make_shared<AnalyticCurve>([z](T t) { return point<T, 3>{std::cos(t), std::sin(t), z}; },
                                                                       std::array<T, 2>{a, b})));
        };
        auto line = [&](T t) {
            return unwrap(make_edge(m, point<T, 3>{std::cos(t), std::sin(t), 0.}, point<T, 3>{std::cos(t), std::sin(t), 1.}));
        };
        auto w = unwrap(make_wire(m, std::vector{arc(0.), line(b), arc(1.), line(a)}));
        auto f = unwrap(make_face(m, srf, w));
        ASSERT_NEAR(uv_signed_area(m, w, 2048), pi, 1e-5);
        check_pcurves(m, f);
        // the edge on the seam got the pcurve on the side of the patch
        bool found = false;
        for (const auto &ce : m.wire(w).coedges)
        {
            const auto &ed = m.edge(ce.edge);
            const auto mid = ce.pcurve->value(0.5 * (ed.u1 + ed.u2));
            if (std::abs(mid[0] - seam_u) < 1e-6)
                found = true;
            ASSERT_GE(mid[0], a - 1e-9);
            ASSERT_LE(mid[0], b + 1e-9);
        }
        ASSERT_TRUE(found);
    }
}

TEST(tests_brep_face_wire, errors)
{
    Model<T> m;
    auto srf = xy_plane();

    // edge across the seam of a cylinder
    {
        auto cyl = cylinder(1., 1.);
        auto arc = [&](T z) {
            return unwrap(make_edge(m, std::make_shared<AnalyticCurve>([z](T t) { return point<T, 3>{std::cos(t), std::sin(t), z}; },
                                                                       std::array<T, 2>{-pi / 2, pi / 2})));
        };
        auto line = [&](T y) { return unwrap(make_edge(m, point<T, 3>{0., y, 0.}, point<T, 3>{0., y, 1.})); };
        auto w = unwrap(make_wire(m, std::vector{arc(0.), line(1.), arc(1.), line(-1.)}));
        auto r = make_face(m, cyl, w);
        ASSERT_TRUE(r.error().code == BuildErrc::CrossesSeam);
        ASSERT_EQ(r.error().shapes.size(), 1);
    }
    // edge off the surface
    {
        auto w = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(1., {0., 0., 1.})))}));
        auto r = make_face(m, srf, w);
        ASSERT_TRUE(r.error().code == BuildErrc::EdgeOffSurface);
    }
    // open wire
    {
        auto e = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}));
        ASSERT_TRUE(make_face(m, srf, unwrap(make_wire(m, std::vector{e}))).error().code == BuildErrc::WireNotClosed);
    }
    // wire already bounding a face, or given twice
    {
        auto w = square_wire(m, 0.5);
        ASSERT_TRUE(make_face(m, srf, w).has_value());
        ASSERT_TRUE(make_face(m, srf, w).error().code == BuildErrc::WireInUse);
        auto w2 = square_wire(m, 1.);
        auto h = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.2, {0., 0., 0.})))}));
        ASSERT_TRUE(make_face(m, srf, w2, std::vector{h, h}).error().code == BuildErrc::WireInUse);
        ASSERT_TRUE(make_face(m, srf, WireId{999}).error().code == BuildErrc::InvalidId);
    }
    ASSERT_TRUE(make_face(m, std::shared_ptr<Surface<T, 3>>{}, WireId{0}).error().code == BuildErrc::NullSurface);
    ASSERT_EQ(m.count<ShapeType::Face>(), 1);
}

TEST(tests_brep_face_wire, plane_pcurves_of_isos)
{
    // a planar rational B-spline surface (z = 0.5, irregular poles and weights, an interior knot
    // in each direction): its iso lines are rational curves of the plane; an edge following one
    // (a curve on the surface whose pcurve is a segment at constant u or v) gets an exact pcurve
    // on the plane, whatever the iso (bound or interior, along u or v) and its sense
    using T = double;
    std::mt19937 rng{7};
    std::uniform_real_distribution<T> jitter(-0.1, 0.1), weight(0.6, 1.8);
    std::vector<std::array<T, 4>> poles;
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
        {
            const T w = weight(rng), x = i / 3. + jitter(rng), y = j / 3. + jitter(rng);
            poles.push_back({w * x, w * y, w * 0.5, w});
        }
    const std::vector<T> k{0., 0., 0., 0.5, 1., 1., 1.};
    auto srf = std::make_shared<BSSurfaceRational<T, 3>>(poles, k, k, 2, 2);
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-9., 9., -9., 9.}};
    // an edge on [2, 5] following the segment p0 -> p1 of the parameters of the surface
    auto check = [&](point<T, 2> p0, point<T, 2> p1) {
        INFO("iso from (" << p0[0] << ", " << p0[1] << ") to (" << p1[0] << ", " << p1[1] << ")");
        auto seg = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{p0, p1}, std::vector<T>{2., 2., 5., 5.}, 1);
        auto cos = std::make_shared<CurveOnSurface<T, 3>>(std::static_pointer_cast<Curve<T, 2>>(seg), std::static_pointer_cast<Surface<T, 3>>(srf));
        Edge<T> e{cos, 2., 5., VertexId{}, VertexId{}};
        const auto fit = plane_pcurve(e, plane, 1e-9);
        REQUIRE(fit.has_value());
        REQUIRE(fit->exact);
        for (int i = 0; i <= 20; ++i)
        {
            const T t = 2. + 3. * i / 20.;
            const auto uv = fit->pcurve->value(t);
            REQUIRE(norm(plane.value(uv[0], uv[1]) - cos->value(t)) <= 1e-12);
        }
    };
    for (T v : {0., 1., 0.37, 0.5})
    {
        check({0., v}, {1., v}); // along u
        check({1., v}, {0., v}); // reversed
        check({0.2, v}, {0.8, v}); // a part of it
    }
    for (T u : {0., 1., 0.61})
    {
        check({u, 0.}, {u, 1.}); // along v
        check({u, 0.9}, {u, 0.1}); // reversed part
    }
    // within rounding of a bound: the row of poles, not a knot inserted at 1e-17
    check({0., 1e-17}, {1., 1e-17});
    check({1. - 1e-16, 0.}, {1. - 1e-16, 1.});
    // beyond the bounds: no exact pcurve, and no exception
    {
        auto seg = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0., 1. + 1e-6}, {1., 1. + 1e-6}}, std::vector<T>{0., 0., 1., 1.}, 1);
        auto cos = std::make_shared<CurveOnSurface<T, 3>>(std::static_pointer_cast<Curve<T, 2>>(seg), std::static_pointer_cast<Surface<T, 3>>(srf));
        Edge<T> e{cos, 0., 1., VertexId{}, VertexId{}};
        ASSERT_FALSE(plane_pcurve(e, plane, 1e-9).has_value());
    }
    // a rational pcurve moved by a period: shifted_pcurve, through gbs::translate
    {
        auto seg = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0., 0.37}, {1., 0.37}}, std::vector<T>{0., 0., 1., 1.}, 1);
        auto cos = std::make_shared<CurveOnSurface<T, 3>>(std::static_pointer_cast<Curve<T, 2>>(seg), std::static_pointer_cast<Surface<T, 3>>(srf));
        const auto fit = plane_pcurve(Edge<T>{cos, 0., 1., VertexId{}, VertexId{}}, plane, 1e-9);
        REQUIRE(fit.has_value());
        REQUIRE(dynamic_cast<const BSCurveRational<T, 2> *>(fit->pcurve.get()) != nullptr);
        const point<T, 2> shift{2. * std::numbers::pi, -1.};
        const auto moved = shifted_pcurve(fit->pcurve, shift);
        REQUIRE(moved != nullptr);
        for (int i = 0; i <= 10; ++i)
        {
            const T t = i / 10.;
            const auto a = fit->pcurve->value(t), b = moved->value(t);
            ASSERT_NEAR(b[0], a[0] + shift[0], 1e-12);
            ASSERT_NEAR(b[1], a[1] + shift[1], 1e-12);
        }
        ASSERT_EQ(shifted_pcurve(fit->pcurve, point<T, 2>{0., 0.}), fit->pcurve);
    }
}
