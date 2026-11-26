#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/bselementary.h>
#include <gbs/curvecomposite.h>
#include <gbs/curvetrimmed.h>
#include <gbs/curvecircle.h>
#include <gbs/bssbuild.h>
#include "tests_brep_geom.h"

#include <cmath>
#include <functional>
#include <map>
#include <string>
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

TEST(tests_brep_face_wire, hole_near_a_curved_outer_boundary)
{
    // a hole reaching 0.99 of the radius of a disk, at the middle of a side of the 15-gon inscribed in
    // the circle (its apothem is cos(pi/15) = 0.978): inside
    {
        Model<T> m;
        const T a = pi / 15.;
        auto outer = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(1., {0., 0., 0.})))}));
        auto hole = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.02, {0.97 * std::cos(a), 0.97 * std::sin(a), 0.})))}));
        ASSERT_TRUE(make_face(m, xy_plane(), outer, std::vector{hole}).has_value());
        // the same hole moved out by 0.03 crosses the circle
        auto across = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.02, {std::cos(a), std::sin(a), 0.})))}));
        auto outer2 = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(1., {0., 0., 0.})))}));
        ASSERT_TRUE(make_face(m, xy_plane(), outer2, std::vector{across}).error().code == BuildErrc::HoleOutsideOuter);
    }
    // a square bitten by a half disk of radius 1.5 centered on its top side: a hole between the bite and
    // the chord of the inscribed polygon (sag 8e-3) is outside the face
    {
        Model<T> m;
        const T h = 1.5;
        auto bite = std::make_shared<AnalyticCurve>([h](T t) { return point<T, 3>{h * std::cos(t), h - h * std::sin(t), 0.}; }, std::array<T, 2>{0., pi});
        auto outer = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, point<T, 3>{-h, -h, 0.}, point<T, 3>{h, -h, 0.})),
                                                     unwrap(make_edge(m, point<T, 3>{h, -h, 0.}, point<T, 3>{h, h, 0.})), unwrap(make_edge(m, bite)),
                                                     unwrap(make_edge(m, point<T, 3>{-h, h, 0.}, point<T, 3>{-h, -h, 0.}))}));
        const T a = pi / 30., d = h - 0.004; // middle of the first chord, 4e-3 into the bite
        auto hole = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.003, {d * std::cos(a), h - d * std::sin(a), 0.})))}));
        auto r = make_face(m, xy_plane(), outer, std::vector{hole});
        ASSERT_FALSE(r.has_value());
        ASSERT_TRUE(r.error().code == BuildErrc::HoleOutsideOuter);
        // the same hole 4e-3 below the bite is in the face
        const T e = h + 0.004;
        auto in = unwrap(make_wire(m, std::vector{unwrap(make_edge(m, circle3d(0.003, {e * std::cos(a), h - e * std::sin(a), 0.})))}));
        ASSERT_TRUE(make_face(m, xy_plane(), outer, std::vector{in}).has_value());
    }
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

TEST(tests_brep_face_wire, projected_pcurve_tolerance)
{
    // the tolerance of an edge with a projected pcurve bounds its deviation over the whole edge:
    // the error of the interpolation peaks off the middle of the end spans, where only the middles
    // were measured (a circle of radius 1 on a B-spline plane: tolerance 2.402e-6 for a deviation
    // of 2.625e-6; radius 0.4: the face tolerance, 1e-6, for 1.05e-6)
    for (T r : {1., 0.4})
    {
        INFO("radius " << r);
        Model<T> m;
        const auto srf = xy_plane();
        const auto crv = circle3d(r, {0., 0., 0.});
        const auto e = unwrap(make_edge(m, crv));
        const auto f = unwrap(make_face(m, srf, unwrap(make_wire(m, std::vector{e}))));
        const auto &ce = m.wire(m.face(f).wires.front()).coedges.front();
        const auto &ed = m.edge(e);
        T deviation{};
        constexpr int n = 20'000;
        for (int i = 0; i <= n; ++i)
        {
            const T t = ed.u1 + (ed.u2 - ed.u1) * T(i) / T(n);
            const auto uv = ce.pcurve->value(t);
            deviation = std::max(deviation, distance(srf->value(uv[0], uv[1]), crv->value(t)));
        }
        ASSERT_LE(deviation, ed.tol * (1. + 1e-3));
        ASSERT_GT(deviation, ed.tol * 0.99); // (not a bound far above it)
    }
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

TEST(tests_brep_face_wire, exact_seam_edge_on_either_side)
{
    // the seam edge of a natural face is exact on u = 0; a patch [pi, 2 pi] takes it on u = 2 pi
    auto srf = cylinder(1., 1.);
    for (auto [a, b] : {std::array<T, 2>{0., pi}, std::array<T, 2>{pi, 2. * pi}})
    {
        Model<T> m;
        const auto natural = unwrap(make_face(m, srf));
        const auto nc = m.wire(m.face(natural).wires[0]).coedges;
        const auto seam = std::ranges::find_if(nc, [&](const CoEdge<T> &ce) { return std::ranges::count(nc, ce.edge, &CoEdge<T>::edge) == 2; })->edge;
        auto seg = [&](point<T, 2> p, point<T, 2> q) { return unwrap(make_edge(m, std::make_shared<CurveOnSurface<T, 3>>(line2d(p, q), srf))); };
        auto w = unwrap(make_wire(m, std::vector{seg({a, 0.}, {b, 0.}), seg({pi, 0.}, {pi, 1.}), seg({b, 1.}, {a, 1.}), seam}));
        auto f = unwrap(make_face(m, srf, w));
        ASSERT_NEAR(std::abs(uv_signed_area(m, w)), pi, 1e-12);
        check_pcurves(m, f);
        const auto &ce = *std::ranges::find(m.wire(w).coedges, seam, &CoEdge<T>::edge);
        ASSERT_NEAR(ce.pcurve->value(0.5)[0], a == 0. ? 0. : 2. * pi, 1e-12);
        ASSERT_DOUBLE_EQ(m.edge(seam).tol, tol); // moved by a period, still exact
        check_pcurves(m, natural);
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

TEST(tests_brep_face_wire, options_out_of_range)
{
    // a pcurve degree of 0 is refused before anything is built (the interpolation of degree 0 crashed)
    Model<T> m;
    auto w = square_wire(m, 0.5);
    auto r = make_face(m, xy_plane(), w, std::vector<WireId>{}, MakeFaceOptions<T>{.pcurve_degree = 0});
    ASSERT_FALSE(r.has_value());
    ASSERT_TRUE(r.error().code == BuildErrc::InvalidOptions);
    ASSERT_TRUE(make_face(m, xy_plane(), w, std::vector<WireId>{}, MakeFaceOptions<T>{.pcurve_tol = -1.}).error().code == BuildErrc::InvalidTolerance);
    ASSERT_EQ(m.count<ShapeType::Face>(), 0);
    ASSERT_TRUE(make_face(m, xy_plane(), w, std::vector<WireId>{}, MakeFaceOptions<T>{.pcurve_degree = 1}).has_value());
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

// Issue #260: the same rule for every curve type when an edge leaves the box of a plane. A
// circle of radius 0.1 centred at (1, 0.3), half beyond the box [0, 1]^2 of the plane, and a square
// around it: as an analytic circle, an ellipse, a NURBS circle, a trimmed circle, a composite of
// two arcs; as line curves or B-spline segments. Inside the box, all are built; beyond it, they
// must all be built or all refused (master built the B-splines only, through the exact pcurve of a
// B-spline on a plane, which ignored the box)
TEST(tests_brep_face_wire, same_rule_for_every_curve_beyond_the_box_of_a_plane)
{
    using C = std::shared_ptr<Curve<T, 3>>;
    auto built = [](std::shared_ptr<Surface<T, 3>> srf, const std::vector<C> &cs) {
        Model<T> m;
        std::vector<EdgeId> es;
        for (const auto &c : cs)
            es.push_back(unwrap(make_edge(m, c)));
        const auto f = make_face(m, srf, unwrap(make_wire(m, es)));
        return f.has_value() && check(m, ShapeId{*f}).ok();
    };
    const ax2<T, 3> xy{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    const auto plane = std::make_shared<Plane<T>>(xy, std::array<T, 4>{0., 1., 0., 1.});
    for (const point<T, 3> o : {point<T, 3>{0.5, 0.5, 0.}, point<T, 3>{1., 0.3, 0.}})
    {
        const bool beyond = o[0] > 0.75;
        const ax2<T, 3> ax{{o, {0., 0., 1.}, {1., 0., 0.}}};
        std::vector<std::pair<std::string, std::vector<C>>> cases;
        cases.push_back({"Circle3d", {std::make_shared<Circle3d<T>>(0.1, ax)}});
        cases.push_back({"EllipseCurve", {std::make_shared<EllipseCurve<T>>(0.1, 0.1, Frame3<T>{ax})}});
        cases.push_back({"BSCurveRational", {std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(0.1, 0., 2. * pi, ax))}});
        cases.push_back({"CurveTrimmed", {std::make_shared<CurveTrimmed<T, 3>>(std::make_shared<Circle3d<T>>(0.1, ax), 0., 2. * pi)}});
        cases.push_back({"CurveComposite", {std::make_shared<CurveComposite<T, 3>>(std::vector<C>{std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(0.1, 0., pi, ax)),
                                                                                                  std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(0.1, pi, 2. * pi, ax))})}});
        const std::vector<point<T, 3>> q{o + point<T, 3>{-0.15, -0.15, 0.}, o + point<T, 3>{0.15, -0.15, 0.}, o + point<T, 3>{0.15, 0.15, 0.}, o + point<T, 3>{-0.15, 0.15, 0.}};
        std::vector<C> lines, segs;
        for (int i = 0; i < 4; ++i)
        {
            lines.push_back(std::make_shared<LineCurve<T>>(q[i], q[(i + 1) % 4] - q[i], 0., 1.));
            segs.push_back(std::make_shared<BSCurve<T, 3>>(build_segment(q[i], q[(i + 1) % 4])));
        }
        cases.push_back({"LineCurve square", lines});
        cases.push_back({"BSCurve square", segs});
        std::map<bool, std::vector<std::string>> by;
        for (const auto &[nm, cs] : cases)
            by[built(plane, cs)].push_back(nm);
        std::string seen;
        for (const auto &[ok, names] : by)
        {
            seen += ok ? " built:" : " refused:";
            for (const auto &n : names)
                seen += " " + n;
        }
        INFO((beyond ? "beyond the box:" : "inside the box:") << seen);
        // (a plane is unbounded: all built)
        ASSERT_EQ(by.size(), 1u);
        ASSERT_TRUE(by.begin()->first);
    }
}

namespace
{
    using Curves = std::vector<std::shared_ptr<Curve<T, 3>>>;

    /// The face on srf bounded by the closed chain of curves cs (each on its whole range).
    BuildResult<FaceId> face_of(Model<T> &m, std::shared_ptr<Surface<T, 3>> srf, const Curves &cs)
    {
        std::vector<EdgeId> es;
        for (const auto &c : cs)
            es.push_back(unwrap(make_edge(m, c)));
        return make_face(m, std::move(srf), unwrap(make_wire(m, es)));
    }

    /// A patch of a surface of revolution about z, from the angle 0 to 1 between the heights z1 and z2,
    /// its meridian at radius rho(z): the arcs analytic (ellipses) or rational B-splines, the meridians
    /// line curves or B-spline segments.
    Curves patch(T z1, T z2, const std::function<T(T)> &rho, bool analytic)
    {
        const auto at = [&](T a, T z) { return point<T, 3>{rho(z) * std::cos(a), rho(z) * std::sin(a), z}; };
        auto arc = [&](T z) -> std::shared_ptr<Curve<T, 3>> {
            const ax2<T, 3> ax{{{0., 0., z}, {0., 0., 1.}, {1., 0., 0.}}};
            if (analytic)
                return std::make_shared<EllipseCurve<T>>(rho(z), rho(z), Frame3<T>{ax}, 0., 1.);
            return std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(rho(z), 0., 1., ax));
        };
        auto meridian = [&](T a, T za, T zb) -> std::shared_ptr<Curve<T, 3>> {
            if (analytic)
                return std::make_shared<LineCurve<T>>(at(a, za), at(a, zb) - at(a, za), 0., 1.);
            return std::make_shared<BSCurve<T, 3>>(build_segment(at(a, za), at(a, zb)));
        };
        return {arc(z1), meridian(1., z1, z2), arc(z2), meridian(0., z1, z2)};
    }
}

// Issue #260: the open directions of elementary surfaces (a plane in u and v, a cylinder and a cone
// in v) are unbounded, the face built on a copy of the primitive with a box enclosing its wire (the
// surface given shared when its box already does). The closed directions (the angle), the apex of a
// cone and the other surfaces (B-spline, revolution) are bounded: an edge beyond is refused,
// "beyond the bounds of the surface", whatever its curve
TEST(tests_brep_face_wire, elementary_surfaces_unbounded_in_their_open_directions)
{
    const ax2<T, 3> xy{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    auto refused_beyond = [](const BuildResult<FaceId> &f) { return !f.has_value() && f.error().message.find("beyond the bounds of the surface") != std::string::npos; };
    // a plane: the surface given shared inside its box, a copy enclosing the wire beyond it
    {
        const auto plane = std::make_shared<Plane<T>>(xy, std::array<T, 4>{0., 1., 0., 1.});
        Model<T> m;
        const auto in = unwrap(face_of(m, plane, {std::make_shared<Circle3d<T>>(0.1, ax2<T, 3>{{{0.5, 0.5, 0.}, {0., 0., 1.}, {1., 0., 0.}}})}));
        ASSERT_EQ(m.face(in).surface.get(), plane.get());
        const auto out = unwrap(face_of(m, plane, {std::make_shared<Circle3d<T>>(0.1, ax2<T, 3>{{{1., 0.3, 0.}, {0., 0., 1.}, {1., 0., 0.}}})}));
        ASSERT_NE(m.face(out).surface.get(), plane.get());
        const auto b = m.face(out).surface->bounds();
        ASSERT_TRUE(b[0] <= 0.9 && b[1] >= 1.1 && b[2] <= 0.2 && b[3] >= 0.4);
        ASSERT_TRUE(check(m, ShapeId{out}).ok());
        // (the area in (u, v) of the inscribed polygon of the pcurves: that of the same circle inside)
        ASSERT_NEAR(uv_signed_area(m, m.face(out).wires.front()), uv_signed_area(m, m.face(in).wires.front()), 1e-12);
        // the plane given unchanged
        ASSERT_EQ(plane->bounds(), (std::array<T, 4>{0., 1., 0., 1.}));
    }
    // a cylinder (v in [0, 1]) and a cone (radius 1 at v = 0, a = 30 degrees, v in [0, 1]): a patch
    // up to v = 1.5, analytic curves or B-splines: built, valid
    const auto cyl = std::make_shared<CylindricalSurface<T>>(1., xy, std::array<T, 4>{0., 2. * pi, 0., 1.});
    const T ta = std::tan(pi / 6.);
    const auto cone = std::make_shared<ConicalSurface<T>>(1., pi / 6., xy, std::array<T, 4>{0., 2. * pi, 0., 1.});
    for (const bool analytic : {true, false})
    {
        Model<T> m;
        const auto fc = face_of(m, cyl, patch(0.2, 1.5, [](T) { return 1.; }, analytic));
        INFO("analytic " << analytic << ": " << (fc ? std::string("built") : fc.error().message));
        ASSERT_TRUE(fc.has_value());
        ASSERT_TRUE(check(m, ShapeId{*fc}).ok());
        const auto fk = face_of(m, cone, patch(0.2, 1.5, [&](T z) { return 1. + z * ta; }, analytic));
        ASSERT_TRUE(fk.has_value());
        ASSERT_TRUE(check(m, ShapeId{*fk}).ok());
        // down to v = -1 (radius 0.42), still on the cone's side of its apex (v = -1.73): built
        Model<T> m2;
        ASSERT_TRUE(face_of(m2, cone, patch(-1., 0.5, [&](T z) { return 1. + z * ta; }, analytic)).has_value());
    }
    // past the apex: the other nappe is not the face's; refused, beyond the bounds
    {
        Model<T> m;
        const T apex = -1. / ta;
        const auto f = face_of(m, cone, patch(apex - 1., apex - 0.5, [&](T z) { return std::abs(1. + z * ta); }, true));
        INFO((f ? std::string("built") : f.error().message));
        ASSERT_TRUE(refused_beyond(f) || (!f.has_value() && f.error().code == BuildErrc::EdgeOffSurface));
    }
    // the angle of a cylinder bounded to [0, pi]: a patch at the angles 4 to 5 is refused (a closed
    // direction is not grown)
    {
        const auto half = std::make_shared<CylindricalSurface<T>>(1., xy, std::array<T, 4>{0., pi, 0., 1.});
        Model<T> m;
        // (the patch of the angles 0 to 1 turned to 4 to 5)
        Curves turned;
        const auto rot = [](const point<T, 3> &p) { return point<T, 3>{std::cos(4.) * p[0] - std::sin(4.) * p[1], std::sin(4.) * p[0] + std::cos(4.) * p[1], p[2]}; };
        const ax2<T, 3> a1{{{0., 0., 0.2}, {0., 0., 1.}, {std::cos(4.), std::sin(4.), 0.}}}, a2{{{0., 0., 0.8}, {0., 0., 1.}, {std::cos(4.), std::sin(4.), 0.}}};
        turned.push_back(std::make_shared<EllipseCurve<T>>(1., 1., Frame3<T>{a1}, 0., 1.));
        turned.push_back(std::make_shared<LineCurve<T>>(rot(point<T, 3>{std::cos(1.), std::sin(1.), 0.2}), point<T, 3>{0., 0., 0.6}, 0., 1.));
        turned.push_back(std::make_shared<EllipseCurve<T>>(1., 1., Frame3<T>{a2}, 0., 1.));
        turned.push_back(std::make_shared<LineCurve<T>>(rot(point<T, 3>{1., 0., 0.2}), point<T, 3>{0., 0., 0.6}, 0., 1.));
        const auto f = face_of(m, half, turned);
        INFO((f ? std::string("built") : f.error().message));
        ASSERT_TRUE(refused_beyond(f));
    }
    // bounded surfaces: a B-spline plane [0, 1]^2 and a rational B-spline cylinder (the NURBS of a
    // cylinder, v in [0, 1]): beyond them, refused for every curve
    {
        const auto bsp = std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}}, std::vector<T>{0., 0., 1., 1.},
                                                           std::vector<T>{0., 0., 1., 1.}, 1, 1);
        const ax2<T, 3> ax{{{1., 0.3, 0.}, {0., 0., 1.}, {1., 0., 0.}}};
        for (const auto &c : Curves{std::make_shared<Circle3d<T>>(0.1, ax), std::make_shared<BSCurveRational<T, 3>>(build_circle_arc<T>(0.1, 0., 2. * pi, ax))})
        {
            Model<T> m;
            const auto f = face_of(m, bsp, {c});
            INFO((f ? std::string("built") : f.error().message));
            ASSERT_TRUE(refused_beyond(f));
        }
        const auto nurbs = std::make_shared<CylindricalSurface<T>>(1., xy, std::array<T, 4>{0., pi, 0., 1.})->to_nurbs();
        for (const bool analytic : {true, false})
        {
            Model<T> m;
            const auto f = face_of(m, nurbs, patch(0.2, 1.5, [](T) { return 1.; }, analytic));
            INFO((f ? std::string("built") : f.error().message));
            ASSERT_TRUE(refused_beyond(f));
        }
    }
}

// Two squares side by side on one small plane (box [0, 1]^2), beyond its box: each face on its own
// copy of the plane; sewn, they still merge into one face (`unify_faces` maps the parameters of two
// primitives of one type by their geometry, not by the pointer)
TEST(tests_brep_face_wire, faces_beyond_a_shared_plane_still_merge)
{
    const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., 1., 0., 1.});
    Model<T> m;
    auto square = [&](T x0) {
        std::vector<EdgeId> es;
        const std::vector<point<T, 3>> q{{x0, 0., 0.}, {x0 + 1., 0., 0.}, {x0 + 1., 1., 0.}, {x0, 1., 0.}};
        for (int i = 0; i < 4; ++i)
            es.push_back(unwrap(make_edge(m, q[i], q[(i + 1) % 4])));
        return unwrap(make_face(m, std::shared_ptr<Surface<T, 3>>{plane}, unwrap(make_wire(m, es))));
    };
    const auto a = square(2.), b = square(3.);
    ASSERT_NE(m.face(a).surface.get(), m.face(b).surface.get());
    const auto sewn = unwrap(sew(m, std::vector<FaceId>{a, b}));
    auto uses = m.shell(sewn.shells.front()).faces;
    const auto r = unify_faces(m, uses, T(1e-6));
    ASSERT_EQ(uses.size(), 1u);
    ASSERT_EQ(r.faces_merged, 1u);
    ASSERT_TRUE(check(m, ShapeId{uses.front().face}).ok());
}
