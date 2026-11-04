#include <doctest_gtest.hpp>
#include <gbs-io/knot_cleanup.h>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include <gbs-io/iges_brep.h>
#include <gbs-io/step/geometry_writer.h>
#include <gbs-brep/builders.h>
#include <gbs/bsstools.h>
#include <gbs/primitives.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

using namespace gbs;
using gbs::nurbs_export::clean_knots;
using gbs::nurbs_export::KnotCleanupStats;
using gbs::nurbs_export::knot_detail::distinct;

namespace
{
    /// the shortest span of a knot vector
    double shortest_span(const std::vector<double> &k)
    {
        const auto x = distinct(k);
        double s = x.back() - x.front();
        for (std::size_t i = 0; i + 1 < x.size(); ++i)
            s = std::min(s, x[i + 1] - x[i]);
        return s;
    }

    // a loft knot (0.85 + 0.9 + 0.95) / 3, one ulp above 0.9
    const double k_near = 0.9000000000000001;

    // the surface of issue #213: v knots 0, 0.5, k_near, 1, degree 3 x 3
    BSSurface<double, 3> issue_surface()
    {
        std::vector<std::array<double, 3>> poles;
        for (std::size_t j{}; j < 6; ++j)
            for (std::size_t i{}; i < 4; ++i)
                poles.push_back({i / 3., j / 5., 0.1 * std::sin(double(i + j))});
        return {poles, {0., 0., 0., 0., 1., 1., 1., 1.}, {0., 0., 0., 0., 0.5, k_near, 1., 1., 1., 1.}, 3, 3};
    }

    // the part of it after v = 0.9 as the split of gbs 0.7.0 gave it: 0.9 x 4 then k_near
    BSSurface<double, 3> issue_part_after()
    {
        auto s = issue_surface();
        const auto ku = s.knotsFlatsU();
        // the rows along v of the surface, cut at 0.9 by knot insertion, the poles before it dropped
        std::vector<std::array<double, 3>> poles;
        std::vector<std::vector<std::array<double, 3>>> cols(4);
        std::vector<double> k;
        for (std::size_t i{}; i < 4; ++i)
        {
            for (std::size_t j{}; j < 6; ++j)
                cols[i].push_back(s.poles()[i + 4 * j]);
            k = s.knotsFlatsV();
            insert_knots(0.9, 3, 3, k, cols[i]);
        }
        // k: 0 x4, 0.5, 0.9 x3, k_near, 1 x4; keep from the first 0.9, with one more 0.9 in front
        const auto first = std::size_t(std::lower_bound(k.begin(), k.end(), 0.9) - k.begin());
        std::vector<double> kv_after(k.begin() + std::ptrdiff_t(first), k.end());
        kv_after.insert(kv_after.begin(), 0.9);
        const std::size_t nv = kv_after.size() - 4;
        for (std::size_t j{}; j < nv; ++j)
            for (std::size_t i{}; i < 4; ++i)
                poles.push_back(cols[i][cols[i].size() - nv + j]);
        return {poles, ku, kv_after, 3, 3};
    }
}

// a span of one ulp at the start of a curve: removed, the end multiplicity back to degree + 1
TEST(tests_knot_cleanup, curve_span_of_an_ulp_at_the_start)
{
    const std::vector<std::array<double, 3>> poles0 = {
        {0., 0., 0.}, {1., 2., 0.}, {2., -1., 1.}, {3., 1., 0.}, {4., 0., 2.}, {5., 1., 1.}};
    const BSCurve<double, 3> ref(poles0, {0., 0., 0., 0., 0.5, k_near, 1., 1., 1., 1.}, 3);
    auto k = ref.knotsFlats();
    auto P = ref.poles();
    insert_knots(0.9, 3, 3, k, P);
    // the old cut: 0.9 x 4, then k_near
    const auto first = std::size_t(std::lower_bound(k.begin(), k.end(), 0.9) - k.begin());
    std::vector<double> kt(k.begin() + std::ptrdiff_t(first), k.end());
    kt.insert(kt.begin(), 0.9);
    std::vector<std::array<double, 3>> Pt(P.end() - std::ptrdiff_t(kt.size() - 4), P.end());
    BSCurve<double, 3> part(Pt, kt, 3);
    ASSERT_EQ(part.knotsFlats(), (std::vector<double>{0.9, 0.9, 0.9, 0.9, k_near, 1., 1., 1., 1.}));

    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(part, 1e-7, st));
    ASSERT_EQ(part.knotsFlats(), (std::vector<double>{0.9, 0.9, 0.9, 0.9, 1., 1., 1., 1.}));
    ASSERT_EQ(part.poles().size(), 4);
    ASSERT_EQ(part.bounds()[0], 0.9);
    ASSERT_EQ(part.bounds()[1], 1.);
    ASSERT_EQ(st.cleaned_curves, 1);
    ASSERT_EQ(st.removed_spans, 1);
    ASSERT_EQ(st.kept_spans, 0);
    ASSERT_LT(st.max_deviation, 1e-12);
    for (int i{}; i <= 20; ++i)
    {
        const double u = 0.9 + 0.1 * i / 20.;
        ASSERT_LT(distance(part(u), ref(u)), 1e-12); // the same parameters
    }
    // cleaned again: nothing left to do
    ASSERT_FALSE(clean_knots(part, 1e-7, st));
}

// two interior knots 1e-13 apart: one knot, the curve within the tolerance at the same parameters
TEST(tests_knot_cleanup, curve_close_interior_knots)
{
    const std::vector<std::array<double, 2>> poles = {{0., 0.}, {1., 2.}, {2., -1.}, {3., 1.}, {4., 0.}, {5., 1.}, {6., 3.}};
    const BSCurve<double, 2> ref(poles, {0., 0., 0., 0., 0.3, 0.5, 0.5 + 1e-13, 1., 1., 1., 1.}, 3);
    auto c = ref;
    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(c, 1e-6, st));
    ASSERT_EQ(distinct(c.knotsFlats()).size(), 4); // 0, 0.3, 0.5, 1
    ASSERT_LE(st.max_deviation, 1e-6);
    for (int i{}; i <= 100; ++i)
    {
        const double u = i / 100.;
        ASSERT_LT(distance(c(u), ref(u)), 1e-6);
    }
}

// spans longer than the tolerance, or short in parameter but long in space, are kept
TEST(tests_knot_cleanup, long_spans_kept)
{
    const std::vector<std::array<double, 2>> poles = {{0., 0.}, {1., 2.}, {2., -1.}, {3., 1.}, {4., 0.}};
    BSCurve<double, 2> c(poles, {0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, 3);
    KnotCleanupStats<double> st;
    ASSERT_FALSE(clean_knots(c, 1e-3, st));
    // a span of 1e-13 over which the curve runs from (0, 0) to near (10, 0): a real piece of the curve
    BSCurve<double, 2> jump({{0., 0.}, {10., 0.}, {10., 1.}, {11., 2.}, {12., 0.}}, {0., 0., 0., 0., 1e-13, 1., 1., 1., 1.}, 3);
    ASSERT_FALSE(clean_knots(jump, 1e-3, st));
    ASSERT_EQ(st.removed_spans, 0);
    // a tolerance of 0 does nothing
    ASSERT_FALSE(clean_knots(c, 0., st));
}

// the surface of issue #213 after the cut: the ulp span in v removed, the surface the same
TEST(tests_knot_cleanup, surface_of_issue_213)
{
    const auto ref = issue_surface();
    auto s = issue_part_after();
    ASSERT_EQ(s.knotsFlatsV(), (std::vector<double>{0.9, 0.9, 0.9, 0.9, k_near, 1., 1., 1., 1.}));
    ASSERT_EQ(s.nPolesV(), 5);

    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(s, 1e-7, st));
    ASSERT_EQ(s.knotsFlatsV(), (std::vector<double>{0.9, 0.9, 0.9, 0.9, 1., 1., 1., 1.}));
    ASSERT_EQ(s.knotsFlatsU(), ref.knotsFlatsU());
    ASSERT_EQ(s.nPolesV(), 4);
    ASSERT_EQ(s.nPolesU(), 4);
    ASSERT_EQ(st.cleaned_surfaces, 1);
    ASSERT_LT(st.max_deviation, 1e-12);
    for (int i{}; i <= 10; ++i)
        for (int j{}; j <= 10; ++j)
        {
            const double u = i / 10., v = 0.9 + 0.1 * j / 10.;
            ASSERT_LT(distance(s(u, v), ref(u, v)), 1e-12);
        }
}

// a rational surface, the close knots along u, inside
TEST(tests_knot_cleanup, rational_surface)
{
    std::vector<std::array<double, 4>> poles;
    for (std::size_t j{}; j < 4; ++j)
        for (std::size_t i{}; i < 6; ++i)
        {
            const double w = 1. + 0.3 * std::cos(double(i * j));
            poles.push_back({w * i, w * j, w * 0.2 * std::sin(double(i + 2 * j)), w});
        }
    const BSSurfaceRational<double, 3> ref(poles, {0., 0., 0., 0.4, 0.4 + 1e-14, 0.7, 1., 1., 1.}, {0., 0., 0., 0., 1., 1., 1., 1.}, 2, 3);
    auto s = ref;
    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(s, 1e-7, st));
    ASSERT_EQ(distinct(s.knotsFlatsU()).size(), 4); // 0, 0.4, 0.7, 1
    ASSERT_EQ(s.knotsFlatsV(), ref.knotsFlatsV());
    for (int i{}; i <= 20; ++i)
        for (int j{}; j <= 5; ++j)
        {
            const double d = distance(s(i / 20., j / 5.), ref(i / 20., j / 5.));
            ASSERT_LT(d, 1e-7);
            ASSERT_LE(d, st.max_deviation * (1. + 1e-12) + 1e-15); // a bound, rational too
        }
}

// a cut within knot_eps of a knot is made at the knot: the pieces of the surface of the issue, and a
// curve trimmed an ulp before a knot as an edge parameter gives it, have no span of an ulp
TEST(tests_knot_cleanup, cut_made_at_the_knot)
{
    const auto [before, after] = split_v(issue_surface(), 0.9);
    ASSERT_EQ(after.knotsFlatsV(), (std::vector<double>{k_near, k_near, k_near, k_near, 1., 1., 1., 1.}));
    ASSERT_EQ(before.knotsFlatsV().back(), k_near);
    ASSERT_GT(shortest_span(before.knotsFlatsV()), 0.1);

    const BSCurve<double, 2> pc({{0., 0.}, {0.2, 0.3}, {0.4, 0.1}, {0.7, 0.5}, {1., 1.}}, {0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, 3);
    nurbs_export::ApproxStats<double> as;
    std::visit([](const auto &b) { ASSERT_EQ(b.knotsFlats(), (std::vector<double>{0.5, 0.5, 0.5, 0.5, 1., 1., 1., 1.})); },
               nurbs_export::to_nurbs<double, 2>(pc, std::nextafter(0.5, 0.), 1., 1e-5, as));

    // at a knot of multiplicity degree + 1 (a corner), the piece after it starts there, degree + 1 times
    const std::vector<std::array<double, 2>> poles{{0., 0.}, {1., 1.}, {2., 0.}, {3., 2.}, {4., 0.}, {5., 1.}};
    BSCurve<double, 2> corner(poles, {0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1.}, 2);
    auto piece = corner;
    piece.trim(std::nextafter(0.5, 0.), 1.);
    ASSERT_EQ(piece.knotsFlats(), (std::vector<double>{0.5, 0.5, 0.5, 1., 1., 1.}));
    for (double u : {0.5, 0.7, 1.})
        ASSERT_LT(distance(piece(u), corner(u)), 1e-15);
}

// a pcurve with a span of an ulp at its start, on a cylinder of radius 5: cleaned, its image on the
// surface within the tolerance
TEST(tests_knot_cleanup, pcurve_on_its_surface)
{
    const CylindricalSurface<double> cyl(5., Frame3<double>{}, {0., 2. * std::numbers::pi, 0., 10.});
    const double a = std::nextafter(0.5, 0.);
    const BSCurve<double, 2> ref({{1., 1.}, {1.2, 2.}, {1.5, 4.}, {2., 5.}, {2.4, 7.}}, {0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, 3);
    auto k = ref.knotsFlats();
    auto P = ref.poles();
    insert_knots(a, 3, 3, k, P);
    // the old cut at a: a x 4, then 0.5
    const auto first = std::size_t(std::lower_bound(k.begin(), k.end(), a) - k.begin());
    std::vector<double> kt(k.begin() + std::ptrdiff_t(first), k.end());
    kt.insert(kt.begin(), a);
    BSCurve<double, 2> pc(std::vector<std::array<double, 2>>(P.end() - std::ptrdiff_t(kt.size() - 4), P.end()), kt, 3);
    ASSERT_LT(shortest_span(pc.knotsFlats()), 1e-15);

    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(pc, cyl, 1e-6, st));
    ASSERT_EQ(st.cleaned_pcurves, 1);
    ASSERT_EQ(st.cleaned_curves, 0);
    ASSERT_GT(shortest_span(pc.knotsFlats()), 0.4);
    for (int i{}; i <= 20; ++i)
    {
        const double t = a + (1. - a) * i / 20.;
        const auto uv = pc(t), uv0 = ref(t);
        ASSERT_LT(distance(cyl(uv[0], uv[1]), cyl(uv0[0], uv0[1])), 1e-6);
    }

    // through the STEP writer, given its surface
    step::P21Writer w;
    step::GeometryWriter g(w, 1e-5, 1e-6);
    auto again = BSCurve<double, 2>(std::vector<std::array<double, 2>>(P.end() - std::ptrdiff_t(kt.size() - 4), P.end()), kt, 3);
    g.curve2d(again, a, 1., &cyl);
    ASSERT_EQ(g.knot_stats().cleaned_pcurves, 1);
    g.curve2d(again, a, 1.); // without its surface: as it is
    ASSERT_EQ(g.knot_stats().cleaned_pcurves, 1);
}

// a knot vector that is not clamped is left as it is (its ends are not at degree + 1)
TEST(tests_knot_cleanup, unclamped_left_as_it_is)
{
    BSCurve<double, 3> c({{0., 0., 0.}, {1e-9, 0., 0.}, {2e-9, 1e-9, 0.}, {3e-9, 0., 0.}}, {0., 1., 2., 3., 4., 5., 6., 7.}, 3);
    const auto k0 = c.knotsFlats();
    KnotCleanupStats<double> st;
    ASSERT_FALSE(clean_knots(c, 1e-6, st));
    ASSERT_EQ(c.knotsFlats(), k0);
}

// many short spans (every other one, 1e-12 long): all removed, within the tolerance
TEST(tests_knot_cleanup, many_short_spans)
{
    std::vector<double> k{0., 0., 0., 0.};
    const std::size_t n = 200;
    for (std::size_t i = 1; i < n; ++i)
        k.insert(k.end(), {double(i) / n, double(i) / n + 1e-12});
    k.insert(k.end(), {1., 1., 1., 1.});
    std::vector<std::array<double, 3>> poles;
    for (std::size_t i = 0; i + 4 < k.size(); ++i)
        poles.push_back({double(i), std::sin(0.1 * double(i)), 0.});
    const BSCurve<double, 3> ref(poles, k, 3);
    auto c = ref;
    KnotCleanupStats<double> st;
    ASSERT_TRUE(clean_knots(c, 1e-6, st));
    ASSERT_EQ(st.removed_spans, n - 1);
    ASSERT_EQ(distinct(c.knotsFlats()).size(), n + 1);
    ASSERT_LE(st.max_deviation, 1e-6);
    for (int i{}; i <= 1000; ++i)
        ASSERT_LE(distance(c(i / 1000.), ref(i / 1000.)), st.max_deviation + 1e-15);
}

// STEP and IGES: the surface of a face cleaned by default, as it is with knot_tol = 0
TEST(tests_knot_cleanup, step_and_iges_export)
{
    brep::Model<double> m;
    const auto srf = std::make_shared<BSSurface<double, 3>>(issue_part_after());
    const auto face = brep::make_face(m, srf);
    REQUIRE(face.has_value());
    const brep::ShapeId shape{*face};

    auto read_back_knots_v = [](const std::string &text) {
        brep::Model<double> m2;
        auto r = step::read_step(m2, text);
        REQUIRE(r.has_value());
        REQUIRE(r->report.failed.empty());
        const auto faces = brep::explore<brep::FaceId>(m2, r->root);
        REQUIRE(faces.size() == 1);
        const auto bs = std::dynamic_pointer_cast<BSSurface<double, 3>>(m2.face(faces.front()).surface);
        REQUIRE(bs);
        return bs->knotsFlatsV();
    };

    // by default: cleaned within the uncertainty written
    const auto by_default = step::write_step(m, shape);
    REQUIRE(by_default.has_value());
    ASSERT_EQ(by_default->second.knot_cleanup.cleaned_surfaces, 1);
    ASSERT_GE(by_default->second.knot_cleanup.removed_spans, 1);
    ASSERT_LT(by_default->second.knot_cleanup.max_deviation, 1e-12);
    ASSERT_EQ(read_back_knots_v(by_default->first), (std::vector<double>{0.9, 0.9, 0.9, 0.9, 1., 1., 1., 1.}));

    // knot_tol = 0: written as it is, the knot an ulp off 0.9 a knot of its own (it was merged with 0.9 into an
    // invalid multiplicity of 5), read back as it was
    const auto plain = step::write_step(m, shape, {.knot_tol = 0.});
    REQUIRE(plain.has_value());
    ASSERT_EQ(plain->second.knot_cleanup.removed_spans, 0);
    ASSERT_NE(plain->first.find("(4,1,4),(0.9,0.9000000000000001,1.)"), std::string::npos);
    ASSERT_EQ(read_back_knots_v(plain->first), srf->knotsFlatsV());

    // a tolerance given
    const auto given = step::write_step(m, shape, {.knot_tol = 1e-9});
    REQUIRE(given.has_value());
    // the surface, and the two edges along v (its iso lines, with the same span)
    ASSERT_EQ(given->second.knot_cleanup.cleaned_surfaces, 1);
    ASSERT_EQ(given->second.knot_cleanup.cleaned_curves, 2);
    ASSERT_EQ(given->second.knot_cleanup.removed_spans, 3);

    IgesWriter<double> w;
    const auto rep = add_brep(w, m, shape);
    ASSERT_EQ(rep.knot_cleanup.cleaned_surfaces, 1);
    ASSERT_GE(rep.knot_cleanup.removed_spans, 1);
    IgesWriter<double> w0;
    ASSERT_EQ(add_brep(w0, m, shape, "", IgesExportOptions<double>{.knot_tol = 0.}).knot_cleanup.removed_spans, 0);
}

// close v knots 2e-8 apart inside a face: cleaned at knot_tol, the uncertainty computed covering an edge
// and its face both moved; with an uncertainty given, knot_tol at most half of it
TEST(tests_knot_cleanup, step_uncertainty_and_knot_tol)
{
    std::vector<std::array<double, 3>> poles;
    for (std::size_t j{}; j < 6; ++j)
        for (std::size_t i{}; i < 4; ++i)
            poles.push_back({i / 3., j / 5., 0.1 * std::sin(double(i + j))});
    brep::Model<double> m;
    const auto srf = std::make_shared<BSSurface<double, 3>>(
        poles, std::vector<double>{0., 0., 0., 0., 1., 1., 1., 1.}, std::vector<double>{0., 0., 0., 0., 0.5, 0.5 + 2e-8, 1., 1., 1., 1.}, 3, 3);
    const auto face = brep::make_face(m, srf);
    REQUIRE(face.has_value());
    const brep::ShapeId shape{*face};

    const auto computed = step::write_step(m, shape, {.knot_tol = 1e-5});
    REQUIRE(computed.has_value());
    const auto &rc = computed->second;
    ASSERT_GE(rc.knot_cleanup.removed_spans, 1);
    ASSERT_GT(rc.knot_cleanup.max_deviation, 0.);
    ASSERT_GE(rc.uncertainty, 2 * rc.knot_cleanup.max_deviation);

    const auto given = step::write_step(m, shape, {.uncertainty = 1e-9, .knot_tol = 1e-5});
    REQUIRE(given.has_value());
    ASSERT_EQ(given->second.uncertainty, 1e-9);
    ASSERT_LE(given->second.knot_cleanup.max_deviation, 0.5e-9);
}
