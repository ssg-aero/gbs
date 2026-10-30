#include <doctest_gtest.hpp>
#include <gbs-io/step/geometry.h>
#include <gbs-io/step/geometry_writer.h>
#include <gbs/surfaceofrevolution.h>
#include <gbs/bscbuild.h>

#include <cmath>
#include <functional>
#include <numbers>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::step;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi;
    const ax2<T, 3> tilted{{{1., -2., 0.5}, {1., 1., 1.}, {1., 0., 0.}}};

    T dist(const point<T, 3> &a, const point<T, 3> &b) { return norm(a - b); }

    struct RoundTrip
    {
        P21Writer w;
        GeometryWriter g{w, 1e-6};
        std::uint64_t ctx = g.context(1., 1e-7);
        std::optional<P21File> f;
        std::optional<GeometryReader> r;

        void read()
        {
            auto parsed = parse_p21(w.str());
            REQUIRE(parsed.has_value());
            f = std::move(*parsed);
            r.emplace(*f, read_units(*f, ctx, 1.));
        }
    };

    // analytic curve C(t) given by a lambda (not a NURBS, not a primitive)
    class Lambda : public Curve<T, 3>
    {
        std::function<point<T, 3>(T)> f_;

    public:
        explicit Lambda(std::function<point<T, 3>(T)> f) : f_{std::move(f)} {}
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            if (d == 0)
                return f_(t);
            const T h = 1e-6;
            return (1. / (2 * h)) * (f_(t + h) - f_(t - h));
        }
        auto bounds() const -> std::array<T, 2> override { return {0., 1.}; }
    };

    // analytic surface S(u, v) on [0, 1]^2 given by a lambda (not a NURBS, not a primitive)
    class LambdaSurface : public Surface<T, 3>
    {
        std::function<point<T, 3>(T, T)> f_;

    public:
        explicit LambdaSurface(std::function<point<T, 3>(T, T)> f) : f_{std::move(f)} {}
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const T h = 1e-6;
            if (du == 0 && dv == 0)
                return f_(u, v);
            if (du == 1 && dv == 0)
                return (1. / (2 * h)) * (f_(u + h, v) - f_(u - h, v));
            return (1. / (2 * h)) * (f_(u, v + h) - f_(u, v - h));
        }
        auto bounds() const -> std::array<T, 4> override { return {0., 1., 0., 1.}; }
        auto has_second_derivatives() const -> bool override { return false; }
    };
}

TEST(tests_step_write_geometry, curves)
{
    RoundTrip rt;
    const LineCurve<T> line{{1., 2., 3.}, {0., 2., 1.}, -1., 3.};
    const EllipseCurve<T> circle{2.5, 2.5, tilted, 5.9, 6.6};
    const EllipseCurve<T> ellipse{3., 1.2, tilted, -1., 2.};
    const BSCurve<T, 3> bs{points_vector<T, 3>{{0., 0., 0.}, {1., 2., 0.}, {3., 2., 1.}, {4., 0., 0.}}, std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, 2};
    const auto rational = build_circle_arc<T>(2., 0.3, 2.5, tilted); // rational quadratic arcs
    const Lambda helix{[](T t) { return point<T, 3>{std::cos(4 * t), std::sin(4 * t), t}; }};

    const auto il = rt.g.curve(line, -1., 3.);
    const auto ic = rt.g.curve(circle, 5.9, 6.6);
    const auto ie = rt.g.curve(ellipse, -1., 2.);
    const auto ib = rt.g.curve(bs, 0.2, 0.9); // restricted to the edge range
    const auto ir = rt.g.curve(rational, 0.3, 2.5);
    const auto ih = rt.g.curve(helix, 0., 1.);
    rt.read();
    auto &r = *rt.r;
    ASSERT_EQ(rt.f->instance(il).type(), "LINE");
    ASSERT_EQ(rt.f->instance(ic).type(), "CIRCLE");
    ASSERT_EQ(rt.f->instance(ie).type(), "ELLIPSE");
    ASSERT_TRUE(rt.f->instance(ir).is_complex());
    // same parametrization: the analytic entities have the parameters of the primitives
    for (T t : {-1., 0.3, 3.})
        ASSERT_LT(dist(r.curve<3>(il)->value(t), line.value(t)), 1e-14);
    for (T t : {5.9, 6.3, 6.6})
        ASSERT_LT(dist(r.curve<3>(ic)->value(t), circle.value(t)), 1e-13);
    for (T t : {-1., 0.5, 2.})
        ASSERT_LT(dist(r.curve<3>(ie)->value(t), ellipse.value(t)), 1e-13);
    const auto b = r.curve<3>(ib);
    ASSERT_NEAR(b->domain()[0], 0.2, 1e-15);
    ASSERT_NEAR(b->domain()[1], 0.9, 1e-15);
    for (T t : {0.2, 0.5, 0.77, 0.9})
        ASSERT_LT(dist(b->value(t), bs.value(t)), 1e-13);
    for (T t : {0.3, 1., 2.5})
        ASSERT_LT(dist(r.curve<3>(ir)->value(t), rational.value(t)), 1e-13);
    // not a NURBS: approximated, parameter preserved, counted
    for (T t : {0., 0.37, 1.})
        ASSERT_LT(dist(r.curve<3>(ih)->value(t), helix.value(t)), 1e-6);
    ASSERT_EQ(rt.g.stats().approximated_curves, 1u);
    ASSERT_LE(rt.g.stats().max_deviation, 1e-6);
}

TEST(tests_step_write_geometry, iso_lines_and_pcurves)
{
    RoundTrip rt;
    // sides of natural faces: curves on surfaces whose pcurves are iso segments parametrized like the side
    auto cyl = std::make_shared<CylindricalSurface<T>>(1.5, tilted, std::array<T, 4>{0., 2 * pi, -1., 2.});
    auto seg = [](point<T, 2> a, point<T, 2> b, T t1, T t2) {
        return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{a, b}, std::vector<T>{t1, t1, t2, t2}, 1);
    };
    const CurveOnSurface<T, 3> generatrix{seg({0.7, -1.}, {0.7, 2.}, -1., 2.), cyl};       // u = 0.7: a line
    const CurveOnSurface<T, 3> parallel{seg({0., 0.5}, {2 * pi, 0.5}, 0., 2 * pi), cyl};    // v = 0.5: a circle
    auto nurbs = std::make_shared<BSSurfaceRational<T, 3>>(build_torus<T>(3., 1., tilted, 0., pi, 0., 2 * pi));
    const CurveOnSurface<T, 3> on_nurbs{seg({0.5, 0.}, {0.5, 2 * pi}, 0., 2 * pi), nurbs}; // iso of a B-spline surface
    const auto ig = rt.g.curve(generatrix, -1., 2.);
    const auto ip = rt.g.curve(parallel, 0., 2 * pi);
    const auto in = rt.g.curve(on_nurbs, 0., 2 * pi);
    // pcurves
    const auto il2 = rt.g.curve2d(*seg({1., 2.}, {3., 2.5}, 0.5, 1.5), 0.5, 1.5);
    const BSCurve<T, 2> bs2{points_vector<T, 2>{{0., 0.}, {1., 1.}, {2., 0.}}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2};
    const auto ib2 = rt.g.curve2d(bs2, 0., 1.);
    rt.read();
    auto &r = *rt.r;
    ASSERT_EQ(rt.f->instance(ig).type(), "LINE");
    ASSERT_EQ(rt.f->instance(ip).type(), "CIRCLE");
    for (T t : {-1., 0., 2.})
        ASSERT_LT(dist(r.curve<3>(ig)->value(t), generatrix.value(t)), 1e-13);
    for (T t : {0., 1., 4., 2 * pi})
        ASSERT_LT(dist(r.curve<3>(ip)->value(t), parallel.value(t)), 1e-13);
    for (T t : {0., 1., 4., 2 * pi})
        ASSERT_LT(dist(r.curve<3>(in)->value(t), on_nurbs.value(t)), 1e-12);
    ASSERT_EQ(rt.g.stats().approximated_curves, 0u); // all exact
    ASSERT_EQ(rt.f->instance(il2).type(), "LINE");
    for (T t : {0.5, 1., 1.5})
    {
        const auto a = r.curve<2>(il2)->value(t);
        ASSERT_NEAR(a[0], 1. + 2. * (t - 0.5), 1e-14);
        ASSERT_NEAR(a[1], 2. + 0.5 * (t - 0.5), 1e-14);
    }
    ASSERT_NEAR(r.curve<2>(ib2)->value(0.5)[1], 0.5, 1e-14);
}

TEST(tests_step_write_geometry, surfaces)
{
    RoundTrip rt;
    std::vector<std::shared_ptr<Surface<T, 3>>> ss{
        std::make_shared<Plane<T>>(tilted, std::array<T, 4>{-2., 3., -1., 4.}),
        std::make_shared<CylindricalSurface<T>>(1.5, tilted, std::array<T, 4>{-0.5, 2., -1., 3.}),
        std::make_shared<ConicalSurface<T>>(1., pi / 6, tilted, std::array<T, 4>{0., 2 * pi, -1.5, 2.}),
        std::make_shared<SphericalSurface<T>>(2., tilted, std::array<T, 4>{-1., 4., -pi / 2, 1.}),
        std::make_shared<ToroidalSurface<T>>(3., 0.8, tilted, std::array<T, 4>{5., 7., -1., 4.}),
        std::make_shared<ToroidalSurface<T>>(0.5, 1., tilted, std::array<T, 4>{0., 2., -1., 1.}), // degenerate torus
        std::make_shared<BSSurfaceRational<T, 3>>(build_sphere<T>(2., tilted, 0., 2., -1., 1.)),
        std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {1., 0., 1.}, {0., 1., 0.}, {1., 1., 2.}, {0., 2., 1.}, {1., 2., 0.}},
                                          std::vector<T>{0., 0., 1., 1.}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 1, 2),
    };
    std::vector<std::uint64_t> ids;
    for (const auto &s : ss)
        ids.push_back(rt.g.surface(*s));
    // a surface that is neither: approximated, parameters preserved
    auto profile = std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>({1., 0.}, {2., 1.}));
    const SurfaceOfRevolution<T> rev{profile, ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, 0., 2.};
    const auto irev = rt.g.surface(rev);
    rt.read();
    auto &r = *rt.r;
    const std::vector<std::string> types{"PLANE", "CYLINDRICAL_SURFACE", "CONICAL_SURFACE", "SPHERICAL_SURFACE", "TOROIDAL_SURFACE",
                                         "DEGENERATE_TOROIDAL_SURFACE", "", "B_SPLINE_SURFACE_WITH_KNOTS"};
    for (std::size_t k = 0; k < ss.size(); ++k)
    {
        INFO(k);
        if (!types[k].empty())
            ASSERT_EQ(rt.f->instance(ids[k]).type(), types[k]);
        const auto d = r.surface(ids[k]);
        const auto [u1, u2, v1, v2] = ss[k]->bounds();
        for (int i = 0; i <= 6; ++i)
            for (int j = 0; j <= 6; ++j)
            {
                const T u = u1 + (u2 - u1) * i / 6., v = v1 + (v2 - v1) * j / 6.;
                ASSERT_LT(dist(d->value(u, v), ss[k]->value(u, v)), 1e-12); // same parametrization
            }
    }
    ASSERT_TRUE(rt.f->instance(ids[6]).is_complex()); // rational
    const auto drev = r.surface(irev);
    for (T u : {0.4, 1.3})
        for (T v : {0.2, 0.9})
            ASSERT_LT(dist(drev->value(u, v), rev.value(u, v)), 1e-6);
    ASSERT_EQ(rt.g.stats().approximated_surfaces, 1u);
}

TEST(tests_step_write_geometry, units)
{
    struct Case
    {
        double mm;
        const char *name;
    };
    for (auto [mm, name] : {Case{1., "MILLI METRE"}, Case{10., "CENTI METRE"}, Case{1000., "METRE"}, Case{1e-3, "MICRO METRE"},
                            Case{25.4, "INCH"}, Case{304.8, "FOOT"}, Case{7., "MODEL UNIT"}})
    {
        INFO(name);
        P21Writer w;
        GeometryWriter g{w};
        const auto ctx = g.context(mm, 2e-6);
        const auto f = parse_p21(w.str());
        REQUIRE(f.has_value());
        const auto u = read_units(*f, ctx, 1.);
        ASSERT_NEAR(u.length_mm, mm, 1e-12 * mm);
        ASSERT_EQ(u.length_name, name);
        ASSERT_EQ(u.angle_name, "RADIAN");
        ASSERT_NEAR(*u.uncertainty, 2e-6 * mm, 1e-18 * mm); // in mm (target 1 mm)
    }
    P21Writer w;
    GeometryWriter g{w};
    ASSERT_THROW(g.context(0., 1e-7), std::invalid_argument);
    ASSERT_THROW(g.context(1., 0.), std::invalid_argument);
    ASSERT_THROW(g.direction({0., 0., 0.}), std::invalid_argument);
}

TEST(tests_step_write_geometry, approximations_within_tolerance)
{
    // a wavy curve and a wavy surface, neither primitive nor B-spline, are refined where they deviate until
    // within the tolerance: a uniform sampling stopped at 1025 points for a curve, 129 x 129 for a surface
    nurbs_export::ApproxStats<T> st;
    const Lambda c{[](T t) {
        const T r = 1. + 0.05 * std::sin(150. * t);
        return point<T, 3>{r * std::cos(t), r * std::sin(t), 0.};
    }};
    T dev{};
    std::visit([&](const auto &b) {
        for (int i = 0; i <= 100000; ++i)
            dev = std::max(dev, dist(b.value(2. * pi * i / 100000.), c.value(2. * pi * i / 100000.)));
    }, nurbs_export::to_nurbs<T, 3>(c, 0., 2. * pi, 1e-5, st));
    ASSERT_LE(dev, 1e-5);
    const LambdaSurface s{[](T u, T v) { return point<T, 3>{u, v, 0.05 * std::sin(20. * u) * std::cos(20. * v)}; }};
    dev = 0.;
    std::visit([&](const auto &b) {
        for (int j = 0; j <= 200; ++j)
            for (int i = 0; i <= 200; ++i)
                dev = std::max(dev, dist(b.value(i / 200., j / 200.), s.value(i / 200., j / 200.)));
    }, nurbs_export::to_nurbs<T>(s, 1e-5, st));
    ASSERT_LE(dev, 1e-5);
    ASSERT_EQ(st.approximated_curves, 1u);
    ASSERT_EQ(st.approximated_surfaces, 1u);
    ASSERT_EQ(st.over_tolerance, 0u);
    ASSERT_LE(st.max_deviation, 1e-5);
    // a step cannot be approximated: counted, its deviation reported
    nurbs_export::ApproxStats<T> st2;
    const Lambda step{[](T t) { return point<T, 3>{t, t < 0.5 ? 0. : 1., 0.}; }};
    (void)nurbs_export::to_nurbs<T, 3>(step, 0., 1., 1e-5, st2);
    ASSERT_EQ(st2.over_tolerance, 1u);
    ASSERT_GT(st2.max_deviation, 1e-5);
}
