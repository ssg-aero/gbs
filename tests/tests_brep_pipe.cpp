#include <doctest_gtest.hpp>
#include "tests_first_derivatives.h"
#include <gbs-brep/pipe.h>
#include <gbs-brep/check.h>
#include <gbs-brep/offset.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscinterp.h>
#include <gbs/extrema.h>
#include "tests_brep_geom.h"

#include <array>
#include <cmath>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using brep_tests::pi;
    using T = double;

    void expect_valid(const Model<T> &m, const ShapeId &s)
    {
        const auto rep = check(m, s);
        std::ostringstream os;
        for (const auto &e : rep.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        INFO(os.str());
        REQUIRE(rep.ok());
    }

    T volume(const Model<T> &m, const ShapeId &s)
    {
        T v{};
        for (auto so : explore<SolidId>(m, s))
            v += signed_volume(m, m.solid(so).outer, 160);
        return v;
    }

    template <typename S>
    std::size_t count(const Model<T> &m, const ShapeId &s)
    {
        std::size_t n{};
        for (auto f : explore<FaceId>(m, s))
            n += dynamic_cast<const S *>(m.face(f).surface.get()) != nullptr;
        return n;
    }

    /// The disk of radius r centered at c, normal to n.
    FaceId disk(Model<T> &m, T r, const point<T, 3> &c, const point<T, 3> &n, const point<T, 3> &x)
    {
        const ax2<T, 3> ax{{c, n, x}};
        const auto e = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(r, r, ax)));
        const auto w = unwrap(make_wire(m, std::vector<EdgeId>{e}));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax, std::array<T, 4>{-2 * r, 2 * r, -2 * r, 2 * r}), w));
    }

    /// The square of side a centered at c, in the plane spanned by x and y.
    FaceId square(Model<T> &m, T a, const point<T, 3> &c, const point<T, 3> &x, const point<T, 3> &y)
    {
        const auto h = a / 2;
        const auto p00 = c - h * x - h * y, p10 = c + h * x - h * y, p11 = c + h * x + h * y, p01 = c - h * x + h * y;
        std::vector<EdgeId> e{unwrap(make_edge(m, p00, p10)), unwrap(make_edge(m, p10, p11)), unwrap(make_edge(m, p11, p01)),
                              unwrap(make_edge(m, p01, p00))};
        const auto w = unwrap(make_wire(m, e));
        return unwrap(make_face(m, std::make_shared<Plane<T>>(ax2<T, 3>{{c, cross(x, y), x}}, std::array<T, 4>{-a, a, -a, a}), w));
    }

    /// Cubic B-spline through the points of f on [0, 1].
    template <typename F>
    std::shared_ptr<BSCurve<T, 3>> spline(F f, std::size_t n = 64)
    {
        std::vector<T> u(n + 1);
        std::vector<constrType<T, 3, 1>> q(n + 1);
        for (std::size_t i = 0; i <= n; ++i)
        {
            u[i] = T(i) / T(n);
            q[i] = {f(u[i])};
        }
        return std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(q, u, std::size_t{3}));
    }

    std::array<std::size_t, 7> counts(const Model<T> &m)
    {
        std::array<std::size_t, 7> n{};
        for (std::size_t t = 0; t < n.size(); ++t)
            n[t] = m.count(ShapeType(t));
        return n;
    }

    /// The faces of a tube of radius r around c (on its bounds), but its caps, are at r from c (within tol).
    void expect_tube(const Model<T> &m, const ShapeId &s, const Curve<T, 3> &c, T r, T tol)
    {
        for (auto f : explore<FaceId>(m, s))
        {
            const auto &srf = *m.face(f).surface;
            if (dynamic_cast<const Plane<T> *>(&srf))
                continue;
            const auto b = srf.bounds();
            for (int i = 0; i <= 12; ++i)
                for (int j = 0; j <= 40; ++j)
                {
                    const auto p = srf.value(b[0] + (b[1] - b[0]) * i / 12., b[2] + (b[3] - b[2]) * j / 40.);
                    ASSERT_NEAR(extrema_curve_point(c, p, 1e-12)[1], r, tol);
                }
        }
    }

}

TEST(tests_brep_pipe, along_a_circle)
{
    // a disk r 1 along the circle R 4 of the plane z = 0: a torus, exact, 2 pi² R r²
    Model<T> m;
    const auto d = disk(m, 1., {4., 0., 0.}, {0., 1., 0.}, {1., 0., 0.});
    const ax2<T, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    const auto circle = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(4., 4., ax)));
    auto t = unwrap(make_pipe(m, d, circle));
    expect_valid(m, t.shape);
    ASSERT_EQ(t.approximated, 0u);
    ASSERT_EQ(explore<FaceId>(m, t.shape).size(), 1u);
    ASSERT_EQ(count<ToroidalSurface<T>>(m, t.shape), 1u);
    ASSERT_NEAR(volume(m, t.shape), 2 * pi * pi * 4., 2e-3 * 2 * pi * pi * 4.);
    // Frenet turns the same way along a circle
    auto fr = unwrap(make_pipe(m, d, circle, PipeOptions<T>{.frame = PipeFrame::Frenet}));
    ASSERT_EQ(count<ToroidalSurface<T>>(m, fr.shape), 1u);
    // a quarter of it: a quarter of the torus, closed by the disk and its turned copy
    const auto arc = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(4., 4., ax, 0., pi / 2)));
    auto q = unwrap(make_pipe(m, d, arc));
    expect_valid(m, q.shape);
    ASSERT_EQ(count<ToroidalSurface<T>>(m, q.shape), 1u);
    ASSERT_EQ(count<Plane<T>>(m, q.shape), 2u);
    ASSERT_NEAR(volume(m, q.shape), pi * pi * 2., 2e-3 * pi * pi * 2.);
}

TEST(tests_brep_pipe, along_segments_and_arcs)
{
    // a square 1 x 1 along a segment of 5: a box
    Model<T> m;
    const auto s = square(m, 1., {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.});
    const auto seg = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 5.}));
    auto b = unwrap(make_pipe(m, s, seg));
    expect_valid(m, b.shape);
    ASSERT_EQ(count<Plane<T>>(m, b.shape), 6u);
    ASSERT_NEAR(volume(m, b.shape), 5., 1e-9);

    // a disk r 0.5 along an elbow: segment, quarter of circle R 2, segment; cylinders and a torus,
    // pi r² times the length (Pappus)
    const auto d = disk(m, 0.5, {0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 0., 3.}, {0., 1., 0.}, {-1., 0., 0.}}}, 0., pi / 2);
    ASSERT_NEAR(norm(arc->value(0.) - point<T, 3>{0., 0., 3.}), 0., 1e-14);
    const auto e2 = unwrap(make_edge(m, arc));
    const auto e3 = unwrap(make_edge(m, arc->value(pi / 2), point<T, 3>{6., 0., 5.}));
    const auto elbow = unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3}));
    auto p = unwrap(make_pipe(m, d, elbow));
    expect_valid(m, p.shape);
    ASSERT_EQ(p.approximated, 0u);
    ASSERT_EQ(count<CylindricalSurface<T>>(m, p.shape), 2u);
    ASSERT_EQ(count<ToroidalSurface<T>>(m, p.shape), 1u);
    ASSERT_EQ(count<Plane<T>>(m, p.shape), 2u);
    const T L = 3. + pi + 4.;
    ASSERT_NEAR(volume(m, p.shape), pi * 0.25 * L, 2e-3 * pi * 0.25 * L);
    // the wire profile gives a shell, the edge a shell of a face per path edge
    auto w = unwrap(make_pipe(m, ShapeId{m.face(d).wires.front()}, elbow));
    ASSERT_TRUE(std::holds_alternative<ShellId>(w.shape));
    ASSERT_EQ(m.shell(std::get<ShellId>(w.shape)).faces.size(), 3u);
    auto e = unwrap(make_pipe(m, explore<EdgeId>(m, ShapeId{d}).front(), elbow));
    ASSERT_TRUE(std::holds_alternative<ShellId>(e.shape));
    // one edge along one segment: a face
    auto f = unwrap(make_pipe(m, explore<EdgeId>(m, ShapeId{d}).front(), seg));
    ASSERT_TRUE(std::holds_alternative<FaceId>(f.shape));
}

TEST(tests_brep_pipe, along_a_helix)
{
    // a square 0.4 x 0.4 normal to a helix (radius 2, pitch 1.5, two turns), rotation-minimizing:
    // B-spline faces within approx_tol, volume = area x length (Pappus: the centroid on the path)
    Model<T> m;
    const T a = 2., b = 1.5 / (2 * pi);
    const auto h = spline([&](T u) { const T t = 4 * pi * u; return point<T, 3>{a * std::cos(t), a * std::sin(t), b * t}; }, 96);
    const auto path = unwrap(make_edge(m, h, 0., 1.));
    const auto t0 = h->value(0., 1);
    const point<T, 3> x{1., 0., 0.};
    const auto y = cross(t0, x) / norm(cross(t0, x));
    const auto s = square(m, 0.4, h->value(0.), x, y);
    PipeOptions<T> o;
    o.approx_tol = 1e-5;
    auto r = unwrap(make_pipe(m, s, path, o));
    expect_valid(m, r.shape);
    ASSERT_EQ(r.approximated, 4u);
    ASSERT_LE(r.deviation, 1e-5);
    ASSERT_EQ(count<Plane<T>>(m, r.shape), 2u);
    const T L = gbs::length(*h, 0., 1.);
    ASSERT_NEAR(volume(m, r.shape), 0.16 * L, 2e-3 * 0.16 * L);
    // a disk normal to the path sweeps a tube, within approx_tol
    const auto d = disk(m, 0.3, h->value(0.), t0 / norm(t0), x);
    auto tube = unwrap(make_pipe(m, d, path, o));
    expect_valid(m, tube.shape);
    expect_tube(m, tube.shape, *h, 0.3, 1.1e-5);
    // a fixed binormal along z keeps the square level: same volume
    auto fb = unwrap(make_pipe(m, s, path, PipeOptions<T>{.frame = PipeFrame::FixedBinormal, .binormal = {0., 0., 1.}}));
    expect_valid(m, fb.shape);
    ASSERT_NEAR(volume(m, fb.shape), 0.16 * L, 2e-3 * 0.16 * L);
}

TEST(tests_brep_pipe, along_a_closed_path)
{
    // a disk r 0.5 around a stadium: two half circles R 3 joined by segments of 4, G1 and closed;
    // two tori and two cylinders, no caps, pi r² L
    Model<T> m;
    const auto right = std::make_shared<EllipseCurve<T>>(3., 3., ax2<T, 3>{{{4., 0., 0.}, {0., 0., 1.}, {0., -1., 0.}}}, 0., pi);
    const auto left = std::make_shared<EllipseCurve<T>>(3., 3., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}}}, 0., pi);
    const auto va = unwrap(make_vertex(m, point<T, 3>{4., -3., 0.})), vb = unwrap(make_vertex(m, point<T, 3>{4., 3., 0.}));
    const auto vc = unwrap(make_vertex(m, point<T, 3>{0., 3., 0.})), vd = unwrap(make_vertex(m, point<T, 3>{0., -3., 0.}));
    std::vector<EdgeId> es{unwrap(make_edge(m, right, 0., pi, va, vb)), unwrap(make_edge(m, vb, vc)), unwrap(make_edge(m, left, 0., pi, vc, vd)),
                           unwrap(make_edge(m, vd, va))};
    const auto stadium = unwrap(make_wire(m, es));
    ASSERT_TRUE(m.wire(stadium).closed);
    const auto d = disk(m, 0.5, {4., -3., 0.}, {1., 0., 0.}, {0., 0., 1.});
    auto st = unwrap(make_pipe(m, d, stadium));
    expect_valid(m, st.shape);
    ASSERT_EQ(st.approximated, 0u);
    ASSERT_EQ(count<ToroidalSurface<T>>(m, st.shape), 2u);
    ASSERT_EQ(count<CylindricalSurface<T>>(m, st.shape), 2u);
    ASSERT_EQ(count<Plane<T>>(m, st.shape), 0u);
    const T L = 2 * 3 * pi + 2 * 4.;
    ASSERT_NEAR(volume(m, st.shape), pi * 0.25 * L, 2e-3 * pi * 0.25 * L);
}

namespace
{
    /// (4 cos t, 3 sin t, sin 2t): a closed space curve
    class Saddle : public Curve<T, 3>
    {
    public:
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            const T h = T(d) * pi / 2;
            return {4 * std::cos(t + h), 3 * std::sin(t + h), (d == 0 ? 1. : std::pow(2., T(d))) * std::sin(2 * t + h)};
        }
        auto bounds() const -> std::array<T, 2> override { return {0., 2 * pi}; }
    };
}

TEST(tests_brep_pipe, closed_space_curve)
{
    // a disk around a closed curve out of any plane: the rotation-minimizing frame comes back
    // turned by the holonomy of the curve, spread along it; B-spline faces, no caps, pi r² L
    Model<T> m;
    const auto c = std::make_shared<Saddle>();
    const auto path = unwrap(make_edge(m, c));
    ASSERT_EQ(m.edge(path).v1, m.edge(path).v2);
    const auto t0 = c->value(0., 1);
    const auto d = disk(m, 0.3, c->value(0.), t0, cross(t0, point<T, 3>{0., 0., 1.}));
    auto r = unwrap(make_pipe(m, d, path));
    expect_valid(m, r.shape);
    ASSERT_EQ(count<Plane<T>>(m, r.shape), 0u);
    ASSERT_GE(r.approximated, 1u);
    const T L = gbs::length(*c, 0., 2 * pi);
    ASSERT_NEAR(volume(m, r.shape), pi * 0.09 * L, 2e-3 * pi * 0.09 * L);
    // the disk is normal to the path: whatever the frame, a tube of radius r around it
    ASSERT_LE(r.deviation, 1e-5);
    expect_tube(m, r.shape, *c, 0.3, 1.1e-5);
    // Frenet along it: also closed (the curve is closed with its derivatives)
    auto fr = unwrap(make_pipe(m, d, path, PipeOptions<T>{.frame = PipeFrame::Frenet}));
    expect_valid(m, fr.shape);
}

TEST(tests_brep_pipe, refusals)
{
    Model<T> m;
    const auto d = disk(m, 0.5, {0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    // a corner
    const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., 3.}, point<T, 3>{3., 0., 3.}));
    const auto corner = unwrap(make_wire(m, std::vector<EdgeId>{e1, e2}));
    auto r = make_pipe(m, d, corner);
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    // a profile wider than the radius of curvature: folds
    const auto tight = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(0.4, 0.4, ax2<T, 3>{{{0.4, 0., 0.}, {0., 1., 0.}, {-1., 0., 0.}}}, 0., pi / 2)));
    auto f = make_pipe(m, d, tight);
    ASSERT_FALSE(f.has_value());
    ASSERT_EQ(f.error().code, BuildErrc::NonGeneric);
    // Frenet along a segment, a fixed binormal along the tangent
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.frame = PipeFrame::Frenet}).error().code, BuildErrc::NonGeneric);
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.frame = PipeFrame::FixedBinormal, .binormal = {0., 0., 1.}}).error().code, BuildErrc::NonGeneric);
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.frame = PipeFrame::FixedBinormal}).error().code, BuildErrc::InvalidBounds);
    // a path that is neither an edge nor a wire, a dead id
    ASSERT_EQ(make_pipe(m, d, ShapeId{d}).error().code, BuildErrc::InvalidId);
    ASSERT_EQ(make_pipe(m, d, ShapeId{EdgeId{9999}}).error().code, BuildErrc::InvalidId);
    // a profile that is neither an edge, a wire nor a face, before any check of the geometry
    const auto box = unwrap(make_pipe(m, d, e1));
    ASSERT_EQ(make_pipe(m, box.shape, tight).error().code, BuildErrc::InvalidId);
    ASSERT_EQ(make_pipe(m, ShapeId{m.edge(e1).v1}, e1).error().code, BuildErrc::InvalidId);
    // tolerances that are not positive; an accuracy that the frame cannot reach
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.approx_tol = 0.}).error().code, BuildErrc::InvalidBounds);
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.approx_tol = std::nan("")}).error().code, BuildErrc::InvalidBounds);
    ASSERT_EQ(make_pipe(m, d, e1, PipeOptions<T>{.tol = -1.}).error().code, BuildErrc::InvalidBounds);
    const auto h = spline([](T u) { return point<T, 3>{std::cos(6 * u), std::sin(6 * u), u}; });
    ASSERT_EQ(make_pipe(m, d, unwrap(make_edge(m, h, 0., 1.)), PipeOptions<T>{.approx_tol = 1e-300}).error().code, BuildErrc::PCurveApproximation);
}

namespace
{
    /// A tent from (1, 0, 0) to (2, 0, 0), its top at (1.5, 0.1, 0): C0 only, so that no cubic B-spline
    /// through its points comes within 1e-12 of it.
    class Tent : public Curve<T, 3>
    {
    public:
        auto value(T u, size_t d = 0) const -> std::array<T, 3> override
        {
            const T sg = u < 1.5 ? 1. : -1.;
            if (d == 0)
                return {u, 0.2 * (0.5 - std::abs(u - 1.5)), 0.};
            return d == 1 ? std::array<T, 3>{1., 0.2 * sg, 0.} : std::array<T, 3>{0., 0., 0.};
        }
        auto bounds() const -> std::array<T, 2> override { return {1., 2.}; }
    };
}

TEST(tests_brep_pipe, failure_leaves_the_model)
{
    // a wire profile: a segment, swept into a plane, then a tent that cannot be approximated within
    // approx_tol: refused once the face of the segment is made, which is removed
    Model<T> m;
    const auto a = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}));
    const auto b = unwrap(make_edge(m, std::make_shared<Tent>()));
    const auto profile = unwrap(make_wire(m, std::vector<EdgeId>{a, b}));
    const auto path = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto before = counts(m);
    auto r = make_pipe(m, profile, path, PipeOptions<T>{.approx_tol = 1e-12});
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::PCurveApproximation);
    ASSERT_EQ(counts(m), before);
    // the same profile within 1e-4: a shell of a plane and an extruded B-spline
    auto ok = unwrap(make_pipe(m, profile, path, PipeOptions<T>{.approx_tol = 1e-4}));
    ASSERT_TRUE(std::holds_alternative<ShellId>(ok.shape));
    ASSERT_EQ(count<Plane<T>>(m, ok.shape), 1u);
}

TEST(tests_brep_pipe, self_intersection)
{
    // a square 0.4 x 0.4 along a helix of pitch 0.2: the turns overlap, refused; of pitch 0.5, it passes
    for (const T pitch : {0.2, 0.5})
    {
        Model<T> m;
        const auto h = spline([&](T u) { const T t = 4 * pi * u; return point<T, 3>{2 * std::cos(t), 2 * std::sin(t), pitch * t / (2 * pi)}; }, 96);
        const auto t0 = h->value(0., 1);
        const point<T, 3> x{1., 0., 0.};
        const auto s = square(m, 0.4, h->value(0.), x, cross(t0, x) / norm(cross(t0, x)));
        const auto path = unwrap(make_edge(m, h, 0., 1.));
        const auto before = counts(m);
        auto r = make_pipe(m, s, path);
        if (pitch < 0.4)
        {
            ASSERT_FALSE(r.has_value());
            ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
            ASSERT_EQ(counts(m), before);
        }
        else
        {
            ASSERT_TRUE(r.has_value());
            expect_valid(m, r->shape);
        }
    }
    // a ring r 0.1 along a figure eight (sin t, sin 2t, 0), through its crossing: refused
    Model<T> m;
    const auto eight = spline([](T u) { const T t = -1. + (pi + 2.) * u; return point<T, 3>{std::sin(t), std::sin(2 * t), 0.}; }, 128);
    const auto t0 = eight->value(0., 1) / norm(eight->value(0., 1));
    const ax2<T, 3> ax{{eight->value(0.), t0, {0., 0., 1.}}};
    const auto ring = unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(0.1, 0.1, ax)));
    auto r = make_pipe(m, ring, unwrap(make_edge(m, eight, 0., 1.)));
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
}

TEST(tests_brep_pipe, junctions)
{
    // a corner of 1e-8 rad, within the tolerance: the planes of both segments stay exact
    {
        Model<T> m;
        const auto s = square(m, 1., {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.});
        const T a = 1e-8;
        const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 5.}));
        const auto e2 = unwrap(make_edge(m, point<T, 3>{0., 0., 5.}, point<T, 3>{5 * std::sin(a), 0., 5. + 5 * std::cos(a)}));
        auto b = unwrap(make_pipe(m, s, unwrap(make_wire(m, std::vector<EdgeId>{e1, e2}))));
        expect_valid(m, b.shape);
        ASSERT_EQ(b.approximated, 0u);
        ASSERT_EQ(count<Plane<T>>(m, b.shape), 10u);
        ASSERT_NEAR(volume(m, b.shape), 10., 1e-6);
    }
    // an S of two quarters of circle R 2 turning opposite ways, in the plane z = 0: rotation-minimizing,
    // two exact tori; Frenet reverses at the junction, refused at once
    {
        Model<T> m;
        const auto a1 = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 0., 0.}, {0., 0., -1.}, {-1., 0., 0.}}}, 0., pi / 2);
        const auto a2 = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 4., 0.}, {0., 0., 1.}, {0., -1., 0.}}}, 0., pi / 2);
        ASSERT_NEAR(norm(a1->value(pi / 2) - a2->value(0.)), 0., 1e-14);
        const auto w = unwrap(make_wire(m, std::vector<EdgeId>{unwrap(make_edge(m, a1)), unwrap(make_edge(m, a2))}));
        const auto d = disk(m, 0.5, {0., 0., 0.}, {0., 1., 0.}, {0., 0., 1.});
        auto r = unwrap(make_pipe(m, d, w));
        expect_valid(m, r.shape);
        ASSERT_EQ(r.approximated, 0u);
        ASSERT_EQ(count<ToroidalSurface<T>>(m, r.shape), 2u);
        ASSERT_NEAR(volume(m, r.shape), pi * 0.25 * 2 * pi, 2e-3 * pi * 0.25 * 2 * pi);
        auto f = make_pipe(m, d, w, PipeOptions<T>{.frame = PipeFrame::Frenet});
        ASSERT_FALSE(f.has_value());
        ASSERT_EQ(f.error().code, BuildErrc::NonGeneric);
    }
    // two quarters of circle tangent at their junction, in planes 30 degrees apart: Frenet turns the
    // profile there by 30 degrees, refused; rotation-minimizing, two exact tori
    {
        Model<T> m;
        const auto a1 = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 0., 0.}, {0., 0., -1.}, {-1., 0., 0.}}}, 0., pi / 2);
        const T c = std::cos(pi / 6), sn = std::sin(pi / 6);
        // the second plane: turned by 30 degrees around the tangent at the junction (along x)
        const point<T, 3> z2{0., sn, -c}, x2{0., c, sn};
        const auto end = a1->value(pi / 2);
        const auto a2 = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{end - 2. * x2, z2, x2}}, 0., pi / 2);
        ASSERT_NEAR(norm(a2->value(0.) - end), 0., 1e-14);
        ASSERT_NEAR(norm(cross(a2->value(0., 1), a1->value(pi / 2, 1))), 0., 1e-12);
        const auto w = unwrap(make_wire(m, std::vector<EdgeId>{unwrap(make_edge(m, a1)), unwrap(make_edge(m, a2))}));
        const auto d = disk(m, 0.5, {0., 0., 0.}, {0., 1., 0.}, {0., 0., 1.});
        auto r = unwrap(make_pipe(m, d, w));
        expect_valid(m, r.shape);
        ASSERT_EQ(count<ToroidalSurface<T>>(m, r.shape), 2u);
        auto f = make_pipe(m, d, w, PipeOptions<T>{.frame = PipeFrame::Frenet});
        ASSERT_FALSE(f.has_value());
        ASSERT_EQ(f.error().code, BuildErrc::NonGeneric);
    }
    // a closed stadium whose last edge ends 5e-7 off its first one, at a vertex of tolerance 1e-6:
    // closed within the tolerance of that vertex
    {
        Model<T> m;
        const auto right = std::make_shared<EllipseCurve<T>>(3., 3., ax2<T, 3>{{{4., 0., 0.}, {0., 0., 1.}, {0., -1., 0.}}}, 0., pi);
        const auto left = std::make_shared<EllipseCurve<T>>(3., 3., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}}}, 0., pi);
        const auto va = unwrap(make_vertex(m, point<T, 3>{4., -3., 0.}, 1e-6)), vb = unwrap(make_vertex(m, point<T, 3>{4., 3., 0.}));
        const auto vc = unwrap(make_vertex(m, point<T, 3>{0., 3., 0.})), vd = unwrap(make_vertex(m, point<T, 3>{0., -3., 0.}));
        const auto last = std::make_shared<BSCurve<T, 3>>(std::vector<std::array<T, 3>>{{0., -3., 0.}, {4., -3., 5e-7}}, std::vector<T>{0., 0., 1., 1.}, 1);
        const auto stadium = unwrap(make_wire(m, std::vector<EdgeId>{unwrap(make_edge(m, right, 0., pi, va, vb)), unwrap(make_edge(m, vb, vc)),
                                                                     unwrap(make_edge(m, left, 0., pi, vc, vd)), unwrap(make_edge(m, last, 0., 1., vd, va))}));
        ASSERT_TRUE(m.wire(stadium).closed);
        const auto d = disk(m, 0.5, {4., -3., 0.}, {1., 0., 0.}, {0., 0., 1.});
        auto st = unwrap(make_pipe(m, d, stadium));
        expect_valid(m, st.shape);
        ASSERT_NEAR(volume(m, st.shape), pi * 0.25 * (6 * pi + 8.), 2e-3 * pi * 0.25 * (6 * pi + 8.));
    }
}

TEST(tests_brep_pipe, any_curve_of_the_path)
{
    // the elbow with its arc as a rational B-spline: recognized on its points, an exact torus
    {
        Model<T> m;
        const auto d = disk(m, 0.5, {0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
        const auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{2., 0., 3.}, {0., 1., 0.}, {-1., 0., 0.}}}, 0., pi / 2);
        const auto nurbs = arc->to_nurbs();
        ASSERT_TRUE((dynamic_cast<const BSCurveRational<T, 3> *>(nurbs.get()) != nullptr));
        const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 3.}));
        const auto e2 = unwrap(make_edge(m, nurbs));
        const auto e3 = unwrap(make_edge(m, arc->value(pi / 2), point<T, 3>{6., 0., 5.}));
        auto p = unwrap(make_pipe(m, d, unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3}))));
        expect_valid(m, p.shape);
        ASSERT_EQ(p.approximated, 0u);
        ASSERT_EQ(count<ToroidalSurface<T>>(m, p.shape), 1u);
        ASSERT_EQ(count<CylindricalSurface<T>>(m, p.shape), 2u);
    }
    // a path without second derivatives (an edge of a face offset from a B-spline seen through its
    // first derivatives only): no exception,
    // the folds measured on differences of its tangent; Frenet needs them, refused
    {
        Model<T> m;
        std::vector<std::array<T, 3>> poles;
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 4; ++i)
                poles.push_back({T(i), T(j), 0.1 * i * j});
        const std::vector<T> k{0., 0., 0., 0., 1., 1., 1., 1.};
        const auto f = unwrap(make_face(m, std::make_shared<gbs_tests::FirstDerivativesOnly<T>>(std::make_shared<BSSurface<T, 3>>(poles, k, k, 3, 3))));
        const auto g = unwrap(offset_face(m, f, 0.1));
        const auto pe = explore<EdgeId>(m, ShapeId{g}).front();
        ASSERT_FALSE(m.edge(pe).curve->has_second_derivatives());
        const auto c = m.edge(pe).curve->value(m.edge(pe).u1), t = m.edge(pe).curve->value(m.edge(pe).u1, 1);
        const auto x = cross(t, point<T, 3>{0., 0., 1.}) / norm(cross(t, point<T, 3>{0., 0., 1.}));
        const auto s = square(m, 0.1, c, x, cross(t, x) / norm(cross(t, x)));
        auto r = make_pipe(m, s, pe);
        ASSERT_TRUE(r.has_value());
        expect_valid(m, r->shape);
        auto fr = make_pipe(m, s, pe, PipeOptions<T>{.frame = PipeFrame::Frenet});
        ASSERT_FALSE(fr.has_value());
        ASSERT_EQ(fr.error().code, BuildErrc::NonGeneric);
    }
}
