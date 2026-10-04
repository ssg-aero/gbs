#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/bscbuild.h>

#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <unordered_set>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;
    const T tol = 1e-6;
    constexpr T pi = std::numbers::pi_v<T>;

    // Analytic surface S(u,v) given by a lambda; first derivatives by central differences.
    class Analytic : public Surface<T, 3>
    {
        std::function<point<T, 3>(T, T)> f_;
        std::array<T, 4> b_;

    public:
        Analytic(std::function<point<T, 3>(T, T)> f, std::array<T, 4> b) : f_{std::move(f)}, b_{b} {}
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const T h = 1e-6;
            if (du == 0 && dv == 0)
                return f_(u, v);
            if (du == 1 && dv == 0)
                return (f_(u + h, v) - f_(u - h, v)) * (T(0.5) / h);
            if (du == 0 && dv == 1)
                return (f_(u, v + h) - f_(u, v - h)) * (T(0.5) / h);
            throw std::runtime_error("not implemented");
        }
        auto bounds() const -> std::array<T, 4> override { return b_; }
    };

    std::shared_ptr<Surface<T, 3>> cylinder(T r = 1., T h = 2.)
    {
        return std::make_shared<Analytic>([r](T u, T v) { return point<T, 3>{r * std::cos(u), r * std::sin(u), v}; },
                                          std::array<T, 4>{0., 2. * pi, 0., h});
    }
    std::shared_ptr<Surface<T, 3>> sphere(T r = 1.)
    {
        return std::make_shared<Analytic>(
            [r](T u, T v) { return point<T, 3>{r * std::cos(v) * std::cos(u), r * std::cos(v) * std::sin(u), r * std::sin(v)}; },
            std::array<T, 4>{0., 2. * pi, -pi / 2., pi / 2.});
    }
    std::shared_ptr<Surface<T, 3>> torus(T R = 3., T r = 1.)
    {
        return std::make_shared<Analytic>(
            [R, r](T u, T v) { return point<T, 3>{(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)}; },
            std::array<T, 4>{0., 2. * pi, 0., 2. * pi});
    }
    std::shared_ptr<Surface<T, 3>> cone()
    {
        return std::make_shared<Analytic>([](T u, T v) { return point<T, 3>{v * std::cos(u), v * std::sin(u), v}; },
                                          std::array<T, 4>{0., 2. * pi, 0., 1.});
    }

    // Exact rational NURBS cylinder: the rational circle extruded along z.
    std::shared_ptr<Surface<T, 3>> nurbs_cylinder(T h = 2.)
    {
        auto circle = build_circle<T, 3>(1.);
        points_vector<T, 4> poles;
        for (T z : {T(0), h})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z; // weighted homogeneous coordinates
                poles.push_back(p);
            }
        return std::make_shared<BSSurfaceRational<T, 3>>(poles, circle.knotsFlats(), std::vector<T>{0., 0., h, h},
                                                         circle.degree(), 1);
    }

    std::shared_ptr<Surface<T, 3>> plane()
    {
        points_vector<T, 3> poles{{0., 0., 0.}, {2., 0., 0.}, {0., 1., 0.}, {2., 1., 0.}};
        std::vector<T> k{0., 0., 1., 1.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);
    }

    std::shared_ptr<Surface<T, 3>> revolution(T angle)
    {
        auto profile = std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>({1., 1.}, {2., 3.}));
        return std::make_shared<SurfaceOfRevolution<T>>(profile, ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, 0., angle);
    }

    struct Counts
    {
        size_t vertices, edges, degenerate, seams, free_edges;
        bool closed_shell;
    };

    // Structural checks common to every natural-bounds face, then the topology counts.
    Counts check_face(const Model<T> &m, FaceId fid)
    {
        const auto &f = m.face(fid);
        REQUIRE(f.natural_bounds);
        REQUIRE(f.wires.size() == 1);
        const auto wid = f.wires.front();
        const auto &w = m.wire(wid);
        REQUIRE(w.coedges.size() == 4);
        REQUIRE(is_closed(m, wid));

        // counter-clockwise in (u,v), area of the parametric rectangle
        auto [u1, u2, v1, v2] = f.surface->bounds();
        CHECK(uv_signed_area(m, wid) == doctest::Approx((u2 - u1) * (v2 - v1)).epsilon(1e-9));

        // pcurve through the surface == 3D point of the co-edge, orientation applied
        for (size_t i = 0; i < 4; ++i)
        {
            const auto &ce = w.coedges[i];
            const auto &e = m.edge(ce.edge);
            REQUIRE(ce.pcurve);
            for (T s : {0., 0.2, 0.5, 0.8, 1.})
            {
                const T t = e.u1 + s * (e.u2 - e.u1);
                const auto uv = coedge_uv(m, ce, t);
                CHECK(distance(f.surface->value(uv[0], uv[1]), coedge_point(m, wid, i, t)) <= e.tol + 1e-12);
            }
            CHECK(m.vertex(e.v1).tol >= e.tol);
            CHECK(m.vertex(e.v2).tol >= e.tol);
        }

        Counts c{};
        auto edges = explore<EdgeId>(m, fid);
        c.vertices = explore<VertexId>(m, fid).size();
        c.edges = edges.size();
        for (auto eid : edges)
            if (m.edge(eid).degenerate)
                ++c.degenerate;
        for (auto eid : edges)
        {
            size_t uses = 0;
            for (const auto &ce : w.coedges)
                uses += ce.edge == eid;
            if (uses == 2)
            {
                ++c.seams;
                // a seam is used once in each sense, with two different pcurves
                std::vector<const CoEdge<T> *> u;
                for (const auto &ce : w.coedges)
                    if (ce.edge == eid)
                        u.push_back(&ce);
                CHECK(u[0]->orient != u[1]->orient);
                CHECK(u[0]->pcurve != u[1]->pcurve);
            }
        }
        Model<T> &mm = const_cast<Model<T> &>(m);
        auto sh = mm.add(Shell{{FaceUse{fid, Orientation::Forward}}});
        c.free_edges = free_edges(m, sh).size();
        c.closed_shell = is_closed(m, sh);
        CHECK(is_orientable(m, sh));
        mm.erase(sh);
        return c;
    }
}

TEST(tests_brep_face, surface_closure)
{
    auto c = surface_closure(*cylinder(), tol);
    ASSERT_TRUE(c.closed_u);
    ASSERT_FALSE(c.closed_v);
    ASSERT_FALSE(c.degenerate());

    c = surface_closure(*sphere(), tol);
    ASSERT_TRUE(c.closed_u);
    ASSERT_TRUE(c.degenerate_v1 && c.degenerate_v2);
    ASSERT_FALSE(c.degenerate_u1 || c.degenerate_u2);
    ASSERT_FALSE(c.closed_v); // degenerate sides are never reported as a seam

    c = surface_closure(*torus(), tol);
    ASSERT_TRUE(c.closed_u && c.closed_v);

    c = surface_closure(*cone(), tol);
    ASSERT_TRUE(c.closed_u);
    ASSERT_TRUE(c.degenerate_v1);
    ASSERT_FALSE(c.degenerate_v2);

    c = surface_closure(*plane(), tol);
    ASSERT_FALSE(c.closed_u || c.closed_v || c.degenerate());

    c = surface_closure(*nurbs_cylinder(), tol);
    ASSERT_TRUE(c.closed_u);
    ASSERT_FALSE(c.closed_v || c.degenerate());
}

TEST(tests_brep_face, plane)
{
    Model<T> m;
    auto f = unwrap(make_face(m, plane()));
    auto c = check_face(m, f);
    ASSERT_EQ(c.vertices, 4);
    ASSERT_EQ(c.edges, 4);
    ASSERT_EQ(c.degenerate, 0);
    ASSERT_EQ(c.seams, 0);
    ASSERT_EQ(c.free_edges, 4);
    ASSERT_FALSE(c.closed_shell);
    // the outer wire starts at (u1, v1) along v = v1
    const auto &w = m.wire(m.face(f).wires.front());
    ASSERT_NEAR(distance(m.vertex(coedge_start(m, w.coedges[0])).pnt, point<T, 3>{0., 0., 0.}), 0., 1e-15);
    ASSERT_NEAR(distance(m.vertex(coedge_end(m, w.coedges[0])).pnt, point<T, 3>{2., 0., 0.}), 0., 1e-15);
}

TEST(tests_brep_face, cylinder_seam)
{
    for (auto srf : {cylinder(), nurbs_cylinder()})
    {
        Model<T> m;
        auto f = unwrap(make_face(m, srf));
        auto c = check_face(m, f);
        ASSERT_EQ(c.vertices, 2);
        ASSERT_EQ(c.edges, 3); // two circles + the seam
        ASSERT_EQ(c.seams, 1);
        ASSERT_EQ(c.degenerate, 0);
        ASSERT_EQ(c.free_edges, 2);
        ASSERT_FALSE(c.closed_shell);
        // the circles are closed edges
        size_t closed_edges = 0;
        for (auto eid : explore<EdgeId>(m, f))
            closed_edges += m.edge(eid).v1 == m.edge(eid).v2;
        ASSERT_EQ(closed_edges, 2);
    }
}

TEST(tests_brep_face, sphere_poles)
{
    Model<T> m;
    auto f = unwrap(make_face(m, sphere()));
    auto c = check_face(m, f);
    ASSERT_EQ(c.vertices, 2);   // the poles
    ASSERT_EQ(c.edges, 3);      // two degenerate edges + the seam
    ASSERT_EQ(c.degenerate, 2);
    ASSERT_EQ(c.seams, 1);
    ASSERT_EQ(c.free_edges, 0);
    ASSERT_TRUE(c.closed_shell); // a sphere is closed by itself
    for (auto eid : explore<EdgeId>(m, f))
        if (m.edge(eid).degenerate)
            ASSERT_FALSE(m.edge(eid).curve);
}

TEST(tests_brep_face, torus_two_seams)
{
    Model<T> m;
    auto f = unwrap(make_face(m, torus()));
    auto c = check_face(m, f);
    ASSERT_EQ(c.vertices, 1);
    ASSERT_EQ(c.edges, 2);
    ASSERT_EQ(c.seams, 2);
    ASSERT_EQ(c.degenerate, 0);
    ASSERT_TRUE(c.closed_shell);
}

TEST(tests_brep_face, cone_apex)
{
    Model<T> m;
    auto f = unwrap(make_face(m, cone()));
    auto c = check_face(m, f);
    ASSERT_EQ(c.vertices, 2); // apex + seam end
    ASSERT_EQ(c.edges, 3);    // apex (degenerate), seam, base circle
    ASSERT_EQ(c.degenerate, 1);
    ASSERT_EQ(c.seams, 1);
    ASSERT_EQ(c.free_edges, 1);
}

TEST(tests_brep_face, revolution)
{
    {
        Model<T> m;
        auto f = unwrap(make_face(m, revolution(2. * pi)));
        auto c = check_face(m, f);
        ASSERT_EQ(c.vertices, 2);
        ASSERT_EQ(c.edges, 3);
        ASSERT_EQ(c.seams, 1);
        ASSERT_EQ(c.free_edges, 2);
    }
    {
        Model<T> m;
        auto f = unwrap(make_face(m, revolution(pi)));
        auto c = check_face(m, f);
        ASSERT_EQ(c.vertices, 4);
        ASSERT_EQ(c.edges, 4);
        ASSERT_EQ(c.seams, 0);
        ASSERT_EQ(c.free_edges, 4);
    }
}

TEST(tests_brep_face, shared_surface_and_errors)
{
    Model<T> m;
    auto srf = plane();
    auto f = unwrap(make_face(m, srf));
    ASSERT_TRUE(m.face(f).surface == srf); // the surface is shared, not copied

    const auto n_vtx = m.count<ShapeType::Vertex>();
    ASSERT_TRUE(make_face(m, std::shared_ptr<Surface<T, 3>>{}).error().code == BuildErrc::NullSurface);
    ASSERT_TRUE(make_face(m, srf, 0.).error().code == BuildErrc::InvalidTolerance);
    auto collapsed = std::make_shared<Analytic>([](T, T) { return point<T, 3>{1., 2., 3.}; }, std::array<T, 4>{0., 1., 0., 1.});
    ASSERT_TRUE(make_face(m, collapsed).error().code == BuildErrc::DegenerateSurface);
    auto infinite = std::make_shared<Analytic>([](T u, T v) { return point<T, 3>{u, v, 0.}; },
                                               std::array<T, 4>{0., std::numeric_limits<T>::infinity(), 0., 1.});
    ASSERT_TRUE(make_face(m, infinite).error().code == BuildErrc::UnboundedSurface);
    auto empty = std::make_shared<Analytic>([](T u, T v) { return point<T, 3>{u, v, 0.}; }, std::array<T, 4>{1., 1., 0., 1.});
    ASSERT_TRUE(make_face(m, empty).error().code == BuildErrc::UnboundedSurface);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), n_vtx); // failures left the model unchanged
    ASSERT_EQ(m.count<ShapeType::Face>(), 1);
}

TEST(tests_brep_face, uv_signed_area_sign)
{
    Model<T> m;
    auto f = unwrap(make_face(m, plane()));
    auto wid = m.face(f).wires.front();
    ASSERT_NEAR(uv_signed_area(m, wid), 1., 1e-12);
    // the same co-edges reversed (order and senses) enclose a negative area: a hole
    auto w = m.wire(wid);
    std::ranges::reverse(w.coedges);
    for (auto &ce : w.coedges)
        ce.orient = reverse(ce.orient);
    auto hole = m.add(std::move(w));
    ASSERT_TRUE(is_closed(m, hole));
    ASSERT_NEAR(uv_signed_area(m, hole), -1., 1e-12);
}
