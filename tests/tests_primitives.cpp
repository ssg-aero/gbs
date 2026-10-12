#include <doctest_gtest.hpp>
#include <gbs/primitives.h>
#include <gbs-brep/brep>

#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi;

    // tilted frame, reference direction not orthogonal to the axis
    const ax2<T, 3> tilted{{{1., -2., 0.5}, {1., 1., 1.}, {1., 0., 0.}}};

    T dist(const point<T, 3> &a, const point<T, 3> &b) { return norm(a - b); }

    // primitives on boxes crossing the angle 0 where they have angles
    std::vector<std::shared_ptr<ElementarySurface<T>>> surfaces()
    {
        return {
            std::make_shared<Plane<T>>(tilted, std::array<T, 4>{-2., 3., -1., 4.}),
            std::make_shared<CylindricalSurface<T>>(1.5, tilted, std::array<T, 4>{-0.5, 2., -1., 3.}),
            std::make_shared<ConicalSurface<T>>(1., pi / 6, tilted, std::array<T, 4>{0., 2 * pi, -1.5, 2.}),
            std::make_shared<SphericalSurface<T>>(2., tilted, std::array<T, 4>{-1., 4., -pi / 2, 1.}),
            std::make_shared<ToroidalSurface<T>>(3., 0.8, tilted, std::array<T, 4>{5., 7., -1., 4.}),
        };
    }

    template <typename F>
    void grid(const std::array<T, 4> &b, int n, F f)
    {
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                f(b[0] + (b[1] - b[0]) * i / n, b[2] + (b[3] - b[2]) * j / n);
    }
}

TEST(tests_primitives, surfaces_are_their_exact_nurbs)
{
    for (const auto &s : surfaces())
    {
        const auto n = s->to_nurbs();
        grid(s->bounds(), 12, [&](T u, T v) {
            const auto [nu, nv] = s->nurbs_parameters(u, v);
            ASSERT_LT(dist(n->value(nu, nv), s->value(u, v)), 1e-11);
        });
    }
}

TEST(tests_primitives, surface_formulas)
{
    const Frame3<T> f{tilted};
    auto rho = [&](T u) { return std::cos(u) * f.x + std::sin(u) * f.y; };
    const auto ss = surfaces();
    grid({0., 2., -1., 1.}, 6, [&](T u, T v) {
        ASSERT_LT(dist(ss[0]->value(u, v), f.o + u * f.x + v * f.y), 1e-14);
        ASSERT_LT(dist(ss[1]->value(u, v), f.o + 1.5 * rho(u) + v * f.z), 1e-14);
        ASSERT_LT(dist(ss[2]->value(u, v), f.o + (1. + v * std::tan(pi / 6)) * rho(u) + v * f.z), 1e-14);
        ASSERT_LT(dist(ss[3]->value(u, v), f.o + (2. * std::cos(v)) * rho(u) + (2. * std::sin(v)) * f.z), 1e-14);
        ASSERT_LT(dist(ss[4]->value(u, v), f.o + (3. + 0.8 * std::cos(v)) * rho(u) + (0.8 * std::sin(v)) * f.z), 1e-14);
    });
    // orthonormal right-handed frame, normal of the cylinder outwards
    ASSERT_NEAR(f.x * f.z, 0., 1e-15);
    ASSERT_NEAR(norm(cross(f.x, f.y) - f.z), 0., 1e-15);
    const auto n = cross(ss[1]->value(0.3, 0.5, 1, 0), ss[1]->value(0.3, 0.5, 0, 1));
    ASSERT_GT(n * rho(0.3), 0.);
}

TEST(tests_primitives, derivatives)
{
    const T h = 1e-4;
    for (const auto &s : surfaces())
        grid({0.2, 1.4, -0.6, 0.7}, 3, [&](T u, T v) {
            for (std::size_t du = 0; du <= 2; ++du)
                for (std::size_t dv = 0; du + dv <= 2; ++dv)
                {
                    // central difference of the derivative of lower order
                    point<T, 3> fd;
                    if (du > 0)
                        fd = (1. / (2 * h)) * (s->value(u + h, v, du - 1, dv) - s->value(u - h, v, du - 1, dv));
                    else if (dv > 0)
                        fd = (1. / (2 * h)) * (s->value(u, v + h, du, dv - 1) - s->value(u, v - h, du, dv - 1));
                    else
                        continue;
                    ASSERT_LT(dist(s->value(u, v, du, dv), fd), 1e-6);
                }
            ASSERT_LT(norm(s->value(u, v, 3, 1) - (1. / (2 * h)) * (s->value(u + h, v, 2, 1) - s->value(u - h, v, 2, 1))), 1e-6);
        });
}

TEST(tests_primitives, inverse_and_iso_lines)
{
    for (const auto &s : surfaces())
    {
        const auto b = s->bounds();
        grid(b, 8, [&](T u, T v) {
            const auto p = s->value(u, v);
            const auto uv = s->parameters(p);
            ASSERT_LT(dist(s->value(uv[0], uv[1]), p), 1e-12);
            if (std::abs(std::cos(v)) > 1e-6 || dynamic_cast<SphericalSurface<T> *>(s.get()) == nullptr)
            {
                T d = uv[0] - u; // an angle is found modulo 2 pi (u = 2 pi comes back as 0 on a full turn)
                if (dynamic_cast<Plane<T> *>(s.get()) == nullptr)
                    d -= 2 * pi * std::round(d / (2 * pi));
                ASSERT_NEAR(d, 0., 1e-9 * (1. + std::abs(u)));
                ASSERT_GE(uv[0], b[0] - 1e-12); // in the turn of the box
                ASSERT_LE(uv[0], std::max(b[1], b[0] + 2 * pi) + 1e-12);
            }
            const auto iso = s->iso_u(u, b[2], b[3]);
            ASSERT_LT(dist(iso->value(v), p), 1e-13);
        });
    }
}

TEST(tests_primitives, curves)
{
    const Frame3<T> f{tilted};
    // ellipse arc across the angle 0, parametrized by the eccentric angle
    const EllipseCurve<T> e{3., 1.2, tilted, 5.9, 6.6};
    for (T t : {5.9, 6.2, 6.5, 6.6})
    {
        ASSERT_LT(dist(e.value(t), f.o + (3. * std::cos(t)) * f.x + (1.2 * std::sin(t)) * f.y), 1e-14);
        ASSERT_NEAR(e.parameter(e.value(t)), t, 1e-12);
        ASSERT_LT(dist(e.value(t, 1), (1. / 2e-6) * (e.value(t + 1e-6) - e.value(t - 1e-6))), 1e-6);
        ASSERT_LT(dist(e.to_nurbs()->value(e.nurbs_parameter(t)), e.value(t)), 1e-12);
    }
    ASSERT_THROW(EllipseCurve<T>(1., 1., tilted, 0., 7.), std::invalid_argument);
    ASSERT_THROW(EllipseCurve<T>(0., 1., tilted), std::invalid_argument);
    // reversed: same points, parameter t -> -t, still right-handed
    const auto r = e.reversed();
    ASSERT_EQ(r->bounds()[0], -6.6);
    for (T t : {5.9, 6.3, 6.6})
        ASSERT_LT(dist(r->value(-t), e.value(t)), 1e-14);
    const auto &g = std::dynamic_pointer_cast<EllipseCurve<T>>(r)->frame();
    ASSERT_NEAR(norm(cross(g.x, g.y) - g.z), 0., 1e-15);
    // line with a vector of magnitude 2: the STEP parametrization
    const LineCurve<T> l{{1., 2., 3.}, {0., 2., 0.}, -1., 4.};
    ASSERT_LT(dist(l.value(1.5), point<T, 3>{1., 5., 3.}), 1e-15);
    ASSERT_NEAR(l.parameter({7., 5., 3.}), 1.5, 1e-15);
    ASSERT_LT(dist(l.reversed()->value(-1.5), l.value(1.5)), 1e-15);
    ASSERT_LT(dist(l.to_nurbs()->value(2.), l.value(2.)), 1e-15);
}

TEST(tests_primitives, rigid_moves)
{
    // a quarter turn about z and a translation
    auto move = [](const point<T, 3> &p) { return point<T, 3>{-p[1] + 10., p[0], p[2] + 1.}; };
    for (const auto &s : surfaces())
    {
        const auto m = s->moved(move);
        ASSERT_EQ(m->bounds(), s->bounds());
        grid(s->bounds(), 5, [&](T u, T v) { ASSERT_LT(dist(m->value(u, v), move(s->value(u, v))), 1e-13); });
    }
    const EllipseCurve<T> e{2., 1., tilted, 0., 3.};
    const auto me = e.moved(move);
    ASSERT_LT(dist(me->value(1.7), move(e.value(1.7))), 1e-14);
}

TEST(tests_primitives, invalid_arguments)
{
    ASSERT_THROW(CylindricalSurface<T>(-1., tilted, {0., 1., 0., 1.}), std::invalid_argument);
    ASSERT_THROW(CylindricalSurface<T>(1., tilted, {0., 1., 1., 1.}), std::invalid_argument);
    ASSERT_THROW(Plane<T>(tilted, {0., std::numeric_limits<T>::infinity(), 0., 1.}), std::invalid_argument);
    ASSERT_THROW(SphericalSurface<T>(1., tilted, {0., 1., -2., 0.}), std::invalid_argument);
    ASSERT_THROW(ConicalSurface<T>(1., pi / 2, tilted, {0., 1., 0., 1.}), std::invalid_argument);
    ASSERT_THROW((Plane<T>(ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {0., 0., 3.}}}, {0., 1., 0., 1.})), std::invalid_argument);
}

TEST(tests_primitives, brep_faces_and_solids)
{
    using namespace gbs::brep;
    const ax2<T, 3> z{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    // natural faces: cylinder (seam), sphere (seam, poles), torus (two seams)
    Model<T> m;
    for (std::shared_ptr<Surface<T, 3>> s : std::vector<std::shared_ptr<Surface<T, 3>>>{
             std::make_shared<CylindricalSurface<T>>(1., z, std::array<T, 4>{0., 2 * pi, 0., 2.}),
             std::make_shared<SphericalSurface<T>>(1., z, std::array<T, 4>{0., 2 * pi, -pi / 2, pi / 2}),
             std::make_shared<ToroidalSurface<T>>(3., 1., z, std::array<T, 4>{0., 2 * pi, 0., 2 * pi})})
    {
        const auto f = unwrap(make_face(m, s));
        ASSERT_TRUE(check(m, f).ok());
    }
    // closed cylinder from ordered loops: circle edges are ellipse primitives, pcurves are segments
    const T R = 1.5, h = 2.;
    auto lateral = std::make_shared<CylindricalSurface<T>>(R, z, std::array<T, 4>{0., 2 * pi, 0., h});
    const auto vb = unwrap(make_vertex(m, {R, 0., 0.})), vt = unwrap(make_vertex(m, {R, 0., h}));
    auto circle = [&](T zc, VertexId v) {
        auto c = std::make_shared<EllipseCurve<T>>(R, R, ax2<T, 3>{{{0., 0., zc}, {0., 0., 1.}, {1., 0., 0.}}});
        return unwrap(make_edge(m, c, 0., 2 * pi, v, v));
    };
    const auto bottom = circle(0., vb), top = circle(h, vt);
    const auto seam = unwrap(make_edge(m, std::make_shared<LineCurve<T>>(point<T, 3>{R, 0., 0.}, point<T, 3>{0., 0., 1.}, 0., h), 0., h, vb, vt));
    constexpr auto Fw = Orientation::Forward, Rv = Orientation::Reversed;
    const auto w = unwrap(make_wire_ordered(m, {{bottom, Fw}, {seam, Fw}, {top, Rv}, {seam, Rv}}));
    const auto use = unwrap(make_face_use(m, lateral, w));
    ASSERT_EQ(use.orient, Fw);
    ASSERT_NEAR(uv_signed_area(m, w), 2 * pi * h, 1e-9); // the parameter is the angle
    auto disk = [&](T zc, EdgeId e, Orientation o) {
        auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., zc}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{-2., 2., -2., 2.});
        return unwrap(make_face_use(m, plane, unwrap(make_wire_ordered(m, {{e, o}}))));
    };
    const auto sh = m.add(Shell{{use, disk(0., bottom, Rv), disk(h, top, Fw)}});
    const auto so = unwrap(make_solid(m, sh));
    ASSERT_NEAR(signed_volume(m, sh, 256), pi * R * R * h, 2e-3);
    ASSERT_TRUE(check(m, so).ok());
    // a placed copy keeps primitives
    const auto c = copy_shape(m, so, Placement<T>::frame({5., 0., 0.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.}));
    for (auto f : explore<FaceId>(m, c))
        ASSERT_TRUE(std::dynamic_pointer_cast<ElementarySurface<T>>(m.face(f).surface) != nullptr);
    for (auto e : explore<EdgeId>(m, c))
        ASSERT_TRUE(std::dynamic_pointer_cast<ElementaryCurve<T>>(m.edge(e).curve) != nullptr);
    ASSERT_TRUE(check(m, c).ok());
    ASSERT_NEAR(signed_volume(m, m.solid(std::get<SolidId>(c)).outer, 256), pi * R * R * h, 2e-3);
}
