#include <doctest_gtest.hpp>
#include <gbs-brep/section.h>
#include <gbs-brep/edit.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

#include <cmath>
#include <numbers>
#include <string>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;

namespace
{
    using brep_tests::pi;
    using T = double;
    constexpr T tol = 1e-6;

    std::string file(const std::string &name) { return std::string(GBS_TESTS_IN_DIR) + "/step/" + name; }

    ShapeId read(Model<T> &m, const std::string &name)
    {
        auto r = step::read_step_file(m, file(name));
        REQUIRE(r.has_value());
        return r->root;
    }

    std::shared_ptr<Plane<T>> plane(const point<T, 3> &o, const point<T, 3> &n, T half = 100.)
    {
        const auto z = intersect_detail::unit<T>(n);
        return std::make_shared<Plane<T>>(intersect_detail::frame_of<T>(o, z, intersect_detail::orthogonal<T>(z)), std::array<T, 4>{-half, half, -half, half});
    }

    T length(const Model<T> &m, EdgeId e)
    {
        const auto &ed = m.edge(e);
        T l{};
        constexpr int n = 2000;
        for (int k = 0; k < n; ++k)
            l += norm(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * (k + 1.) / n) - ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(k) / n));
        return l;
    }

    // every edge on both surfaces through its pcurves; wires closed or not, total length
    T check_section(const Model<T> &m, const SectionResult<T> &r, const Surface<T, 3> *tool = nullptr)
    {
        REQUIRE(r.complete);
        T total{};
        for (const auto &se : r.edges)
        {
            const auto &e = m.edge(se.edge);
            const auto &sa = *m.face(se.face_a).surface;
            const auto &sb = tool ? *tool : *m.face(se.face_b).surface;
            for (int k = 0; k <= 10; ++k)
            {
                const T t = e.u1 + (e.u2 - e.u1) * k / 10.;
                const auto p = e.curve->value(t);
                const auto qa = se.pcurve_a->value(t), qb = se.pcurve_b->value(t);
                REQUIRE(norm(sa.value(qa[0], qa[1]) - p) <= e.tol + 1e-12);
                REQUIRE(norm(sb.value(qb[0], qb[1]) - p) <= e.tol + 1e-12);
            }
            REQUIRE(norm(m.vertex(e.v1).pnt - e.curve->value(e.u1)) <= m.vertex(e.v1).tol);
            REQUIRE(norm(m.vertex(e.v2).pnt - e.curve->value(e.u2)) <= m.vertex(e.v2).tol);
            total += length(m, se.edge);
        }
        return total;
    }
}

TEST(tests_brep_section, box_by_planes)
{
    Model<T> m;
    const auto box = read(m, "box.stp"); // 2 x 3 x 4
    const auto vertices = explore<VertexId>(m, box).size();
    const auto cut = plane({0., 0., 1.}, {0., 0., 1.});
    auto r = unwrap(section(m, box, cut));
    ASSERT_EQ(r.edges.size(), 4u);
    ASSERT_NEAR(check_section(m, r, cut.get()), 10., 1e-9);
    ASSERT_EQ(r.wires.size(), 1u);
    ASSERT_TRUE(m.wire(r.wires[0]).closed);
    ASSERT_EQ(explore<VertexId>(m, box).size(), vertices); // the box is not modified

    // through two vertical edges: the diagonals of the top and the bottom, ending on the
    // vertices of the box, and the two vertical edges themselves: a closed rectangle
    auto d = unwrap(section(m, box, plane({0., 0., 0.}, {3., -2., 0.})));
    ASSERT_EQ(d.edges.size(), 2u);
    ASSERT_EQ(d.boundary_edges.size(), 2u);
    ASSERT_EQ(d.wires.size(), 1u);
    ASSERT_TRUE(m.wire(d.wires[0]).closed);
    const auto box_vertices = explore<VertexId>(m, box);
    for (const auto &se : d.edges)
    {
        const auto &e = m.edge(se.edge);
        ASSERT_TRUE(std::ranges::find(box_vertices, e.v1) != box_vertices.end());
        ASSERT_TRUE(std::ranges::find(box_vertices, e.v2) != box_vertices.end());
    }
    // on the top face: coincident; the section is the boundary of the top face
    auto top = unwrap(section(m, box, plane({0., 0., 4.}, {0., 0., 1.})));
    ASSERT_TRUE(top.edges.empty());
    ASSERT_EQ(top.coincident.size(), 1u);
    ASSERT_EQ(top.boundary_edges.size(), 4u);
    ASSERT_EQ(top.wires.size(), 1u);
    // away from the box: empty
    auto far = unwrap(section(m, box, plane({0., 0., 9.}, {0., 0., 1.})));
    ASSERT_TRUE(far.edges.empty());
}

TEST(tests_brep_section, cylinder_and_torus)
{
    Model<T> m;
    const auto cyl = read(m, "cylinder.stp"); // R 5, h 10
    // across: a closed circle
    const auto p1 = plane({0., 0., 5.}, {0., 0., 1.});
    auto c = unwrap(section(m, cyl, p1));

    ASSERT_EQ(c.edges.size(), 1u);
    ASSERT_NEAR(check_section(m, c, p1.get()), 10 * pi, 1e-4); // lengths measured by 2000 chords
    ASSERT_EQ(c.wires.size(), 1u);
    ASSERT_TRUE(m.wire(c.wires[0]).closed);
    // oblique: an ellipse, cut at the seam of the cylinder, one closed wire
    const T a = pi / 6;
    const auto p2 = plane({0., 0., 5.}, {0., std::sin(a), std::cos(a)});
    auto e = unwrap(section(m, cyl, p2));
    const T A = 5. / std::cos(a), B = 5., h = (A - B) * (A - B) / ((A + B) * (A + B));
    ASSERT_NEAR(check_section(m, e, p2.get()), pi * (A + B) * (1 + 3 * h / (10 + std::sqrt(4 - 3 * h))), 1e-4);
    ASSERT_EQ(e.wires.size(), 1u);
    ASSERT_TRUE(m.wire(e.wires[0]).closed);
    // along the axis: a generatrix, the two diameters of the caps, and the seam edge of the
    // cylinder: a closed rectangle
    const auto p3 = plane({0., 0., 0.}, {0., 1., 0.});
    auto g = unwrap(section(m, cyl, p3));
    ASSERT_NEAR(check_section(m, g, p3.get()), 10. + 2 * 10., 1e-9);
    ASSERT_EQ(g.boundary_edges.size(), 1u);
    ASSERT_EQ(g.wires.size(), 1u);
    ASSERT_TRUE(m.wire(g.wires[0]).closed);

    // the partial torus of OpenCascade, by an oblique plane (marching)
    Model<T> mt;
    const auto torus = read(mt, "occt_partial_torus.stp");
    auto bb = bounding_box(mt, torus);
    const auto mid = (bb.min + bb.max) / 2.;
    const auto p4 = plane(mid, {0.2, 0.1, 1.});
    auto t = unwrap(section(mt, torus, p4));
    ASSERT_FALSE(t.edges.empty());
    check_section(mt, t, p4.get());
    for (auto w : t.wires)
        ASSERT_TRUE(mt.wire(w).closed);
}

TEST(tests_brep_section, two_solids)
{
    using namespace brep_tests;
    // boxes 0..2 and 1..3: the section is the skew hexagon of the edges of the cube 1..2 that join
    // a face of each box, length 6
    Model<T> m;
    auto a = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front()));
    auto b = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {1., 1., 1.}, 2.))).shells.front()));
    auto r = unwrap(section(m, a, b));
    ASSERT_EQ(r.edges.size(), 6u);
    ASSERT_NEAR(check_section(m, r), 6., 1e-9);
    ASSERT_EQ(r.wires.size(), 1u);
    ASSERT_TRUE(m.wire(r.wires[0]).closed);

    // a sphere R 3 centred on the side of the cylinder R 5, turned so that its poles are off the
    // cylinder: one closed loop on the lateral face, across its seam (marching)
    Model<T> mc;
    const auto cyl = read(mc, "cylinder.stp");
    const auto sph = read(mc, "sphere.stp");
    const auto moved = copy_shape(mc, sph, Placement<T>{.r = {{{1., 0., 0.}, {0., 0., -1.}, {0., 1., 0.}}}, .t = {5., 0., 5.}});
    auto s = unwrap(section(mc, cyl, moved));
    ASSERT_FALSE(s.edges.empty());
    check_section(mc, s);
    ASSERT_EQ(s.wires.size(), 1u);
    ASSERT_TRUE(mc.wire(s.wires[0]).closed);

    // errors
    ASSERT_FALSE(section(mc, cyl, ShapeId{SolidId{999}}).has_value());
    ASSERT_FALSE(section(mc, cyl, std::shared_ptr<Surface<T, 3>>{}).has_value());
}
