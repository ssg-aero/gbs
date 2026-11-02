// Faces of one domain unified (gbs-brep/unify.h), outside the booleans: a pair that cannot be
// put on its plane, primitives on boxes shorter than a turn, pieces whose boundary starts at a
// pole, a shared edge that another face uses, the tolerances of a plane in B-spline form, the
// entities left unused.
#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include <gbs/curveline.h>

#include <algorithm>
#include <map>
#include <memory>
#include <numbers>
#include <vector>

using namespace gbs;
using namespace gbs::brep;
using T = double;

namespace
{
    const T pi = std::numbers::pi;
    const Frame3<T> xyz{ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}};

    FaceId natural(Model<T> &m, std::shared_ptr<Surface<T, 3>> s, T tol = brep_default_tolerance<T>)
    {
        auto f = make_face(m, std::move(s), tol);
        REQUIRE(f.has_value());
        return *f;
    }

    /// The faces sewn into one shell, then unified: the faces left.
    std::vector<FaceUse> sewn_unified(Model<T> &m, const std::vector<FaceId> &faces, UnifyReport *report = nullptr)
    {
        const auto sewn = sew(m, faces);
        REQUIRE(sewn.has_value());
        REQUIRE(sewn->shells.size() == 1u);
        auto uses = m.shell(sewn->shells.front()).faces;
        const auto r = unify_faces(m, uses, T(1e-6));
        if (report)
            *report = r;
        return uses;
    }

    /// Whether the boundary of a face runs an edge twice (a seam).
    bool has_seam(const Model<T> &m, FaceId f)
    {
        std::map<std::uint32_t, int> runs;
        for (auto w : m.face(f).wires)
            for (const auto &ce : m.wire(w).coedges)
                if (++runs[ce.edge.index] == 2)
                    return true;
        return false;
    }

    auto plane(const std::array<T, 4> &box) { return std::make_shared<Plane<T>>(xyz, box); }
}

// A plane in B-spline form next to a coplanar face, one of its pcurves a 2D line (not a
// B-spline): it cannot be put on its plane, the pair is left as it is and the unification ends.
TEST(tests_brep_unify, pair_that_cannot_be_put_on_its_plane)
{
    Model<T> m;
    // z = 0, x, y in [0, 2], its parameters in [0, 4] (the map to the plane is not the identity)
    auto s = std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.}},
                                               std::vector<T>{0., 0., 4., 4.}, std::vector<T>{0., 0., 4., 4.}, 1, 1);
    const auto a = natural(m, s);
    // the bottom side's pcurve (u, 0) as a 2D line: the same values
    m.wire(m.face(a).wires[0]).coedges[0].pcurve = std::make_shared<Line<T, 2>>(point<T, 2>{0., 0.}, point<T, 2>{1., 0.});
    const auto b = natural(m, plane({2., 3., 0., 2.}));
    UnifyReport r;
    const auto faces = sewn_unified(m, {a, b}, &r);
    ASSERT_EQ(faces.size(), 2u);
    ASSERT_EQ(r.faces_merged, 0u);
}

// Faces on cylinders whose boxes are shorter than a turn: the angle is periodic all the same.
// Two pieces covering 1.9 pi become one face; two halves closing a ring become one face with a
// seam, however tall (one tolerance per parameter: an angle is not compared with a length).
TEST(tests_brep_unify, primitives_on_partial_boxes)
{
    auto cylinder = [](T h) { return CylindricalSurface<T>(1., xyz, {0., 2. * pi, 0., h}); };
    {
        Model<T> m;
        const auto c = cylinder(2.);
        const auto faces = sewn_unified(m, {natural(m, c.with_box({0., pi, 0., 2.})), natural(m, c.with_box({pi, 1.9 * pi, 0., 2.}))});
        ASSERT_EQ(faces.size(), 1u);
        ASSERT_FALSE(has_seam(m, faces.front().face));
        ASSERT_TRUE(check(m, ShapeId{faces.front().face}).ok());
    }
    for (T h : {10., 1e5})
    {
        Model<T> m;
        const auto c = cylinder(h);
        const auto faces = sewn_unified(m, {natural(m, c.with_box({0., pi, 0., h})), natural(m, c.with_box({pi, 2. * pi, 0., h}))});
        ASSERT_EQ(faces.size(), 1u);
        ASSERT_TRUE(has_seam(m, faces.front().face));
        ASSERT_TRUE(check(m, ShapeId{faces.front().face}).ok());
    }
}

// Two pieces of a cone (of a sphere), the boundary of the second starting at its apex (at its
// south pole), where the surface has no normal and no metric: the orientation and the first map
// are taken at a regular point, the samples there are only checked; they merge.
TEST(tests_brep_unify, boundary_starting_at_an_apex_or_a_pole)
{
    const ConicalSurface<T> cone(0., pi / 6., xyz, {0., 2. * pi, 0., 1.}); // (apex at v = 0)
    const SphericalSurface<T> sphere(1., xyz, {0., 2. * pi, -pi / 2., pi / 2.});
    for (const ElementarySurface<T> *s : {static_cast<const ElementarySurface<T> *>(&cone), static_cast<const ElementarySurface<T> *>(&sphere)})
    {
        Model<T> m;
        const auto [u1, u2, v1, v2] = s->bounds();
        const auto a = natural(m, s->with_box({0., pi, v1, v2}));
        const auto b = natural(m, s->with_box({pi, 1.5 * pi, v1, v2}));
        REQUIRE(m.edge(m.wire(m.face(b).wires[0]).coedges[0].edge).degenerate);
        const auto faces = sewn_unified(m, {a, b});
        ASSERT_EQ(faces.size(), 1u);
        ASSERT_TRUE(check(m, ShapeId{faces.front().face}).ok());
    }
}

// Two coplanar faces whose shared edge a third face also uses: they are not merged (the third
// face would keep an edge that no longer bounds them).
TEST(tests_brep_unify, shared_edge_another_face_uses)
{
    Model<T> m;
    const auto a = natural(m, plane({0., 1., 0., 1.})), b = natural(m, plane({1., 2., 0., 1.}));
    const auto sewn = sew(m, std::vector<FaceId>{a, b});
    REQUIRE(sewn.has_value());
    REQUIRE(sewn->merged.size() == 1u);
    const auto shared = sewn->merged.front().first;
    // the face x = 1, y and z in [0, 1], its bottom side made that edge
    const auto c = natural(m, std::make_shared<Plane<T>>(Frame3<T>{ax2<T, 3>{{{1., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}}}, std::array<T, 4>{0., 1., 0., 1.}));
    m.wire(m.face(c).wires[0]).coedges[0].edge = shared;
    ASSERT_FALSE(merge_faces(m, FaceUse{a, Orientation::Forward}, FaceUse{b, Orientation::Forward}, T(1e-6)).has_value());
    ASSERT_TRUE(m.alive(a));
    ASSERT_TRUE(m.alive(b));
}

// The planar part of a B-spline put on its plane: the face and the edge that bulge off the plane
// (within tol) get tolerances that cover the distance, their vertices too.
TEST(tests_brep_unify, plane_in_bspline_form_tolerances)
{
    Model<T> m;
    // z = 0 over [0, 2] x [0, 2], the middle pole of the far side lifted: that side bulges by 5e-7
    points_vector<T, 3> p;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            p.push_back({T(i), T(j), i == 1 && j == 2 ? 1e-6 : 0.});
    auto s = std::make_shared<BSSurface<T, 3>>(p, std::vector<T>{0., 0., 0., 1., 1., 1.}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2, 2);
    UnifyReport r;
    const auto faces = sewn_unified(m, {natural(m, s, 1e-7), natural(m, plane({2., 3., 0., 2.}), 1e-7)}, &r);
    ASSERT_EQ(r.faces_planed, 1u);
    ASSERT_EQ(faces.size(), 1u);
    const auto f = faces.front().face;
    ASSERT_GE(m.face(f).tol, 4.9e-7);
    T worst{};
    for (auto e : explore<EdgeId>(m, ShapeId{f}))
    {
        const auto &ed = m.edge(e);
        worst = std::max(worst, ed.tol);
        ASSERT_GE(m.vertex(ed.v1).tol, ed.tol);
        ASSERT_GE(m.vertex(ed.v2).tol, ed.tol);
    }
    ASSERT_GE(worst, 4.9e-7);
    ASSERT_TRUE(check(m, ShapeId{f}).ok());
}

// The edge two merged faces shared, and nothing else, is erased with the unification.
TEST(tests_brep_unify, unused_entities_erased)
{
    Model<T> m;
    const auto faces = sewn_unified(m, {natural(m, plane({0., 1., 0., 1.})), natural(m, plane({1., 2., 0., 1.}))});
    ASSERT_EQ(faces.size(), 1u);
    const auto f = faces.front().face;
    ASSERT_EQ(m.count(ShapeType::Edge), explore<EdgeId>(m, ShapeId{f}).size());
    ASSERT_EQ(m.count(ShapeType::Vertex), explore<VertexId>(m, ShapeId{f}).size());
}
