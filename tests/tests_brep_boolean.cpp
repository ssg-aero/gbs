#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/check.h>
#include <gbs-brep/sweep.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <sstream>
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

    std::string file(const std::string &name) { return std::string(GBS_TESTS_IN_DIR) + "/step/" + name; }

    ShapeId read(Model<T> &m, const std::string &name)
    {
        auto r = step::read_step_file(m, file(name));
        REQUIRE(r.has_value());
        return r->root;
    }

    void expect_valid(const Model<T> &m, const ShapeId &s)
    {
        const auto rep = check(m, s);
        std::ostringstream os;
        for (const auto &e : rep.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        INFO(os.str());
        REQUIRE(rep.ok());
    }

    // volume by the sampling of the faces (signed_volume), a few 1e-4 relative on trimmed faces
    T volume(const Model<T> &m, const ShapeId &s)
    {
        T v{};
        for (auto so : explore<SolidId>(m, s))
        {
            v += signed_volume(m, m.solid(so).outer, 160);
            for (auto c : m.solid(so).voids)
                v += signed_volume(m, c, 160);
        }
        return v;
    }

    SolidId box(Model<T> &m, const point<T, 3> &o, T size)
    {
        using namespace brep_tests;
        return unwrap(make_solid(m, unwrap(sew(m, box_faces(m, o, size))).shells.front()));
    }

    void expect_volume(const Model<T> &m, const BooleanResult<T> &r, T v, std::size_t solids = 1)
    {
        expect_valid(m, r.shape);
        REQUIRE(r.solids.size() == solids);
        REQUIRE(std::abs(volume(m, r.shape) - v) <= 2e-3 * std::max(T(1), std::abs(v)));
    }
}

TEST(tests_brep_boolean, boxes)
{
    // boxes 0..2 and 1..3: overlap 1..2
    Model<T> m;
    const auto a = box(m, {0., 0., 0.}, 2.), b = box(m, {1., 1., 1.}, 2.);
    const auto faces_a = explore<FaceId>(m, ShapeId{a}).size();
    auto f = unwrap(fuse(m, a, b));
    expect_volume(m, f, 8. + 8. - 1.);
    ASSERT_EQ(explore<FaceId>(m, f.shape).size(), 12u);
    auto c = unwrap(common(m, a, b));
    expect_volume(m, c, 1.);
    ASSERT_EQ(explore<FaceId>(m, c.shape).size(), 6u);
    auto d = unwrap(cut(m, a, b));
    expect_volume(m, d, 7.);
    ASSERT_EQ(explore<FaceId>(m, d.shape).size(), 9u);
    // the operands are not modified
    ASSERT_EQ(explore<FaceId>(m, ShapeId{a}).size(), faces_a);
    expect_valid(m, ShapeId{a});
    expect_valid(m, ShapeId{b});
    // disjoint boxes: two solids for the fuse, nothing in common
    const auto far = box(m, {5., 5., 5.}, 1.);
    expect_volume(m, unwrap(fuse(m, a, far)), 9., 2);
    ASSERT_TRUE(unwrap(common(m, a, far)).solids.empty());
    // a box inside another: cut makes a cavity
    const auto inner = box(m, {0.5, 0.5, 0.5}, 1.);
    auto hole = unwrap(cut(m, a, inner));
    expect_volume(m, hole, 7.);
    ASSERT_EQ(m.solid(hole.solids.front()).voids.size(), 1u);
}

TEST(tests_brep_boolean, cylinders_and_spheres)
{
    // a box 0..4 (bilinear B-spline faces) and the sphere R 3 centred at its corner (4, 4, 4),
    // axis along (1, 1, 1) so that its poles and seam stay off the edges of the box (marching)
    Model<T> m;
    const auto bx = ShapeId{box(m, {0., 0., 0.}, 4.)};
    const auto sph = read(m, "sphere.stp"); // R 3 at the origin
    const point<T, 3> z{1. / std::sqrt(3.), 1. / std::sqrt(3.), 1. / std::sqrt(3.)}, x{1. / std::sqrt(2.), -1. / std::sqrt(2.), 0.};
    const auto s = copy_shape(m, sph, Placement<T>::frame({4., 4., 4.}, x, cross(z, x), z));
    // the eighth of the sphere at the corner lies in the box
    const T eighth = 4. / 3. * pi * 27. / 8.;
    auto c = unwrap(common(m, bx, s));
    expect_volume(m, c, eighth);
    auto d = unwrap(cut(m, bx, s));
    expect_volume(m, d, 64. - eighth);
    auto f = unwrap(fuse(m, bx, s));
    expect_volume(m, f, 64. + 36. * pi - eighth);

    // the cylinder R 5 h 10 minus a box inside it: a cavity
    const auto cyl = read(m, "cylinder.stp");
    const auto slab = copy_shape(m, read(m, "box.stp"), Placement<T>{.t = {-1., -1.5, 3.}}); // -1..1 x -1.5..1.5 x 3..7
    auto r = unwrap(cut(m, cyl, slab));
    expect_volume(m, r, 250. * pi - 24.);
    ASSERT_EQ(m.solid(r.solids.front()).voids.size(), 1u);
}

TEST(tests_brep_boolean, coincident_faces)
{
    Model<T> m;
    const auto a = read(m, "box.stp"); // 0..2 x 0..3 x 0..4 (planes)
    // shifted in x and y: the bottoms and the tops are coplanar
    const auto b = copy_shape(m, a, Placement<T>{.t = {1., 1., 0.}});
    expect_volume(m, unwrap(fuse(m, a, b)), 24. + 24. - 8.);
    expect_volume(m, unwrap(common(m, a, b)), 8.);
    expect_volume(m, unwrap(cut(m, a, b)), 16.);
    expect_valid(m, a);
    expect_valid(m, b);
    // stacked, face against face: one solid, nothing in common, the cut keeps the first box
    const auto up = copy_shape(m, a, Placement<T>{.t = {0., 0., 4.}});
    auto f = unwrap(fuse(m, a, up));
    expect_volume(m, f, 48.);
    ASSERT_EQ(explore<FaceId>(m, f.shape).size(), 10u); // the touching faces are gone, the sides split in two
    ASSERT_TRUE(unwrap(common(m, a, up)).solids.empty());
    expect_volume(m, unwrap(cut(m, a, up)), 24.);
    // stacked and shifted: an L-shaped step
    const auto step = copy_shape(m, a, Placement<T>{.t = {1., 0., 4.}});
    expect_volume(m, unwrap(fuse(m, a, step)), 48.);
    // the same box twice
    const auto same = copy_shape(m, a, Placement<T>{});
    auto s = unwrap(fuse(m, a, same));
    expect_volume(m, s, 24.);
    ASSERT_EQ(explore<FaceId>(m, s.shape).size(), 6u);
    expect_volume(m, unwrap(common(m, a, same)), 24.);
    ASSERT_TRUE(unwrap(cut(m, a, same)).solids.empty());

    // coaxial cylinders R 5, z 0..10 and 5..15: the lateral faces coincide on 5..10
    const auto c1 = read(m, "cylinder.stp");
    const auto c2 = copy_shape(m, c1, Placement<T>{.t = {0., 0., 5.}});
    expect_volume(m, unwrap(fuse(m, c1, c2)), 25. * pi * 15.);
    expect_volume(m, unwrap(common(m, c1, c2)), 25. * pi * 5.);
    expect_volume(m, unwrap(cut(m, c1, c2)), 25. * pi * 5.);
}

TEST(tests_brep_boolean, non_generic)
{
    // a box turned by 45 degrees about y, resting on an edge on the top of a box with B-spline
    // faces: the edge lies on a face that is not a primitive
    Model<T> m;
    const auto bx = read(m, "box.stp");
    const auto cube = ShapeId{box(m, {-1., -1., 0.}, 4.)};
    const T c = std::sqrt(0.5);
    const auto s = copy_shape(m, bx, Placement<T>::frame({1., 2.5, 4.}, {c, 0., c}, {-c, 0., c}, {0., -1., 0.}));
    const auto edges = m.count<ShapeType::Edge>(), vertices = m.count<ShapeType::Vertex>();
    auto r = fuse(m, cube, s);
    if (r)
        FAIL("expected NonGeneric");
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    // nothing is left of the attempt
    ASSERT_EQ(m.count<ShapeType::Edge>(), edges);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), vertices);
    // the operands are kept
    expect_valid(m, cube);
    expect_valid(m, s);
    // errors
    ASSERT_FALSE(fuse(m, bx, ShapeId{SolidId{999}}).has_value());
    const auto face = explore<FaceId>(m, bx).front();
    ASSERT_FALSE(fuse(m, bx, ShapeId{face}).has_value());
}

TEST(tests_brep_boolean, tangent_contacts)
{
    Model<T> m;
    const auto bx = read(m, "box.stp"); // 0..2 x 0..3 x 0..4
    const auto sph = read(m, "sphere.stp"); // R 3
    const auto cyl = read(m, "cylinder.stp"); // R 5, z 0..10
    const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};

    // a sphere touching the top of the box at one point
    {
        const auto s = copy_shape(m, sph, Placement<T>::frame({1., 1.5, 7.}, X, cross(Y, X), Y));
        const auto edges = m.count<ShapeType::Edge>();
        const auto f = unwrap(fuse(m, bx, s));
        expect_volume(m, f, 24. + 36. * pi, 2);
        // the operation leaves only the edges of its result
        ASSERT_EQ(m.count<ShapeType::Edge>(), edges + explore<EdgeId>(m, f.shape).size());
        ASSERT_TRUE(unwrap(common(m, bx, s)).solids.empty());
        expect_volume(m, unwrap(cut(m, bx, s)), 24.);
    }
    // the cylinder lying on the top of the box, tangent along y at x = 1, z = 4
    {
        const auto c = copy_shape(m, cyl, Placement<T>::frame({1., -3.5, 9.}, Z, cross(Y, Z), Y));
        expect_volume(m, unwrap(fuse(m, bx, c)), 24. + 250. * pi, 2);
        ASSERT_TRUE(unwrap(common(m, bx, c)).solids.empty());
        expect_volume(m, unwrap(cut(m, bx, c)), 24.);
    }
    // the box 3..5 x -1.5..1.5 x 3..7 against the cylinder: its face x = 5 is tangent to the
    // cylinder (turned so that the line of contact is off the seam); the cut leaves two slivers
    // touching along the line of contact
    {
        const T a = 1.;
        const auto c = copy_shape(m, cyl, Placement<T>::frame({0., 0., 0.}, {std::cos(a), std::sin(a), 0.}, {-std::sin(a), std::cos(a), 0.}, Z));
        const auto b = copy_shape(m, bx, Placement<T>{.t = {3., -1.5, 3.}});
        const T y = 1.5, sliver = 4. * (5. * y - (y / 2. * std::sqrt(25. - y * y) + 12.5 * std::asin(y / 5.)));
        const auto d = unwrap(cut(m, b, c));
        expect_volume(m, d, 2. * sliver, 2);
        // the slivers touch along the line of contact, each with its own edge there
        const auto e0 = explore<EdgeId>(m, ShapeId{d.solids[0]}), e1 = explore<EdgeId>(m, ShapeId{d.solids[1]});
        ASSERT_TRUE(std::ranges::none_of(e0, [&](EdgeId e) { return std::ranges::find(e1, e) != e1.end(); }));
        expect_volume(m, unwrap(common(m, b, c)), 24. - 2. * sliver);
        expect_volume(m, unwrap(fuse(m, b, c)), 250. * pi + 2. * sliver);
    }
}

TEST(tests_brep_boolean, edges_on_faces)
{
    Model<T> m;
    const auto bx = read(m, "box.stp"); // 0..2 x 0..3 x 0..4
    const auto cyl = read(m, "cylinder.stp"); // R 5, z 0..10
    const point<T, 3> Y{0., 1., 0.}, Z{0., 0., 1.};
    const T c = std::sqrt(0.5);

    // the box turned by 45 degrees about y, resting on its edge of length 4 along the top of the
    // other at x = 1, from y = -0.5 to 3.5: the edge crosses the boundary of the top face
    {
        const auto s = copy_shape(m, bx, Placement<T>::frame({1., 3.5, 4.}, {c, 0., c}, {-c, 0., c}, {0., -1., 0.}));
        const auto f = unwrap(fuse(m, bx, s));
        expect_volume(m, f, 48., 2);
        // the top is split in two along the edge
        ASSERT_EQ(explore<FaceId>(m, f.shape).size(), 13u);
        ASSERT_TRUE(unwrap(common(m, bx, s)).solids.empty());
        const auto d = unwrap(cut(m, bx, s));
        expect_volume(m, d, 24.);
        ASSERT_EQ(explore<FaceId>(m, d.shape).size(), 7u);
    }
    // boxes touching along a part of an edge: x = 2, y = 3, z in 1..4
    {
        const auto s = copy_shape(m, bx, Placement<T>{.t = {2., 3., 1.}});
        expect_volume(m, unwrap(fuse(m, bx, s)), 48., 2);
        ASSERT_TRUE(unwrap(common(m, bx, s)).solids.empty());
        expect_volume(m, unwrap(cut(m, bx, s)), 24.);
    }
    // the cylinder lying on the box with its seam on the top: the seam lies on the face
    {
        const auto s = copy_shape(m, cyl, Placement<T>::frame({1., -3.5, 9.}, T(-1) * Z, {-1., 0., 0.}, Y));
        const auto f = unwrap(fuse(m, bx, s));
        expect_volume(m, f, 24. + 250. * pi, 2);
        ASSERT_TRUE(unwrap(common(m, bx, s)).solids.empty());
        expect_volume(m, unwrap(cut(m, bx, s)), 24.);
    }
}

TEST(tests_brep_boolean, branch_points)
{
    // two cylinders R 5, length 30, axes z and x crossing at the origin (Steinmetz): their
    // intersection curves cross at the branch points (0, +-5, 0), where the surfaces are tangent
    Model<T> m;
    const auto cyl = read(m, "cylinder.stp"); // R 5, z 0..10
    FaceId disc;
    for (auto f : explore<FaceId>(m, cyl))
        if (dynamic_cast<const Plane<T> *>(m.face(f).surface.get()) && classify(m, f, point<T, 2>{0., 0.}) != Location::Out &&
            std::abs(m.face(f).surface->value(0., 0.)[2]) < 1e-9)
            disc = f;
    REQUIRE(disc.valid());
    const auto bar = unwrap(make_prism(m, disc, point<T, 3>{0., 0., 30.})).shape;
    const T a = 1., ca = std::cos(a), sa = std::sin(a);
    const auto c1 = copy_shape(m, bar, Placement<T>::frame({0., 0., -15.}, {ca, sa, 0.}, {-sa, ca, 0.}, {0., 0., 1.}));
    const auto c2 = copy_shape(m, bar, Placement<T>::frame({-15., 0., 0.}, {0., ca, sa}, {0., -sa, ca}, {1., 0., 0.}));
    const T steinmetz = 16. * 125. / 3., v = 25. * pi * 30.;
    expect_volume(m, unwrap(common(m, c1, c2)), steinmetz);
    // the cut leaves two halves of c1, touching at the branch points
    expect_volume(m, unwrap(cut(m, c1, c2)), v - steinmetz, 2);
    expect_volume(m, unwrap(fuse(m, c1, c2)), 2. * v - steinmetz);
}

TEST(tests_brep_boolean, poles)
{
    // the sphere R 3 centred on the side of the cylinder R 5, axes parallel: both poles of the
    // sphere lie on the cylinder, and the section goes through them
    Model<T> m;
    const auto sph = read(m, "sphere.stp");
    const auto cyl = read(m, "cylinder.stp");
    const T a = 1.;
    const auto s = copy_shape(m, sph, Placement<T>{.t = {5. * std::cos(a), 5. * std::sin(a), 5.}});
    // reference: slices of the sphere, the lens of two discs (radii 5 and r(z), centres 5 apart)
    auto lens = [](T R, T r, T d) {
        const T x = r * r * std::acos((d * d + r * r - R * R) / (2. * d * r)) + R * R * std::acos((d * d + R * R - r * r) / (2. * d * R));
        return x - 0.5 * std::sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R));
    };
    T inside{};
    constexpr int n = 20000;
    for (int k = 0; k < n; ++k)
    {
        const T z = -3. + 6. * (k + 0.5) / n;
        inside += lens(5., std::sqrt(9. - z * z), 5.) * 6. / n;
    }
    const T ball = 36. * pi;
    expect_volume(m, unwrap(common(m, cyl, s)), inside);
    expect_volume(m, unwrap(cut(m, cyl, s)), 250. * pi - inside);
    expect_volume(m, unwrap(fuse(m, cyl, s)), 250. * pi + ball - inside);
    expect_volume(m, unwrap(cut(m, s, cyl)), ball - inside);
}
