#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/check.h>
#include <gbs-brep/sweep.h>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>
#include "tests_brep_geom.h"

#include <algorithm>
#include <array>
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

    /// The number of shapes of each kind in the model: a refused operation leaves it unchanged.
    std::array<std::size_t, 7> counts(const Model<T> &m)
    {
        std::array<std::size_t, 7> n{};
        for (std::size_t t = 0; t < n.size(); ++t)
            n[t] = m.count(ShapeType(t));
        return n;
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

TEST(tests_brep_boolean, bspline_tangency)
{
    // the cylinder lying on the top of a cube with B-spline faces: a tangency along a line found
    // by marching, traced as such (refused as `NonGeneric` before)
    Model<T> m;
    const auto cyl = read(m, "cylinder.stp"); // R 5, z 0..10
    const auto cube = ShapeId{box(m, {-1., -1., 0.}, 4.)};
    const point<T, 3> Y{0., 1., 0.}, Z{0., 0., 1.};
    const auto s = copy_shape(m, cyl, Placement<T>::frame({1., -3.5, 9.}, Z, cross(Y, Z), Y));
    expect_volume(m, unwrap(fuse(m, cube, s)), 64. + 250. * pi, 2);
    ASSERT_TRUE(unwrap(common(m, cube, s)).solids.empty());
    expect_volume(m, unwrap(cut(m, cube, s)), 64.);
    // an incomplete intersection (too few pairs of patches allowed) is refused, and nothing is
    // left of the attempt
    const auto before = counts(m);
    BooleanOptions<T> few;
    few.intersect.max_pairs = 4;
    auto r = fuse(m, cube, s, few);
    if (r)
        FAIL("expected NonGeneric");
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    ASSERT_EQ(counts(m), before);
    // the operands are kept
    expect_valid(m, cube);
    expect_valid(m, s);
    // errors
    ASSERT_FALSE(fuse(m, cube, ShapeId{SolidId{999}}).has_value());
    const auto face = explore<FaceId>(m, cube).front();
    ASSERT_FALSE(fuse(m, cube, ShapeId{face}).has_value());
}

TEST(tests_brep_boolean, tangency_ending_in_a_face)
{
    // a solid -1..1 x -1..1 under the B-spline face z = -bend max(x, 0)^3 (a cubic in x with a
    // knot at 0, exact), down to z = -2; the cylinder of radius 5 along x lies on it along y = 0,
    // over x in -5..5 (clearance 5 - sqrt(24) = 0.1 above the sides y = +-1). With bend = 1 it is
    // tangent for x <= 0 and the face bends away beyond: the tangency runs from the edge with the
    // side x = -1 to x = 0.01, where the gap x^3 reaches the tolerance, inside the face
    using brep_tests::quad;
    auto operands = [](Model<T> &m, T bend) {
        const std::vector<T> kx{-1., -1., -1., -1., 0., 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const std::array<T, 5> gx{-1., -2. / 3., 0., 2. / 3., 1.}, zx{0., 0., 0., 0., -bend};
        auto bspline = [&](const points_vector<T, 3> &poles) -> std::shared_ptr<Surface<T, 3>> {
            return std::make_shared<BSSurface<T, 3>>(poles, kx, k1, 3, 1);
        };
        // the top, between y = -1 and 1; a side, ruled between z = -2 and the curve of the top
        auto top = [&] {
            points_vector<T, 3> p;
            for (T y : {T(-1), T(1)})
                for (int i = 0; i < 5; ++i)
                    p.push_back({gx[i], y, zx[i]});
            return p;
        };
        auto side = [&](T y) {
            points_vector<T, 3> p;
            for (int i = 0; i < 5; ++i)
                p.push_back({gx[i], y, -2.});
            for (int i = 0; i < 5; ++i)
                p.push_back({gx[i], y, zx[i]});
            return p;
        };
        std::vector<FaceId> faces;
        for (const auto &srf : std::array{bspline(top()), bspline(side(-1.)), bspline(side(1.)),
                                          quad({-1., -1., -2.}, {1., -1., -2.}, {-1., 1., -2.}, {1., 1., -2.}),
                                          quad({-1., -1., -2.}, {-1., 1., -2.}, {-1., -1., 0.}, {-1., 1., 0.}),
                                          quad({1., -1., -2.}, {1., 1., -2.}, {1., -1., -bend}, {1., 1., -bend})})
            faces.push_back(unwrap(make_face(m, srf)));
        const auto block = ShapeId{unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()))};
        const point<T, 3> X{1., 0., 0.}, Y{0., 1., 0.}, Z{0., 0., 1.};
        const auto cyl = copy_shape(m, read(m, "cylinder.stp"), Placement<T>::frame({-5., 0., 5.}, Y, Z, X)); // R 5, along x
        return std::pair{block, cyl};
    };
    // control: the flat top (bend = 0), tangent all along, is fused like the B-spline cube
    {
        Model<T> m;
        const auto [block, cyl] = operands(m, 0.);
        expect_volume(m, unwrap(fuse(m, block, cyl)), 8. + 250. * pi, 2);
        ASSERT_TRUE(unwrap(common(m, block, cyl)).solids.empty());
        expect_volume(m, unwrap(cut(m, block, cyl)), 8.);
    }
    // the section ends inside the face: one tangent edge, from (-1, 0, 0) to x = 0.01; singular,
    // and complete (the refusal below does not come from a cap). On its own model: the section
    // adds its edges to the model
    {
        Model<T> m;
        const auto [block, cyl] = operands(m, 1.);
        const auto sec = unwrap(section(m, block, cyl));
        CHECK(sec.singular);
        CHECK(sec.complete);
        REQUIRE_EQ(sec.edges.size(), 1u);
        CHECK(sec.edges[0].tangent);
        const auto &e = m.edge(sec.edges[0].edge);
        auto x1 = m.vertex(e.v1).pnt, x2 = m.vertex(e.v2).pnt;
        if (x1[0] > x2[0])
            std::swap(x1, x2);
        CHECK(norm(x1 - point<T, 3>{-1., 0., 0.}) <= 1e-6);
        CHECK(x2[0] > 0.);
        CHECK(x2[0] < 0.05);
    }
    // refused as `NonGeneric`, nothing left of the attempt, the operands kept
    Model<T> m;
    const auto [block, cyl] = operands(m, 1.);
    expect_valid(m, block);
    const T v_block = volume(m, block);
    // 8 under z = 0, less the integral of x^3 over 0..1 times the width 2
    REQUIRE(v_block == doctest::Approx(7.5).epsilon(1e-3));
    const auto before = counts(m);
    for (auto op : {BooleanOp::Fuse, BooleanOp::Common, BooleanOp::Cut})
    {
        auto r = boolean(m, op, block, cyl);
        if (r)
            FAIL("expected NonGeneric");
        CHECK_EQ(r.error().code, BuildErrc::NonGeneric);
        CHECK(r.error().message.find("singular or incomplete intersection") != std::string::npos);
        CHECK_EQ(counts(m), before);
    }
    expect_valid(m, block);
    expect_valid(m, cyl);
    CHECK_EQ(volume(m, block), v_block);
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

TEST(tests_brep_boolean, coincident_bsplines)
{
    // cubes with bilinear B-spline faces, and a box with planes, whose faces coincide
    Model<T> m;
    const auto a = ShapeId{box(m, {0., 0., 0.}, 2.)};
    // shifted in x and y: bottoms and tops coincide on 1..2 x 1..2
    const auto b = ShapeId{box(m, {1., 1., 0.}, 2.)};
    expect_volume(m, unwrap(fuse(m, a, b)), 14.);
    expect_volume(m, unwrap(common(m, a, b)), 2.);
    expect_volume(m, unwrap(cut(m, a, b)), 6.);
    // stacked, face against face: one solid, nothing in common
    const auto up = ShapeId{box(m, {0., 0., 2.}, 2.)};
    const auto f = unwrap(fuse(m, a, up));
    expect_volume(m, f, 16.);
    ASSERT_EQ(explore<FaceId>(m, f.shape).size(), 10u);
    ASSERT_TRUE(unwrap(common(m, a, up)).solids.empty());
    expect_volume(m, unwrap(cut(m, a, up)), 8.);
    // a box with planes (0..2 x 0..3 x 0..4) and a B-spline cube 1..3 x 1..3 x 0..2: bottoms on
    // the same plane, one a primitive, the other a B-spline
    const auto planes = read(m, "box.stp");
    const auto cube = ShapeId{box(m, {1., 1., 0.}, 2.)};
    expect_volume(m, unwrap(fuse(m, planes, cube)), 28.);
    expect_volume(m, unwrap(common(m, planes, cube)), 4.);
    expect_volume(m, unwrap(cut(m, planes, cube)), 20.);
    expect_volume(m, unwrap(cut(m, cube, planes)), 4.);
}

namespace
{
    /// A cylinder of radius r and height h on the z axis whose side is a closed rational B-spline
    /// (the disc bounded by the rational B-spline circle, extruded): seam at x = r.
    ShapeId nurbs_cylinder(Model<T> &m, T r, T h)
    {
        auto circle = std::make_shared<BSCurveRational<T, 3>>(build_circle<T, 3>(r));
        const auto e = unwrap(make_edge(m, std::static_pointer_cast<Curve<T, 3>>(circle)));
        const auto w = unwrap(make_wire(m, std::vector<EdgeId>{e}));
        auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{point<T, 3>{0., 0., 0.}, point<T, 3>{0., 0., 1.}, point<T, 3>{1., 0., 0.}},
                                                std::array<T, 4>{-2. * r, 2. * r, -2. * r, 2. * r});
        const auto disc = unwrap(make_face(m, std::static_pointer_cast<Surface<T, 3>>(plane), w));
        return unwrap(make_prism(m, ShapeId{disc}, point<T, 3>{0., 0., h})).shape;
    }
}

TEST(tests_brep_boolean, coincident_closed_bsplines)
{
    Model<T> m;
    const auto a = nurbs_cylinder(m, 1., 2.);
    // the side is a closed B-spline, not a primitive
    std::size_t bsplines{};
    for (auto f : explore<FaceId>(m, a))
        bsplines += dynamic_cast<const ElementarySurface<T> *>(m.face(f).surface.get()) == nullptr;
    ASSERT_EQ(bsplines, 1u);
    expect_volume(m, unwrap(common(m, a, a)), pi * 2.);

    // coaxial, shifted by 1 and turned by 1 rad: the sides coincide on 1..2, the seams differ and
    // the circles of each cross the seam of the other
    const T c = std::cos(1.), s = std::sin(1.);
    const auto b = copy_shape(m, a, Placement<T>::frame({0., 0., 1.}, {c, s, 0.}, {-s, c, 0.}, {0., 0., 1.}));
    expect_volume(m, unwrap(fuse(m, a, b)), pi * 3.);
    expect_volume(m, unwrap(common(m, a, b)), pi);
    expect_volume(m, unwrap(cut(m, a, b)), pi);

    // stacked: the sides touch along a circle only (no point of either faces the other)
    const auto up = copy_shape(m, a, Placement<T>{.t = {0., 0., 2.}});
    const auto f = unwrap(fuse(m, a, up));
    expect_volume(m, f, pi * 4.);
    ASSERT_TRUE(unwrap(common(m, a, up)).solids.empty());

    // a primitive cylinder (R 5, z 0..10) and a B-spline one of the same radius on 5..15
    const auto cyl = read(m, "cylinder.stp");
    const auto big = copy_shape(m, nurbs_cylinder(m, 5., 10.), Placement<T>::frame({0., 0., 5.}, {c, s, 0.}, {-s, c, 0.}, {0., 0., 1.}));
    expect_volume(m, unwrap(fuse(m, cyl, big)), 25. * pi * 15.);
    expect_volume(m, unwrap(common(m, cyl, big)), 25. * pi * 5.);
    expect_volume(m, unwrap(cut(m, big, cyl)), 25. * pi * 5.);
    // a box resting on the top of the B-spline cylinder (-1..1 x -1.5..1.5 on z = 15)
    const auto lid = copy_shape(m, read(m, "box.stp"), Placement<T>{.t = {-1., -1.5, 15.}});
    expect_volume(m, unwrap(fuse(m, big, lid)), 25. * pi * 10. + 24.);
}

TEST(tests_brep_boolean, exact_plane_pcurves)
{
    // the caps of a cylinder whose side is a rational B-spline: their edges are iso lines of the
    // side (curves on it), their pcurves on the planes exact rational B-splines of the circle's 9
    // poles, instead of interpolations of hundreds of samples
    Model<T> m;
    const auto a = nurbs_cylinder(m, 1., 2.);
    std::size_t caps{};
    for (auto f : explore<FaceId>(m, a))
    {
        const auto &srf = *m.face(f).surface;
        if (!dynamic_cast<const Plane<T> *>(&srf))
            continue;
        ++caps;
        for (auto w : m.face(f).wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                const auto *pc = dynamic_cast<const BSCurveRational<T, 2> *>(ce.pcurve.get());
                REQUIRE(pc != nullptr);
                ASSERT_EQ(pc->poles().size(), 9u);
                const auto &e = m.edge(ce.edge);
                for (int k = 0; k <= 50; ++k)
                {
                    const T t = e.u1 + (e.u2 - e.u1) * k / 50.;
                    const auto uv = pc->value(t);
                    ASSERT_LE(norm(srf.value(uv[0], uv[1]) - e.curve->value(t)), 1e-12);
                }
            }
    }
    ASSERT_EQ(caps, 2u);
    // written to STEP and read back exactly (the rational pcurves as rational B-splines)
    {
        auto w = step::write_step(m, a);
        if (!w)
            FAIL(w.error().message);
        ASSERT_EQ(w->second.approximated_curves, 0u);
        Model<T> back;
        auto r = step::read_step(back, w->first);
        REQUIRE(r.has_value());
        REQUIRE(r->report.check.ok());
        ASSERT_NEAR(volume(back, r->root), pi * 2., 2e-3 * pi * 2.);
    }
    // a plane face bounded by B-spline segments: segments in (u, v)
    const auto plane = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 1.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{-5., 5., -5., 5.});
    std::vector<EdgeId> edges;
    const std::vector<point<T, 3>> pts{{0., 0., 1.}, {2., 0., 1.}, {2., 1., 1.}, {0., 1., 1.}};
    for (std::size_t k = 0; k < pts.size(); ++k)
        edges.push_back(unwrap(make_edge(m, pts[k], pts[(k + 1) % pts.size()])));
    const auto f = unwrap(make_face(m, std::static_pointer_cast<Surface<T, 3>>(plane), unwrap(make_wire(m, edges))));
    for (const auto &ce : m.wire(m.face(f).wires.front()).coedges)
    {
        const auto *pc = dynamic_cast<const BSCurve<T, 2> *>(ce.pcurve.get());
        REQUIRE(pc != nullptr);
        ASSERT_EQ(pc->poles().size(), 2u);
    }
}
