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

    /// The box [a, b]: the prism of a bilinear face.
    SolidId block(Model<T> &m, const point<T, 3> &a, const point<T, 3> &b)
    {
        using namespace brep_tests;
        const auto f = unwrap(make_face(m, quad({a[0], a[1], a[2]}, {b[0], a[1], a[2]}, {a[0], b[1], a[2]}, {b[0], b[1], a[2]})));
        return explore<SolidId>(m, unwrap(make_prism(m, ShapeId{f}, point<T, 3>{0., 0., b[2] - a[2]})).shape).front();
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
    // over x in -5..5 (clearance 5 - sqrt(24) = 0.1 above the sides y = +-1). It is tangent for
    // x <= 0. With bend = 1 the face bends away beyond: the tangency runs from the edge with the
    // side x = -1 to x = 0.01, where the gap x^3 reaches the tolerance, and ends there inside the
    // face, which it does not split. With bend = -1 the face rises into the cylinder beyond: the
    // surfaces cross, two transversal curves start at the end of the tangency, refused
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
    // the section ends inside the face: one tangent edge, from (-1, 0, 0) to x = 0.01, neither
    // singular nor incomplete. On its own model: the section adds its edges to the model
    {
        Model<T> m;
        const auto [block, cyl] = operands(m, 1.);
        const auto sec = unwrap(section(m, block, cyl));
        CHECK_FALSE(sec.singular);
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
    // the tangency ending inside the face is accepted: an edge with a free end, which splits
    // nothing; volumes 7.5 (8 under z = 0, less the integral of x^3 over 0..1 times the width 2)
    // and 250 pi
    {
        Model<T> m;
        const auto [block, cyl] = operands(m, 1.);
        REQUIRE(volume(m, block) == doctest::Approx(7.5).epsilon(1e-3));
        expect_volume(m, unwrap(fuse(m, block, cyl)), 7.5 + 250. * pi, 2);
        ASSERT_TRUE(unwrap(common(m, block, cyl)).solids.empty());
        expect_volume(m, unwrap(cut(m, block, cyl)), 7.5);
    }
    // the face rising into the cylinder (bend = -1): the surfaces cross beyond the end of the
    // tangency, refused as `NonGeneric`, nothing left of the attempt, the operands kept
    Model<T> m;
    const auto [block, cyl] = operands(m, -1.);
    expect_valid(m, block);
    const T v_block = volume(m, block);
    REQUIRE(v_block == doctest::Approx(8.5).epsilon(1e-3));
    {
        Model<T> ms;
        const auto [b2, c2] = operands(ms, -1.);
        const auto sec = unwrap(section(ms, b2, c2));
        CHECK(sec.singular);
        CHECK(sec.complete);
    }
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

TEST(tests_brep_boolean, partial_coincidence)
{
    // the block -1..1 x -1..1 under the B-spline face z = -bend max(x, 0)^3 (exact, a knot at 0),
    // down to z = -2, and the box -1..0.5 x -1..1 x 0..1 on it: the bottom of the box lies on the
    // top of the block over its patch x <= 0. With bend = 1 the top bends away beyond (a gap
    // x^3 under the box); with bend = -1 it rises into the box. Either way the surfaces do not
    // cross beyond the knot line x = 0, the border of the common region
    using brep_tests::quad;
    // the top as a B-spline in (x, y), or transposed (y, x: its knot line along u), or rational
    // (unit weights); the box over x0..0.5
    enum class Top
    {
        plain,
        transposed,
        rational
    };
    auto operands = [](Model<T> &m, T bend, Top kind = Top::plain, T x0 = -1.) {
        const std::vector<T> kx{-1., -1., -1., -1., 0., 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const std::array<T, 5> gx{-1., -2. / 3., 0., 2. / 3., 1.}, zx{0., 0., 0., 0., -bend};
        auto bspline = [&](const points_vector<T, 3> &poles) -> std::shared_ptr<Surface<T, 3>> {
            return std::make_shared<BSSurface<T, 3>>(poles, kx, k1, 3, 1);
        };
        auto top_surface = [&]() -> std::shared_ptr<Surface<T, 3>> {
            points_vector<T, 3> p;
            if (kind == Top::transposed)
            {
                for (int i = 0; i < 5; ++i)
                    for (T y : {T(-1), T(1)})
                        p.push_back({gx[i], y, zx[i]});
                return std::make_shared<BSSurface<T, 3>>(p, k1, kx, 1, 3);
            }
            for (T y : {T(-1), T(1)})
                for (int i = 0; i < 5; ++i)
                    p.push_back({gx[i], y, zx[i]});
            if (kind == Top::plain)
                return bspline(p);
            points_vector<T, 4> w;
            for (const auto &x : p)
                w.push_back({x[0], x[1], x[2], 1.});
            return std::make_shared<BSSurfaceRational<T, 3>>(w, kx, k1, 3, 1);
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
        for (const auto &srf : std::array{top_surface(), bspline(side(-1.)), bspline(side(1.)),
                                          quad({-1., -1., -2.}, {1., -1., -2.}, {-1., 1., -2.}, {1., 1., -2.}),
                                          quad({-1., -1., -2.}, {-1., 1., -2.}, {-1., -1., 0.}, {-1., 1., 0.}),
                                          quad({1., -1., -2.}, {1., 1., -2.}, {1., -1., -bend}, {1., 1., -bend})})
            faces.push_back(unwrap(make_face(m, srf)));
        const auto block = ShapeId{unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()))};
        const point<T, 3> lo{x0, -1., 0.}, hi{0.5, 1., 1.};
        std::vector<FaceId> box_faces;
        for (const auto &srf : std::array{quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}),
                                          quad({lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}),
                                          quad({lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}),
                                          quad({hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]}, {hi[0], lo[1], hi[2]}, {hi[0], hi[1], hi[2]})})
            box_faces.push_back(unwrap(make_face(m, srf)));
        const auto box = ShapeId{unwrap(make_solid(m, unwrap(sew(m, box_faces)).shells.front()))};
        return std::tuple{block, box, faces.front(), box_faces};
    };
    // the detection: the top of the block is the carrier, its border the knot line x = 0 across
    // the face, and the gap keeps one sign beyond
    {
        Model<T> m;
        const auto [block, box, top, box_faces] = operands(m, 1.);
        const auto bottom = box_faces[0];
        REQUIRE(m.alive(top));
        REQUIRE(m.alive(bottom));
        const auto p = partial_coincidence(m, top, bottom);
        REQUIRE(p);
        CHECK_EQ(p->carrier, top);
        CHECK_EQ(p->other, bottom);
        CHECK(p->parts());
        REQUIRE_EQ(p->borders.size(), 1u);
        const auto &b = p->borders.front();
        CHECK_FALSE(b.along_u); // u = 0 along v
        CHECK_EQ(b.fixed, 0.);
        CHECK_EQ(b.lo, 0.); // (v in 0..1 across y)
        CHECK_EQ(b.hi, 1.);
        // faces wholly on each other, or apart, are not partial
        CHECK_FALSE(partial_coincidence(m, box_faces[1], top)); // the top of the box, apart
    }
    // the transposed top: its border is the knot line v = 0, along u
    {
        Model<T> m;
        const auto [block, box, top, box_faces] = operands(m, 1., Top::transposed);
        const auto p = partial_coincidence(m, top, box_faces[0]);
        REQUIRE(p);
        REQUIRE_EQ(p->borders.size(), 1u);
        CHECK(p->borders.front().along_u);
        CHECK_EQ(p->borders.front().fixed, 0.);
        CHECK_EQ(p->borders.front().lo, 0.);
        CHECK_EQ(p->borders.front().hi, 1.);
    }
    // the same booleans with the top transposed, rational, and under a box from x0 = -0.3: the
    // patch -1..0 of the top then holds points under the box (on its bottom) and points beyond it
    // (not facing it), on one piece, which coincides with the bottom of the box
    for (auto [kind, x0] : {std::pair{Top::transposed, -1.}, std::pair{Top::rational, -1.}, std::pair{Top::plain, -0.3}})
    {
        INFO("kind " << int(kind) << ", x0 " << x0);
        Model<T> m;
        const auto [block, box, top, box_faces] = operands(m, 1., kind, x0);
        const T v_box = (0.5 - x0) * 2.;
        REQUIRE(volume(m, block) == doctest::Approx(7.5).epsilon(1e-3));
        expect_volume(m, unwrap(fuse(m, block, box)), 7.5 + v_box);
        CHECK(unwrap(common(m, block, box)).solids.empty());
        expect_volume(m, unwrap(cut(m, block, box)), 7.5);
    }
    // a crossing narrower than the sampling: the top z = -f, f = x ((x - 1/2)^2 - e) on 0..1 (a
    // free Bézier piece after a knot of multiplicity 3 at 0), e = 1/200: under the plane z = 0
    // but on the strip 1/2 +- sqrt(e), 0.14 wide, between the samples of the patch (step 1/6).
    // The Bernstein coefficients of the gap, the heights of the poles above the plane, have both
    // signs: the patch is intersected, and the common part is the sliver over the strip
    {
        Model<T> m;
        const T e = 1. / 200.;
        const std::vector<T> kx{-1., -1., -1., -1., 0., 0., 0., 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const std::array<T, 7> gx{-1., -2. / 3., -1. / 3., 0., 1. / 3., 2. / 3., 1.};
        const std::array<T, 7> zx{0., 0., 0., 0., -(0.25 - e) / 3., (0.5 + 2. * e) / 3., -(0.25 - e)};
        auto bspline = [&](const points_vector<T, 3> &poles) -> std::shared_ptr<Surface<T, 3>> {
            return std::make_shared<BSSurface<T, 3>>(poles, kx, k1, 3, 1);
        };
        auto rows = [&](T y, bool from_bottom) {
            points_vector<T, 3> p;
            if (from_bottom)
                for (int i = 0; i < 7; ++i)
                    p.push_back({gx[i], y, -2.});
            for (int i = 0; i < 7; ++i)
                p.push_back({gx[i], y, zx[i]});
            return p;
        };
        points_vector<T, 3> top = rows(-1., false), top_hi = rows(1., false);
        top.insert(top.end(), top_hi.begin(), top_hi.end());
        const auto top_srf = bspline(top);
        auto f = [&](T x) { return x * ((x - 0.5) * (x - 0.5) - e); };
        for (T x : {0.1, 0.5, 0.9})
            REQUIRE(top_srf->value(x, 0.5)[2] == doctest::Approx(-f(x)).epsilon(1e-12));
        std::vector<FaceId> faces;
        for (const auto &srf : std::array{top_srf, bspline(rows(-1., true)), bspline(rows(1., true)),
                                          quad({-1., -1., -2.}, {1., -1., -2.}, {-1., 1., -2.}, {1., 1., -2.}),
                                          quad({-1., -1., -2.}, {-1., 1., -2.}, {-1., -1., 0.}, {-1., 1., 0.}),
                                          quad({1., -1., -2.}, {1., 1., -2.}, {1., -1., -f(1.)}, {1., 1., -f(1.)})})
            faces.push_back(unwrap(make_face(m, srf)));
        const auto block = ShapeId{unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()))};
        const point<T, 3> lo{-0.8, -0.8, 0.}, hi{0.9, 0.8, 1.};
        std::vector<FaceId> box_faces;
        for (const auto &srf : std::array{quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}),
                                          quad({lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}),
                                          quad({lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}),
                                          quad({hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]}, {hi[0], lo[1], hi[2]}, {hi[0], hi[1], hi[2]})})
            box_faces.push_back(unwrap(make_face(m, srf)));
        const auto box = ShapeId{unwrap(make_solid(m, unwrap(sew(m, box_faces)).shells.front()))};
        const auto p = partial_coincidence(m, faces.front(), box_faces.front());
        REQUIRE(p);
        REQUIRE_EQ(p->crossing.size(), 1u);
        // (drawn back from the border x = 0, where the surfaces meet at a crease, an edge of both)
        CHECK(p->crossing.front()[0] > 0.);
        CHECK(p->crossing.front()[0] < 1e-3);
        CHECK_EQ(p->crossing.front()[1], 1.);
        // the sliver: 1.6 times the integral of -f over the strip, F(x) = x^4 / 4 - x^3 / 3 + (1/4 - e) x^2 / 2
        auto F = [&](T x) { return std::pow(x, 4) / 4. - std::pow(x, 3) / 3. + (0.25 - e) * x * x / 2.; };
        const T sliver = -1.6 * (F(0.5 + std::sqrt(e)) - F(0.5 - std::sqrt(e)));
        const auto c = unwrap(common(m, block, box));
        REQUIRE_EQ(c.solids.size(), 1u);
        expect_valid(m, c.shape);
        CHECK(volume(m, c.shape) == doctest::Approx(sliver).epsilon(0.05));
    }
    // an L-shaped solid (two boxes of heights 1 and 0.5 side by side in y, fused) on the top: its
    // two bottom faces lie on the same patch of the top and part from it beyond x = 0; the top is
    // analysed against both before it is split
    {
        Model<T> m;
        const auto [block, box, top, box_faces] = operands(m, 1.);
        auto cuboid = [&](const point<T, 3> &lo, const point<T, 3> &hi) {
            std::vector<FaceId> f;
            for (const auto &srf : std::array{quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}),
                                              quad({lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                              quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}),
                                              quad({lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                              quad({lo[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}),
                                              quad({hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]}, {hi[0], lo[1], hi[2]}, {hi[0], hi[1], hi[2]})})
                f.push_back(unwrap(make_face(m, srf)));
            return ShapeId{unwrap(make_solid(m, unwrap(sew(m, f)).shells.front()))};
        };
        const auto l = unwrap(fuse(m, cuboid({-0.8, -0.8, 0.}, {0.5, 0., 1.}), cuboid({-0.8, 0., 0.}, {0.5, 0.8, 0.5})));
        REQUIRE_EQ(l.solids.size(), 1u);
        const T v_l = 1.3 * 0.8 * 1.5;
        expect_volume(m, l, v_l);
        expect_volume(m, unwrap(fuse(m, block, l.shape)), 7.5 + v_l);
        CHECK(unwrap(common(m, block, l.shape)).solids.empty());
        expect_volume(m, unwrap(cut(m, block, l.shape)), 7.5);
    }
    // volumes: the block 8 - bend / 2 (the integral of x^3 over 0..1 times the width 2), the box
    // 3, their overlap where the top rises into the box: 2 times the integral of x^3 over 0..0.5
    for (T bend : {T(1), T(-1)})
    {
        INFO("bend " << bend);
        Model<T> m;
        const auto [block, box, top, box_faces] = operands(m, bend);
        const T v_block = 8. - bend / 2., overlap = bend < 0 ? 2. * std::pow(0.5, 4) / 4. : 0.;
        REQUIRE(volume(m, block) == doctest::Approx(v_block).epsilon(1e-3));
        REQUIRE(volume(m, box) == doctest::Approx(3.).epsilon(1e-3));
        expect_volume(m, unwrap(fuse(m, block, box)), v_block + 3. - overlap);
        const auto c = unwrap(common(m, block, box));
        if (bend > 0)
            CHECK(c.solids.empty());
        else
            expect_volume(m, c, overlap);
        expect_volume(m, unwrap(cut(m, block, box)), v_block - overlap);
        expect_valid(m, block);
        expect_valid(m, box);
    }
    // the top z = -f(x), f = max(x, 0)^3 - 16 max(x - 0.5, 0)^3 (knots at 0 and 0.5, coefficients
    // by blossoming), under a box: on its bottom for x <= 0, parting from
    // it over 0..0.5, then crossing it at x0 = 0.5 c / (c - 1), c = 16^(1/3), and rising into the
    // box. The patch 0..0.5 is not intersected, the patch 0.5..1 is, where the surfaces cross
    {
        Model<T> m;
        const std::vector<T> kx{-1., -1., -1., -1., 0., 0.5, 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
        const std::array<T, 6> gx{-1., -2. / 3., -1. / 6., 0.5, 5. / 6., 1.}, zx{0., 0., 0., 0., -0.5, 1.};
        auto bspline = [&](const points_vector<T, 3> &poles) -> std::shared_ptr<Surface<T, 3>> {
            return std::make_shared<BSSurface<T, 3>>(poles, kx, k1, 3, 1);
        };
        auto rows = [&](T y, bool from_bottom) {
            points_vector<T, 3> p;
            if (from_bottom)
                for (int i = 0; i < 6; ++i)
                    p.push_back({gx[i], y, -2.});
            for (int i = 0; i < 6; ++i)
                p.push_back({gx[i], y, zx[i]});
            return p;
        };
        points_vector<T, 3> top = rows(-1., false), top_hi = rows(1., false);
        top.insert(top.end(), top_hi.begin(), top_hi.end());
        std::vector<FaceId> faces;
        for (const auto &srf : std::array{bspline(top), bspline(rows(-1., true)), bspline(rows(1., true)),
                                          quad({-1., -1., -2.}, {1., -1., -2.}, {-1., 1., -2.}, {1., 1., -2.}),
                                          quad({-1., -1., -2.}, {-1., 1., -2.}, {-1., -1., 0.}, {-1., 1., 0.}),
                                          quad({1., -1., -2.}, {1., 1., -2.}, {1., -1., 1.}, {1., 1., 1.})})
            faces.push_back(unwrap(make_face(m, srf)));
        const auto block = ShapeId{unwrap(make_solid(m, unwrap(sew(m, faces)).shells.front()))};
        // the box -0.8..0.9 x -0.8..0.8 x 0..1, inside the block in x and y
        const point<T, 3> lo{-0.8, -0.8, 0.}, hi{0.9, 0.8, 1.};
        std::vector<FaceId> box_faces;
        for (const auto &srf : std::array{quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}),
                                          quad({lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {hi[0], lo[1], lo[2]}, {lo[0], lo[1], hi[2]}, {hi[0], lo[1], hi[2]}),
                                          quad({lo[0], hi[1], lo[2]}, {hi[0], hi[1], lo[2]}, {lo[0], hi[1], hi[2]}, {hi[0], hi[1], hi[2]}),
                                          quad({lo[0], lo[1], lo[2]}, {lo[0], hi[1], lo[2]}, {lo[0], lo[1], hi[2]}, {lo[0], hi[1], hi[2]}),
                                          quad({hi[0], lo[1], lo[2]}, {hi[0], hi[1], lo[2]}, {hi[0], lo[1], hi[2]}, {hi[0], hi[1], hi[2]})})
            box_faces.push_back(unwrap(make_face(m, srf)));
        const auto box = ShapeId{unwrap(make_solid(m, unwrap(sew(m, box_faces)).shells.front()))};
        // the block: 8 (f integrates to 0 over 0..1); the overlap: 1.6 times the integral of -f
        // over x0..0.9, F(x) = 4 (x - 0.5)^4 - x^4 / 4 its primitive
        const T c = std::cbrt(16.), x0 = 0.5 * c / (c - 1.), v_box = 1.7 * 1.6;
        auto F = [](T x) { return 4. * std::pow(x - 0.5, 4) - std::pow(x, 4) / 4.; };
        const T overlap = 1.6 * (F(0.9) - F(x0));
        REQUIRE(volume(m, block) == doctest::Approx(8.).epsilon(1e-3));
        const auto p = partial_coincidence(m, faces.front(), box_faces.front());
        REQUIRE(p);
        CHECK_FALSE(p->parts());
        REQUIRE_EQ(p->crossing.size(), 1u);
        CHECK_EQ(p->crossing.front()[0], 0.5);
        expect_volume(m, unwrap(fuse(m, block, box)), 8. + v_box - overlap);
        expect_volume(m, unwrap(common(m, block, box)), overlap);
        expect_volume(m, unwrap(cut(m, block, box)), 8. - overlap);
    }
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

TEST(tests_brep_boolean, section_edges_meeting_at_a_vertex)
{
    // the cone of cone.stp moved by (6.7, 0, 0.3): its base disk overlaps the bottom of the
    // cylinder by 0.3 near its side. Section edges of a face meeting at a vertex cross there a
    // little inside one of them: split_face took the point for a crossing inside, and failed
    // (`Disconnected`) with the cuts split at a vertex of their own
    Model<T> m;
    const auto cyl = read(m, "cylinder.stp"); // R 5, z 0..10
    const auto cone = copy_shape(m, read(m, "cone.stp"), Placement<T>{.t = {6.7, 0., 0.3}});
    const T vc = 250. * pi, vk = volume(m, cone);
    const auto c = unwrap(common(m, cyl, cone));
    expect_valid(m, c.shape);
    const T inter = volume(m, c.shape);
    ASSERT_GT(inter, 0.05);
    ASSERT_LT(inter, 0.1);
    expect_volume(m, unwrap(cut(m, cyl, cone)), vc - inter);
    expect_volume(m, unwrap(fuse(m, cyl, cone)), vc + vk - inter);
}

TEST(tests_brep_boolean, cavity_inside_a_ring)
{
    // a cube with a cavity, inside a thin square ring around it: the cavity was given to the
    // smallest solid whose box holds it, the ring (cube 8, ring 0.488, both valid)
    Model<T> m;
    const auto cube = unwrap(cut(m, block(m, {0., 0., 0.}, {2., 2., 2.}), block(m, {0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}))).shape;
    const auto ring = unwrap(cut(m, block(m, {-0.6, -0.6, 0.4}, {2.6, 2.6, 1.6}), block(m, {-0.5, -0.5, 0.3}, {2.5, 2.5, 1.7}))).shape;
    const T v_ring = (3.2 * 3.2 - 9.) * 1.2;
    const auto r = unwrap(fuse(m, cube, ring));
    expect_volume(m, r, 7. + v_ring, 2);
    for (auto s : r.solids)
    {
        const T v = volume(m, ShapeId{s});
        if (m.solid(s).voids.empty())
            ASSERT_NEAR(v, v_ring, 2e-3);
        else
            ASSERT_NEAR(v, 7., 2e-3);
    }
}

TEST(tests_brep_boolean, failure_leaves_the_model)
{
    // an operand with an edge on a curve that is neither a primitive nor a B-spline: copy_shape
    // threw, now it moves it exactly (`MovedCurve`), and the fuse goes through
    Model<T> m;
    const auto odd = read(m, "box.stp"); // 0..2 x 0..3 x 0..4
    const auto e = explore<EdgeId>(m, odd).front();
    const auto line = m.edge(e).curve;
    m.edge(e).curve = std::make_shared<brep_tests::AnalyticCurve>([line](T t) { return line->value(t); }, line->bounds());
    const auto b = block(m, {1., 1., 1.}, {3., 4., 5.});
    expect_volume(m, unwrap(fuse(m, odd, b)), 24. + 24. - 6.);
    // a failure after the copies of both operands: the cylinder lying on a cube with B-spline
    // faces (its edge on such a curve too), too few pairs allowed to their intersection; whatever
    // the order of the operands, everything the operation added is erased
    const auto cube = ShapeId{box(m, {-1., -1., 0.}, 4.)};
    const point<T, 3> Y{0., 1., 0.}, Z{0., 0., 1.};
    const auto cyl = copy_shape(m, read(m, "cylinder.stp"), Placement<T>::frame({1., -3.5, 9.}, Z, cross(Y, Z), Y));
    const auto ce = explore<EdgeId>(m, cyl).front();
    const auto c0 = m.edge(ce).curve;
    m.edge(ce).curve = std::make_shared<brep_tests::AnalyticCurve>([c0](T t) { return c0->value(t); }, c0->bounds());
    BooleanOptions<T> few;
    few.intersect.max_pairs = 4;
    const auto before = counts(m);
    for (const auto &[x, y] : {std::pair{cube, cyl}, std::pair{cyl, cube}})
    {
        auto r = fuse(m, x, y, few);
        if (r)
            FAIL("expected NonGeneric");
        ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
        ASSERT_EQ(counts(m), before);
    }
}

TEST(tests_brep_boolean, faces_a_few_tolerances_apart)
{
    // faces within the confusion radius of each other (brep_confusion_factor tol, as their
    // vertices) coincide. box.stp and its copy moved by (-1, 0.5, -d) have their tops, and their
    // bottoms, d apart: the faces coincided within tol only, and the pieces of the section along a
    // part of an edge at d were dropped (1.5 tol: "a section edge does not split a face") or cut
    // the faces at a sliver from it (2.5 and 5 tol: "the cuts and the boundary do not form closed
    // loops"). As at d = 0 now, the edges imprinted within d of the faces they lie on
    constexpr T tol = brep_default_tolerance<T>;
    for (T d : {1.5 * tol, 2.5 * tol, 5. * tol})
    {
        INFO("d = " << d / tol << " tol");
        Model<T> m;
        const auto a = read(m, "box.stp"); // 0..2 x 0..3 x 0..4
        const auto b = copy_shape(m, a, Placement<T>{.t = {-1., 0.5, -d}});
        const T both = 2.5 * (4. - d);
        const auto f = unwrap(fuse(m, a, b));
        expect_volume(m, f, 48. - both);
        ASSERT_EQ(explore<FaceId>(m, f.shape).size(), 14u);
        const auto c = unwrap(common(m, a, b));
        expect_volume(m, c, both);
        ASSERT_EQ(explore<FaceId>(m, c.shape).size(), 6u);
        expect_volume(m, unwrap(cut(m, a, b)), 24. - both);
        // stacked, overlapping by d: face against face, one solid, nothing in common (refused: a
        // shell of the result open at 1.5 tol, the cuts not closing beyond)
        const auto up = copy_shape(m, a, Placement<T>{.t = {0.5, 0.5, 4. - d}});
        const auto s = unwrap(fuse(m, a, up));
        expect_volume(m, s, 48.);
        ASSERT_EQ(explore<FaceId>(m, s.shape).size(), 12u);
        ASSERT_TRUE(unwrap(common(m, a, up)).solids.empty());
        expect_volume(m, unwrap(cut(m, a, up)), 24.);
    }
}

TEST(tests_brep_boolean, coincident_within_the_confusion_radius)
{
    // faces a few tolerances apart, of every kind of coincidence test (refused on the previous
    // version from 1.5 tol)
    constexpr T tol = brep_default_tolerance<T>, d = 5. * tol;
    {
        // bilinear B-spline cubes 0..2 and 1..3 x 1..3 x d..2 + d (points of the faces projected)
        Model<T> m;
        const auto a = ShapeId{box(m, {0., 0., 0.}, 2.)}, b = ShapeId{box(m, {1., 1., d}, 2.)};
        expect_volume(m, unwrap(fuse(m, a, b)), 14.);
        expect_volume(m, unwrap(common(m, a, b)), 2.);
        expect_volume(m, unwrap(cut(m, a, b)), 6.);
    }
    {
        // box.stp and its copy moved by (-1, 0.5, 0), turned by d rad about the line y of (0.5, ., 4):
        // the tops cross at x = 0.5, d / 2 apart at the ends of their overlap, but the planes are
        // not parallel within the confusion radius across their region (planes by their largest gap)
        Model<T> m;
        const auto a = read(m, "box.stp");
        const point<T, 3> X{std::cos(d), 0., -std::sin(d)}, Y{0., 1., 0.}, Z{std::sin(d), 0., std::cos(d)}, P{0.5, 0.5, 4.};
        const auto b = copy_shape(m, a, Placement<T>::frame(P - 1.5 * X - 4. * Z, X, Y, Z));
        expect_volume(m, unwrap(fuse(m, a, b)), 38.);
        expect_volume(m, unwrap(common(m, a, b)), 10.);
        expect_volume(m, unwrap(cut(m, a, b)), 14.);
    }
    {
        // cylinders R 5, z 0..10 and 5..15, axes d apart (closed forms, confirmed by points); the
        // circles of each imprinted on the side of the other, within d of it: the deviation of
        // their pcurves, at their seam, was measured between samples, under the gap there
        Model<T> m;
        const auto a = read(m, "cylinder.stp");
        const auto b = copy_shape(m, a, Placement<T>{.t = {d, 0., 5.}});
        expect_volume(m, unwrap(fuse(m, a, b)), 25. * pi * 15.);
        expect_volume(m, unwrap(common(m, a, b)), 25. * pi * 5.);
        expect_volume(m, unwrap(cut(m, a, b)), 25. * pi * 5.);
    }
    {
        // the same with closed rational B-spline sides, R 1, z 0..2 and 1..3
        Model<T> m;
        const auto a = nurbs_cylinder(m, 1., 2.);
        const auto b = copy_shape(m, a, Placement<T>{.t = {d, 0., 1.}});
        expect_volume(m, unwrap(fuse(m, a, b)), 3. * pi);
        expect_volume(m, unwrap(common(m, a, b)), pi);
        expect_volume(m, unwrap(cut(m, a, b)), pi);
    }
}

TEST(tests_brep_boolean, dead_reference_keeps_its_reason)
{
    // an operand referring to an erased vertex: InvalidId with the vertex (it was NonGeneric, the
    // reason lost), and the model left as it was
    Model<T> m;
    auto a = box(m, {0., 0., 0.}, 1.);
    auto b = box(m, {0.5, 0.5, 0.5}, 1.);
    const auto v = explore<VertexId>(m, b).front();
    m.erase(v);
    const auto faces = m.count<ShapeType::Face>();
    auto r = boolean(m, BooleanOp::Fuse, a, b);
    ASSERT_FALSE(r.has_value());
    ASSERT_TRUE(r.error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(r.error().shapes.front() == ShapeId{v});
    ASSERT_EQ(r.error().message.rfind("boolean: ", 0), 0u);
    ASSERT_EQ(m.count<ShapeType::Face>(), faces);
}
