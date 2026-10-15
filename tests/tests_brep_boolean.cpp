#include <doctest_gtest.hpp>
#include <gbs-brep/boolean.h>
#include <gbs-brep/check.h>
#include <gbs-io/step/reader.h>
#include "tests_brep_geom.h"

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

TEST(tests_brep_boolean, non_generic)
{
    // boxes sharing a face plane: coincident faces are not handled yet
    Model<T> m;
    const auto a = box(m, {0., 0., 0.}, 2.), b = box(m, {1., 1., 0.}, 2.);
    auto r = fuse(m, a, b);
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::NonGeneric);
    // the operands are kept
    expect_valid(m, ShapeId{a});
    expect_valid(m, ShapeId{b});
    // errors
    ASSERT_FALSE(fuse(m, ShapeId{a}, ShapeId{SolidId{999}}).has_value());
    const auto face = explore<FaceId>(m, ShapeId{a}).front();
    ASSERT_FALSE(fuse(m, ShapeId{a}, ShapeId{face}).has_value());
}
