#include <doctest_gtest.hpp>
#include <gbs-brep/split.h>
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

    std::shared_ptr<Plane<T>> plane(const point<T, 3> &o, const point<T, 3> &n, T half = 100.)
    {
        const auto z = intersect_detail::unit<T>(n);
        return std::make_shared<Plane<T>>(intersect_detail::frame_of<T>(o, z, intersect_detail::orthogonal<T>(z)), std::array<T, 4>{-half, half, -half, half});
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

    // volume by the sampling of the faces (signed_volume): about 1e-4 relative on trimmed faces
    T volume(const Model<T> &m, const ShapeId &s)
    {
        T v{};
        for (auto so : explore<SolidId>(m, s))
            v += signed_volume(m, m.solid(so).outer, 256);
        return v;
    }

    // area of the faces by the surface element on their interior points, summed over a grid
    T area(const Model<T> &m, FaceId f)
    {
        const auto &s = *m.face(f).surface;
        const auto b = classify_detail::face_uv_box(m, f);
        constexpr int n = 60;
        T a{};
        const T du = (b[1] - b[0]) / n, dv = (b[3] - b[2]) / n;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
            {
                const point<T, 2> uv{b[0] + (i + 0.5) * du, b[2] + (j + 0.5) * dv};
                if (classify(m, f, uv) == Location::In)
                    a += norm(cross(s.value(uv[0], uv[1], 1, 0), s.value(uv[0], uv[1], 0, 1))) * du * dv;
            }
        return a;
    }
}

TEST(tests_brep_split, plate_through_its_hole)
{
    // the plate 4 x 2 in with a hole of radius 0.5 in at (1, 1), cut by the plane x = 1 in: two
    // faces, the hole cut in two halves
    Model<T> m;
    const auto plate = read(m, "plate_inch.stp");
    const auto f0 = explore<FaceId>(m, plate).front();
    const T a0 = area(m, f0);
    auto r = unwrap(split(m, plate, plane({25.4, 0., 0.}, {1., 0., 0.})));
    ASSERT_EQ(r.faces.size(), 1u);
    ASSERT_EQ(r.faces[0].second.size(), 2u);
    ASSERT_TRUE(r.dangling.empty());
    const auto faces = explore<FaceId>(m, plate);
    ASSERT_EQ(faces.size(), 2u);
    T a{};
    for (auto f : faces)
    {
        ASSERT_EQ(m.face(f).wires.size(), 1u); // no hole left in either piece
        a += area(m, f);
    }
    ASSERT_NEAR(a, a0, 1e-2 * a0);
    expect_valid(m, plate);
}

TEST(tests_brep_split, solids_stay_closed)
{
    // a box cut by an oblique plane: the crossed faces are split, the shell stays closed
    Model<T> m;
    const auto box = read(m, "box.stp"); // 2 x 3 x 4
    const T v0 = volume(m, box);
    auto r = unwrap(split(m, box, plane({1., 1.5, 2.}, {0.3, 0.2, 1.})));
    ASSERT_EQ(r.faces.size(), 4u); // the four sides
    ASSERT_EQ(explore<FaceId>(m, box).size(), 10u);
    expect_valid(m, box);
    ASSERT_NEAR(volume(m, box), v0, 1e-3 * v0);

    // the cylinder cut across: the lateral face in two, across its seam
    const auto cyl = read(m, "cylinder.stp");
    const T vc = volume(m, cyl);
    auto c = unwrap(split(m, cyl, plane({0., 0., 5.}, {0., 0., 1.})));
    ASSERT_EQ(c.faces.size(), 1u);
    ASSERT_EQ(c.faces[0].second.size(), 2u);
    ASSERT_EQ(explore<FaceId>(m, cyl).size(), 4u);
    expect_valid(m, cyl);
    ASSERT_NEAR(volume(m, cyl), vc, 1e-3 * vc);
    // obliquely: an ellipse cut at the seam into two edges
    const auto cyl2 = read(m, "cylinder.stp");
    unwrap(split(m, cyl2, plane({0., 0., 5.}, {0., 0.5, 1.})));
    ASSERT_EQ(explore<FaceId>(m, cyl2).size(), 4u);
    expect_valid(m, cyl2);
    ASSERT_NEAR(volume(m, cyl2), vc, 1e-3 * vc);
}

TEST(tests_brep_split, holes_and_dangling_cuts)
{
    // a box pierced by a thin cylinder: the top and the bottom get a disc and a holed face
    Model<T> m;
    const auto box = read(m, "box.stp");
    const T v0 = volume(m, box);
    const auto tool = std::make_shared<CylindricalSurface<T>>(0.5, ax2<T, 3>{{{1., 1.5, 0.}, {0., 0., 1.}, {1., 0., 0.}}}, std::array<T, 4>{0., 2 * pi, -1., 5.});
    auto r = unwrap(split(m, box, tool));
    ASSERT_EQ(r.faces.size(), 2u);
    for (const auto &[f, pieces] : r.faces)
    {
        ASSERT_EQ(pieces.size(), 2u);
        std::size_t holed{};
        for (auto p : pieces)
            holed += m.face(p).wires.size() == 2;
        ASSERT_EQ(holed, 1u);
    }
    ASSERT_EQ(explore<FaceId>(m, box).size(), 8u);
    expect_valid(m, box);
    ASSERT_NEAR(volume(m, box), v0, 1e-3 * v0);

    // a cut from the boundary into the middle of a face: dangling, the face is kept
    using namespace brep_tests;
    Model<T> mp;
    const auto f = unwrap(make_face(mp, quad({0., 0., 0.}, {2., 0., 0.}, {0., 2., 0.}, {2., 2., 0.})));
    const auto e = unwrap(make_edge(mp, point<T, 3>{1., 0., 0.}, point<T, 3>{1., 1., 0.}));
    const auto &ed = mp.edge(e);
    auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0.5, 0.}, {0.5, 0.5}}, std::vector<T>{ed.u1, ed.u1, ed.u2, ed.u2}, 1);
    auto d = unwrap(split_face(mp, f, std::vector<FaceCut<T>>{{e, pc}}));
    ASSERT_EQ(d.faces.size(), 1u);
    ASSERT_EQ(d.faces[0], f);
    ASSERT_EQ(d.dangling.size(), 1u);
}

TEST(tests_brep_split, by_another_shape)
{
    // boxes 0..2 and 1..3: the faces of the first one crossed by the second one are split
    using namespace brep_tests;
    Model<T> m;
    auto a = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {0., 0., 0.}, 2.))).shells.front()));
    auto b = unwrap(make_solid(m, unwrap(sew(m, box_faces(m, {1., 1., 1.}, 2.))).shells.front()));
    const T va = volume(m, a);
    auto r = unwrap(split(m, a, b));
    ASSERT_EQ(r.faces.size(), 3u); // the faces x = 2, y = 2, z = 2
    ASSERT_EQ(explore<FaceId>(m, a).size(), 9u);
    expect_valid(m, a);
    expect_valid(m, b);
    ASSERT_NEAR(volume(m, a), va, 1e-3 * va);
    // errors
    ASSERT_FALSE(split_face(m, FaceId{999}, std::vector<FaceCut<T>>{}).has_value());
}
