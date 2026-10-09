#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include "tests_brep_geom.h"
#include "tests_brep_helpers.h"

#include <sstream>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using namespace brep_tests;

    std::string dump(const CheckReport &r)
    {
        std::ostringstream os;
        for (const auto &e : r.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        return os.str();
    }

    ShellId sewn_box(Model<T> &m, point<T, 3> origin = {}, T size = 1.)
    {
        auto r = unwrap(sew(m, box_faces(m, origin, size)));
        REQUIRE(r.shells.size() == 1);
        return r.shells.front();
    }

    std::shared_ptr<Surface<T, 3>> plane_z(T z)
    {
        return quad({-2., -2., z}, {2., -2., z}, {-2., 2., z}, {2., 2., z});
    }

    std::vector<Orientation> uses(const Model<T> &m, ShellId s)
    {
        std::vector<Orientation> o;
        for (const auto &fu : m.shell(s).faces)
            o.push_back(fu.orient);
        return o;
    }
}

TEST(tests_brep_solid, box_from_faces_to_solid)
{
    Model<T> m;
    auto sh = sewn_box(m);
    const T v0 = signed_volume(m, sh);
    ASSERT_NEAR(std::abs(v0), 1., 1e-12); // planar faces: the midpoint rule is exact

    auto so = unwrap(make_solid(m, sh));
    ASSERT_TRUE(m.solid(so).outer == sh);
    ASSERT_TRUE(m.solid(so).voids.empty());
    ASSERT_NEAR(signed_volume(m, sh), 1., 1e-12); // turned outwards when needed
    auto c = check(m, so);
    INFO(dump(c));
    ASSERT_TRUE(c.ok());

    // the hand-made box of the model tests is already outward
    auto b = brep_tests::make_box();
    ASSERT_NEAR(signed_volume(b.m, b.shell), 1., 1e-12);
}

TEST(tests_brep_solid, inward_shell_is_turned_and_checked)
{
    Model<T> m;
    auto sh = sewn_box(m);
    auto so = unwrap(make_solid(m, sh));
    // reverse every face use: still closed and orientable, but inward
    for (auto &fu : m.shell(sh).faces)
        fu.orient = reverse(fu.orient);
    ASSERT_NEAR(signed_volume(m, sh), -1., 1e-12);
    ASSERT_TRUE(check(m, so).has(Issue::SolidNotOutward, so));
    ASSERT_TRUE(check(m, so, CheckOptions{.geometry = false}).ok()); // a geometric issue

    auto so2 = unwrap(make_solid(m, sh)); // make_solid turns it back
    ASSERT_NEAR(signed_volume(m, sh), 1., 1e-12);
    ASSERT_TRUE(check(m, so2).ok());
    ASSERT_TRUE(check(m, so).ok()); // same shell
}

TEST(tests_brep_solid, curved_solids_volumes)
{
    // sphere and torus: closed by a single natural face
    {
        Model<T> m;
        auto sh = m.add(Shell{{FaceUse{unwrap(make_face(m, sphere(2.)))}}});
        auto so = unwrap(make_solid(m, sh));
        ASSERT_NEAR(signed_volume(m, sh), 4. / 3. * pi * 8., 1e-3 * 4. / 3. * pi * 8.);
        ASSERT_TRUE(check(m, so).ok());
    }
    {
        Model<T> m;
        auto sh = m.add(Shell{{FaceUse{unwrap(make_face(m, torus(3., 1.)))}}});
        auto so = unwrap(make_solid(m, sh));
        ASSERT_NEAR(signed_volume(m, sh), 2. * pi * pi * 3., 1e-3 * 2. * pi * pi * 3.);
        ASSERT_TRUE(check(m, so).ok());
    }
    // cylinder closed by two trimmed disks
    {
        Model<T> m;
        const T h = 2.;
        auto lateral = unwrap(make_face(m, cylinder(1., h)));
        auto disk = [&](T z) {
            auto c = std::make_shared<AnalyticCurve>([z](T t) { return point<T, 3>{std::cos(t), std::sin(t), z}; }, std::array<T, 2>{0., 2. * pi});
            return unwrap(make_face(m, plane_z(z), unwrap(make_wire(m, std::vector{unwrap(make_edge(m, c))}))));
        };
        auto r = unwrap(sew(m, std::vector{lateral, disk(0.), disk(h)}, {.tol = 1e-5}));
        auto so = unwrap(make_solid(m, r.shells.front()));
        ASSERT_NEAR(signed_volume(m, r.shells.front()), pi * h, 3e-3 * pi * h); // trimmed disks: cell-center membership, ~1e-3 at n = 64
        auto c = check(m, so);
        INFO(dump(c));
        ASSERT_TRUE(c.ok());
    }
}

TEST(tests_brep_solid, hollow_box_with_a_cavity)
{
    Model<T> m;
    auto outer = sewn_box(m, {0., 0., 0.}, 3.);
    auto cavity = sewn_box(m, {1., 1., 1.}, 1.);
    auto so = unwrap(make_solid(m, outer, std::vector{cavity}));
    ASSERT_EQ(m.solid(so).voids.size(), 1);
    ASSERT_NEAR(signed_volume(m, outer), 27., 1e-9);
    ASSERT_NEAR(signed_volume(m, cavity), -1., 1e-12); // normals into the cavity
    ASSERT_NEAR(signed_volume(m, outer) + signed_volume(m, cavity), 26., 1e-9);
    ASSERT_TRUE(check(m, so).ok());

    for (auto &fu : m.shell(cavity).faces)
        fu.orient = reverse(fu.orient);
    ASSERT_TRUE(check(m, so).has(Issue::VoidNotInward, cavity));
}

TEST(tests_brep_solid, errors_leave_model_unchanged)
{
    Model<T> m;
    auto sh = sewn_box(m);
    const auto before = uses(m, sh);

    // open shell
    auto faces = box_faces(m, {5., 0., 0.});
    faces.pop_back();
    auto open = unwrap(sew(m, faces)).shells.front();
    ASSERT_TRUE(make_solid(m, open).error().code == BuildErrc::ShellNotClosed);

    // closed but one face use flipped: not orientable
    Shell bad = m.shell(sh);
    bad.faces[0].orient = reverse(bad.faces[0].orient);
    auto bad_id = m.add(bad);
    ASSERT_TRUE(make_solid(m, bad_id).error().code == BuildErrc::ShellNotOrientable);

    // dead shell, shell given twice, a good outer with a bad cavity
    ASSERT_TRUE(make_solid(m, ShellId{99}).error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(make_solid(m, sh, std::vector{sh}).error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(make_solid(m, sh, std::vector{open}).error().code == BuildErrc::ShellNotClosed);

    ASSERT_TRUE(uses(m, sh) == before); // the outer shell was not turned
    ASSERT_EQ(m.count<ShapeType::Solid>(), 0);
}

TEST(tests_brep_solid, compound)
{
    Model<T> m;
    auto so = unwrap(make_solid(m, sewn_box(m)));
    auto free_face = unwrap(make_face(m, plane_z(5.)));
    auto inner = unwrap(make_compound(m, std::vector<ShapeId>{so}));
    auto c = unwrap(make_compound(m, std::vector<ShapeId>{inner, free_face}));
    ASSERT_EQ(explore<FaceId>(m, c).size(), 7);
    ASSERT_EQ(explore<SolidId>(m, c).size(), 1);
    ASSERT_EQ(explore<CompoundId>(m, c).size(), 2);
    ASSERT_TRUE(check(m, c).ok());
    ASSERT_TRUE(make_compound(m, std::vector<ShapeId>{so, FaceId{999}}).error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(make_compound(m, std::vector<ShapeId>{}).has_value()); // an empty compound is allowed
}
