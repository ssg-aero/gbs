#include <doctest_gtest.hpp>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <gbs-brep/brep>
#include "tests_brep_geom.h"

#include <algorithm>
#include <random>
#include <sstream>

using namespace gbs;
using namespace gbs::brep;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using namespace brep_tests;
    const T tol = 1e-6;

    std::string dump(const CheckReport &r)
    {
        std::ostringstream os;
        for (const auto &e : r.entries)
            os << to_string(e.issue) << " on #" << shape_index(e.shape) << " " << e.detail << "\n";
        return os.str();
    }

    // Bilinear patch through four corners, S(0,0) = a, S(1,0) = b, S(0,1) = c, S(1,1) = d.
    std::shared_ptr<Surface<T, 3>> quad(point<T, 3> a, point<T, 3> b, point<T, 3> c, point<T, 3> d)
    {
        std::vector<T> k{0., 0., 1., 1.};
        return std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{a, b, c, d}, k, k, 1, 1);
    }

    point<T, 3> corner(unsigned i, point<T, 3> origin = {})
    {
        return origin + point<T, 3>{T(i & 1), T((i >> 1) & 1), T((i >> 2) & 1)};
    }

    // Six independent natural faces of a unit box; every face has the same (u, v) axes
    // order on both sides of an axis, so half of the normals point inwards.
    std::vector<FaceId> box_faces(Model<T> &m, point<T, 3> origin = {})
    {
        std::vector<FaceId> f;
        for (unsigned axis = 0; axis < 3; ++axis)
            for (unsigned side = 0; side < 2; ++side)
            {
                const unsigned a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
                auto idx = [&](unsigned u, unsigned v) { return (side << axis) | (u << a1) | (v << a2); };
                f.push_back(unwrap(make_face(m, quad(corner(idx(0, 0), origin), corner(idx(1, 0), origin),
                                                     corner(idx(0, 1), origin), corner(idx(1, 1), origin)))));
            }
        return f;
    }

    std::shared_ptr<Surface<T, 3>> plane_z(T z)
    {
        return quad({-2., -2., z}, {2., -2., z}, {-2., 2., z}, {2., 2., z});
    }

    std::shared_ptr<Curve<T, 2>> seg2(point<T, 2> a, point<T, 2> b)
    {
        return std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>(a, b));
    }
}

TEST(tests_brep_sew, box_from_six_independent_faces)
{
    Model<T> m;
    auto faces = box_faces(m);
    std::mt19937 rng{42};
    std::ranges::shuffle(faces, rng);
    ASSERT_EQ(m.count<ShapeType::Edge>(), 24);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 24);

    auto r = unwrap(sew(m, faces, {.tol = tol}));
    ASSERT_EQ(r.shells.size(), 1);
    ASSERT_EQ(r.merged.size(), 12);
    ASSERT_TRUE(r.free_edges.empty());
    ASSERT_TRUE(r.rejected.empty());
    ASSERT_TRUE(r.orientable);
    ASSERT_EQ(m.count<ShapeType::Edge>(), 12);
    ASSERT_EQ(m.count<ShapeType::Vertex>(), 8);

    const auto sh = r.shells.front();
    ASSERT_TRUE(m.shell(sh).closed);
    ASSERT_TRUE(is_closed(m, sh));
    ASSERT_TRUE(is_manifold(m, sh));
    ASSERT_TRUE(is_orientable(m, sh));
    // half of the face uses had to be reversed
    size_t reversed = 0;
    for (const auto &fu : m.shell(sh).faces)
        reversed += fu.orient == Orientation::Reversed;
    ASSERT_EQ(reversed, 3); // three normals inwards, relative to the first face
    auto c = check(m, sh);
    INFO(dump(c));
    ASSERT_TRUE(c.ok());
}

TEST(tests_brep_sew, cylinder_closed_by_two_disks)
{
    Model<T> m;
    const T h = 2.;
    auto lateral = unwrap(make_face(m, cylinder(1., h)));
    // bottom: angular circle, same parametrization as the cylinder's iso
    auto bottom_circle = std::make_shared<AnalyticCurve>([](T t) { return point<T, 3>{std::cos(t), std::sin(t), 0.}; },
                                                         std::array<T, 2>{0., 2. * pi});
    auto bottom = unwrap(make_face(m, plane_z(0.), unwrap(make_wire(m, std::vector{unwrap(make_edge(m, bottom_circle))}))));
    // top: rational circle, parameter in [0, 1] not proportional to the angle
    auto rc = build_circle<T, 3>(1.);
    auto poles = rc.poles();
    for (auto &p : poles)
        p[2] += p[3] * h;
    auto top_circle = std::make_shared<BSCurveRational<T, 3>>(poles, rc.knotsFlats(), rc.degree());
    auto top = unwrap(make_face(m, plane_z(h), unwrap(make_wire(m, std::vector{unwrap(make_edge(m, top_circle))}))));

    auto r = unwrap(sew(m, std::vector{top, lateral, bottom}, {.tol = 1e-5}));
    ASSERT_EQ(r.shells.size(), 1);
    ASSERT_EQ(r.merged.size(), 2);
    ASSERT_TRUE(r.free_edges.empty());
    ASSERT_TRUE(r.orientable);
    const auto sh = r.shells.front();
    ASSERT_TRUE(m.shell(sh).closed);
    ASSERT_EQ(explore<EdgeId>(m, sh).size(), 3);   // seam + two circles
    ASSERT_EQ(explore<VertexId>(m, sh).size(), 2);
    auto c = check(m, sh);
    INFO(dump(c));
    ASSERT_TRUE(c.ok()); // includes SameParameter of the rebuilt pcurves
    // the rebuilt pcurves of the disks follow the circle once, without folding back at the
    // closing point: the (u,v) area is the disk's (plane_z maps [-2,2]^2 onto [0,1]^2)
    for (auto f : {top, bottom})
        ASSERT_NEAR(uv_signed_area(m, m.face(f).wires.front(), 4096), pi / 16., 1e-6);
}

TEST(tests_brep_sew, components_and_open_shell)
{
    Model<T> m;
    auto a = box_faces(m);
    auto b = box_faces(m, {3., 0., 0.});
    b.pop_back(); // second box open
    std::vector<FaceId> all = a;
    all.insert(all.end(), b.begin(), b.end());
    std::mt19937 rng{7};
    std::ranges::shuffle(all, rng);

    auto r = unwrap(sew(m, all));
    ASSERT_EQ(r.shells.size(), 2);
    size_t closed = 0, open = 0;
    for (auto sh : r.shells)
        (m.shell(sh).closed ? closed : open) += 1;
    ASSERT_EQ(closed, 1);
    ASSERT_EQ(open, 1);
    ASSERT_EQ(r.free_edges.size(), 4); // rim of the open box
    for (auto sh : r.shells)
        ASSERT_TRUE(check(m, sh).ok());
}

TEST(tests_brep_sew, t_junction_stays_free)
{
    Model<T> m;
    // A spans x in [0, 2]; B1 and B2 sit on top of it, each over half of A's top edge
    auto fa = unwrap(make_face(m, quad({0., 0., 0.}, {2., 0., 0.}, {0., 1., 0.}, {2., 1., 0.})));
    auto fb1 = unwrap(make_face(m, quad({0., 1., 0.}, {1., 1., 0.}, {0., 2., 0.}, {1., 2., 0.})));
    auto fb2 = unwrap(make_face(m, quad({1., 1., 0.}, {2., 1., 0.}, {1., 2., 0.}, {2., 2., 0.})));
    auto r = unwrap(sew(m, std::vector{fa, fb1, fb2}));
    ASSERT_EQ(r.merged.size(), 1);       // only B1 | B2
    ASSERT_EQ(r.free_edges.size(), 10);  // 4 + 3 + 3: A's top edge is not cut
    ASSERT_EQ(r.shells.size(), 2);       // A alone, B1 + B2
}

TEST(tests_brep_sew, gap_absorbed_in_tolerances)
{
    Model<T> m;
    const T gap = 0.4e-3;
    auto f1 = unwrap(make_face(m, quad({0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.})));
    auto f2 = unwrap(make_face(m, quad({1., 0., gap}, {2., 0., 0.}, {1., 1., gap}, {2., 1., 0.})));
    ASSERT_EQ(unwrap(sew(m, std::vector{f1, f2}, {.tol = 1e-4})).merged.size(), 0); // gap above tol
    auto r = unwrap(sew(m, std::vector{f1, f2}, {.tol = 1e-3}));
    ASSERT_EQ(r.merged.size(), 1);
    const auto e = r.merged.front().first;
    ASSERT_GE(m.edge(e).tol, gap * (1. - 1e-9));
    for (auto v : {m.edge(e).v1, m.edge(e).v2})
    {
        ASSERT_GE(m.vertex(v).tol, m.edge(e).tol);
        ASSERT_GE(m.vertex(v).tol, gap * (1. - 1e-9)); // the absorbed corner lies in the ball
    }
    ASSERT_TRUE(check(m, r.shells.front()).ok());
}

TEST(tests_brep_sew, same_ends_different_curves_rejected)
{
    Model<T> m;
    // plane z = 0 whose parameters are x and y on [-2, 2]^2
    std::vector<T> k{-2., -2., 2., 2.};
    auto srf = std::make_shared<BSSurface<T, 3>>(points_vector<T, 3>{{-2., -2., 0.}, {2., -2., 0.}, {-2., 2., 0.}, {2., 2., 0.}}, k, k, 1, 1);
    // A: segment + arc bulging up; B: segment + arc bulging down; both share the segment and the ends of the arcs
    auto arc = [](T bulge) {
        return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{1., 0.}, {0.5, bulge}, {0., 0.}}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2);
    };
    auto fa = unwrap(make_face(m, srf, std::vector<std::shared_ptr<Curve<T, 2>>>{seg2({0., 0.}, {1., 0.}), arc(0.5)}));
    auto fb = unwrap(make_face(m, srf, std::vector<std::shared_ptr<Curve<T, 2>>>{seg2({0., 0.}, {1., 0.}), arc(-0.5)}));
    auto r = unwrap(sew(m, std::vector{fa, fb}));
    ASSERT_EQ(r.merged.size(), 1); // the segments
    ASSERT_EQ(r.free_edges.size(), 2); // the arcs
    // all four edges join (0,0) and (1,0): the five other pairs, including the two arcs and the
    // segment / arc pairs inside a face, have matching ends but curves apart
    ASSERT_EQ(r.rejected.size(), 5);
    const std::pair<EdgeId, EdgeId> arcs{r.free_edges[0], r.free_edges[1]};
    ASSERT_TRUE(std::ranges::find(r.rejected, arcs) != r.rejected.end());
}

TEST(tests_brep_sew, moebius_strip_is_not_orientable)
{
    Model<T> m;
    const T R = 3., w = 0.5;
    auto P = [&](T th, T s) {
        return point<T, 3>{(R + s * std::cos(th / 2)) * std::cos(th), (R + s * std::cos(th / 2)) * std::sin(th), s * std::sin(th / 2)};
    };
    std::vector<FaceId> faces;
    const int n = 6;
    for (int i = 0; i < n; ++i)
    {
        const T t0 = 2. * pi * i / n, t1 = 2. * pi * (i + 1) / n;
        faces.push_back(unwrap(make_face(m, quad(P(t0, -w), P(t1, -w), P(t0, w), P(t1, w)))));
    }
    auto r = unwrap(sew(m, faces, {.tol = 1e-9}));
    ASSERT_EQ(r.shells.size(), 1);
    ASSERT_EQ(r.merged.size(), n); // n - 1 inner seams + the twisted closing one
    ASSERT_FALSE(r.orientable);
    ASSERT_FALSE(is_orientable(m, r.shells.front()));
    ASSERT_TRUE(check(m, r.shells.front()).has(Issue::ShellNotOrientable, r.shells.front()));
}

TEST(tests_brep_sew, errors_leave_model_unchanged)
{
    Model<T> m;
    auto faces = box_faces(m);
    const auto ne = m.count<ShapeType::Edge>();
    ASSERT_TRUE(sew(m, std::vector<FaceId>{}).error().code == BuildErrc::EmptyInput);
    ASSERT_TRUE(sew(m, faces, {.tol = 0.}).error().code == BuildErrc::InvalidTolerance);
    ASSERT_TRUE(sew(m, std::vector{faces[0], FaceId{99}}).error().code == BuildErrc::InvalidId);
    ASSERT_TRUE(sew(m, std::vector{faces[0], faces[0]}).error().code == BuildErrc::InvalidId);
    ASSERT_EQ(m.count<ShapeType::Edge>(), ne);
    ASSERT_EQ(m.count<ShapeType::Shell>(), 0);
}
