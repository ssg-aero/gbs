// Exact writing and reading of a BREP model in JSON (gbs-io/brep_json.h), exact equality of two
// models (gbs-brep/equal.h). Design: docs/sources/design/io_brep_json.md.
#include <doctest_gtest.hpp>
#include <gbs-io/brep_json.h>
#include <gbs-io/step/reader.h>
#include <gbs-brep/brep>
#include <gbs-brep/equal.h>
#include <blade_fixture.h>
#include "tests_brep_geom.h"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <filesystem>
#include <format>
#include <future>
#include <limits>
#include <numbers>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using namespace gbs::brep;
using namespace gbs::io;
using rapidjson::SizeType;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    /// b is a, entity by entity, bit for bit, sharing included (brep::equal), with the difference if not.
    template <std::floating_point R>
    void expect_equal(const Model<R> &a, const Model<R> &b)
    {
        const auto d = first_difference(a, b);
        INFO(d.value_or(""));
        REQUIRE_FALSE(d.has_value());
    }

    /// Written, read back, equal; written again: the same text; the indented text, the same model.
    template <std::floating_point R>
    auto round_trip(const Model<R> &m) -> std::string
    {
        const auto text = unwrap(write_brep(m));
        const auto back = unwrap(read_brep<R>(text));
        expect_equal(m, back);
        REQUIRE(unwrap(write_brep(back)) == text);
        const auto pretty = unwrap(write_brep(m, {.pretty = true}));
        REQUIRE(pretty != text);
        REQUIRE(unwrap(write_brep(unwrap(read_brep<R>(pretty)))) == text);
        return text;
    }

    /// The number of geometries of this type in the file.
    std::size_t count_type(const std::string &text, std::string_view type)
    {
        rapidjson::Document d;
        d.Parse(text.data(), text.size());
        std::size_t n{};
        for (const auto &g : d["geometry"].GetArray())
            n += std::string_view{g["type"].GetString()} == type;
        return n;
    }

    std::string step_dir() { return std::string(GBS_TESTS_IN_DIR) + "/step"; }

    /// A directory of its own in the temporary directory, removed with its files when the guard
    /// goes, after an assertion that fails too.
    struct TempDir
    {
        std::filesystem::path path;

        explicit TempDir(std::string_view prefix)
            : path{std::filesystem::temp_directory_path() / (std::string{prefix} + std::to_string(std::random_device{}()))}
        {
            std::filesystem::create_directories(path);
        }
        TempDir(const TempDir &) = delete;
        TempDir &operator=(const TempDir &) = delete;
        ~TempDir()
        {
            std::error_code ignored;
            std::filesystem::remove_all(path, ignored);
        }
    };

    /// The text of a model whose geometry is a line and n reversed curves, each of the one before
    /// it (geometry i at depth i + 1), the last one the curve of an edge.
    std::string chain_text(std::size_t n)
    {
        std::string t = R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":1.0,"next_part":1,"geometry":[)"
                        R"({"type":"line","dim":3,"origin":[0.0,0.0,0.0],"vector":[1.0,0.0,0.0],"range":[0.0,1.0]})";
        for (std::size_t i = 1; i <= n; ++i)
            t += std::format(R"(,{{"type":"reversed_curve","dim":3,"basis":{}}})", i - 1);
        t += std::format(R"(],"vertices":[{{"p":[0.0,0.0,0.0],"tol":1e-7}}],"edges":[{{"curve":{},"u":[0.0,1.0],"v":[0,0],"tol":1e-7,)"
                         R"("degenerate":false,"same_parameter":true}}],"wires":[],"faces":[],"shells":[],"solids":[],"compounds":[]}})",
                         n);
        return t;
    }

    // ---- models ----------------------------------------------------------------

    std::shared_ptr<BSCurve<T, 2>> bs2(T y = 0.)
    {
        return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{{0., y}, {0.3, y + 0.2}, {0.7, y - 0.1}, {1., y}}, std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, 2);
    }
    std::shared_ptr<BSCurve<T, 3>> bs3()
    {
        return std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {1., 0.3, 0.1}, {2., -0.2, 0.3}, {3., 0., 0.}},
                                               std::vector<T>{0., 0., 0., 0., 1., 1., 1., 1.}, 3);
    }
    BSCfunction<T> law() { return BSCfunction<T>{std::vector<T>{0.05, 0.1, 0.07}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2}; }
    std::shared_ptr<BSSurface<T, 3>> bss3()
    {
        points_vector<T, 3> p;
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 4; ++i)
                p.push_back({T(i), T(j), 0.1 * T(i * j) - 0.05 * T(i)});
        return std::make_shared<BSSurface<T, 3>>(p, std::vector<T>{0., 0., 0., 0.5, 1., 1., 1.}, std::vector<T>{0., 0., 1., 2., 2.}, 2, 1);
    }

    /// One object of every class the format holds, in edges, co-edges and faces, some shared,
    /// some nested: the entities hold them, the topology does not need to make sense.
    Model<T> every_geometry()
    {
        Model<T> m;
        const auto v0 = m.add(Vertex<T>{{0., 0., 0.}, 1e-7});
        const auto v1 = m.add(Vertex<T>{{1., 0., 0.}, 2e-7});
        auto edge = [&](std::shared_ptr<Curve<T, 3>> c) {
            const auto [u1, u2] = c->bounds();
            const bool finite = std::isfinite(u2 - u1);
            return m.add(Edge<T>{c, finite ? u1 : 0., finite ? u2 : 1., v0, v1, 1e-7});
        };
        // 2D curves
        auto b2 = bs2();
        std::vector<std::shared_ptr<Curve<T, 2>>> c2{
            b2,
            std::make_shared<BSCurveRational<T, 2>>(points_vector<T, 2>{{0., 0.}, {1., 1.}, {2., 0.}}, std::vector<T>{0., 0., 0., 1., 1., 1.},
                                                    std::vector<T>{1., std::sqrt(0.5), 1.}, 2),
            std::make_shared<CurveTrimmed<T, 2>>(b2, 0.2, 0.7),
            std::make_shared<CurveReversed<T, 2>>(b2),
            std::make_shared<CurveReparametrized<T, 2>>(b2, BSCfunction<T>{std::vector<T>{0., 0.4, 1.}, std::vector<T>{0., 0., 0., 2., 2., 2.}, 2}),
            std::make_shared<CurveExtended<T, 2>>(b2, true),
            std::make_shared<CurveOffset2D<T, BSCfunction<T>>>(b2, law()),
        };
        c2.push_back(std::make_shared<CurveComposite<T, 2>>(std::vector<std::shared_ptr<Curve<T, 2>>>{c2[2], c2[3]}));
        {
            points_vector<T, 2> p{{0., 0.}, {1., 0.}, {0., 1.}, {1.2, 1.1}};
            auto s2 = std::make_shared<BSSurface<T, 2>>(p, std::vector<T>{0., 0., 1., 1.}, std::vector<T>{0., 0., 1., 1.}, 1, 1);
            c2.push_back(std::make_shared<CurveOnSurface<T, 2>>(bs2(0.5), s2));
        }
        // 3D surfaces
        Frame3<T> fr{ax2<T, 3>{{{1., 2., 3.}, {0.1, 0.2, 1.}, {1., 0., 0.}}}};
        auto plane = std::make_shared<Plane<T>>(fr, std::array<T, 4>{-1., 2., -3., 4.});
        auto b3s = bss3();
        auto trimmed_u = bss3(); // U bounds an ulp inside the end knots, kept by trimU
        trimmed_u->trimU(std::nextafter(0., 1.), std::nextafter(1., 0.));
        REQUIRE(trimmed_u->bounds()[0] > 0.);
        auto nurbs_s = std::make_shared<BSSurfaceRational<T, 3>>(
            points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 1.}}, std::vector<T>{1., 0.8, 0.9, 1.2}, std::vector<T>{0., 0., 1., 1.},
            std::vector<T>{0., 0., 1., 1.}, 1, 1);
        auto revolution = std::make_shared<SurfaceOfRevolution<T>>(bs2(1.), ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, 0.25, 2.5);
        Placement<T> turn = Placement<T>::frame({1., -2., 0.5}, {std::cos(0.3), std::sin(0.3), 0.}, {-std::sin(0.3), std::cos(0.3), 0.}, {0., 0., 1.});
        auto moved_revolution = std::make_shared<MovedSurface<T>>(revolution, turn);
        std::vector<std::shared_ptr<Surface<T, 3>>> s3{
            plane,
            std::make_shared<CylindricalSurface<T>>(2., fr, std::array<T, 4>{-pi, 3. * pi, 0., 5.})->with_box({0.5, 7., -1., 6.}),
            std::make_shared<ConicalSurface<T>>(1.5, 0.3, fr, std::array<T, 4>{0., 2. * pi, -1., 2.}),
            // a sphere a little past its poles, as with_box leaves it (its constructor refuses it)
            SphericalSurface<T>(3., fr, std::array<T, 4>{0., pi, -1., 1.}).with_box({0., 2. * pi, -0.5 * pi * (1. + 1e-11), 0.5 * pi}),
            std::make_shared<ToroidalSurface<T>>(4., 1., fr, std::array<T, 4>{0., pi, -pi, pi}),
            b3s,
            trimmed_u,
            nurbs_s,
            revolution,
            std::make_shared<OffsetSurface<T>>(b3s, 0.2),
            std::make_shared<OffsetSurface<T>>(moved_revolution, -0.1),
            moved_revolution,
            std::make_shared<SeamlessSurface<T>>(b3s, std::array<T, 2>{1., 0.}),
        };
        // 3D curves
        auto b3 = bs3();
        auto unclamped_bounds = bs3();
        unclamped_bounds->trim(0.1, 0.9, false); // the bounds only
        Frame3<T> moved_frame = fr.moved([&](const point<T, 3> &p) { return turn.apply(p); }); // not quite orthonormal
        auto ellipse = std::make_shared<EllipseCurve<T>>(2., 1., moved_frame, 5.9, 6.6);
        auto composite3 = std::make_shared<CurveComposite<T, 3>>(
            std::vector<std::shared_ptr<Curve<T, 3>>>{std::make_shared<CurveTrimmed<T, 3>>(b3, 0., 0.5), std::make_shared<CurveTrimmed<T, 3>>(b3, 0.5, 1.)});
        std::vector<std::shared_ptr<Curve<T, 3>>> c3{
            b3,
            unclamped_bounds,
            std::make_shared<BSCurveRational<T, 3>>(build_ellipse_arc<T, 3>(1., 2., 0.2, 2.9, {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.})),
            std::make_shared<LineCurve<T>>(point<T, 3>{1., 2., 3.}, point<T, 3>{0.1, -0.7, 0.3}, -2., 5.),
            ellipse,
            composite3,
            std::make_shared<MovedCurve<T>>(composite3, turn),
            std::make_shared<CurveReversed<T, 3>>(ellipse),
            std::make_shared<CurveExtended<T, 3>>(b3, false),
            std::make_shared<CurveReparametrized<T, 3>>(b3, law()),
            std::make_shared<CurveOffset3D<T, BSCfunction<T>>>(b3, law(), std::array<T, 3>{0., 0., 1.}),
            std::make_shared<CurveOnSurface<T, 3>>(c2[0], plane),     // the surface of a face, the pcurve of a co-edge
            std::make_shared<CurveOnSurface<T, 3>>(bs2(0.3), s3.back()), // a fillet's contact on a seamless view
        };
        for (const auto &c : c3)
            edge(c);
        edge(b3); // shared by two edges
        // a wire holding every 2D curve, one twice
        Wire<T> w;
        for (std::size_t k = 0; k < c2.size(); ++k)
            w.coedges.push_back({EdgeId{std::uint32_t(k % 4)}, k % 2 ? Orientation::Reversed : Orientation::Forward, c2[k]});
        w.coedges.push_back({EdgeId{0}, Orientation::Forward, c2[0]});
        w.coedges.push_back({EdgeId{1}, Orientation::Forward, nullptr}); // a free co-edge
        const auto wid = m.add(w);
        for (const auto &s : s3)
            m.add(Face<T>{s, {wid}, 1e-6, false});
        m.add(Face<T>{plane, {}, 1e-6, true}); // a surface shared by two faces
        return m;
    }

    /// The edge of the solid whose middle is nearest to p.
    EdgeId edge_near(const Model<T> &m, SolidId s, point<T, 3> p)
    {
        EdgeId best;
        T d = std::numeric_limits<T>::infinity();
        for (auto e : explore<EdgeId>(m, ShapeId{s}))
        {
            const auto &ed = m.edge(e);
            if (!ed.curve)
                continue;
            if (const T x = norm(ed.curve->value(0.5 * (ed.u1 + ed.u2)) - p); x < d)
                d = x, best = e;
        }
        return best;
    }

    /// The face of the solid whose surface, at the middle of its bounds, is nearest to p.
    FaceId face_near(const Model<T> &m, SolidId s, point<T, 3> p)
    {
        FaceId best;
        T d = std::numeric_limits<T>::infinity();
        for (auto f : explore<FaceId>(m, ShapeId{s}))
        {
            const auto b = m.face(f).surface->bounds();
            if (const T x = norm(m.face(f).surface->value(0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])) - p); x < d)
                d = x, best = f;
        }
        return best;
    }

    SolidId block(Model<T> &m, const point<T, 3> &a, const point<T, 3> &b)
    {
        const auto f = unwrap(make_face(m, brep_tests::quad({a[0], a[1], a[2]}, {b[0], a[1], a[2]}, {a[0], b[1], a[2]}, {b[0], b[1], a[2]})));
        return explore<SolidId>(m, unwrap(make_prism(m, ShapeId{f}, point<T, 3>{0., 0., b[2] - a[2]})).shape).front();
    }

    /// A text of a small model, with geometry of several kinds, for the corrupt files.
    std::string small_text()
    {
        Model<T> m;
        const auto s = block(m, {0., 0., 0.}, {2., 1., 1.});
        m.setName(ShapeId{s}, "block");
        return unwrap(write_brep(m));
    }

    /// The text with one change, made on its JSON document.
    template <typename F>
    std::string edited(const std::string &text, F &&f)
    {
        rapidjson::Document d;
        d.Parse<rapidjson::kParseFullPrecisionFlag>(text.data(), text.size());
        REQUIRE(!d.HasParseError());
        f(d, d.GetAllocator());
        rapidjson::StringBuffer b;
        rapidjson::Writer<rapidjson::StringBuffer> w(b);
        d.Accept(w);
        return {b.GetString(), b.GetSize()};
    }

    /// read_brep refuses the text with this code, the message at this place with this reason.
    void expect_refused(const std::string &text, BuildErrc code, std::string_view where, std::string_view reason_part = {})
    {
        const auto r = read_brep<T>(text, "file.brep");
        REQUIRE_FALSE(r.has_value());
        const auto &msg = r.error().message;
        INFO(msg);
        REQUIRE(r.error().code == code);
        REQUIRE(msg.starts_with(where.empty() ? std::string{"file.brep: "} : "file.brep: " + std::string{where} + ": "));
        if (where.empty())
            REQUIRE_FALSE(msg.starts_with("file.brep: /"));
        if (!reason_part.empty())
            REQUIRE(msg.find(reason_part) != std::string::npos);
    }

    /// write_brep refuses the model with this code, naming this entity.
    void expect_write_refused(const Model<T> &m, BuildErrc code, const ShapeId &entity, std::string_view message_part)
    {
        const auto r = write_brep(m);
        REQUIRE_FALSE(r.has_value());
        INFO(r.error().message);
        REQUIRE(r.error().code == code);
        REQUIRE(r.error().shapes == std::vector<ShapeId>{entity});
        REQUIRE(r.error().message.find(message_part) != std::string::npos);
    }
}

TEST(tests_io_brep_json, empty_model)
{
    Model<T> m;
    const auto text = unwrap(write_brep(m));
    ASSERT_EQ(text, R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":1.0,"next_part":1,"geometry":[],"vertices":[],"edges":[],)"
                    R"("wires":[],"faces":[],"shells":[],"solids":[],"compounds":[]})");
    round_trip(m);
}

TEST(tests_io_brep_json, every_geometry_class)
{
    const auto m = every_geometry();
    const auto text = round_trip(m);
    const auto back = unwrap(read_brep<T>(text));
    // shared again once read: the curve of the first and last edges, the plane of the first and
    // last faces (acceptance test 4) and of the curve on surface of edge 11, the pcurve of the
    // first and last co-edges but one
    const auto last_edge = EdgeId{static_cast<std::uint32_t>(back.capacity<ShapeType::Edge>() - 1)};
    const auto last_face = FaceId{static_cast<std::uint32_t>(back.capacity<ShapeType::Face>() - 1)};
    ASSERT_EQ(back.edge(EdgeId{0}).curve.get(), back.edge(last_edge).curve.get());
    ASSERT_EQ(back.face(FaceId{0}).surface.get(), back.face(last_face).surface.get());
    const auto *cos = dynamic_cast<const CurveOnSurface<T, 3> *>(back.edge(EdgeId{11}).curve.get());
    ASSERT_TRUE(cos != nullptr);
    ASSERT_EQ(cos->p_basisSurface().get(), back.face(FaceId{0}).surface.get());
    const auto &w = back.wire(WireId{0});
    ASSERT_EQ(w.coedges.front().pcurve.get(), w.coedges[w.coedges.size() - 2].pcurve.get());
    ASSERT_EQ(cos->p_basisCurve().get(), w.coedges.front().pcurve.get());
    // each object once in the file: the plane (two faces and a curve on surface), the B-spline
    // surfaces (two in 3D, one the basis of an offset and of a seamless view too; one in 2D)
    ASSERT_EQ(count_type(text, "plane"), 1u);
    ASSERT_EQ(count_type(text, "bspline_surface"), 3u);
}

TEST(tests_io_brep_json, equality_sees_every_difference)
{
    const auto a = every_geometry();
    // a copy shares the geometry: equal
    ASSERT_TRUE(equal(a, a));
    auto b = a;
    ASSERT_TRUE(equal(a, b));
    // a model built again: other objects, the same bits and the same sharing
    ASSERT_TRUE(equal(a, every_geometry()));
    // one bit of a pole
    {
        auto c = a;
        auto bs = std::make_shared<BSCurve<T, 3>>(*std::dynamic_pointer_cast<BSCurve<T, 3>>(c.edge(EdgeId{0}).curve));
        auto p = bs->poles();
        p[1][0] = std::nextafter(p[1][0], 2.);
        bs->movePoles(p);
        c.edge(EdgeId{0}).curve = bs;
        ASSERT_EQ(first_difference(a, c).value_or(""), "edge 0: B-spline poles differ");
    }
    // -0 against 0
    {
        auto c = a;
        c.vertex(VertexId{0}).pnt[0] = -0.;
        ASSERT_EQ(first_difference(a, c).value_or(""), "vertex 0: points or tolerances differ");
    }
    // the same values, but the curve of the last edge no longer shared with the first
    {
        auto c = a;
        const auto last = EdgeId{static_cast<std::uint32_t>(c.capacity<ShapeType::Edge>() - 1)};
        c.edge(last).curve = std::make_shared<BSCurve<T, 3>>(*std::dynamic_pointer_cast<BSCurve<T, 3>>(c.edge(last).curve));
        ASSERT_NE(first_difference(a, c).value_or("").find("shared on one side"), std::string::npos);
    }
    // a dead entity, a name, an instance, the unit, the next part key
    {
        auto c = a;
        c.erase(ShapeId{FaceId{2}});
        ASSERT_EQ(first_difference(a, c).value_or(""), "face 2: alive on one side only");
        auto d = a;
        d.setName(ShapeId{EdgeId{3}}, "x");
        ASSERT_EQ(first_difference(a, d).value_or(""), "edge 3: names differ");
        auto e = a;
        e.setUnitScale(25.4);
        ASSERT_EQ(first_difference(a, e).value_or(""), "units differ");
        auto f = a;
        (void)f.newPart();
        ASSERT_EQ(first_difference(a, f).value_or(""), "next part keys differ");
    }
}

TEST(tests_io_brep_json, dead_entities_ids_and_attributes)
{
    Model<T> m;
    const auto s = block(m, {0., 0., 0.}, {2., 1., 1.});
    const auto t = block(m, {3., 0., 0.}, {4., 1., 2.});
    const auto c = m.add(Compound{{ShapeId{s}, ShapeId{t}, ShapeId{VertexId{2}}}});
    // attributes of every kind, names in UTF-8 with what JSON escapes
    m.setName(ShapeId{s}, "pièce \"A\"\\ \n\t\x01 é");
    m.setName(ShapeId{c}, std::string("nul\0inside", 10));
    m.setName(ShapeId{FaceId{3}}, "a face");
    m.setExternalId(ShapeId{s}, std::numeric_limits<std::int64_t>::min());
    m.setExternalId(ShapeId{EdgeId{5}}, std::numeric_limits<std::int64_t>::max());
    const auto part = m.newPart();
    (void)m.newPart(); // a key given, not used
    m.setInstance(ShapeId{s}, InstanceOf<T>{part, Placement<T>::frame({1., 2., 3.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.})});
    m.setUnitScale(25.4);
    // dead entities everywhere, some still referred to (erase does not cascade), invalid ids
    m.erase(ShapeId{t});
    m.erase(ShapeId{VertexId{2}});
    m.erase(ShapeId{EdgeId{7}});
    m.erase(ShapeId{FaceId{1}});
    m.add(Edge<T>{nullptr, 0., 1., VertexId{}, VertexId{0}, 1e-6, true}); // degenerate, an invalid id
    m.add(Shell{});
    m.erase(ShapeId{ShellId{static_cast<std::uint32_t>(m.capacity<ShapeType::Shell>() - 1)}});
    m.add(Solid{ShellId{}, {}});
    const auto text = round_trip(m);
    const auto back = unwrap(read_brep<T>(text));
    ASSERT_FALSE(back.alive(ShapeId{t}));
    ASSERT_TRUE(back.alive(ShapeId{s}));
    ASSERT_EQ(back.compound(c).shapes[1], ShapeId{t}); // still referred to, dead
    ASSERT_EQ(back.name(ShapeId{c}), std::string("nul\0inside", 10));
    ASSERT_EQ(back.unitScale(), 25.4);
    ASSERT_TRUE(back.instance(ShapeId{s}).has_value());
    // the part keys go on from where they were
    auto copy = back;
    ASSERT_EQ(copy.newPart(), 3u);
}

TEST(tests_io_brep_json, empty_and_absent_names)
{
    // (acceptance test 3) a model holds no empty name: setName(s, "") removes it, name() is then
    // empty as for a shape never named; written and read, the same
    Model<T> m;
    const auto a = m.add(Vertex<T>{});
    const auto b = m.add(Vertex<T>{});
    m.setName(ShapeId{a}, "a");
    m.setName(ShapeId{a}, "");
    const auto text = round_trip(m);
    ASSERT_EQ(text.find("\"name\""), std::string::npos);
    const auto back = unwrap(read_brep<T>(text));
    ASSERT_EQ(back.name(ShapeId{a}), "");
    ASSERT_EQ(back.name(ShapeId{b}), "");
    // an empty name in a file is refused: no model could have written it
    expect_refused(R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":1.0,"next_part":1,"geometry":[],)"
                   R"("vertices":[{"p":[0.0,0.0,0.0],"tol":1e-06,"name":""}],"edges":[],"wires":[],"faces":[],"shells":[],"solids":[],"compounds":[]})",
                   BuildErrc::InvalidFile, "/vertices/0/name");
}

TEST(tests_io_brep_json, compacted)
{
    // (acceptance test 2) a model with dead entities, written compacted: the remap is that of
    // compact(), the model read is the compacted model; the model written is not modified
    Model<T> m;
    const auto s = block(m, {0., 0., 0.}, {2., 1., 1.});
    const auto t = block(m, {3., 0., 0.}, {4., 1., 2.});
    m.setName(ShapeId{t}, "second");
    m.erase(ShapeId{s});
    m.erase(ShapeId{FaceId{2}});
    m.erase(ShapeId{VertexId{0}});
    const auto before = m;
    const auto [text, remap] = unwrap(write_brep_compact(m));
    expect_equal(before, m);
    auto compacted = m;
    const auto expected = compacted.compact();
    ASSERT_TRUE(remap.tables == expected.tables);
    const auto back = unwrap(read_brep<T>(text));
    expect_equal(compacted, back);
    ASSERT_EQ(back.capacity<ShapeType::Solid>(), 1u);
    ASSERT_EQ(back.name(remap.map(ShapeId{t})), "second");
    ASSERT_FALSE(valid(remap.map(ShapeId{s})));
}

TEST(tests_io_brep_json, numbers_bit_for_bit)
{
    // reals of every kind: random bits (subnormals included), -0, the extremes
    Model<T> m;
    std::mt19937_64 g{12345};
    std::vector<T> x{-0., 0., std::numeric_limits<T>::denorm_min(), -std::numeric_limits<T>::denorm_min(), std::numeric_limits<T>::max(),
                     std::numeric_limits<T>::lowest(), std::numeric_limits<T>::min(), 0.1, 1. / 3.};
    while (x.size() < 30000)
    {
        if (const auto d = std::bit_cast<T>(g()); std::isfinite(d))
            x.push_back(d);
        x.push_back(std::bit_cast<T>(g() & 0x800FFFFFFFFFFFFFull)); // subnormal, either sign
    }
    for (std::size_t i = 0; i + 3 < x.size(); i += 4)
        m.add(Vertex<T>{{x[i], x[i + 1], x[i + 2]}, x[i + 3]});
    round_trip(m);
}

TEST(tests_io_brep_json, not_finite_and_not_utf8_refused)
{
    // the text is standard JSON: a number not finite, a name not UTF-8 are refused, naming the entity
    Model<T> m;
    m.add(Vertex<T>{});
    const auto v = m.add(Vertex<T>{{0., std::numeric_limits<T>::infinity(), 0.}, 1e-7});
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{v}, "vertex 1: a number not finite");
    m.vertex(v).pnt[1] = std::numeric_limits<T>::quiet_NaN();
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{v}, "vertex 1: a number not finite");
    m.vertex(v).pnt[1] = 0.;
    // in a geometry: named by the first entity that uses it
    const auto e = m.add(Edge<T>{std::make_shared<LineCurve<T>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}, 0., 1.), 0., 1., v, v});
    auto bs = std::make_shared<BSCurve<T, 3>>(points_vector<T, 3>{{0., 0., 0.}, {std::numeric_limits<T>::infinity(), 0., 0.}}, std::vector<T>{0., 0., 1., 1.}, 1);
    m.edge(e).curve = bs;
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "geometry 0 (of edge 0): a number not finite");
    m.edge(e).curve = nullptr;
    for (const char *bad : {"\xC3", "\xC0\x80", "\xED\xA0\x80", "\xF4\x90\x80\x80", "a\xFF"})
    {
        m.setName(ShapeId{e}, bad);
        expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "edge 0: a name that is not UTF-8");
    }
    m.setName(ShapeId{e}, "\xF0\x9F\x99\x82 \xE2\x82\xAC \xC3\xA9");
    round_trip(m);
    // not JSON numbers in a text either
    expect_refused(R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":NaN})", BuildErrc::InvalidFile, "");
    expect_refused(R"({"format":"gbs-brep","version":1,"real":"double","unit_scale":1e400})", BuildErrc::InvalidFile, "");
    expect_refused("{\"format\":\"gbs-brep\",\"version\":1,\"real\":\"\xFF\"}", BuildErrc::InvalidFile, "");
}

TEST(tests_io_brep_json, float_model)
{
    Model<float> m;
    const auto v = m.add(Vertex<float>{{0.1f, std::nextafter(1.f, 2.f), -0.f}, 1e-5f});
    auto c = std::make_shared<BSCurve<float, 3>>(points_vector<float, 3>{{0.f, 0.f, 0.f}, {1.f / 3.f, 0.2f, 0.f}, {1.f, 0.f, 0.f}},
                                                 std::vector<float>{0.f, 0.f, 0.f, 1.f, 1.f, 1.f}, 2);
    m.add(Edge<float>{c, 0.f, 1.f, v, v, 1e-5f});
    const auto text = round_trip(m);
    ASSERT_NE(text.find(R"("real":"float")"), std::string::npos);
    // a float model is not read as a double one, nor the reverse
    expect_refused(text, BuildErrc::InvalidFile, "/real");
    const auto r = read_brep<float>(unwrap(write_brep(Model<T>{})));
    ASSERT_FALSE(r.has_value());
    ASSERT_NE(r.error().message.find(": /real: "), std::string::npos);
    // a number beyond the range of a float: refused before its conversion, undefined there
    for (const double x : {1e300, -1e39})
    {
        const auto huge = read_brep<float>(edited(text, [x](auto &d, auto &) { d["vertices"][0]["tol"].SetDouble(x); }), "file.brep");
        ASSERT_FALSE(huge.has_value());
        ASSERT_TRUE(huge.error().message.starts_with("file.brep: /vertices/0/tol: a number of a float expected"));
    }
}

TEST(tests_io_brep_json, every_operation)
{
    // (acceptance test 1) the result of each operation of gbs-brep stays in a model that is
    // written and read back equal: no class of computation (FieldCurve, PathCurve, LoopCurve)
    // remains in it, the write would refuse it
    using brep_tests::quad;
    Model<T> m;
    const ax1<T, 3> z_axis{{{0., 0., 0.}, {0., 0., 1.}}};
    std::vector<ShapeId> results;
    // primitives and sweeps: a prism of a B-spline face, revolutions (cylinders, cone, sphere, torus)
    const auto blk = block(m, {0., 0., 0.}, {2., 1.5, 1.});
    results.push_back(ShapeId{blk});
    const auto rect = unwrap(make_face(m, quad({2., 0., 0.}, {3., 0., 0.}, {2., 0., 1.}, {3., 0., 1.})));
    results.push_back(unwrap(make_revolution(m, rect, z_axis)).shape);
    results.push_back(unwrap(make_revolution(m, rect, z_axis, pi)).shape);
    const auto xz = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, {1., 0., 0.}}}, std::array<T, 4>{-5., 5., -5., 5.});
    const auto e1 = unwrap(make_edge(m, point<T, 3>{0., 0., 0.}, point<T, 3>{2., 0., 0.}));
    const auto e2 = unwrap(make_edge(m, point<T, 3>{2., 0., 0.}, point<T, 3>{0., 0., 3.}));
    const auto e3 = unwrap(make_edge(m, point<T, 3>{0., 0., 3.}, point<T, 3>{0., 0., 0.}));
    results.push_back(unwrap(make_revolution(m, unwrap(make_face(m, xz, unwrap(make_wire(m, std::vector<EdgeId>{e1, e2, e3})))), z_axis)).shape);
    const auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, {1., 0., 0.}}}, -pi / 2, pi / 2);
    const auto hd = unwrap(make_face(m, xz, unwrap(make_wire(m, std::vector<EdgeId>{unwrap(make_edge(m, arc)), unwrap(make_edge(m, point<T, 3>{0., 0., 2.}, point<T, 3>{0., 0., -2.}))}))));
    results.push_back(unwrap(make_revolution(m, hd, z_axis)).shape);
    // STEP reading, and a torus by revolution of the cap of the cylinder
    const auto cyl = unwrap(step::read_step_file(m, step_dir() + "/cylinder.stp")).root; // R 5, h 10
    FaceId cap;
    for (auto f : explore<FaceId>(m, cyl))
        if (dynamic_cast<const Plane<T> *>(m.face(f).surface.get()) && std::abs(m.face(f).surface->value(0., 0.)[2]) < 1e-9)
            cap = f;
    results.push_back(unwrap(make_revolution(m, cap, ax1<T, 3>{{{10., 0., 0.}, {0., 1., 0.}}})).shape);
    // a pipe: a disk along a quarter circle
    const ax2<T, 3> dax{{{4., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}}};
    const auto disk = unwrap(make_face(m, std::make_shared<Plane<T>>(dax, std::array<T, 4>{-2., 2., -2., 2.}),
                                       unwrap(make_wire(m, std::vector<EdgeId>{unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(1., 1., dax)))}))));
    const ax2<T, 3> pax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};
    results.push_back(unwrap(make_pipe(m, disk, unwrap(make_edge(m, std::make_shared<EllipseCurve<T>>(4., 4., pax, 0., pi / 2))))).shape);
    // booleans (faces unified on their surfaces), section, split
    const auto plate = ShapeId{block(m, {-6., -6., 1.}, {6., 6., 4.})};
    const auto holed = unwrap(cut(m, plate, cyl));
    results.push_back(holed.shape);
    results.push_back(unwrap(fuse(m, ShapeId{block(m, {20., 0., 0.}, {22., 2., 2.})}, ShapeId{block(m, {21., 1., 1.}, {23., 3., 3.})})).shape);
    results.push_back(unwrap(common(m, ShapeId{block(m, {30., 0., 0.}, {32., 2., 2.})}, ShapeId{block(m, {31., 1., 1.}, {33., 3., 3.})})).shape);
    const auto cutting = std::make_shared<Plane<T>>(ax2<T, 3>{{{0., 0., 2.}, {0.1, 0.2, 1.}, {1., 0., 0.}}}, std::array<T, 4>{-20., 20., -20., 20.});
    (void)unwrap(section(m, holed.shape, cutting));
    (void)unwrap(split(m, holed.shape, cutting));
    // fillets and chamfers of a block, a chordal fillet, an offset solid, a thick solid, an offset face
    const auto b1 = block(m, {40., 0., 0.}, {42., 1.5, 1.});
    results.push_back(ShapeId{unwrap(fillet(m, b1, std::vector<EdgeId>{edge_near(m, b1, {42., 0.75, 1.})}, 0.3))});
    const auto b2 = block(m, {50., 0., 0.}, {52., 1.5, 1.});
    results.push_back(ShapeId{unwrap(chamfer(m, b2, std::vector<EdgeId>{edge_near(m, b2, {52., 0.75, 1.})}, 0.2))});
    const auto b3 = block(m, {60., 0., 0.}, {62., 1.5, 1.});
    results.push_back(ShapeId{unwrap(fillet_chordal(m, b3, std::vector<EdgeId>{edge_near(m, b3, {62., 0.75, 1.})}, 0.3))});
    const auto b4 = block(m, {70., 0., 0.}, {72., 1.5, 1.});
    results.push_back(ShapeId{unwrap(offset_solid(m, b4, 0.2))});
    const auto b5 = block(m, {80., 0., 0.}, {82., 1.5, 1.});
    FaceId top;
    for (auto f : explore<FaceId>(m, ShapeId{b5}))
        if (std::abs(m.face(f).surface->value(0.5, 0.5)[2] - 1.) < 1e-9)
            top = f;
    results.push_back(ShapeId{unwrap(make_thick(m, b5, std::vector<FaceId>{top}, 0.1))});
    // (review of PR 252) the blends that make the most new geometry: the corner of three fillets of
    // three radii (a vertex blend, B-spline), of one radius (a sphere), a conic section, a fillet
    // between two faces over a chamfer
    for (const bool same : {false, true})
    {
        const T x0 = same ? 100. : 90.;
        const auto c = block(m, {x0, 0., 0.}, {x0 + 3., 3., 3.});
        const auto ex = edge_near(m, c, {x0 + 1.5, 3., 3.}), ey = edge_near(m, c, {x0 + 3., 1.5, 3.}), ez = edge_near(m, c, {x0 + 3., 3., 1.5});
        const RadiusLaw<T> radii{{{ex, 0.5, 0.5}, {ey, 0.5, same ? 0.5 : 0.7}, {ez, 0.5, same ? 0.5 : 0.9}}};
        results.push_back(ShapeId{unwrap(fillet(m, c, std::vector{ex, ey, ez}, radii))});
    }
    {
        const auto b6 = block(m, {110., 0., 0.}, {112., 1.5, 1.});
        FilletOptions<T> conic;
        conic.section = BlendSection<T>::conic(0.6);
        results.push_back(ShapeId{unwrap(fillet(m, b6, std::vector<EdgeId>{edge_near(m, b6, {112., 0.75, 1.})}, 0.3, conic))});
        const auto b7 = block(m, {120., 0., 0.}, {122., 1., 1.});
        const auto ch = unwrap(chamfer(m, b7, std::vector<EdgeId>{edge_near(m, b7, {121., 1., 1.})}, 0.1));
        const auto top_face = face_near(m, ch, {121., 0.5, 1.}), side_face = face_near(m, ch, {121., 1., 0.5});
        results.push_back(ShapeId{unwrap(fillet_faces(m, ch, std::vector<FaceId>{top_face}, std::vector<FaceId>{side_face}, 0.3))});
    }
    const auto bump = unwrap(make_face(m, quad({0., 0., 10.}, {2., 0., 10.5}, {0., 2., 10.3}, {2., 2., 11.})));
    const auto off = unwrap(offset_face(m, bump, 0.2)); // an OffsetSurface
    results.push_back(ShapeId{off});
    // copies: identity (geometry shared) and placed (an offset surface moved: MovedSurface)
    const Placement<T> turn = Placement<T>::frame({1., 2., 3.}, {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.});
    results.push_back(copy_shape(m, holed.shape, Placement<T>{}));
    results.push_back(copy_shape(m, ShapeId{off}, turn));
    results.push_back(copy_shape(m, results[1], turn));
    // the root of a blade on its hub filleted: contacts on a seamless view (SeamlessSurface)
    const auto [blade, root] = blade_fixture::blade_on_hub(m);
    results.push_back(ShapeId{unwrap(fillet(m, blade, root, 0.03))});

    // the classes met, by their names
    std::set<std::string> classes;
    for (auto f : m.ids<FaceId>())
    {
        const auto &s = *m.face(f).surface; // typeid of a call would evaluate it
        classes.insert(typeid(s).name());
    }
    for (auto e : m.ids<EdgeId>())
        if (const auto &c = m.edge(e).curve)
        {
            const auto &curve = *c;
            classes.insert(typeid(curve).name());
            if (auto x = dynamic_cast<const CurveOnSurface<T, 3> *>(c.get()))
            {
                const auto &s = *x->p_basisSurface();
                classes.insert(typeid(s).name());
            }
        }
    for (const auto &name : {typeid(Plane<T>).name(), typeid(CylindricalSurface<T>).name(), typeid(ConicalSurface<T>).name(),
                             typeid(SphericalSurface<T>).name(), typeid(ToroidalSurface<T>).name(), typeid(BSSurface<T, 3>).name(),
                             typeid(OffsetSurface<T>).name(), typeid(MovedSurface<T>).name(), typeid(SeamlessSurface<T>).name(),
                             typeid(LineCurve<T>).name(), typeid(EllipseCurve<T>).name(), typeid(CurveOnSurface<T, 3>).name()})
    {
        INFO(name);
        ASSERT_TRUE(classes.contains(name));
    }
    const auto text = round_trip(m);
    // the results read back: the same check and volumes, bit for bit
    const auto back = unwrap(read_brep<T>(text));
    for (const auto &r : results)
    {
        ASSERT_EQ(check(back, r).entries.size(), check(m, r).entries.size());
        for (auto so : explore<SolidId>(back, r))
            ASSERT_EQ(std::bit_cast<std::uint64_t>(signed_volume(back, back.solid(so).outer, 32)),
                      std::bit_cast<std::uint64_t>(signed_volume(m, m.solid(so).outer, 32)));
    }
}

TEST(tests_io_brep_json, step_files)
{
    // every STEP file of the tests, each in a model of its own, then all in one model with
    // some of them erased
    Model<T> all;
    std::vector<ShapeId> roots;
    std::size_t files{};
    for (const auto &entry : std::filesystem::directory_iterator(step_dir()))
    {
        if (entry.path().extension() != ".stp")
            continue;
        INFO(entry.path().filename().string());
        Model<T> m;
        auto r = step::read_step_file(m, entry.path(), step::StepReadOptions{.target_unit_mm = 1.});
        if (!r)
            continue; // a file the reader refuses holds no model
        ++files;
        round_trip(m);
        roots.push_back(all.append(m).map(r->root));
    }
    ASSERT_GE(files, 20u);
    for (std::size_t k = 0; k < roots.size(); k += 3)
        all.erase(roots[k]);
    round_trip(all);
}

TEST(tests_io_brep_json, deterministic_and_pure)
{
    // (acceptance test 6) the same model gives the same bytes: written twice, built twice (other
    // addresses), written by several threads at once from the same const model
    const auto a = every_geometry();
    const auto text = unwrap(write_brep(a));
    ASSERT_EQ(unwrap(write_brep(a)), text);
    ASSERT_EQ(unwrap(write_brep(every_geometry())), text);
    std::vector<std::future<std::string>> writers;
    for (int k = 0; k < 8; ++k)
        writers.push_back(std::async(std::launch::async, [&a] { return unwrap(write_brep(a)); }));
    for (auto &w : writers)
        ASSERT_EQ(w.get(), text);
}

TEST(tests_io_brep_json, files_streams_and_values)
{
    const auto m = every_geometry();
    const TempDir tmp{"gbs_brep_json_"};
    const auto &dir = tmp.path;
    const auto path = dir / "modèle.brep";
    unwrap(write_brep_file(m, path));
    expect_equal(m, unwrap(read_brep_file<T>(path)));
    unwrap(write_brep_file(m, path, {.pretty = true})); // replaced whole
    expect_equal(m, unwrap(read_brep_file<T>(path)));
    std::stringstream ss;
    unwrap(write_brep(m, ss));
    expect_equal(m, unwrap(read_brep<T>(ss)));
    // a file missing, a file cut short: InvalidFile, with its name
    const auto missing = read_brep_file<T>(dir / "missing.brep");
    ASSERT_FALSE(missing.has_value());
    ASSERT_EQ(missing.error().code, BuildErrc::InvalidFile);
    write_text_atomic(dir / "cut.brep", unwrap(write_brep(m)).substr(0, 100));
    const auto cut_short = read_brep_file<T>(dir / "cut.brep");
    ASSERT_FALSE(cut_short.has_value());
    ASSERT_NE(cut_short.error().message.find("cut.brep"), std::string::npos);
    // inside an application's own document: a value of it
    rapidjson::StringBuffer b;
    rapidjson::Writer<rapidjson::StringBuffer> w(b);
    w.StartObject();
    w.Key("application");
    w.StartObject();
    w.Key("roots");
    w.StartArray();
    w.Uint(3);
    w.EndArray();
    w.EndObject();
    w.Key("model");
    unwrap(write_brep_value(m, w));
    w.EndObject();
    rapidjson::Document d;
    d.Parse<brep_json_parse_flags>(b.GetString(), b.GetSize());
    ASSERT_FALSE(d.HasParseError());
    expect_equal(m, unwrap(read_brep_value<T>(d["model"], "document.app: model")));
    const auto refused = read_brep_value<T>(d["application"], "document.app: application");
    ASSERT_FALSE(refused.has_value());
    ASSERT_TRUE(refused.error().message.starts_with("document.app: application: /format: "));
    // a document parsed with the application's flags, the encoding not checked: a name that is not
    // UTF-8 refused by the reader all the same (it could not be written back)
    Model<T> named;
    named.setName(ShapeId{named.add(Vertex<T>{})}, "ok");
    auto bad = unwrap(write_brep(named));
    bad.replace(bad.find("\"ok\"") + 1, 2, "\xFF\xFE");
    rapidjson::Document loose;
    loose.Parse<rapidjson::kParseFullPrecisionFlag>(bad.data(), bad.size());
    ASSERT_FALSE(loose.HasParseError());
    const auto not_utf8 = read_brep_value<T>(loose, "document.app: model");
    ASSERT_FALSE(not_utf8.has_value());
    ASSERT_TRUE(not_utf8.error().message.starts_with("document.app: model: /vertices/0/name: "));
}

TEST(tests_io_brep_json, version_and_format)
{
    const auto text = small_text();
    const auto newer = edited(text, [](auto &d, auto &) { d["version"].SetUint(2); });
    expect_refused(newer, BuildErrc::UnsupportedVersion, "/version", "version 2 of the format, newer than this reader (version 1)");
    // a newer version is reported as such even if what follows is unknown
    expect_refused(edited(text,
                          [](auto &d, auto &a) {
                              d["version"].SetUint(7);
                              d.AddMember("new_table", rapidjson::Value{rapidjson::kArrayType}, a);
                          }),
                   BuildErrc::UnsupportedVersion, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d["version"].SetUint(0); }), BuildErrc::InvalidFile, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d.RemoveMember("version"); }), BuildErrc::InvalidFile, "/version");
    expect_refused(edited(text, [](auto &d, auto &) { d["format"].SetString("step"); }), BuildErrc::InvalidFile, "/format", "not a BREP model of gbs");
    expect_refused(edited(text, [](auto &d, auto &) { d.RemoveMember("format"); }), BuildErrc::InvalidFile, "/format");
    expect_refused("[1, 2]", BuildErrc::InvalidFile, "", "not a BREP model of gbs");
    expect_refused(R"({"a": 1})", BuildErrc::InvalidFile, "/format");
}

TEST(tests_io_brep_json, truncated_and_malformed)
{
    const auto text = small_text();
    // every text cut short is refused as such, never read
    for (std::size_t n = 0; n < text.size(); ++n)
        expect_refused(text.substr(0, n), BuildErrc::InvalidFile, "");
    expect_refused(text + "x", BuildErrc::InvalidFile, "", "line 1, column");
    expect_refused(text.substr(0, 40) + '\0' + text.substr(40), BuildErrc::InvalidFile, "");
    // a deep nesting: refused without exhausting the stack (the parser is not recursive)
    expect_refused(std::string(1'000'000, '['), BuildErrc::InvalidFile, "");
    auto deep = text;
    const std::string real = R"("real":"double")";
    deep.replace(deep.find(real), real.size(), R"("real":)" + std::string(1'000'000, '[') + std::string(1'000'000, ']'));
    expect_refused(deep, BuildErrc::InvalidFile, "/real");
}

TEST(tests_io_brep_json, corrupt_values)
{
    // (acceptance test 5) indices out of bounds, wrong types, unknown keys: each refused at its place
    const auto text = small_text();
    using V = rapidjson::Value;
    rapidjson::Document doc;
    doc.Parse(text.data(), text.size());
    const auto n_geometry = doc["geometry"].Size(), n_compounds = doc["compounds"].Size();
    const auto edge_curve = doc["edges"][0]["curve"].GetUint(), face_surface = doc["faces"][0]["surface"].GetUint();
    SizeType bspline = n_geometry; // the first B-spline curve
    for (SizeType i = 0; i < n_geometry && bspline == n_geometry; ++i)
        if (std::string_view{doc["geometry"][i]["type"].GetString()} == "bspline_curve")
            bspline = i;
    REQUIRE(bspline < n_geometry);
    const auto at_bspline = "/geometry/" + std::to_string(bspline);
    auto refused = [&](auto &&edit, std::string_view where, std::string_view reason = {}) {
        expect_refused(edited(text, edit), BuildErrc::InvalidFile, where, reason);
    };
    refused([](auto &d, auto &) { d["edges"][0]["v"][0].SetUint(999); }, "/edges/0/v/0", "vertex index below");
    refused([](auto &d, auto &) { d["edges"][0]["v"][1].SetInt(-1); }, "/edges/0/v/1");
    refused([](auto &d, auto &) { d["wires"][0]["coedges"][0][0].SetUint(10'000); }, "/wires/0/coedges/0/0");
    refused([](auto &d, auto &) { d["wires"][0]["coedges"][0][1].SetString("X"); }, "/wires/0/coedges/0/1");
    refused([](auto &d, auto &a) { d["compounds"].PushBack(V{rapidjson::kObjectType}, a); }, "/compounds/" + std::to_string(n_compounds),
            "member \"shapes\" missing");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("colour", V{1}, a); }, "/vertices/0", "unknown member \"colour\"");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("tol", V{1.}, a); }, "/vertices/0", "given twice");
    refused([](auto &d, auto &) { d["vertices"][0]["tol"].SetString("small"); }, "/vertices/0/tol", "a number");
    refused([](auto &d, auto &) { d["vertices"][0]["p"].PopBack(); }, "/vertices/0/p", "an array of 3");
    refused([](auto &d, auto &) { d["faces"][0]["natural_bounds"].SetInt(1); }, "/faces/0/natural_bounds");
    refused([](auto &d, auto &) { d["unit_scale"].SetDouble(-1.); }, "/unit_scale");
    refused([](auto &d, auto &a) { d["vertices"][0].AddMember("name", V{"", 0}, a); }, "/vertices/0/name");
    refused([](auto &d, auto &a) {
        V inst{rapidjson::kObjectType};
        V place{rapidjson::kObjectType};
        V r{rapidjson::kArrayType}, t{rapidjson::kArrayType};
        for (int i = 0; i < 3; ++i)
        {
            V row{rapidjson::kArrayType};
            for (int j = 0; j < 3; ++j)
                row.PushBack(V{i == j ? 1. : 0.}, a);
            r.PushBack(row, a);
            t.PushBack(V{0.}, a);
        }
        place.AddMember("r", r, a).AddMember("t", t, a);
        inst.AddMember("part", V{1u}, a).AddMember("placement", place, a);
        d["solids"][0].AddMember("instance", inst, a);
    }, "/solids/0/instance/part", "below next_part");
    // geometry: of the wrong kind, referring to itself or to one after it, out of its range
    refused([&](auto &d, auto &) { d["edges"][0]["curve"].SetUint(face_surface); }, "/edges/0/curve", "a 3D curve expected");
    refused([&](auto &d, auto &) { d["faces"][0]["surface"].SetUint(edge_curve); }, "/faces/0/surface", "a 3D surface expected");
    refused([&](auto &d, auto &) { d["edges"][0]["curve"].SetUint(d["geometry"].Size()); }, "/edges/0/curve", "below");
    refused(
        [&](auto &d, auto &a) {
            V c{rapidjson::kObjectType};
            V list{rapidjson::kArrayType};
            list.PushBack(V{d["geometry"].Size()}, a); // itself
            c.AddMember("type", "composite_curve", a).AddMember("dim", V{3}, a).AddMember("curves", list, a);
            d["geometry"].PushBack(c, a);
        },
        "/geometry/" + std::to_string(n_geometry) + "/curves/0", "below " + std::to_string(n_geometry));
    refused([&](auto &d, auto &) { d["geometry"][bspline]["type"].SetString("spiral"); }, at_bspline + "/type", "curve of type \"spiral\"");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["dim"].SetUint(4); }, at_bspline + "/dim");
    refused([&](auto &d, auto &) {
        auto &k = d["geometry"][bspline]["knots"];
        std::swap(k[0], k[k.Size() - 1]); // out of order
    }, at_bspline + "/knots", "in order");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["degree"].SetUint(10); }, at_bspline + "/knots");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["poles"].PopBack(); }, at_bspline + "/knots");
    // (review of PR 252) what the format bounds: a degree above 24, evaluated at an exponential
    // cost; a composite curve of a curve whose bounds decrease, where CurveComposite::value read
    // out of its list; a next part key above 2^63, past which newPart() would run out
    refused([&](auto &d, auto &) { d["geometry"][bspline]["degree"].SetUint(25); }, at_bspline + "/degree", "a degree up to 24");
    refused([&](auto &d, auto &) { d["geometry"][bspline]["degree"].SetUint(40); }, at_bspline + "/degree");
    const auto dim_b = doc["geometry"][bspline]["dim"].GetUint();
    refused(
        [&](auto &d, auto &a) {
            V range{rapidjson::kArrayType};
            range.PushBack(V{1.}, a).PushBack(V{0.}, a);
            V t{rapidjson::kObjectType};
            t.AddMember("type", "trimmed_curve", a).AddMember("dim", V{dim_b}, a).AddMember("basis", V{bspline}, a).AddMember("range", range, a);
            d["geometry"].PushBack(t, a);
            V list{rapidjson::kArrayType};
            list.PushBack(V{n_geometry}, a);
            V c{rapidjson::kObjectType};
            c.AddMember("type", "composite_curve", a).AddMember("dim", V{dim_b}, a).AddMember("curves", list, a);
            d["geometry"].PushBack(c, a);
        },
        "/geometry/" + std::to_string(n_geometry + 1) + "/curves/0", "bounds increase");
    refused([](auto &d, auto &) { d["next_part"].SetUint64((std::uint64_t{1} << 63) + 1); }, "/next_part", "from 1 to 2^63");
}

TEST(tests_io_brep_json, constructors_refusals_reported)
{
    // values a constructor of gbs refuses: reported at the geometry
    Model<T> m;
    const auto v = m.add(Vertex<T>{});
    m.add(Edge<T>{std::make_shared<EllipseCurve<T>>(1., 2., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}), 0., 1., v, v});
    const auto text = edited(unwrap(write_brep(m)), [](auto &d, auto &) { d["geometry"][0]["radii"][0].SetDouble(-1.); });
    expect_refused(text, BuildErrc::InvalidFile, "/geometry/0", "radii must be positive");
}

TEST(tests_io_brep_json, random_damage)
{
    // bytes changed at random: the text is read or refused with a BuildError, never an exception
    const auto text = small_text();
    std::mt19937 g{7};
    const std::string alphabet = "0123456789-+.eE,:[]{}\"nultrfaseNIy x";
    std::size_t read{}, refused{};
    for (int k = 0; k < 3000; ++k)
    {
        auto t = text;
        for (int j = 0; j < 1 + k % 3; ++j)
            t[g() % t.size()] = alphabet[g() % alphabet.size()];
        try
        {
            (read_brep<T>(t).has_value() ? read : refused)++;
        }
        catch (const std::exception &e)
        {
            INFO(e.what());
            FAIL("an exception out of read_brep");
        }
    }
    ASSERT_GT(read, 0u);
    ASSERT_GT(refused, 0u);
}

TEST(tests_io_brep_json, refused_writes)
{
    // a geometry that cannot be read back exactly: an old Line (its constructor normalizes), a
    // class of an application (even derived from a B-spline): refused, naming the entity
    struct MyCurve : BSCurve<T, 3>
    {
        MyCurve(const points_vector<T, 3> &p, const std::vector<T> &k, std::size_t d) : BSCurve<T, 3>(p, k, d) {}
    };
    Model<T> m;
    const auto v = m.add(Vertex<T>{});
    const auto e = m.add(Edge<T>{std::make_shared<Line<T, 3>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 1., 0.}), 0., 1., v, v});
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "edge 0: curve of class");
    m.edge(e).curve = std::make_shared<MyCurve>(points_vector<T, 3>{{0., 0., 0.}, {1., 0., 0.}}, std::vector<T>{0., 0., 1., 1.}, 1);
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "cannot be written");
    // inside another geometry too
    m.edge(e).curve = std::make_shared<CurveReversed<T, 3>>(std::make_shared<Circle3d<T>>(1., ax1<T, 3>{{{0., 0., 0.}, {0., 0., 1.}}}));
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "cannot be written");
    // (review of PR 252) a composite curve of a curve whose bounds decrease, a B-spline of degree
    // 25 (evaluated at an exponential cost)
    const auto line = std::make_shared<LineCurve<T>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}, -10., 10.);
    m.edge(e).curve = std::make_shared<CurveComposite<T, 3>>(std::vector<std::shared_ptr<Curve<T, 3>>>{std::make_shared<CurveTrimmed<T, 3>>(line, 1., 0.)});
    expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "edge 0: composite curve of a curve whose bounds decrease");
    {
        points_vector<T, 3> p(26);
        for (std::size_t i = 0; i < p.size(); ++i)
            p[i] = {T(i), T(i % 2), 0.};
        std::vector<T> k(52, 0.);
        std::fill(k.begin() + 26, k.end(), 1.);
        m.edge(e).curve = std::make_shared<BSCurve<T, 3>>(p, k, 25);
        expect_write_refused(m, BuildErrc::UnsupportedEntity, ShapeId{e}, "B-spline of degree 25, above 24");
    }
    // an id beyond its table: nothing written, the file in place kept; compacted, refused as well
    // (compact() mapped it to an invalid id in silence)
    m.edge(e).curve = nullptr;
    m.edge(e).v1 = VertexId{3};
    expect_write_refused(m, BuildErrc::InvalidId, ShapeId{e}, "edge 0: vertex id 3 beyond its table (1 vertices)");
    const auto compacted = write_brep_compact(m);
    ASSERT_FALSE(compacted.has_value());
    ASSERT_EQ(compacted.error().code, BuildErrc::InvalidId);
    ASSERT_EQ(compacted.error().shapes, std::vector<ShapeId>{ShapeId{e}});
    ASSERT_EQ(compacted.error().message, "edge 0: vertex id 3 beyond its table (1 vertices)");
    const TempDir tmp{"gbs_brep_json_w_"};
    const auto &dir = tmp.path;
    write_text_atomic(dir / "kept.brep", "old");
    ASSERT_FALSE(write_brep_file(m, dir / "kept.brep").has_value());
    ASSERT_EQ(read_text(dir / "kept.brep").value(), "old");
}

TEST(tests_io_brep_json, geometry_depth_bounded)
{
    // (review of PR 252) the geometry nested brep_json_max_depth levels at most: what goes through
    // a geometry (the writer, brep::equal, an evaluation, its destruction) recurses once per level,
    // and a file of 93 KB nesting 2,000 levels exhausted the stack of a thread
    ASSERT_EQ(brep_json_max_depth, 256u);
    const auto m = unwrap(read_brep<T>(chain_text(brep_json_max_depth - 1))); // its last geometry the deepest held
    round_trip(m);
    expect_refused(chain_text(brep_json_max_depth), BuildErrc::InvalidFile, "/geometry/256", "a geometry nested at most 256 levels deep expected");
    expect_refused(chain_text(20'000), BuildErrc::InvalidFile, "/geometry/256");
    // written: one level more refused, naming the edge, the depth counted below a geometry met before
    Model<T> w = m;
    const auto e = w.add(Edge<T>{std::make_shared<CurveReversed<T, 3>>(w.edge(EdgeId{0}).curve), 0., 1., VertexId{0}, VertexId{0}, 1e-7});
    expect_write_refused(w, BuildErrc::UnsupportedEntity, ShapeId{e}, "edge 1: geometry nested deeper than 256 levels");
    // a chain built in memory, deeper still: refused without going down to its end
    std::shared_ptr<Curve<T, 3>> c = std::make_shared<LineCurve<T>>(point<T, 3>{0., 0., 0.}, point<T, 3>{1., 0., 0.}, 0., 1.);
    for (int k = 0; k < 1'000; ++k)
        c = std::make_shared<CurveReversed<T, 3>>(c);
    Model<T> d;
    const auto dv = d.add(Vertex<T>{});
    const auto de = d.add(Edge<T>{c, 0., 1., dv, dv, 1e-7});
    expect_write_refused(d, BuildErrc::UnsupportedEntity, ShapeId{de}, "edge 0: geometry nested deeper than 256 levels");
}

TEST(tests_io_brep_json, offset_curve_directions)
{
    // (review of PR 252) the direction of an offset curve given back as it is: normalized again by
    // the constructor, a third of the directions changed by an ulp ((0, 1, 1) among them), and the
    // writer refused them
    Model<T> m;
    const auto v = m.add(Vertex<T>{});
    std::vector<std::array<T, 3>> directions{{0., 1., 1.}, {1., 0., 1.}, {1., 2., 3.}, {0., 0., -1.}};
    std::mt19937_64 g{252};
    std::uniform_real_distribution<T> u{-1., 1.};
    for (int k = 0; k < 200; ++k)
        directions.push_back({u(g), u(g), u(g)});
    for (const auto &d : directions)
        m.add(Edge<T>{std::make_shared<CurveOffset3D<T, BSCfunction<T>>>(bs3(), law(), d), 0., 1., v, v, 1e-7});
    const auto text = round_trip(m); // brep::equal compares the directions bit for bit
    // taken as it is only when unit, in memory and in a file
    using C = CurveOffset3D<T, BSCfunction<T>>;
    ASSERT_THROW((C{bs3(), law(), std::array<T, 3>{0., 1., 1.}, unit_direction}), std::invalid_argument);
    const auto unit = normalized(std::array<T, 3>{0., 1., 1.});
    ASSERT_TRUE((C{bs3(), law(), unit, unit_direction}.direction() == unit));
    rapidjson::Document doc;
    doc.Parse(text.data(), text.size());
    SizeType off = 0;
    while (std::string_view{doc["geometry"][off]["type"].GetString()} != "offset_curve")
        ++off;
    expect_refused(edited(text, [off](auto &d, auto &) { d["geometry"][off]["direction"][0].SetDouble(0.5); }), BuildErrc::InvalidFile,
                   "/geometry/" + std::to_string(off) + "/direction", "a unit direction expected");
}

TEST(tests_io_brep_json, next_part_bounded)
{
    // (review of PR 252) the next part key at most 2^63, so that newPart() never runs out: read and
    // written up to it, refused above
    const auto top = std::uint64_t{1} << 63;
    const auto text = edited(small_text(), [top](auto &d, auto &) { d["next_part"].SetUint64(top); });
    const auto back = unwrap(read_brep<T>(text));
    ASSERT_EQ(back.nextPart(), top);
    ASSERT_EQ(unwrap(write_brep(back)), text);
    Model<T> m;
    m.reservePartKeys(top + 1);
    const auto r = write_brep(m);
    ASSERT_FALSE(r.has_value());
    ASSERT_EQ(r.error().code, BuildErrc::UnsupportedEntity);
    ASSERT_EQ(r.error().message, "next part key 9223372036854775809, above 2^63, the most the format holds");
}
