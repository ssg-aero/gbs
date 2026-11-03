#include <doctest_gtest.hpp>
#include <gbs/subd/limit.h>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Cages of reference (scripts/subd_reference/cages.py) and the values OpenSubdiv gives for them
// (scripts/subd_reference/reference.cpp), tests/in/subd/. Design: docs/sources/design/subd_pr01_cage.md.

using namespace gbs;
using subd::Cage;
using subd::CageIssue;
using subd::LimitSurface;
using subd::Topology;
using subd::VertexRule;

namespace
{
    using T = double;
    using P3 = point<T, 3>;

    const std::vector<std::string> names{"box", "torus", "ngons", "grid", "fans", "creases", "hexprism"};

    auto read_json(const std::string &name) -> rapidjson::Document
    {
        std::ifstream in{std::string(GBS_TESTS_IN_DIR) + "/subd/" + name};
        rapidjson::IStreamWrapper isw{in};
        rapidjson::Document d;
        d.ParseStream(isw);
        if (d.HasParseError())
            throw std::runtime_error{"cannot read " + name};
        return d;
    }

    auto to_point(const rapidjson::Value &v) -> P3 { return {v[0].GetDouble(), v[1].GetDouble(), v[2].GetDouble()}; }

    auto read_cage(const std::string &name) -> Cage<T>
    {
        const auto d = read_json(name + ".cage.json");
        Cage<T> c;
        for (const auto &v : d["vertices"].GetArray())
            c.vertices.push_back(to_point(v));
        for (const auto &f : d["faces"].GetArray())
        {
            std::vector<std::size_t> fv;
            for (const auto &i : f.GetArray())
                fv.push_back(i.GetUint());
            c.faces.push_back(fv);
        }
        for (const auto &e : d["creases"].GetArray())
            c.creases.push_back({e[0].GetUint(), e[1].GetUint()});
        for (const auto &v : d["corners"].GetArray())
            c.corners.push_back(v.GetUint());
        return c;
    }

    auto size_of(const Cage<T> &c) -> T
    {
        P3 lo = c.vertices[0], hi = c.vertices[0];
        for (const auto &p : c.vertices)
            for (std::size_t i{}; i < 3; ++i)
                lo[i] = std::min(lo[i], p[i]), hi[i] = std::max(hi[i], p[i]);
        return norm(hi - lo);
    }

    auto angle(const P3 &a, const P3 &b) -> T { return std::atan2(norm(a ^ b), a * b); }

    /// For each point of ref, the index of ours within tol of it; a bijection.
    auto match(const points_vector<T, 3> &ours, const rapidjson::Value &ref, T tol) -> std::vector<std::size_t>
    {
        std::vector<std::size_t> map;
        std::vector<bool> used(ours.size(), false);
        for (const auto &r : ref.GetArray())
        {
            const auto p = to_point(r);
            std::size_t best{};
            T d = std::numeric_limits<T>::infinity();
            for (std::size_t i{}; i < ours.size(); ++i)
                if (const auto di = distance(ours[i], p); di < d)
                    d = di, best = i;
            EXPECT_LE(d, tol);
            EXPECT_FALSE(used[best]);
            used[best] = true;
            map.push_back(best);
        }
        return map;
    }

    /// The patch and its (u, v) at parameter t along half-edge h (from its origin).
    auto on_edge(const LimitSurface<T> &s, std::size_t h, T t) -> std::tuple<std::size_t, T, T>
    {
        const auto &topo = s.topology();
        const auto f = topo.face(h), n = topo.face_size(f), i = h - topo.face_begin(f);
        if (n == 4)
        {
            switch (i)
            {
            case 0: return {s.patch(f), t, T(0)};
            case 1: return {s.patch(f), T(1), t};
            case 2: return {s.patch(f), 1 - t, T(1)};
            default: return {s.patch(f), T(0), 1 - t};
            }
        }
        if (t <= T(0.5))
            return {s.patch(f, i), 2 * t, T(0)};
        return {s.patch(f, (i + 1) % n), T(0), 2 * (1 - t)};
    }

    /// Patch and (u, v) of the corner of the first face around vertex v.
    auto at_vertex(const LimitSurface<T> &s, std::size_t h) -> std::tuple<std::size_t, T, T>
    {
        const auto &topo = s.topology();
        const auto f = topo.face(h), n = topo.face_size(f), i = h - topo.face_begin(f);
        if (n != 4)
            return {s.patch(f, i), T(0), T(0)};
        constexpr std::array<std::array<T, 2>, 4> uv{{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
        return {s.patch(f), uv[i][0], uv[i][1]};
    }

    void expect_issue(CageIssue issue, std::size_t index, const Cage<T> &cage)
    {
        try
        {
            [[maybe_unused]] const Topology t{cage};
            FAIL("cage accepted");
        }
        catch (const subd::CageError &e)
        {
            EXPECT_EQ(int(e.issue()), int(issue));
            EXPECT_EQ(e.index(), index);
        }
    }
}

TEST(tests_subd, refusals)
{
    const points_vector<T, 3> v{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}};
    const auto nan = std::numeric_limits<T>::quiet_NaN();
    expect_issue(CageIssue::NoFace, subd::none, {v, {}, {}, {}});
    expect_issue(CageIssue::FaceTooSmall, 1, {v, {{0, 1, 2}, {0, 1}}, {}, {}});
    expect_issue(CageIssue::VertexOutOfRange, 0, {v, {{0, 1, 9}}, {}, {}});
    expect_issue(CageIssue::RepeatedVertex, 0, {v, {{0, 1, 1, 2}}, {}, {}});
    expect_issue(CageIssue::NonFiniteVertex, 2, {{{0, 0, 0}, {1, 0, 0}, {nan, 1, 0}}, {{0, 1, 2}}, {}, {}});
    expect_issue(CageIssue::NonManifoldEdge, 2, {v, {{0, 1, 2}, {1, 0, 3}, {0, 1, 4}}, {}, {}});
    expect_issue(CageIssue::InconsistentOrientation, 1, {v, {{0, 1, 2}, {0, 1, 3}, {1, 2, 3, 4}}, {}, {}});
    expect_issue(CageIssue::NonManifoldVertex, 0, {v, {{0, 1, 2}, {0, 3, 4}}, {}, {}});
    expect_issue(CageIssue::IsolatedVertex, 4, {v, {{0, 1, 2, 3}}, {}, {}});
    expect_issue(CageIssue::InteriorValenceTwo, 0, {v, {{0, 1, 2, 3}, {0, 3, 2, 1}}, {}, {}}); // a pillow
    expect_issue(CageIssue::CreaseNotAnEdge, 1, {v, {{0, 1, 2, 3}, {0, 3, 4}}, {{0, 1}, {0, 2}}, {}});
    expect_issue(CageIssue::CornerOutOfRange, 0, {v, {{0, 1, 2, 3}, {0, 3, 4}}, {}, {9}});
    // the reference cages are all accepted
    for (const auto &name : names)
        EXPECT_NO_THROW(Topology{read_cage(name)});
}

TEST(tests_subd, rules)
{
    const auto count = [](const Topology &t, VertexRule r) {
        std::size_t n{};
        for (std::size_t v{}; v < t.n_vertices(); ++v)
            n += t.rule(v) == r;
        return n;
    };
    // creases: darts at the end of the crease from a corner and at the ends of the two creases of
    // the top; corners where that crease meets the equator and the tagged one; creases along the
    // equator and at the meeting of the two creases of the top
    const Topology c{read_cage("creases")};
    EXPECT_EQ(count(c, VertexRule::Dart), 3u);
    EXPECT_EQ(count(c, VertexRule::Corner), 2u);
    EXPECT_EQ(count(c, VertexRule::Crease), 8u);
    EXPECT_EQ(count(c, VertexRule::Smooth), 13u);
    // grid: its four corners belong to a single face
    const Topology g{read_cage("grid")};
    EXPECT_EQ(count(g, VertexRule::Corner), 4u);
    EXPECT_EQ(count(g, VertexRule::Crease), 10u);
    // hexprism: twelve vertices of valence 3, closed (Euler characteristic 2); refined, the centres of
    // its two hexagons are extraordinary vertices of valence 6
    const Topology h{read_cage("hexprism")};
    for (std::size_t v{}; v < h.n_vertices(); ++v)
        EXPECT_EQ(h.valence(v), 3u);
    EXPECT_EQ(h.n_vertices() + h.n_faces(), h.n_edges() + 2);
    const Topology r{subd::refine(read_cage("hexprism"))};
    std::size_t six{};
    for (std::size_t v{}; v < r.n_vertices(); ++v)
        six += r.valence(v) == 6;
    EXPECT_EQ(six, 2u);
}

TEST(tests_subd, refine_as_opensubdiv)
{
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const auto ref = read_json(name + ".osd.json");
        const auto tol = 1e-13 * size_of(cage);
        for (std::size_t level = 1; level <= 2; ++level)
        {
            const auto r = subd::refine(cage, level);
            const auto &osd = ref[level == 1 ? "level1" : "level2"];
            ASSERT_EQ(r.vertices.size(), osd["vertices"].Size());
            ASSERT_EQ(r.faces.size(), osd["faces"].Size());
            const auto map = match(r.vertices, osd["vertices"], tol);
            for (std::size_t f{}; f < r.faces.size(); ++f)
            {
                const auto &of = osd["faces"][rapidjson::SizeType(f)];
                ASSERT_EQ(r.faces[f].size(), of.Size());
                for (std::size_t i{}; i < r.faces[f].size(); ++i)
                    EXPECT_EQ(r.faces[f][i], map[of[rapidjson::SizeType(i)].GetUint()]);
            }
        }
    }
}

TEST(tests_subd, limit_as_opensubdiv)
{
    // OpenSubdiv's limit masks at the vertices of the cage refined twice: positions everywhere but
    // at darts, normals where its tangent masks are exact (smooth vertices; boundary vertices of
    // two faces and corners of one, on their single side)
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const auto ref = read_json(name + ".osd.json");
        const auto r = subd::refine(cage, 2);
        const auto map = match(r.vertices, ref["level2"]["vertices"], 1e-13 * size_of(cage));
        const LimitSurface<T> s{r};
        const auto &t = s.topology();
        const auto &lim = ref["limit2"];
        std::size_t normals{};
        for (std::size_t i{}; i < map.size(); ++i)
        {
            const auto v = map[i];
            if (t.rule(v) != VertexRule::Dart)
                EXPECT_LE(distance(s.vertex_limit(v), to_point(lim["positions"][rapidjson::SizeType(i)])), 1e-13 * size_of(cage));
            const bool smooth = t.rule(v) == VertexRule::Smooth;
            const bool boundary = t.is_boundary(v) && ((t.rule(v) == VertexRule::Crease && t.face_count(v) == 2) || (t.rule(v) == VertexRule::Corner && t.face_count(v) == 1));
            if (!smooth && !boundary)
                continue;
            const auto [p, u, w] = at_vertex(s, t.outgoing(v));
            const auto n = to_point(lim["tangent1"][rapidjson::SizeType(i)]) ^ to_point(lim["tangent2"][rapidjson::SizeType(i)]);
            EXPECT_LE(angle(s.normal(p, u, w), n), 1e-11);
            ++normals;
        }
        EXPECT_GT(normals, 0u);
    }
}

TEST(tests_subd, patches_as_opensubdiv)
{
    // the limit surface and its derivatives at random points of OpenSubdiv's regular patches; on the
    // patches of an n-gon, OpenSubdiv measures the derivatives in the parameters of the n-gon, where
    // a patch spans half of them: twice ours per order
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const auto ref = read_json(name + ".osd.json");
        const LimitSurface<T> s{cage};
        ASSERT_EQ(s.patch_count(), ref["patches"]["count"].GetUint());
        const auto size = size_of(cage);
        constexpr std::array<std::array<std::size_t, 2>, 6> orders{{{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}}};
        constexpr std::array<const char *, 6> keys{"p", "du", "dv", "duu", "duv", "dvv"};
        for (const auto &sample : ref["patches"]["samples"].GetArray())
        {
            const auto p = sample["patch"].GetUint();
            const auto u = sample["u"].GetDouble(), v = sample["v"].GetDouble();
            const bool ngon = s.topology().face_size(s.patch_face(p)[0]) != 4;
            for (std::size_t k{}; k < orders.size(); ++k)
            {
                const auto expected = to_point(sample[keys[k]]);
                const auto scale = ngon ? std::ldexp(1., int(orders[k][0] + orders[k][1])) : 1.;
                EXPECT_LE(distance(scale * s.value(p, u, v, orders[k][0], orders[k][1]), expected), 1e-11 * (size + norm(expected)));
            }
        }
    }
}

TEST(tests_subd, vertex_limits_by_subdivision)
{
    // every vertex: its limit is where its vertex points go, subdividing its ring 80 times (the
    // slowest of these rings shrinks by 0.64 per step); darts included, where OpenSubdiv is off
    using namespace subd::limit_detail;
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const LimitSurface<T> s{cage};
        const auto &t = s.topology();
        for (std::size_t v{}; v < cage.vertices.size(); ++v)
        {
            auto piece = vertex_ring(cage, t, v);
            for (int level{}; level < 80; ++level)
            {
                const Topology tp{piece, false};
                const auto r = subd::refine(piece, tp);
                piece = vertex_ring(r, Topology{r, false}, 0);
            }
            EXPECT_LE(distance(s.vertex_limit(v), piece.vertices[0]), 1e-12 * size_of(cage));
        }
    }
}

TEST(tests_subd, normals_at_irregular_vertices)
{
    // the normal given at an irregular vertex is the limit of the normals of the patch around it: the
    // gap at 2^-depth shrinks like (λ₃ / λ₂)^depth, 0.82^depth at worst (a crease of three faces
    // on a side). Not given where λ₂ is tied to λ₃: in fans, the crease of four faces (on its
    // four patches); in creases, the crease of one face on a side, the corner of two faces below
    // the equator, the corner tagged without sharp edges (3 faces)
    const std::map<std::string, std::size_t> not_given{{"fans", 4}, {"creases", 6}};
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const LimitSurface<T> s{cage};
        const auto &t = s.topology();
        std::size_t none{};
        for (std::size_t h{}; h < t.n_half_edges(); ++h)
        {
            const auto v = t.origin(h);
            if (t.rule(v) == VertexRule::Smooth && t.valence(v) == 4)
                continue;
            const auto [p, u, w] = at_vertex(s, h);
            const auto toward = [&](T d) { return s.normal(p, u == 0 ? d : 1 - d, w == 0 ? d : 1 - d); };
            try
            {
                const auto n = s.normal(p, u, w);
                const auto far = angle(n, toward(std::ldexp(1., -20))), near = angle(n, toward(std::ldexp(1., -50)));
                EXPECT_LE(near, 1e-4);
                EXPECT_LE(near, std::max(far * std::pow(0.83, 30.), 1e-14));
            }
            catch (const std::domain_error &)
            {
                ++none;
            }
        }
        EXPECT_EQ(none, not_given.contains(name) ? not_given.at(name) : 0u);
    }
}

TEST(tests_subd, deep_near_extraordinary_vertices)
{
    // at 2^-600 and 2^-1070 from smooth vertices of valence 3 and 6 (along the diagonal of a patch
    // starting there), the derivatives of order k stay finite and scale like (2^k λ)^depth, λ the
    // subdominant eigenvalue; the normals are the one at the vertex. Each quarter is evaluated at its
    // own scale: at the scale of the cage, S_u x S_v underflowed from a depth of 500
    struct Case
    {
        std::string name;
        std::size_t levels, valence;
    };
    for (const auto &[name, levels, valence] : {Case{"box", 0, 3}, Case{"hexprism", 1, 6}})
    {
        INFO(name);
        auto cage = subd::refine(read_cage(name), levels);
        const Topology t{cage};
        std::size_t vertex = subd::none;
        for (std::size_t v{}; v < t.n_vertices() && vertex == subd::none; ++v)
            if (t.valence(v) == valence && t.rule(v) == VertexRule::Smooth)
                vertex = v;
        ASSERT_NE(vertex, subd::none);
        const auto f = t.face(t.outgoing(vertex));
        std::ranges::rotate(cage.faces[f], std::ranges::find(cage.faces[f], vertex)); // the vertex at (0, 0) of the patch
        const LimitSurface<T> s{cage};
        const auto p = s.patch(f);
        const T n = T(valence), theta = 2 * std::numbers::pi / n;
        const T lambda = (5 + std::cos(theta) + std::cos(theta / 2) * std::sqrt(2 * (9 + std::cos(theta)))) / 16;
        const auto normal = s.normal(p, 0., 0.);
        for (const int depth : {600, 1070})
        {
            const T d = std::ldexp(1., -depth), d1 = std::ldexp(1., -depth - 1);
            for (const auto [du, dv] : {std::array<std::size_t, 2>{1, 0}, {0, 1}, {2, 0}, {1, 1}})
            {
                if (du + dv == 2 && std::pow(4 * lambda, depth) > 1e300)
                    continue; // beyond the largest double: (4 λ)^1070 for valence 6
                const auto length = [](const P3 &x) { return std::hypot(x[0], x[1], x[2]); }; // no overflow of the squares
                const auto a = s.value(p, d, d, du, dv), b = s.value(p, d1, d1, du, dv);
                ASSERT_TRUE(std::isfinite(length(a)) && length(a) > 0);
                EXPECT_NEAR(length(b) / length(a), std::ldexp(lambda, int(du + dv)), 1e-12); // 2^order λ
            }
            EXPECT_LE(angle(s.normal(p, d, d), normal), 1e-12);
            EXPECT_LE(distance(s.value(p, d, d), s.vertex_limit(vertex)), 1e-14 * size_of(cage));
        }
    }
}

TEST(tests_subd, continuity)
{
    // the patches on both sides of an edge share its points; their normals agree across a smooth edge
    for (const auto &name : names)
    {
        INFO(name);
        const auto cage = read_cage(name);
        const LimitSurface<T> s{cage};
        const auto &t = s.topology();
        const auto size = size_of(cage);
        for (std::size_t e{}; e < t.n_edges(); ++e)
        {
            const auto h = t.edge_half_edge(e), g = t.twin(h);
            if (g == subd::none)
                continue;
            for (const T x : {0.1, 0.37, 0.5, 0.83})
            {
                const auto [p, u, v] = on_edge(s, h, x);
                const auto [q, a, b] = on_edge(s, g, 1 - x);
                EXPECT_LE(distance(s.value(p, u, v), s.value(q, a, b)), 1e-13 * size);
                if (!t.is_crease(e))
                    EXPECT_LE(angle(s.normal(p, u, v), s.normal(q, a, b)), 1e-9);
            }
        }
    }
}

TEST(tests_subd, errors)
{
    const LimitSurface<T> s{read_cage("box")};
    EXPECT_THROW(s.value(0, 0., 0., 1, 0), std::domain_error); // no derivative at an extraordinary vertex
    EXPECT_NO_THROW(s.value(0, 1e-300, 0., 1, 0));
    EXPECT_THROW(s.value(0, -0.1, 0.5), std::domain_error);
    EXPECT_THROW(s.value(0, 0.5, 1.1), std::domain_error);
    EXPECT_THROW(s.value(6, 0.5, 0.5), std::out_of_range);
    EXPECT_THROW(s.patch(0, 1), std::out_of_range);
    EXPECT_THROW(s.vertex_limit(8), std::out_of_range);
    // derivatives of every order beyond 3 vanish on a regular patch
    const LimitSurface<T> torus{read_cage("torus")};
    EXPECT_EQ(norm(torus.value(3, 0.3, 0.6, 4, 0)), 0.);
}
