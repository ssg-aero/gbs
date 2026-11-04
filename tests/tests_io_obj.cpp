// Wavefront OBJ files of polygon meshes (gbs-io/obj.h): faces of any number of sides, lines and
// points, the forms of the references, exact round trips, the faults named by their line, a byte
// order mark and the line ends of every system, the flat polygons.
#include <doctest_gtest.hpp>
#include <gbs-io/obj.h>

#include <algorithm>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <vector>

using namespace gbs;
using namespace gbs::io;
namespace fs = std::filesystem;

namespace
{
    /// The ObjError of a text, or a failure if it reads.
    template <typename T = double>
    auto refusal(std::string_view text) -> ObjError
    {
        try
        {
            (void)read_obj<T>(text, "test.obj");
        }
        catch (const ObjError &e)
        {
            return e;
        }
        FAIL("read, not refused");
        return ObjError{"", 0, ""};
    }

    /// Whether a polygon has these vertices.
    auto has(std::span<const std::uint32_t> p, std::vector<std::uint32_t> v) -> bool { return std::ranges::equal(p, v); }

    /// The locale of the numbers of the process, put back on leaving.
    struct NumericLocale
    {
        std::string saved = std::setlocale(LC_NUMERIC, nullptr);
        ~NumericLocale() { std::setlocale(LC_NUMERIC, saved.c_str()); }
    };

    /// A temporary directory of its own, removed on leaving (a failure included).
    struct TempDir
    {
        fs::path path;
        TempDir()
        {
            std::random_device r;
            path = fs::temp_directory_path() / ("gbs_tests_io_obj_" + std::to_string(r()) + std::to_string(r()));
            fs::create_directories(path);
        }
        ~TempDir()
        {
            std::error_code e;
            fs::remove_all(path, e);
        }
    };

    /// A cube of quadrangles, with an edge as a line and a corner as a point.
    auto cube() -> PolygonMesh<double>
    {
        PolygonMesh<double> m;
        for (int i = 0; i < 8; ++i)
            m.vertices.push_back({double(i & 1), double((i >> 1) & 1), double((i >> 2) & 1)});
        m.faces = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
        m.lines = {{0, 1}, {1, 3, 7}};
        m.points = {0, 7};
        return m;
    }
}

TEST(io_obj, a_modeller_file)
{
    // as Blender writes it: an object, groups, normals and texture coordinates, smoothing, materials
    const auto r = read_obj(R"(# Blender 4.2
mtllib part.mtl
o Part
v 1.000000 1.000000 -1.000000
v 1.000000 -1.000000 -1.000000
v 1.0 1.0 1.0
v 1.0 -1.0 1.0 1.0
v -1 1 -1   0.5 0.5 0.5
vn -0.0000 1.0000 -0.0000
vt 0.625000 0.500000
s 0
usemtl Material
f 1/1/1 5/1/1 3/1/1
f 2//1 4//1 3//1 1//1
f -1 -2 -4
l 1 2
p 3 4
)");
    const auto &m = r.mesh;
    ASSERT_EQ(m.vertices.size(), 5u);
    EXPECT_EQ(m.vertices[3], (point<double, 3>{1.0, -1.0, 1.0})); // the weight not kept
    EXPECT_EQ(m.vertices[4], (point<double, 3>{-1.0, 1.0, -1.0})); // nor a colour
    ASSERT_EQ(m.faces.size(), 3u);
    EXPECT_TRUE(has(m.faces[0], {0, 4, 2}));
    EXPECT_TRUE(has(m.faces[1], {1, 3, 2, 0})); // a quadrangle
    EXPECT_TRUE(has(m.faces[2], {4, 3, 1}));    // relative references
    EXPECT_EQ(m.lines, (Polygons{{0, 1}}));
    EXPECT_EQ(m.points, (std::vector<std::uint32_t>{2, 3}));
    // what is skipped is counted, in the order met
    ASSERT_EQ(r.skipped.size(), 6u);
    EXPECT_EQ(r.skipped[0], (std::pair<std::string, std::size_t>{"mtllib", 1}));
    EXPECT_EQ(r.skipped[1].first, "o");
    EXPECT_EQ(r.skipped[2].first, "vn");
    EXPECT_EQ(r.skipped[3].first, "vt");
}

TEST(io_obj, round_trip_exact)
{
    // the shortest decimal forms read back bit for bit, n-gons, lines and points kept
    auto m = cube();
    std::mt19937_64 rng{3};
    std::uniform_real_distribution<double> u{-1e3, 1e3};
    for (auto &p : m.vertices)
        p = {u(rng), u(rng) * 1e-9, std::nextafter(u(rng), 0.0)};
    m.vertices.push_back({5e-324, -0.0, 1.0 / 3});
    m.faces.push_back({0, 1, 8, 3, 2}); // a pentagon
    const auto text = write_obj(m, "two lines\nof comment");
    EXPECT_TRUE(text.starts_with("# two lines\n# of comment\nv "));
    const auto back = read_obj(text).mesh;
    EXPECT_EQ(back, m);
    EXPECT_TRUE(std::signbit(back.vertices.back()[1]));
    EXPECT_EQ(write_obj(back, "two lines\nof comment"), text);
    // float meshes too, in the shortest form of a float (read as a float: one rounding)
    PolygonMesh<float> f{{{0.1f, 0.2f, 0.3f}, {1, 0, 0}, {0, 1, 0}}, {{0, 1, 2}}, {}, {}};
    EXPECT_EQ(write_obj(f), "v 0.1 0.2 0.3\nv 1 0 0\nv 0 1 0\nf 1 2 3\n");
    EXPECT_EQ(read_obj<float>(write_obj(f)).mesh, f);
    for (int k = 0; k < 2000; ++k)
    {
        const auto x = static_cast<float>(u(rng));
        PolygonMesh<float> g{{{x, std::nextafter(x, 0.f), x * 1e-30f}}, {}, {}, {0}};
        ASSERT_EQ(read_obj<float>(write_obj(g)).mesh, g);
    }
    // a trailing new line of the comment adds no empty comment line
    EXPECT_TRUE(write_obj(f, "one\n").starts_with("# one\nv "));
}

TEST(io_obj, whatever_the_locale)
{
    // a decimal comma in the locale of the process changes neither the text nor its reading
    const NumericLocale saved;
    const bool comma = std::setlocale(LC_NUMERIC, "fr_FR.UTF-8") || std::setlocale(LC_NUMERIC, "de_DE.UTF-8") ||
                       std::setlocale(LC_NUMERIC, "French_France.1252");
    const auto m = cube();
    PolygonMesh<double> half = m;
    for (auto &p : half.vertices)
        p = {p[0] / 2, p[1] / 3, p[2] + 0.25};
    const auto text = write_obj(half);
    EXPECT_EQ(text.find(','), std::string::npos);
    EXPECT_EQ(read_obj(text).mesh, half);
    if (!comma)
        WARN_MESSAGE(false, "no locale with a decimal comma here: written and read in the C locale only, not checked");
}

TEST(io_obj, continuations_comments_and_line_ends)
{
    const auto m = read_obj("v 0 0 0 # a comment\r\nv 1 0 0\r\nv 1 1 0\nv 0 1 0\n\n   \t\nf 1 2 \\\n 3 4\n# the end").mesh;
    EXPECT_EQ(m.vertices.size(), 4u);
    EXPECT_EQ(m.faces, (Polygons{{0, 1, 2, 3}}));
    // a reference to a vertex defined further on
    EXPECT_EQ(read_obj("f 1 2 3\nv 0 0 0\nv 1 0 0\nv 0 1 0").mesh.faces.size(), 1u);
    EXPECT_TRUE(read_obj("").mesh.vertices.empty());
}

// A UTF-8 byte order mark (written by some Windows tools) is skipped, not read as a statement
// "\xEF\xBB\xBFv" that would drop the first vertex and shift every face; the line ends of the old
// Mac OS (a lone \r) end lines as \n and \r\n do.
TEST(io_obj, byte_order_mark_and_line_ends)
{
    const auto b = read_obj("\xEF\xBB\xBFv 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 3\n");
    EXPECT_EQ(b.mesh.vertices.size(), 4u);
    EXPECT_EQ(b.mesh.vertices[0], (point<double, 3>{0, 0, 0}));
    EXPECT_TRUE(b.skipped.empty());
    const auto r = read_obj("v 0 0 0\rv 1 0 0\rv 0 1 0\rf 1 2 3\r");
    EXPECT_EQ(r.mesh.vertices.size(), 3u);
    EXPECT_EQ(r.mesh.faces, (Polygons{{0, 1, 2}}));
    EXPECT_EQ(refusal("v 0 0 0\rv 1 0 0\rf 1 2\r").line, 3u);
    EXPECT_EQ(refusal("v 0 0 0\r\nv 1 0 0\nf 1 2\r").line, 3u);
}

// Polygons, flat: built from lists, indexed and run through as spans, compared.
TEST(io_obj, polygons)
{
    Polygons p{{0, 1, 2}, {2, 3, 4, 5}};
    p.push_back(std::vector<std::size_t>{6, 7, 8});
    ASSERT_EQ(p.size(), 3u);
    EXPECT_FALSE(p.empty());
    EXPECT_EQ(p.indices.size(), 10u);
    EXPECT_EQ(p.offsets, (std::vector<std::size_t>{0, 3, 7, 10}));
    EXPECT_TRUE(has(p[1], {2, 3, 4, 5}));
    std::size_t sides = 0;
    for (const auto f : p)
        sides += f.size();
    EXPECT_EQ(sides, 10u);
    EXPECT_EQ(std::ranges::distance(p), 3);
    EXPECT_EQ(p, (Polygons{std::vector<std::vector<int>>{{0, 1, 2}, {2, 3, 4, 5}, {6, 7, 8}}}));
    p.clear();
    EXPECT_TRUE(p.empty());
    EXPECT_EQ(p, Polygons{});
}

TEST(io_obj, faults_named_by_their_line)
{
    const std::string v3 = "v 0 0 0\nv 1 0 0\nv 0 1 0\n";
    auto e = refusal(v3 + "f 1 2\n");
    EXPECT_EQ(e.line, 4u);
    EXPECT_EQ(e.reason(), "a face needs 3 vertices or more");
    EXPECT_EQ(std::string{e.what()}, "test.obj: line 4: a face needs 3 vertices or more");
    EXPECT_EQ(refusal(v3 + "f 1 2 0").reason(), "vertex reference 0 in \"0\" (they count from 1)");
    EXPECT_EQ(refusal(v3 + "f 1 2 -4").reason(), "vertex reference -4 before the first vertex (3 read so far)");
    EXPECT_EQ(refusal(v3 + "f 1 2 x").reason(), "malformed vertex reference \"x\"");
    EXPECT_EQ(refusal(v3 + "f 1 2 /3").reason(), "malformed vertex reference \"/3\"");
    EXPECT_EQ(refusal(v3 + "\n\nf 1 2 9").line, 6u); // never defined: the line that refers to it
    EXPECT_EQ(refusal(v3 + "f 1 2 9").reason(), "vertex 9 referred to, the file has 3");
    EXPECT_EQ(refusal("v 0 0").reason(), "a vertex needs three coordinates");
    EXPECT_EQ(refusal("v 0 0 1,5").reason(), "coordinate \"1,5\" is not a finite number");
    EXPECT_EQ(refusal("v 0 0 nan").reason(), "coordinate \"nan\" is not a finite number");
    EXPECT_EQ(refusal("v 0 0 1e999").reason(), "coordinate \"1e999\" is not a finite number");
    EXPECT_EQ(refusal(v3 + "l 1").reason(), "a line needs 2 vertices or more");
    EXPECT_EQ(refusal(v3 + "p").reason(), "a point element needs a vertex");
    EXPECT_EQ(refusal("v 0 0 \\\n0\nv 1 \\\n 2").line, 3u); // a statement on two lines: its first
    // after the coordinates, numbers only (a weight, a colour), at most 4
    EXPECT_EQ(refusal("v 0 0 0 f 1 2 3").reason(), "coordinate \"f\" is not a finite number");
    EXPECT_EQ(refusal("v 1 2 3 oops").reason(), "coordinate \"oops\" is not a finite number");
    EXPECT_EQ(refusal("v 1 2 3 1 1 1 1 1").reason(), "\"1\": a vertex has at most 4 numbers after its coordinates");
    EXPECT_EQ(read_obj("v 1 2 3 0.5 0.5 0.5 1").mesh.vertices.size(), 1u);
    // a hexadecimal real is not a decimal one, on every platform (strtod reads it)
    EXPECT_EQ(refusal("v 0x10 0 0").reason(), "coordinate \"0x10\" is not a finite number");
    EXPECT_EQ(refusal("v -0X1p-2 0 0").reason(), "coordinate \"-0X1p-2\" is not a finite number");
    // finite as a double, not as a float
    EXPECT_EQ(refusal<float>("v 1e39 0 0").reason(), "coordinate \"1e39\" is not a finite number (as a float)");
    EXPECT_EQ(read_obj<float>("v 1e-50 0 0").mesh.vertices[0][0], 0.f); // an underflow is 0
    // references at the ends of the integers
    EXPECT_EQ(refusal(v3 + "f 1 2 -9223372036854775808").reason(), "vertex reference -9223372036854775808 before the first vertex (3 read so far)");
    EXPECT_EQ(refusal(v3 + "f 1 2 9223372036854775807").reason(), "vertex reference 9223372036854775807 beyond the 4294967295 vertices a mesh holds");
    EXPECT_EQ(refusal(v3 + "f 1 2 99999999999999999999").reason(), "malformed vertex reference \"99999999999999999999\"");
    EXPECT_EQ(refusal(v3 + "f 1 2 4294967295").reason(), "vertex 4294967295 referred to, the file has 3");
    // what read_obj refuses is not written
    auto m = cube();
    m.faces.push_back({0, 1});
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
    m = cube();
    m.lines.push_back({0, 8});
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
    m = cube();
    m.vertices[3][1] = std::numeric_limits<double>::infinity();
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
}

TEST(io_obj, files)
{
    const TempDir temp;
    const auto &dir = temp.path;
    const auto path = dir / "cube.obj";
    write_obj_file(cube(), path, "a cube");
    const auto back = read_obj_file(path);
    EXPECT_EQ(back.mesh, cube());
    EXPECT_TRUE(back.skipped.empty());
    EXPECT_THROW((void)read_obj_file(dir / "none.obj"), fs::filesystem_error);
    EXPECT_THROW((void)read_obj_file(dir), fs::filesystem_error);
    // the name of the file in the message of a fault
    {
        std::ofstream{dir / "bad.obj"} << "v 0 0 0\nf 1 1\n";
    }
    try
    {
        (void)read_obj_file(dir / "bad.obj");
        FAIL("read, not refused");
    }
    catch (const ObjError &e)
    {
        EXPECT_NE(std::string{e.what()}.find("bad.obj: line 2"), std::string::npos);
    }
}
