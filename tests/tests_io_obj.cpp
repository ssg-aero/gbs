// Wavefront OBJ files of polygon meshes (gbs-io/obj.h): faces of any number of sides, lines and
// points, the forms of the references, exact round trips, the faults named by their line.
#include <doctest_gtest.hpp>
#include <gbs-io/obj.h>

#include <clocale>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <string_view>

using namespace gbs;
using namespace gbs::io;
namespace fs = std::filesystem;

namespace
{
    /// The ObjError of a text, or a failure if it reads.
    auto refusal(std::string_view text) -> ObjError
    {
        try
        {
            (void)read_obj(text, "test.obj");
        }
        catch (const ObjError &e)
        {
            return e;
        }
        FAIL("read, not refused");
        return ObjError{"", 0, ""};
    }

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
    EXPECT_EQ(m.faces[0], (std::vector<std::size_t>{0, 4, 2}));
    EXPECT_EQ(m.faces[1], (std::vector<std::size_t>{1, 3, 2, 0})); // a quadrangle
    EXPECT_EQ(m.faces[2], (std::vector<std::size_t>{4, 3, 1}));    // relative references
    EXPECT_EQ(m.lines, (std::vector<std::vector<std::size_t>>{{0, 1}}));
    EXPECT_EQ(m.points, (std::vector<std::size_t>{2, 3}));
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
    // float meshes too
    PolygonMesh<float> f{{{0.1f, 0.2f, 0.3f}, {1, 0, 0}, {0, 1, 0}}, {{0, 1, 2}}, {}, {}};
    EXPECT_EQ(read_obj<float>(write_obj(f)).mesh, f);
}

TEST(io_obj, whatever_the_locale)
{
    // a decimal comma in the locale of the process changes neither the text nor its reading
    const std::string saved = std::setlocale(LC_NUMERIC, nullptr);
    const bool comma = std::setlocale(LC_NUMERIC, "fr_FR.UTF-8") || std::setlocale(LC_NUMERIC, "de_DE.UTF-8") ||
                       std::setlocale(LC_NUMERIC, "French_France.1252");
    const auto m = cube();
    PolygonMesh<double> half = m;
    for (auto &p : half.vertices)
        p = {p[0] / 2, p[1] / 3, p[2] + 0.25};
    const auto text = write_obj(half);
    EXPECT_EQ(text.find(','), std::string::npos);
    EXPECT_EQ(read_obj(text).mesh, half);
    std::setlocale(LC_NUMERIC, saved.c_str());
    if (!comma)
        MESSAGE("no locale with a decimal comma here: written and read in the C locale only");
}

TEST(io_obj, continuations_comments_and_line_ends)
{
    const auto m = read_obj("v 0 0 0 # a comment\r\nv 1 0 0\r\nv 1 1 0\nv 0 1 0\n\n   \t\nf 1 2 \\\n 3 4\n# the end").mesh;
    EXPECT_EQ(m.vertices.size(), 4u);
    EXPECT_EQ(m.faces, (std::vector<std::vector<std::size_t>>{{0, 1, 2, 3}}));
    // a reference to a vertex defined further on
    EXPECT_EQ(read_obj("f 1 2 3\nv 0 0 0\nv 1 0 0\nv 0 1 0").mesh.faces.size(), 1u);
    EXPECT_TRUE(read_obj("").mesh.vertices.empty());
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
    // what read_obj refuses is not written
    auto m = cube();
    m.faces[2] = {0, 1};
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
    m = cube();
    m.lines[0] = {0, 8};
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
    m = cube();
    m.vertices[3][1] = std::numeric_limits<double>::infinity();
    EXPECT_THROW((void)write_obj(m), std::invalid_argument);
}

TEST(io_obj, files)
{
    std::random_device r;
    const auto dir = fs::temp_directory_path() / ("gbs_tests_io_obj_" + std::to_string(r()) + std::to_string(r()));
    fs::create_directories(dir);
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
    fs::remove_all(dir);
}
