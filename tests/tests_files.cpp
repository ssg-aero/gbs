#include <doctest_gtest.hpp>
#include <gbs-io/files.h>
#include <gbs-io/fromjson.h>
#include <gbs-io/json.h>
#include <gbs-io/tojson.h>

#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <thread>
#include <vector>

using namespace gbs;
using namespace std::string_literals;
namespace fs = std::filesystem;

namespace
{
    /// A directory of its own under the temporary directory, removed at the end of the test.
    struct ScratchDir
    {
        fs::path path;
        ScratchDir()
        {
            std::random_device r;
            path = fs::temp_directory_path() / ("gbs_tests_files_" + std::to_string(r()) + std::to_string(r()));
            fs::create_directories(path);
        }
        ~ScratchDir()
        {
            std::error_code ignored;
            fs::remove_all(path, ignored);
        }
        auto entries() const -> std::size_t { return static_cast<std::size_t>(std::distance(fs::directory_iterator{path}, fs::directory_iterator{})); }
    };

    /// Doubles of every magnitude and digit count, the ends of the range among them.
    auto some_doubles(std::size_t n) -> std::vector<double>
    {
        std::mt19937_64 g{42};
        std::uniform_int_distribution<int> e{-300, 300};
        std::uniform_real_distribution<double> m{-1., 1.};
        std::vector<double> x{0., -0., 1., 0.1, 1. / 3., 5e-324, 2.2250738585072014e-308, 1.7976931348623157e308};
        while (x.size() < n)
            x.push_back(std::ldexp(m(g), e(g)));
        return x;
    }

    auto read(const fs::path &p) -> std::string
    {
        std::ifstream in(p, std::ios::binary);
        return {std::istreambuf_iterator<char>{in}, std::istreambuf_iterator<char>{}};
    }
}

// What write_js_doc writes, parse_file reads back bit for bit (the default parser of RapidJSON
// gives a neighbour of some of these doubles).
TEST(tests_files, json_file_round_trip_is_exact)
{
    ScratchDir dir;
    const auto x = some_doubles(20000);
    rapidjson::Document d;
    d.SetArray();
    for (const double v : x)
        d.PushBack(v, d.GetAllocator());
    const auto file = dir.path / "doubles.json";
    write_js_doc(d, file);

    rapidjson::Document back;
    parse_file(file, back);
    ASSERT_EQ(back.Size(), x.size());
    std::size_t exact{}, default_exact{};
    rapidjson::Document fast;
    fast.Parse(read(file).c_str());
    for (rapidjson::SizeType i{}; i < back.Size(); ++i)
    {
        exact += std::bit_cast<std::uint64_t>(back[i].GetDouble()) == std::bit_cast<std::uint64_t>(x[i]);
        default_exact += fast[i].GetDouble() == x[i];
    }
    ASSERT_EQ(exact, x.size());
    ASSERT_LT(default_exact, x.size()); // why kParseFullPrecisionFlag
}

// A curve written to a file and read back has the same poles and knots, bit for bit.
TEST(tests_files, curve_file_round_trip_is_exact)
{
    ScratchDir dir;
    std::mt19937_64 g{7};
    std::uniform_real_distribution<double> m{-1., 1.};
    points_vector<double, 3> poles;
    for (int i = 0; i < 10; ++i)
        poles.push_back({1e3 * m(g), m(g), 1e-3 * m(g)});
    const BSCurve<double, 3> crv{poles, {0., 0., 0., 0., 1. / 7., 2. / 7., 3. / 7., 4. / 7., 5. / 7., 6. / 7., 1., 1., 1., 1.}, 3};

    rapidjson::Document d;
    d.SetObject();
    auto v = make_json(crv, d.GetAllocator());
    d.AddMember("curve", v, d.GetAllocator());
    write_js_doc(d, dir.path / "curve.json");

    rapidjson::Document back;
    parse_file(dir.path / "curve.json", back);
    const auto crv_back = make_bscurve<double, 3>(back["curve"]);
    ASSERT_EQ(crv_back.poles(), crv.poles());
    ASSERT_EQ(crv_back.knotsFlats(), crv.knotsFlats());
}

TEST(tests_files, json_syntax_error_says_where_and_why)
{
    rapidjson::Document d;
    try
    {
        parse_json("{\"a\": [1, 2,\n  3 4]}", d, "text.json");
        FAIL("no JsonError");
    }
    catch (const JsonError &e)
    {
        ASSERT_EQ(e.line, 2u);
        ASSERT_EQ(e.column, 5u);
        ASSERT_EQ(e.offset, 17u); // the 4
        ASSERT_FALSE(e.reason.empty());
        ASSERT_EQ(std::string{e.what()}, "text.json: line 2, column 5: " + e.reason);
    }

    ScratchDir dir;
    std::ofstream{dir.path / "bad.json", std::ios::binary} << "[1,\r\n2,]";
    ASSERT_THROW(parse_file(dir.path / "bad.json", d), JsonError);
    ASSERT_THROW(parse_file(dir.path / "missing.json", d), std::runtime_error);
    ASSERT_THROW(parse_json("", d), JsonError); // empty text
}

TEST(tests_files, read_text)
{
    ScratchDir dir;
    const std::string text = "line 1\r\nline 2\n\0binary"s;
    std::ofstream{dir.path / "a.txt", std::ios::binary} << text;
    ASSERT_EQ(read_text(dir.path / "a.txt"), text); // bytes as they are
    ASSERT_FALSE(read_text(dir.path / "missing.txt"));
    ASSERT_FALSE(read_text(dir.path)); // a directory
}

TEST(tests_files, write_text_atomic_replaces_whole)
{
    ScratchDir dir;
    const auto file = dir.path / "a.txt";
    write_text_atomic(file, "first");
    ASSERT_EQ(read(file), "first");
    write_text_atomic(file, std::string(1 << 20, 'x'));
    ASSERT_EQ(read(file), std::string(1 << 20, 'x'));
    write_text_atomic(file, "");
    ASSERT_EQ(read(file), "");
    ASSERT_EQ(dir.entries(), 1u); // no temporary file left
}

TEST(tests_files, write_text_atomic_failure_leaves_everything_as_it_was)
{
    ScratchDir dir;
    // the renaming fails: a non-empty directory where the file goes
    const auto target = dir.path / "target";
    fs::create_directories(target / "inside");
    ASSERT_THROW(write_text_atomic(target, "text"), fs::filesystem_error);
    ASSERT_TRUE(fs::is_directory(target / "inside"));
    ASSERT_EQ(dir.entries(), 1u); // the temporary file removed

    // the temporary file cannot be created: no such directory
    try
    {
        write_text_atomic(dir.path / "no" / "such.txt", "text");
        FAIL("no filesystem_error");
    }
    catch (const fs::filesystem_error &e)
    {
        ASSERT_EQ(e.path1(), dir.path / "no" / "such.txt");
        ASSERT_TRUE(e.code());
    }
    ASSERT_EQ(dir.entries(), 1u);
}

// Concurrent writers: each renaming replaces the file whole, the last one wins. (On Windows a
// renaming may fail while another one replaces the file: the file stays whole.)
TEST(tests_files, write_text_atomic_concurrent_writers)
{
    ScratchDir dir;
    const auto file = dir.path / "a.txt";
    std::vector<std::thread> writers;
    for (char c : {'a', 'b', 'c', 'd'})
        writers.emplace_back([&, c] {
            for (int k = 0; k < 20; ++k)
                try
                {
                    write_text_atomic(file, std::string(1 << 16, c));
                }
                catch (const fs::filesystem_error &)
                {
                }
        });
    for (auto &w : writers)
        w.join();
    const auto text = read(file);
    ASSERT_EQ(text.size(), std::size_t{1 << 16});
    ASSERT_EQ(text, std::string(text.size(), text.front()));
    ASSERT_EQ(dir.entries(), 1u);
}
