#include <doctest_gtest.hpp>
#include <gbs-io/step/p21.h>

#include <atomic>
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <new>
#include <optional>
#include <ranges>
#include <string>
#include <thread>

using namespace gbs::step;

// Heap allocations of the test binary, to check that a type lookup does not allocate.
namespace
{
    std::atomic<std::size_t> allocations{0};
}

void *operator new(std::size_t n)
{
    ++allocations;
    if (void *p = std::malloc(n ? n : 1))
        return p;
    throw std::bad_alloc{};
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }

namespace
{
    bool valid_utf8(std::string_view s)
    {
        for (std::size_t i = 0; i < s.size();)
        {
            const auto c = static_cast<unsigned char>(s[i]);
            const std::size_t n = c < 0x80 ? 0 : (c >> 5) == 6 ? 1 : (c >> 4) == 14 ? 2 : (c >> 3) == 30 ? 3 : 4;
            if (n == 4 || i + n >= s.size() + (n == 0 ? 1 : 0))
                return false;
            std::uint32_t cp = n == 0 ? c : c & (0x3Fu >> n);
            for (std::size_t k = 1; k <= n; ++k)
            {
                const auto d = static_cast<unsigned char>(s[i + k]);
                if ((d & 0xC0) != 0x80)
                    return false;
                cp = cp << 6 | (d & 0x3Fu);
            }
            constexpr std::uint32_t min_cp[] = {0, 0x80, 0x800, 0x10000};
            if (cp < min_cp[n] || (cp >= 0xD800 && cp < 0xE000) || cp > 0x10FFFF)
                return false;
            i += n + 1;
        }
        return true;
    }

    std::string wrap(const std::string &data, const std::string &schema = "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }")
    {
        return "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n"
               "FILE_NAME('part.stp','2026-10-09T10:00:00',('me'),(''),'','','');\n"
               "FILE_SCHEMA(('" + schema + "'));\nENDSEC;\nDATA;\n" + data + "ENDSEC;\nEND-ISO-10303-21;\n";
    }

    P21File parse_ok(const std::string &text)
    {
        auto r = parse_p21(text);
        if (!r)
            FAIL("line " << r.error().line << " col " << r.error().column << ": " << r.error().message);
        return std::move(*r);
    }

    P21Error parse_err(const std::string &text)
    {
        auto r = parse_p21(text);
        REQUIRE_FALSE(r.has_value());
        return r.error();
    }
}

TEST(tests_step_p21, header_and_simple_instances)
{
    auto f = parse_ok(wrap(
        "#10=CARTESIAN_POINT('origin',(0.,1.5,-2.E-3));\n"
        "#11=DIRECTION('',(0.,0.,1.));\n"
        "#12=AXIS2_PLACEMENT_3D('',#10,#11,$);\n"));
    ASSERT_EQ(f.size(), 3);
    ASSERT_EQ(f.schemas().size(), 1);
    ASSERT_EQ(f.schemas()[0], "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }");
    auto header = f.header();
    ASSERT_EQ(header.size(), 3);
    ASSERT_EQ(header[1].type(), "FILE_NAME");
    ASSERT_EQ(header[1].args()[0].as_string(), "part.stp");

    auto p = f.instance(10);
    ASSERT_EQ(p.id(), 10);
    ASSERT_EQ(p.type(), "CARTESIAN_POINT");
    ASSERT_FALSE(p.is_complex());
    ASSERT_EQ(p.line(), 8);
    ASSERT_EQ(p.args()[0].as_string(), "origin");
    auto xyz = p.args()[1].as_list();
    ASSERT_EQ(xyz.size(), 3);
    ASSERT_DOUBLE_EQ(xyz[0].as_real(), 0.);
    ASSERT_DOUBLE_EQ(xyz[1].as_real(), 1.5);
    ASSERT_DOUBLE_EQ(xyz[2].as_real(), -2e-3);

    auto ax = f.instance(12);
    ASSERT_EQ(ax.args()[1].as_ref(), 10);
    ASSERT_EQ(ax.args()[2].as_ref(), 11);
    ASSERT_TRUE(ax.args()[3].is_null());
    ASSERT_TRUE(f.instance(ax.args()[1].as_ref()).has_type("cartesian_point")); // case-insensitive lookup

    ASSERT_FALSE(f.contains(99));
    ASSERT_FALSE(f.find(99).has_value());
    ASSERT_THROW((void)f.instance(99), P21AccessError);
}

TEST(tests_step_p21, values)
{
    auto f = parse_ok(wrap(
        "#1=THING(42,-7,+3,1.,1.E-07,-2.5E+3,.5,.T.,.F.,.U.,.UNSPECIFIED.,$,*,\"0A1F\",(1,(2,3),()),LENGTH_MEASURE(1.E-07),POSITIVE_LENGTH_MEASURE(2));\n"));
    auto a = f.instance(1).args();
    ASSERT_EQ(a.size(), 17);
    ASSERT_EQ(a[0].as_int(), 42);
    ASSERT_EQ(a[1].as_int(), -7);
    ASSERT_EQ(a[2].as_int(), 3);
    ASSERT_DOUBLE_EQ(a[0].as_real(), 42.); // an integer reads as a real
    ASSERT_TRUE(a[3].is(ValueKind::Real));
    ASSERT_DOUBLE_EQ(a[3].as_real(), 1.);
    ASSERT_DOUBLE_EQ(a[4].as_real(), 1e-7);
    ASSERT_DOUBLE_EQ(a[5].as_real(), -2500.);
    ASSERT_DOUBLE_EQ(a[6].as_real(), 0.5);
    ASSERT_TRUE(a[7].as_bool());
    ASSERT_FALSE(a[8].as_bool());
    ASSERT_FALSE(a[9].as_logical().has_value());
    ASSERT_EQ(a[10].as_enum(), "UNSPECIFIED");
    ASSERT_TRUE(a[11].is(ValueKind::Omitted));
    ASSERT_TRUE(a[12].is(ValueKind::Derived));
    ASSERT_EQ(a[13].as_binary(), "0A1F");
    auto l = a[14].as_list();
    ASSERT_EQ(l.size(), 3);
    ASSERT_EQ(l[1].as_list()[1].as_int(), 3);
    ASSERT_TRUE(l[2].as_list().empty());
    size_t n = 0;
    for (auto v : l)
        n += v.is(ValueKind::List);
    ASSERT_EQ(n, 2);
    ASSERT_EQ(a[15].typed_name(), "LENGTH_MEASURE");
    ASSERT_DOUBLE_EQ(a[15].typed_value().as_real(), 1e-7);
    ASSERT_DOUBLE_EQ(a[15].as_measure(), 1e-7);
    ASSERT_DOUBLE_EQ(a[16].as_measure(), 2.);
    // wrong kind
    ASSERT_THROW((void)a[0].as_string(), P21AccessError);
    ASSERT_THROW((void)a[10].as_bool(), P21AccessError);
    ASSERT_THROW((void)a.at(17), P21AccessError);
}

TEST(tests_step_p21, strings)
{
    auto f = parse_ok(wrap(
        "#1=S('it''s');\n"
        "#2=S('back\\\\slash');\n"
        "#3=S('caf\\X2\\00E9\\X0\\');\n"
        "#4=S('caf\\X\\E9');\n"
        "#5=S('\\S\\i');\n"
        "#6=S('\\X4\\0001F600\\X0\\');\n"
        "#7=S('\\X2\\D83DDE00\\X0\\');\n"
        "#8=S('long\nline');\n"
        "#9=S('\\PA\\x');\n"
        "#10=S('');\n"));
    ASSERT_EQ(f.instance(1).args()[0].as_string(), "it's");
    ASSERT_EQ(f.instance(2).args()[0].as_string(), "back\\slash");
    ASSERT_EQ(f.instance(3).args()[0].as_string(), "caf\xC3\xA9");     // é
    ASSERT_EQ(f.instance(4).args()[0].as_string(), "caf\xC3\xA9");
    ASSERT_EQ(f.instance(5).args()[0].as_string(), "\xC3\xA9");        // 'i' + 128 = é
    ASSERT_EQ(f.instance(6).args()[0].as_string(), "\xF0\x9F\x98\x80"); // U+1F600
    ASSERT_EQ(f.instance(7).args()[0].as_string(), "\xF0\x9F\x98\x80"); // same, as a UTF-16 surrogate pair
    ASSERT_EQ(f.instance(8).args()[0].as_string(), "longline");        // end of line inside a string ignored
    ASSERT_EQ(f.instance(9).args()[0].as_string(), "x");
    ASSERT_EQ(f.instance(10).args()[0].as_string(), "");
}

TEST(tests_step_p21, complex_instances_and_type_index)
{
    auto f = parse_ok(wrap(
        "#20=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#30,#31,#32),.UNSPECIFIED.,.F.,.F.)\n"
        "B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)CURVE()\n"
        "GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.707106781186548,1.))\n"
        "REPRESENTATION_ITEM(''));\n"
        "#30=CARTESIAN_POINT('',(1.,0.,0.));\n"
        "#31=CARTESIAN_POINT('',(1.,1.,0.));\n"
        "#32=CARTESIAN_POINT('',(0.,1.,0.));\n"
        "#40=B_SPLINE_CURVE_WITH_KNOTS('',1,(#30,#31),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);\n"));
    auto c = f.instance(20);
    ASSERT_TRUE(c.is_complex());
    ASSERT_EQ(c.part_count(), 7);
    ASSERT_EQ(c.type(), "BOUNDED_CURVE");
    ASSERT_TRUE(c.has_type("RATIONAL_B_SPLINE_CURVE"));
    ASSERT_FALSE(c.has_type("LINE"));
    ASSERT_EQ(c.args("B_SPLINE_CURVE")[0].as_int(), 2);
    ASSERT_EQ(c.args("B_SPLINE_CURVE")[1].as_list().size(), 3);
    ASSERT_NEAR(c.args("RATIONAL_B_SPLINE_CURVE")[0].as_list()[1].as_real(), std::sqrt(0.5), 1e-15);
    ASSERT_EQ(c.args("REPRESENTATION_ITEM")[0].as_string(), "");
    ASSERT_THROW((void)c.args("LINE"), P21AccessError);

    // the type index lists a complex instance under each of its partial types
    ASSERT_EQ(f.instances_of("CARTESIAN_POINT").size(), 3);
    auto bs = f.instances_of("B_SPLINE_CURVE_WITH_KNOTS");
    ASSERT_EQ(bs.size(), 2);
    ASSERT_EQ(bs[0].id(), 20);
    ASSERT_EQ(bs[1].id(), 40);
    ASSERT_EQ(f.instances_of("RATIONAL_B_SPLINE_CURVE").size(), 1);
    ASSERT_TRUE(f.instances_of("PLANE").empty());

    // file order iteration
    std::vector<std::uint64_t> ids;
    for (auto i : f.instances())
        ids.push_back(i.id());
    ASSERT_TRUE((ids == std::vector<std::uint64_t>{20, 30, 31, 32, 40}));
}

TEST(tests_step_p21, comments_sections_and_sparse_ids)
{
    const std::string text =
        "ISO-10303-21;\n/* header comment */\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n"
        "FILE_NAME('a','',(''),(''),'','','');\nFILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\n"
        "ANCHOR;\n<a> = #1;\nENDSEC;\n"
        "DATA;\n#1 = /* inline */ POINT ( 1 , 2 ) ;\n#1000000000=POINT(3,4);\nENDSEC;\n"
        "DATA('second',('CONFIG_CONTROL_DESIGN'));\n#7=POINT(5,6);\nENDSEC;\n"
        "END-ISO-10303-21;\n";
    auto r = parse_p21(text);
    REQUIRE(r.has_value());
    ASSERT_EQ(r->size(), 3);
    ASSERT_EQ(r->schemas()[0], "CONFIG_CONTROL_DESIGN");
    ASSERT_EQ(r->instance(1000000000).args()[1].as_int(), 4); // sparse ids: hash index
    ASSERT_EQ(r->instance(7).args()[0].as_int(), 5);           // second DATA section
}

TEST(tests_step_p21, syntax_errors)
{
    auto e = parse_err(wrap("#1=POINT(1,2)\n#2=POINT(3,4);\n")); // missing ';'
    ASSERT_EQ(e.line, 9);
    ASSERT_NE(e.message.find("';'"), std::string::npos);

    e = parse_err(wrap("#1=S('unterminated);\n"));
    ASSERT_NE(e.message.find("unterminated string"), std::string::npos);

    e = parse_err(wrap("#1=POINT(1,,2);\n"));
    ASSERT_EQ(e.line, 8);
    ASSERT_EQ(e.column, 12);

    e = parse_err(wrap("#1=POINT(1);\n#1=POINT(2);\n"));
    ASSERT_NE(e.message.find("defined twice"), std::string::npos);

    e = parse_err("not a step file");
    ASSERT_NE(e.message.find("ISO-10303-21"), std::string::npos);

    e = parse_err(wrap("#1=POINT(1);\n").substr(0, 120)); // truncated
    ASSERT_FALSE(e.message.empty());

    e = parse_err(wrap("#1=();\n"));
    ASSERT_NE(e.message.find("empty complex"), std::string::npos);

    e = parse_err(wrap("#1=POINT(1.2.3);\n"));
    ASSERT_NE(e.message.find("invalid real"), std::string::npos);
}

TEST(tests_step_p21, read_file_and_dangling_reference)
{
    const auto path = std::filesystem::temp_directory_path() / "gbs_p21_test.stp";
    {
        std::ofstream out(path);
        out << wrap("#1=LINE('',#2,#3);\n#2=CARTESIAN_POINT('',(0.,0.,0.));\n");
    }
    auto f = read_p21(path);
    REQUIRE(f.has_value());
    auto line = f->instance(1);
    ASSERT_TRUE(f->contains(line.args()[1].as_ref()));
    ASSERT_THROW((void)f->instance(line.args()[2].as_ref()), P21AccessError); // #3 is missing
    ASSERT_FALSE(read_p21(path.string() + ".missing").has_value());
}

TEST(tests_step_p21, large_file_performance)
{
    // ~200 000 instances, about 15 MB of text
    std::string data;
    data.reserve(16u << 20);
    const int n = 50000;
    for (int i = 0; i < n; ++i)
    {
        const int b = 4 * i + 1;
        data += "#" + std::to_string(b) + "=CARTESIAN_POINT('',(" + std::to_string(i) + ".5,-1.25E-3," + std::to_string(i % 97) + ".));\n";
        data += "#" + std::to_string(b + 1) + "=DIRECTION('',(0.,0.,1.));\n";
        data += "#" + std::to_string(b + 2) + "=AXIS2_PLACEMENT_3D('',#" + std::to_string(b) + ",#" + std::to_string(b + 1) + ",$);\n";
        data += "#" + std::to_string(b + 3) + "=(GEOMETRIC_REPRESENTATION_ITEM()REPRESENTATION_ITEM('item " + std::to_string(i) + "'));\n";
    }
    const auto text = wrap(data);
    const auto t0 = std::chrono::steady_clock::now();
    auto f = parse_p21(text);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t0).count();
    REQUIRE(f.has_value());
    MESSAGE("parsed " << text.size() / (1 << 20) << " MB, " << f->size() << " instances in " << ms << " ms");
    ASSERT_EQ(f->size(), 4 * n);
    ASSERT_EQ(f->instances_of("AXIS2_PLACEMENT_3D").size(), n);
    ASSERT_DOUBLE_EQ(f->instance(4 * 12345 + 1).args()[1].as_list()[0].as_real(), 12345.5);
    ASSERT_LT(ms, 10000); // generous bound for slow CI runners
}

// ---- code review, lot 8 ----------------------------------------------------------

TEST(tests_step_p21, deep_nesting_is_an_error)
{
    // a few kilobytes of '(' or of 'T(' used to overflow the stack (a thread has 512 KiB on macOS)
    const std::size_t depth = 100000;
    const auto lists = wrap("#1=A(" + std::string(depth, '(') + std::string(depth, ')') + ");\n");
    std::string typed = "#1=A(";
    for (std::size_t i = 0; i < depth; ++i)
        typed += "T(";
    typed = wrap(typed + "1" + std::string(depth, ')') + ");\n");
    std::optional<P21Error> e_lists, e_typed;
    std::thread([&] {
        if (auto r = parse_p21(lists); !r)
            e_lists = r.error();
        if (auto r = parse_p21(typed); !r)
            e_typed = r.error();
    }).join();
    REQUIRE(e_lists.has_value());
    ASSERT_NE(e_lists->message.find("nesting too deep"), std::string::npos);
    ASSERT_EQ(e_lists->line, 8);
    ASSERT_EQ(e_lists->column, 5 + 256); // the 257th '(': the argument list at column 5 is the first level
    REQUIRE(e_typed.has_value());
    ASSERT_NE(e_typed->message.find("nesting too deep"), std::string::npos);

    // within the limit
    auto f = parse_ok(wrap("#1=A(" + std::string(255, '(') + "1" + std::string(255, ')') + ");\n"));
    auto v = f.instance(1).args()[0];
    for (int i = 0; i < 254; ++i)
        v = v.as_list()[0];
    ASSERT_EQ(v.as_list()[0].as_int(), 1);
}

TEST(tests_step_p21, list_index_out_of_range_throws)
{
    auto f = parse_ok(wrap("#1=PRODUCT_DEFINITION('id');\n#2=FOO('a',#7,#8);\n#9=CARTESIAN_POINT('');\n"));
    const auto a = f.instance(1).args();
    ASSERT_EQ(a.size(), 1);
    ASSERT_EQ(a[0].as_string(), "id");
    ASSERT_THROW((void)a[1], P21AccessError); // used to read the arguments of #2
    ASSERT_THROW((void)f.instance(9).args()[1], P21AccessError); // used to read past the item table
    std::string what;
    try
    {
        (void)a[2];
    }
    catch (const P21AccessError &e)
    {
        what = e.what();
    }
    ASSERT_NE(what.find("index 2"), std::string::npos);
    ASSERT_NE(what.find("size 1"), std::string::npos);
}

TEST(tests_step_p21, malformed_escapes_kept_as_text)
{
    // a Windows path in FILE_NAME used to reject the whole file
    const std::string text =
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\n"
        "FILE_NAME('C:\\Users\\me\\part.stp','',(''),(''),'','','');\n"
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n"
        "#1=S('D:\\Projects\\part.stp');\n"
        "#2=S('a\\b');\n"
        "#3=S('\\X2\\00E\\X0\\');\n" // 3 hexadecimal digits
        "#4=S('\\X\\G1');\n"
        "#5=S('ends with \\');\n"
        "#6=S('\\X2\\00E9');\n" // no terminator
        "#7=S('\\S\\');\n"
        "#8=S('\\PZ\\');\n"
        "#9=S('C:\\\\dir\\X\\E9t\\X2\\00E9\\X0\\\\S\\i');\n" // valid escapes still decoded
        "#10=S('\\S\\''x');\n"                                   // \S\ applied to an apostrophe
        "ENDSEC;\nEND-ISO-10303-21;\n";
    auto f = parse_ok(text);
    auto s = [&](std::uint64_t id) { return std::string(f.instance(id).args()[0].as_string()); };
    ASSERT_EQ(f.header()[1].args()[0].as_string(), "C:\\Users\\me\\part.stp");
    ASSERT_EQ(s(1), "D:\\Projects\\part.stp");
    ASSERT_EQ(s(2), "a\\b");
    ASSERT_EQ(s(3), "\\X2\\00E\\X0\\");
    ASSERT_EQ(s(4), "\\X\\G1");
    ASSERT_EQ(s(5), "ends with \\");
    ASSERT_EQ(s(6), "\\X2\\00E9");
    ASSERT_EQ(s(7), "\\S\\");
    ASSERT_EQ(s(8), "\\PZ\\");
    ASSERT_EQ(s(9), "C:\\dir\xC3\xA9t\xC3\xA9\xC3\xA9");
    ASSERT_EQ(s(10), "\xC2\xA7x");
}

TEST(tests_step_p21, utf8_byte_order_mark_skipped)
{
    auto f = parse_ok("\xEF\xBB\xBF" + wrap("#1=A(1);\n"));
    ASSERT_EQ(f.instance(1).line(), 8);
    auto e = parse_err("\xEF\xBB\xBF" + wrap("#1=A(1)\n#2=A(2);\n"));
    ASSERT_EQ(e.line, 9);
    ASSERT_EQ(e.column, 1);
}

TEST(tests_step_p21, reals_read_alike_on_every_platform)
{
    auto real = [](const std::string &token) { return parse_p21(wrap("#1=A(" + token + ");\n")); };
    auto value = [&](const std::string &token) {
        auto r = real(token);
        if (!r)
            FAIL(token << ": " << r.error().message);
        return r->instance(1).args()[0].as_real();
    };
    ASSERT_DOUBLE_EQ(value("+1.5"), 1.5); // used to be rejected by std::from_chars (libstdc++, MSVC)
    ASSERT_DOUBLE_EQ(value("+1.E2"), 100.);
    ASSERT_DOUBLE_EQ(value("-1.5E+1"), -15.);
    // overflow is an error (used to be read as an infinity with libc++)
    for (const auto &t : {std::string("1.E400"), std::string("-1.E400"), std::string("1.E99999999999999999999"), "1" + std::string(400, '0') + "."})
    {
        auto r = real(t);
        REQUIRE_FALSE(r.has_value());
        ASSERT_NE(r.error().message.find("out of range"), std::string::npos);
    }
    // underflow gives 0 or a denormal (used to be rejected by std::from_chars)
    ASSERT_EQ(value("1.E-400"), 0.);
    ASSERT_EQ(value("-1.E-400"), 0.);
    ASSERT_EQ(value("0." + std::string(400, '0') + "1"), 0.);
    ASSERT_EQ(value("1.E-99999999999999999999"), 0.);
    const double tiny = value("4.94065645841247E-324");
    ASSERT_TRUE(tiny >= 0. && tiny < 1e-300);
    // overflow or underflow, told from the digits when std::from_chars leaves the value unset
    ASSERT_TRUE(detail::real_underflows("1.E-400"));
    ASSERT_TRUE(detail::real_underflows("-123.E-400"));
    ASSERT_TRUE(detail::real_underflows(".5E-400"));
    ASSERT_TRUE(detail::real_underflows("0." + std::string(400, '0') + "1"));
    ASSERT_FALSE(detail::real_underflows("1.E400"));
    ASSERT_FALSE(detail::real_underflows("+1.E+400"));
    ASSERT_FALSE(detail::real_underflows("1" + std::string(400, '0') + ".E-5"));
    ASSERT_FALSE(detail::real_underflows("0.001E99999999999999999999"));
    // no length limit (used to fail from 128 characters with libc++)
    ASSERT_DOUBLE_EQ(value("1." + std::string(130, '0')), 1.);
    ASSERT_DOUBLE_EQ(value("0." + std::string(300, '3')), 1. / 3.);
}

TEST(tests_step_p21, instance_number_too_large)
{
    // 2^64 used to be read as #0
    auto e = parse_err(wrap("#18446744073709551616=CARTESIAN_POINT('',(0.,0.,0.));\n"));
    ASSERT_NE(e.message.find("instance number too large"), std::string::npos);
    ASSERT_EQ(e.line, 8);
    ASSERT_EQ(e.column, 1);
    e = parse_err(wrap("#5=LINE('',#99999999999999999999,$);\n"));
    ASSERT_NE(e.message.find("instance number too large"), std::string::npos);
    ASSERT_EQ(e.column, 12);
    auto f = parse_ok(wrap("#18446744073709551615=A(#18446744073709551615);\n"));
    ASSERT_EQ(f.instance(18446744073709551615ull).args()[0].as_ref(), 18446744073709551615ull);
}

TEST(tests_step_p21, enumerations_and_binaries_validated)
{
    auto e = parse_err(wrap("#1=B_SPLINE_CURVE_WITH_KNOTS('',1,(#2,#3),.UNSPEC\nIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);\n"));
    ASSERT_NE(e.message.find("enumeration"), std::string::npos);
    ASSERT_EQ(e.line, 8);
    e = parse_err(wrap("#1=A(.T,1.5);\n"));
    ASSERT_NE(e.message.find("enumeration"), std::string::npos);
    for (const char *b : {"\"hello world\"", "\"4F\"", "\"\"", "\"0f\"", "\"0F F\""})
    {
        e = parse_err(wrap(std::string("#1=A(") + b + ");\n"));
        ASSERT_NE(e.message.find("binary"), std::string::npos);
    }

    // underscores and digits; lowercase letters are uppercased, as in keywords; a binary may span lines
    auto f = parse_ok(wrap("#1=A(.A_1.,._X.,.unspecified.,\"0FF\nFF\r\nFF\",\"3\");\n#2=B(1);\n"));
    const auto a = f.instance(1).args();
    ASSERT_EQ(a[0].as_enum(), "A_1");
    ASSERT_EQ(a[1].as_enum(), "_X");
    ASSERT_EQ(a[2].as_enum(), "UNSPECIFIED");
    ASSERT_EQ(a[3].as_binary(), "0FFFFFF");
    ASSERT_EQ(a[4].as_binary(), "3");
    ASSERT_EQ(f.instance(2).line(), 11);

    // line numbers after a binary over three lines (were 2 lines short)
    e = parse_err(wrap("#1=A(\"0FF\nFF\nFF\");\n#2=B(1);\n#3=C(1)\n#4=D(1);\n"));
    ASSERT_EQ(e.line, 13);
    ASSERT_EQ(e.column, 1);
}

TEST(tests_step_p21, strings_always_valid_utf8)
{
    auto f = parse_ok(wrap(
        "#1=S('\\X2\\DC00\\X0\\');\n"                  // lone low surrogate
        "#2=S('\\X2\\D83D\\X0\\\\X2\\DE00\\X0\\');\n" // pair split over two adjacent blocks
        "#3=S('\\X2\\D83D0041\\X0\\');\n"              // high surrogate, then a BMP character
        "#4=S('\\X4\\7FFFFFFF\\X0\\');\n"              // beyond U+10FFFF
        "#5=S('caf\xE9');\n"                           // raw ISO 8859-1, the Part 21 default
        "#6=S('caf\xC3\xA9');\n"                       // raw UTF-8: kept
        "#7=S('\\X2\\D83D\\X0\\x');\n"                 // high surrogate, then text
        "#8=S('\\X4\\0000D800\\X0\\');\n"              // surrogate code point in a 32-bit block
        "#9=S('\xC0\x80\xED\xA0\x80\xF4\x90\x80\x80');\n" // overlong, encoded surrogate, beyond U+10FFFF
        "#10=S('\xE2\x82');\n"                         // truncated sequence
        "#11=S('\\X2\\D83D\\X0\\\n\\X2\\DE00\\X0\\');\n")); // pair split by a line break
    auto s = [&](std::uint64_t id) { return std::string(f.instance(id).args()[0].as_string()); };
    const std::string fffd = "\xEF\xBF\xBD";
    ASSERT_EQ(s(1), fffd);
    ASSERT_EQ(s(2), "\xF0\x9F\x98\x80");
    ASSERT_EQ(s(3), fffd + "A");
    ASSERT_EQ(s(4), fffd);
    ASSERT_EQ(s(5), "caf\xC3\xA9");
    ASSERT_EQ(s(6), "caf\xC3\xA9");
    ASSERT_EQ(s(7), fffd + "x");
    ASSERT_EQ(s(8), fffd);
    ASSERT_EQ(s(9), "\xC3\x80\xC2\x80"
                    "\xC3\xAD\xC2\xA0\xC2\x80"
                    "\xC3\xB4\xC2\x90\xC2\x80\xC2\x80");
    ASSERT_EQ(s(10), "\xC3\xA2\xC2\x82");
    ASSERT_EQ(s(11), "\xF0\x9F\x98\x80");
    for (std::uint64_t id = 1; id <= 11; ++id)
        ASSERT_TRUE(valid_utf8(s(id)));
}

TEST(tests_step_p21, lexer_ignores_the_global_locale)
{
    struct Restore
    {
        std::string saved;
        ~Restore() { std::setlocale(LC_ALL, saved.c_str()); }
    } restore{std::setlocale(LC_ALL, nullptr)};
    const auto text = wrap("#1=cartesian_point('',(1.5,2.5E-3,0.));\n#2=A(.unspecified.,.t.);\n");
    const auto latin1 = wrap("#1=CAF\xC9(1);\n#2=A(1)\xA0;\n"); // keyword with E acute, no-break space
    int tried = 0;
    for (const char *l : {"C", "fr_FR.ISO8859-1", "fr_FR.UTF-8", "tr_TR.ISO8859-9", "tr_TR.UTF-8", "de_DE.UTF-8", "French_France.1252", "Turkish_Turkey.1254"})
    {
        if (!std::setlocale(LC_ALL, l))
            continue;
        ++tried;
        INFO("locale " << l);
        auto f = parse_ok(text);
        ASSERT_TRUE(f.instance(1).has_type("CARTESIAN_POINT"));
        ASSERT_DOUBLE_EQ(f.instance(1).args()[1].as_list()[0].as_real(), 1.5);
        ASSERT_DOUBLE_EQ(f.instance(1).args()[1].as_list()[1].as_real(), 2.5e-3);
        ASSERT_EQ(f.instance(2).args()[0].as_enum(), "UNSPECIFIED");
        ASSERT_TRUE(f.instance(2).args()[1].as_bool());
        auto e = parse_err(latin1);
        ASSERT_NE(e.message.find("unexpected character"), std::string::npos);
    }
    MESSAGE(tried << " locales available");
}

TEST(tests_step_p21, iterators_outlive_their_view)
{
    auto f = parse_ok(wrap("#1=CARTESIAN_POINT('',(1.,2.,3.));\n#2=A(1);\n"));
    ASSERT_TRUE(std::ranges::borrowed_range<ListView>);
    ASSERT_TRUE(std::ranges::forward_range<ListView>);
    ASSERT_TRUE(std::ranges::borrowed_range<decltype(f.instances())>);
    auto it = f.instance(1).args()[1].as_list().begin(); // the ListView temporary dies here
    ++it;
    ASSERT_DOUBLE_EQ((*it).as_real(), 2.);
    auto in = f.instances().begin(); // so does the range of instances
    ++in;
    ASSERT_EQ((*in).id(), 2);
    auto found = std::ranges::find_if(f.instance(1).args()[1].as_list(), [](ValueView v) { return v.as_real() > 2.5; });
    ASSERT_DOUBLE_EQ((*found).as_real(), 3.); // a borrowed range: not std::ranges::dangling
    ASSERT_DOUBLE_EQ(f.instance(1).args()[1].as_list().back().as_real(), 3.);
    ASSERT_EQ(f.instances()[1].id(), 2);
}

TEST(tests_step_p21, read_text_whole_file)
{
    const auto dir = std::filesystem::temp_directory_path();
    const auto path = dir / "gbs_p21_read_text.stp";
    const auto content = wrap("#1=A('x');\r\n") + std::string("\0\xFF\r\n", 4);
    {
        std::ofstream out(path, std::ios::binary);
        out << content;
    }
    const auto text = detail::read_text(path);
    REQUIRE(text.has_value());
    ASSERT_EQ(*text, content); // binary read: bytes unchanged
    ASSERT_FALSE(detail::read_text(path.string() + ".missing").has_value());
    ASSERT_FALSE(detail::read_text(dir).has_value());
    std::ofstream(path, std::ios::binary | std::ios::trunc).close();
    ASSERT_EQ(detail::read_text(path), std::optional<std::string>(""));
    std::filesystem::remove(path);
}

TEST(tests_step_p21, value_nodes_are_compact)
{
    // 24 bytes per value made a parsed file 5 to 8 times the size of its text
    ASSERT_LE(sizeof(detail::Node), 16u);
}

TEST(tests_step_p21, type_lookup_does_not_allocate)
{
    auto f = parse_ok(wrap("#1=B_SPLINE_CURVE_WITH_KNOTS('',1,(#2,#3),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);\n"));
    const auto in = f.instance(1);
    const auto before = allocations.load();
    const bool ok = in.has_type("B_SPLINE_CURVE_WITH_KNOTS") && !in.has_type("RATIONAL_B_SPLINE_SURFACE_WITH_KNOTS") &&
                    in.part("B_SPLINE_CURVE_WITH_KNOTS").has_value();
    const auto after = allocations.load();
    ASSERT_TRUE(ok);
    ASSERT_EQ(after - before, 0u); // an uppercase copy of the name used to be made on every call
    // the lookup stays case-insensitive
    ASSERT_TRUE(in.has_type("b_spline_curve_with_knots"));
    ASSERT_TRUE(in.has_type("B_Spline_Curve_With_Knots"));
    ASSERT_FALSE(in.has_type("b_spline_surface"));
    ASSERT_EQ(f.instances_of("b_spline_curve_with_knots").size(), 1);
}
