#include <doctest_gtest.hpp>
#include <gbs-io/step/p21.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

using namespace gbs::step;

namespace
{
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

    e = parse_err(wrap("#1=S('\\X2\\00E\\X0\\');\n"));
    ASSERT_NE(e.message.find("escape"), std::string::npos);

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
