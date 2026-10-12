#include <doctest_gtest.hpp>
#include <gbs-io/step/p21.h>
#include <gbs-io/step/p21_writer.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <random>
#include <string>

using namespace gbs::step;
using V = P21Value;

namespace
{
    P21File parse(const std::string &s)
    {
        auto f = parse_p21(s);
        if (!f)
            throw std::runtime_error(f.error().message + " at line " + std::to_string(f.error().line));
        return std::move(*f);
    }

    bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
}

TEST(tests_step_p21_writer, values_read_back)
{
    P21Writer w;
    const auto p = w.add("cartesian_point", {V::string("origin"), V::reals({0., 1.5, -2.25})});
    const auto q = w.add("MIXED", {V::ref(p), V::integer(-42), V::enumeration("unspecified"), V::boolean(true), V::boolean(false),
                                   V::logical(std::nullopt), V::omitted(), V::derived(), V::list({}),
                                   V::list({V::integer(1), V::list({V::real(2.), V::string("x")})}),
                                   V::typed("length_measure", V::real(1e-7)), V::refs({p, p}), V::integers({3, 3})});
    const auto f = parse(w.str());
    ASSERT_EQ(f.size(), 2u);
    const auto pa = f.instance(p).args();
    ASSERT_EQ(f.instance(p).type(), "CARTESIAN_POINT");
    ASSERT_EQ(pa[0].as_string(), "origin");
    ASSERT_EQ(pa[1].as_list()[2].as_real(), -2.25);
    const auto a = f.instance(q).args();
    ASSERT_EQ(a.size(), 13u);
    ASSERT_EQ(a[0].as_ref(), p);
    ASSERT_EQ(a[1].as_int(), -42);
    ASSERT_EQ(a[2].as_enum(), "UNSPECIFIED");
    ASSERT_TRUE(a[3].as_bool());
    ASSERT_FALSE(a[4].as_bool());
    ASSERT_FALSE(a[5].as_logical().has_value());
    ASSERT_TRUE(a[6].is(ValueKind::Omitted));
    ASSERT_TRUE(a[7].is(ValueKind::Derived));
    ASSERT_TRUE(a[8].as_list().empty());
    ASSERT_EQ(a[9].as_list()[1].as_list()[1].as_string(), "x");
    ASSERT_EQ(a[10].typed_name(), "LENGTH_MEASURE");
    ASSERT_EQ(a[10].as_measure(), 1e-7);
    ASSERT_EQ(a[11].as_list().size(), 2u);
    ASSERT_EQ(a[12].as_list()[1].as_int(), 3);
}

TEST(tests_step_p21_writer, reals)
{
    // exact round trip of every double, a decimal point always present
    std::vector<double> xs{0., -0., 1., -1., 0.1, 1. / 3., 2. / 3., 1e-7, 1.5e-300, 4.9e-324, 1e300,
                           std::numeric_limits<double>::max(), std::numeric_limits<double>::denorm_min(),
                           std::numeric_limits<double>::epsilon(), 123456789012345678., 25.4, 0.707106781186548};
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> mant(-1., 1.);
    std::uniform_int_distribution<int> expo(-300, 300);
    for (int i = 0; i < 2000; ++i)
        xs.push_back(std::ldexp(mant(rng), expo(rng)));
    P21Writer w;
    const auto id = w.add("REALS", {V::reals(xs)});
    const auto f = parse(w.str());
    const auto l = f.instance(id).args()[0].as_list();
    ASSERT_EQ(l.size(), xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i)
        ASSERT_TRUE(same_bits(l[i].as_real(), xs[i]));
    // Part 21 forms
    ASSERT_EQ(V::real(1.).text(), "1.");
    ASSERT_EQ(V::real(-0.).text(), "-0.");
    ASSERT_EQ(V::real(0.5).text(), "0.5");
    ASSERT_EQ(V::real(1e-7).text(), "1.E-7");
    ASSERT_EQ(V::real(1e300).text(), "1.E300");
    ASSERT_THROW(V::real(std::numeric_limits<double>::infinity()), std::invalid_argument);
    ASSERT_THROW(V::real(std::nan("")), std::invalid_argument);
}

TEST(tests_step_p21_writer, strings)
{
    const std::vector<std::string> ss{
        "", "it's", "back\\slash", "caf\xC3\xA9",                  // é
        "\xE6\x97\xA5\xE6\x9C\xAC",                               // 日本
        "smile \xF0\x9F\x98\x80 ok",                               // 😀, beyond the BMP
        "line\nbreak", "mix \xC3\xA9\xF0\x9F\x98\x80\xC3\xA8 end"};
    P21Writer w;
    std::vector<std::uint64_t> ids;
    for (const auto &s : ss)
        ids.push_back(w.add("TEXT", {V::string(s)}));
    const auto text = w.str();
    ASSERT_EQ(text.find('\xC3'), std::string::npos); // the file itself is ASCII
    const auto f = parse(text);
    for (std::size_t i = 0; i < ss.size(); ++i)
        ASSERT_EQ(f.instance(ids[i]).args()[0].as_string(), ss[i]);
    ASSERT_EQ(V::string("it's").text(), "'it''s'");
    ASSERT_EQ(V::string("caf\xC3\xA9").text(), "'caf\\X2\\00E9\\X0\\'");
    ASSERT_THROW(V::string("bad \xC3"), std::invalid_argument);
    ASSERT_THROW(V::string("bad \xFF"), std::invalid_argument);
}

TEST(tests_step_p21_writer, complex_instances_header_and_lines)
{
    P21Header h;
    h.description = {"gbs test"};
    h.name = "part.stp";
    h.time_stamp = "2026-10-12T10:00:00";
    h.author = {"S\xC3\xA9" "bastien"};
    h.organization = {"ssg-aero"};
    h.schemas = {std::string(schema::ap242)};
    P21Writer w{h};
    // parts given in any order are written in alphabetical (ASCII) order
    const auto c = w.add_complex({{"rational_b_spline_curve", {V::reals({1., 0.5, 1.})}},
                                  {"b_spline_curve", {V::integer(2), V::list({}), V::enumeration("UNSPECIFIED"), V::boolean(false), V::boolean(false)}},
                                  {"bounded_curve", {}},
                                  {"representation_item", {V::string("")}}});
    // a long instance is wrapped after commas, outside strings
    std::vector<double> many(200, 1.25);
    const auto l = w.add("LONG", {V::string("a, string, with, commas, longer than a line of eighty characters for sure, yes"), V::reals(many)});
    const auto text = w.str();
    for (std::size_t pos = 0, next; pos < text.size(); pos = next + 1)
    {
        next = text.find('\n', pos);
        if (next == std::string::npos)
            break;
        ASSERT_LE(next - pos, 120u); // no endless line (a long string may exceed 80)
    }
    // ASCII order, as OCCT writes it: BOUNDED_CURVE before B_SPLINE_CURVE ('O' < '_')
    ASSERT_NE(text.find("#1=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(),.UNSPECIFIED.,.F.,.F.)RATIONAL_B_SPLINE_CURVE("), std::string::npos);
    const auto f = parse(text);
    ASSERT_EQ(f.schemas()[0], std::string(schema::ap242));
    const auto hdr = f.header();
    ASSERT_EQ(hdr[1].args()[0].as_string(), "part.stp");
    ASSERT_EQ(hdr[1].args()[1].as_string(), "2026-10-12T10:00:00");
    ASSERT_EQ(hdr[1].args()[2].as_list()[0].as_string(), "S\xC3\xA9" "bastien");
    ASSERT_TRUE(f.instance(c).is_complex());
    ASSERT_EQ(f.instance(c).args("RATIONAL_B_SPLINE_CURVE")[0].as_list()[1].as_real(), 0.5);
    ASSERT_EQ(f.instance(l).args()[0].as_string(), "a, string, with, commas, longer than a line of eighty characters for sure, yes");
    ASSERT_EQ(f.instance(l).args()[1].as_list().size(), 200u);
    // the time stamp defaults to now
    P21Writer w2;
    ASSERT_EQ(parse(w2.str()).header()[1].args()[1].as_string().size(), 19u);
}

TEST(tests_step_p21_writer, errors_and_file)
{
    P21Writer w;
    ASSERT_THROW(w.add("", {}), std::invalid_argument);
    ASSERT_THROW(w.add("BAD NAME", {}), std::invalid_argument);
    ASSERT_THROW(w.add_complex({}), std::invalid_argument);
    ASSERT_THROW(V::ref(0), std::invalid_argument);
    ASSERT_THROW(V::enumeration("not ok"), std::invalid_argument);
    ASSERT_THROW(V::typed("1bad", V::integer(1)), std::invalid_argument);
    const auto id = w.add("CARTESIAN_POINT", {V::string(""), V::reals({1., 2., 3.})});
    const auto path = (std::filesystem::temp_directory_path() / "gbs_p21_writer.stp").string();
    w.write(path);
    auto f = read_p21(path);
    ASSERT_TRUE(f.has_value());
    ASSERT_EQ(f->instance(id).args()[1].as_list()[2].as_real(), 3.);
    ASSERT_THROW(w.write("/no/such/dir/x.stp"), std::runtime_error);
}
