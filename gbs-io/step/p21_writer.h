#pragma once

/**
 * @file p21_writer.h
 * @brief Writer of the ISO 10303-21 ("Part 21") exchange structure, the text
 * format of STEP files: counterpart of the parser of p21.h.
 *
 * Design: docs/sources/design/step_writer.md, section 2; architecture note
 * docs/sources/design/step_wpr01_p21_writer.md.
 *
 * Values are formatted when they are made (`P21Value`), instances are numbered
 * in the order they are added (#1, #2…), so that leaves are written first and
 * no forward reference is needed. Reals are written in their shortest form
 * that reads back to the same double, always with a decimal point and
 * independently of the locale; strings are encoded as Part 21 requires (non
 * ASCII characters as \X2\…\X0\ or \X4\…\X0\).
 */

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if defined(_LIBCPP_VERSION)
#include <clocale>
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif

namespace gbs::step
{
    /// Schema names for FILE_SCHEMA.
    namespace schema
    {
        inline constexpr std::string_view ap214 = "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }";
        inline constexpr std::string_view ap242 = "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }";
        inline constexpr std::string_view ap203 = "CONFIG_CONTROL_DESIGN";
    }

    /// A parameter value, already formatted in Part 21 syntax.
    class P21Value
    {
        std::string text_;
        explicit P21Value(std::string t) : text_{std::move(t)} {}

        static std::string format_real(double x)
        {
            if (!std::isfinite(x))
                throw std::invalid_argument("P21 real must be finite");
            char buf[64];
            std::string s;
#if defined(_LIBCPP_VERSION)
            // std::to_chars on floating point is unavailable with some libc++ (as from_chars, see p21.h):
            // the shortest of %.15g to %.17g that reads back exactly, in the "C" locale.
            static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t(0));
            for (int prec = 15; prec <= 17; ++prec)
            {
                snprintf_l(buf, sizeof buf, c_locale, "%.*g", prec, x);
                if (strtod_l(buf, nullptr, c_locale) == x)
                    break;
            }
            s = buf;
#else
            auto [end, ec] = std::to_chars(buf, buf + sizeof buf, x); // shortest exact form
            s.assign(buf, end);
#endif
            // Part 21: a decimal point is required, the exponent mark is E
            const auto e = s.find_first_of("eE");
            std::string mant = s.substr(0, e), expo = e == std::string::npos ? "" : s.substr(e + 1);
            if (mant.find('.') == std::string::npos)
                mant += '.';
            if (!expo.empty())
            {
                // drop a '+' and leading zeros of the exponent: 1.E-07 is fine, 1.E+300 -> 1.E300
                std::string sign;
                if (expo[0] == '+' || expo[0] == '-')
                {
                    if (expo[0] == '-')
                        sign = "-";
                    expo.erase(0, 1);
                }
                expo.erase(0, std::min(expo.find_first_not_of('0'), expo.size() - 1));
                return mant + "E" + sign + expo;
            }
            return mant;
        }

        static void append_utf8_escaped(std::string &out, std::string_view s)
        {
            // decode UTF-8, then: printable ASCII as is ('' and \\ doubled), others in \X2\ (BMP) or \X4\ runs
            std::vector<std::uint32_t> cps;
            for (std::size_t i = 0; i < s.size();)
            {
                const auto c = static_cast<unsigned char>(s[i]);
                std::uint32_t cp;
                std::size_t n;
                if (c < 0x80)
                    cp = c, n = 1;
                else if ((c >> 5) == 0x6)
                    cp = c & 0x1F, n = 2;
                else if ((c >> 4) == 0xE)
                    cp = c & 0x0F, n = 3;
                else if ((c >> 3) == 0x1E)
                    cp = c & 0x07, n = 4;
                else
                    throw std::invalid_argument("P21 string: invalid UTF-8");
                if (i + n > s.size())
                    throw std::invalid_argument("P21 string: truncated UTF-8");
                for (std::size_t k = 1; k < n; ++k)
                {
                    const auto cc = static_cast<unsigned char>(s[i + k]);
                    if ((cc >> 6) != 0x2)
                        throw std::invalid_argument("P21 string: invalid UTF-8");
                    cp = (cp << 6) | (cc & 0x3F);
                }
                cps.push_back(cp);
                i += n;
            }
            auto hex = [&out](std::uint32_t v, int digits) {
                static constexpr char h[] = "0123456789ABCDEF";
                for (int d = digits - 1; d >= 0; --d)
                    out += h[(v >> (4 * d)) & 0xF];
            };
            for (std::size_t i = 0; i < cps.size();)
            {
                const auto cp = cps[i];
                if (cp >= 0x20 && cp <= 0x7E)
                {
                    if (cp == '\'')
                        out += "''";
                    else if (cp == '\\')
                        out += "\\\\";
                    else
                        out += static_cast<char>(cp);
                    ++i;
                    continue;
                }
                const bool wide = cp > 0xFFFF;
                out += wide ? "\\X4\\" : "\\X2\\";
                while (i < cps.size() && !(cps[i] >= 0x20 && cps[i] <= 0x7E) && (cps[i] > 0xFFFF) == wide)
                    hex(cps[i++], wide ? 8 : 4);
                out += "\\X0\\";
            }
        }

        static bool valid_keyword(std::string_view k)
        {
            return !k.empty() && std::isalpha(static_cast<unsigned char>(k[0])) &&
                   std::ranges::all_of(k, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
        }

    public:
        [[nodiscard]] const std::string &text() const noexcept { return text_; }

        /// #id
        static P21Value ref(std::uint64_t id)
        {
            if (id == 0)
                throw std::invalid_argument("P21 reference #0");
            return P21Value{"#" + std::to_string(id)};
        }
        static P21Value integer(std::int64_t i) { return P21Value{std::to_string(i)}; }
        static P21Value real(double x) { return P21Value{format_real(x)}; }
        /// 'text' (UTF-8 in, Part 21 encoding out)
        static P21Value string(std::string_view s)
        {
            std::string out{"'"};
            append_utf8_escaped(out, s);
            out += '\'';
            return P21Value{std::move(out)};
        }
        /// .NAME.
        static P21Value enumeration(std::string_view name)
        {
            if (!valid_keyword(name))
                throw std::invalid_argument("P21 enumeration name: " + std::string(name));
            std::string u{name};
            std::ranges::transform(u, u.begin(), [](unsigned char c) { return char(std::toupper(c)); });
            return P21Value{"." + u + "."};
        }
        static P21Value boolean(bool b) { return P21Value{b ? ".T." : ".F."}; }
        /// .T. / .F. / .U. (nullopt)
        static P21Value logical(std::optional<bool> b) { return P21Value{!b ? ".U." : (*b ? ".T." : ".F.")}; }
        /// $
        static P21Value omitted() { return P21Value{"$"}; }
        /// *
        static P21Value derived() { return P21Value{"*"}; }
        /// ( a, b, … )
        static P21Value list(const std::vector<P21Value> &items)
        {
            std::string out{"("};
            for (std::size_t i = 0; i < items.size(); ++i)
            {
                if (i)
                    out += ',';
                out += items[i].text_;
            }
            out += ')';
            return P21Value{std::move(out)};
        }
        static P21Value list(std::initializer_list<P21Value> items) { return list(std::vector<P21Value>(items)); }
        /// List of references.
        static P21Value refs(const std::vector<std::uint64_t> &ids)
        {
            std::vector<P21Value> v;
            v.reserve(ids.size());
            for (auto id : ids)
                v.push_back(ref(id));
            return list(v);
        }
        /// List of reals.
        static P21Value reals(const std::vector<double> &xs)
        {
            std::vector<P21Value> v;
            v.reserve(xs.size());
            for (auto x : xs)
                v.push_back(real(x));
            return list(v);
        }
        /// List of integers.
        static P21Value integers(const std::vector<std::int64_t> &xs)
        {
            std::vector<P21Value> v;
            v.reserve(xs.size());
            for (auto x : xs)
                v.push_back(integer(x));
            return list(v);
        }
        /// TYPE(value), e.g. LENGTH_MEASURE(1.E-07)
        static P21Value typed(std::string_view type, const P21Value &value)
        {
            if (!valid_keyword(type))
                throw std::invalid_argument("P21 type name: " + std::string(type));
            std::string u{type};
            std::ranges::transform(u, u.begin(), [](unsigned char c) { return char(std::toupper(c)); });
            return P21Value{u + "(" + value.text_ + ")"};
        }
    };

    /// Header of a Part 21 file.
    struct P21Header
    {
        std::vector<std::string> description{""};
        std::string implementation_level{"2;1"};
        std::string name{};
        std::string time_stamp{}; ///< ISO 8601; empty: the current UTC time
        std::vector<std::string> author{""};
        std::vector<std::string> organization{""};
        std::string preprocessor_version{"gbs"};
        std::string originating_system{"gbs"};
        std::string authorization{};
        std::vector<std::string> schemas{std::string(schema::ap214)};
    };

    /**
     * @brief Writes a Part 21 file: header, then the data section built by add() / add_complex().
     */
    class P21Writer
    {
        std::vector<std::string> data_; // instance texts, #i at index i - 1
        P21Header header_;

        static std::string upper_keyword(std::string_view type)
        {
            std::string u{type};
            std::ranges::transform(u, u.begin(), [](unsigned char c) { return char(std::toupper(c)); });
            if (u.empty() || !std::isalpha(static_cast<unsigned char>(u[0])) ||
                !std::ranges::all_of(u, [](char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; }))
                throw std::invalid_argument("P21 entity name: " + std::string(type));
            return u;
        }

        static std::string args_text(const std::vector<P21Value> &args)
        {
            std::string s{"("};
            for (std::size_t i = 0; i < args.size(); ++i)
            {
                if (i)
                    s += ',';
                s += args[i].text();
            }
            return s + ")";
        }

        // Breaks a long instance after commas outside strings, so that lines stay readable.
        static void append_wrapped(std::string &out, const std::string &line, std::size_t width = 80)
        {
            std::size_t col = 0;
            bool in_string = false;
            for (std::size_t i = 0; i < line.size(); ++i)
            {
                const char c = line[i];
                out += c;
                ++col;
                if (c == '\'')
                    in_string = !in_string; // '' toggles twice
                if (!in_string && c == ',' && col > width - 8 && i + 1 < line.size())
                {
                    out += "\n  ";
                    col = 2;
                }
            }
            out += '\n';
        }

        static std::string now_utc()
        {
            const std::time_t t = std::time(nullptr);
            std::tm tm{};
#if defined(_WIN32)
            gmtime_s(&tm, &t);
#else
            gmtime_r(&t, &tm);
#endif
            char buf[32];
            std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S", &tm);
            return buf;
        }

    public:
        explicit P21Writer(P21Header h = {}) : header_{std::move(h)} {}

        [[nodiscard]] P21Header &header() noexcept { return header_; }
        [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

        /// Adds a simple instance TYPE(args); returns its #id.
        std::uint64_t add(std::string_view type, const std::vector<P21Value> &args)
        {
            data_.push_back(upper_keyword(type) + args_text(args));
            return data_.size();
        }
        std::uint64_t add(std::string_view type, std::initializer_list<P21Value> args) { return add(type, std::vector<P21Value>(args)); }

        /**
         * @brief Adds a complex instance (A(…)B(…)…); the partial entities are written in
         * alphabetical order, as Part 21 requires.
         */
        std::uint64_t add_complex(std::vector<std::pair<std::string, std::vector<P21Value>>> parts)
        {
            if (parts.empty())
                throw std::invalid_argument("P21 complex instance without part");
            for (auto &[t, a] : parts)
                t = upper_keyword(t);
            std::ranges::sort(parts, {}, [](const auto &p) { return p.first; });
            std::string s{"("};
            for (const auto &[t, a] : parts)
                s += t + args_text(a);
            data_.push_back(s + ")");
            return data_.size();
        }

        /// Text of the whole exchange structure.
        [[nodiscard]] std::string str() const
        {
            auto strings = [](const std::vector<std::string> &v) {
                std::vector<P21Value> l;
                for (const auto &x : v)
                    l.push_back(P21Value::string(x));
                return P21Value::list(l).text();
            };
            const auto &h = header_;
            std::string out = "ISO-10303-21;\nHEADER;\n";
            append_wrapped(out, "FILE_DESCRIPTION(" + strings(h.description) + "," + P21Value::string(h.implementation_level).text() + ");");
            append_wrapped(out, "FILE_NAME(" + P21Value::string(h.name).text() + "," +
                                    P21Value::string(h.time_stamp.empty() ? now_utc() : h.time_stamp).text() + "," + strings(h.author) + "," +
                                    strings(h.organization) + "," + P21Value::string(h.preprocessor_version).text() + "," +
                                    P21Value::string(h.originating_system).text() + "," + P21Value::string(h.authorization).text() + ");");
            append_wrapped(out, "FILE_SCHEMA(" + strings(h.schemas) + ");");
            out += "ENDSEC;\nDATA;\n";
            for (std::size_t i = 0; i < data_.size(); ++i)
                append_wrapped(out, "#" + std::to_string(i + 1) + "=" + data_[i] + ";");
            out += "ENDSEC;\nEND-ISO-10303-21;\n";
            return out;
        }

        /// Writes the file; throws std::runtime_error if it cannot be written.
        void write(const std::filesystem::path &path) const
        {
            std::ofstream f(path, std::ios::binary);
            if (!f)
                throw std::runtime_error("cannot open " + path.string() + " for writing");
            const auto s = str();
            f.write(s.data(), static_cast<std::streamsize>(s.size()));
            if (!f)
                throw std::runtime_error("cannot write " + path.string());
        }
    };
}
