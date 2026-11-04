#pragma once
/**
 * @file real_format.h
 * @brief The shortest decimal form of a double that reads back exactly, with a point as decimal
 * separator whatever the locale of the process (shared by the STEP, STL and OBJ writers), and
 * the reading of a decimal real, alike on every platform (the STEP and OBJ readers).
 */
#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>
#include <system_error>

#if defined(_LIBCPP_VERSION)
#include <clocale>
#include <cstdio>
#include <cstdlib>
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif

namespace gbs::io
{
    /// Appends the shortest form of x that reads back exactly (`%g` style: 1.5, 1e-07, 3).
    inline void append_real(std::string &out, double x)
    {
        char buf[64];
#if defined(_LIBCPP_VERSION)
        // std::to_chars on floating point is unavailable with some libc++ (as from_chars, see
        // step/p21.h): the shortest of %.15g to %.17g that reads back exactly, in the "C" locale
        static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t(0));
        for (int prec = 15; prec <= 17; ++prec)
        {
            snprintf_l(buf, sizeof buf, c_locale, "%.*g", prec, x);
            if (strtod_l(buf, nullptr, c_locale) == x)
                break;
        }
        out += buf;
#else
        const auto [end, ec] = std::to_chars(buf, buf + sizeof buf, x);
        out.append(buf, end);
#endif
    }

    /// The shortest form of x that reads back exactly (see `append_real`).
    [[nodiscard]] inline std::string format_real(double x)
    {
        std::string s;
        append_real(s, x);
        return s;
    }

    namespace real_format_detail
    {
        /**
         * For a real token whose value is beyond the range of double: true if it is tiny (underflow),
         * false if it is huge (overflow). Read from the digits, as std::from_chars leaves the value unset:
         * the sign of the decimal exponent of the first significant digit.
         */
        inline bool real_underflows(std::string_view s)
        {
            if (s.starts_with('+') || s.starts_with('-'))
                s.remove_prefix(1);
            const auto e = std::min(s.find_first_of("Ee"), s.size());
            const auto mantissa = s.substr(0, e);
            std::int64_t exponent = 0; // saturated: beyond 10^9, only its sign matters
            if (e < s.size())
            {
                auto x = s.substr(e + 1);
                const bool negative = x.starts_with('-');
                if (negative || x.starts_with('+'))
                    x.remove_prefix(1);
                for (const char c : x)
                    if (c >= '0' && c <= '9' && exponent < 1'000'000'000)
                        exponent = exponent * 10 + (c - '0');
                if (negative)
                    exponent = -exponent;
            }
            const auto first = mantissa.find_first_of("123456789");
            if (first == std::string_view::npos)
                return true; // zero
            const auto dot = std::min(mantissa.find('.'), mantissa.size());
            const auto magnitude = first < dot ? std::int64_t(dot - first) - 1 : -std::int64_t(first - dot);
            return magnitude + exponent < 0;
        }
    }

    /**
     * Locale-independent conversion of a decimal real, alike on every platform: a leading '+' is
     * accepted, an overflow gives result_out_of_range, an underflow gives 0 (or the denormal), the
     * length is not limited. Returns invalid_argument for a malformed token.
     */
    inline std::errc parse_real(std::string_view s, double &out)
    {
        if (s.starts_with('+'))
        {
            s.remove_prefix(1);
            if (s.starts_with('+') || s.starts_with('-'))
                return std::errc::invalid_argument;
        }
#if defined(_LIBCPP_VERSION)
        // libc++ (conda-forge on macOS) marks floating-point from_chars unavailable:
        // use strtod_l with the "C" locale, as exact and locale independent.
        static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t(0));
        char small[64];
        std::string large;
        const char *text = small;
        if (s.size() < sizeof small)
        {
            std::ranges::copy(s, small);
            small[s.size()] = '\0';
        }
        else
            text = (large = s).c_str();
        char *end = nullptr;
        errno = 0;
        out = strtod_l(text, &end, c_locale);
        if (end != text + s.size())
            return std::errc::invalid_argument;
        return errno == ERANGE && std::isinf(out) ? std::errc::result_out_of_range : std::errc{};
#else
        const auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
        if (p != s.data() + s.size())
            return std::errc::invalid_argument;
        if (ec == std::errc::result_out_of_range && real_format_detail::real_underflows(s))
        {
            out = s.starts_with('-') ? -0. : 0.;
            return std::errc{};
        }
        return ec;
#endif
    }
}
