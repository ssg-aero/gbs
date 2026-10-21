#pragma once
/**
 * @file real_format.h
 * @brief The shortest decimal form of a double that reads back exactly, with a point as decimal
 * separator whatever the locale of the process (shared by the STEP, STL and OBJ writers).
 */
#include <charconv>
#include <string>

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
}
