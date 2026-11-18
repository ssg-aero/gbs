// The C library defines minor(dev), major(dev) and makedev(maj, min) as function-like macros:
// <sys/sysmacros.h> with glibc (included by <sys/types.h> before glibc 2.28, as in the 2.17 sysroot
// of conda-forge), <sys/types.h> on macOS. A function or a lambda of gbs called by one of these
// names stops compiling after them (gbs-brep/fillet.h had a lambda `minor`). They are defined here
// before any header, as the C library has them where it can, as macros of the same shape elsewhere.
#if __has_include(<sys/sysmacros.h>)
#include <sys/sysmacros.h>
#endif
#if __has_include(<sys/types.h>)
#include <sys/types.h>
#endif
#ifndef minor
#define minor(dev) ((dev) & 0xffu)
#endif
#ifndef major
#define major(dev) (((dev) >> 8) & 0xffu)
#endif
#ifndef makedev
#define makedev(maj, min) (((maj) << 8) | (min))
#endif

#include <doctest_gtest.hpp>
#include <gbs-brep/brep>
#include <gbs-io/iges_brep.h>
#include <gbs-io/step/reader.h>
#include <gbs-io/step/writer.h>

TEST(tests_system_macros, headers_compile_after_them)
{
    // still the macros: the headers compiled above with them in force
    ASSERT_EQ(static_cast<unsigned>(minor(makedev(3u, 5u))), 5u);
    ASSERT_EQ(static_cast<unsigned>(major(makedev(3u, 5u))), 3u);
}
