// C++20 module unit for gbs/vecop.ixx (GBS_USE_MODULES=ON).
// Module directives must not be inside preprocessor conditionals (P1857), so
// they live here, unconditionally; the declarations come from the plain header.
module;
#include <array>
#include <algorithm>
#include <gbs/execution.h>
#include <cmath>

export module vecop;

#define GBS_MODULE_EXPORT export
#include "../vecop.ixx"
