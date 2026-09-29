// C++20 module unit for gbs/math.ixx (GBS_USE_MODULES=ON), see vecop.cppm.
module;
#include <numbers>
#include <vector>
#include <algorithm>
#include <stdexcept>

export module math;
export import vecop; // make_range can use overloaded -

#define GBS_MODULE_EXPORT export
#include "../math.ixx"
