#pragma once
// Precompiled-header target for GBS_USE_PCH builds.
// Pulls in the four core gbs headers (vecop/math/basisfunctions/knotsfunctions)
// plus the STL headers they depend on, so every test TU avoids re-parsing them.
// This header must NOT be included by consumer code; it is set via
// target_precompile_headers() in CMake when GBS_USE_PCH=ON.
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <list>
#include <Eigen/Dense>
#include "gbs/vecop.ixx"
#include "gbs/math.ixx"
#include "gbs/basisfunctions.ixx"
#include "gbs/knotsfunctions.ixx"
