// C++20 module unit for gbs/knotsfunctions.ixx (GBS_USE_MODULES=ON), see vecop.cppm.
module;
#include <cassert>
#include <gbs/gbslib.h>
#include <list>
#include <stdexcept>
#include <utility>
#include <Eigen/Dense>

export module knots_functions;

import vecop;
import basis_functions;
import math;

#define GBS_MODULE_EXPORT export
#include "../knotsfunctions.ixx"
