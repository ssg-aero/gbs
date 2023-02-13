#include "gbsBindShaping.h"

void gbs_bind_shaping(py::module &m)
{
    gbs_bind_shaping<double, 1>(m);
    gbs_bind_shaping<double, 2>(m);
    gbs_bind_shaping<double, 3>(m);
}