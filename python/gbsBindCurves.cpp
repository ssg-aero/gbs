#include "gbsBindCurves.h"

void gbs_bind_curves(py::module &m)
{
        gbs_bind_curves<double,1>(m);
        gbs_bind_curves<double,2>(m);
        gbs_bind_curves<double,3>(m);
}