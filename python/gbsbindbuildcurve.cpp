#include "gbsbindbuildcurve.h"


void gbs_bind_build_curve(py::module &m)
{
        gbs_bind_build_curve<double,1>(m);
        gbs_bind_build_curve<double,2>(m);
        gbs_bind_build_curve<double,3>(m);
}