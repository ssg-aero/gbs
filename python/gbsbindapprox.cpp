#include "gbsbindapprox.h"


void gbs_bind_approx(py::module &m)
{
        gbs_bind_approx<double,1>(m);
        gbs_bind_approx<double,2>(m);
        gbs_bind_approx<double,3>(m);
}