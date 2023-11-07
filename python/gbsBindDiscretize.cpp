#include "gbsBindDiscretize.h"

void gbs_bind_discretize(py::module &m)
{
    gbs_bind_discretize<double,1>(m);
    gbs_bind_discretize<double,2>(m);
    gbs_bind_discretize<double,3>(m);
}