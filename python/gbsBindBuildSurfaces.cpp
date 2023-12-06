#include "gbsBindBuildSurfaces.h"


void gbs_bind_build_surface(py::module &m)
{
        gbs_bind_build_surface<double,1>(m);
        gbs_bind_build_surface<double,2>(m);
        gbs_bind_build_surface<double,3>(m);
}