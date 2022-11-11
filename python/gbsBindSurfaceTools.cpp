#include "gbsBindSurfaceTools.h"
#include <gbs/bsstools.h>

void gbs_bind_surfaceTools(py::module &m)
{

    py::enum_<gbs::SurfaceBound>(m, "SurfaceBound", py::arithmetic())
        .value("U_START", gbs::SurfaceBound::U_START)
        .value("V_START", gbs::SurfaceBound::V_START)
        .value("V_END", gbs::SurfaceBound::V_END)
        .value("U_END", gbs::SurfaceBound::U_END)
    ;

    gbs_bind_surfaceTools<double,1,false>(m);
    gbs_bind_surfaceTools<double,2,false>(m);
    gbs_bind_surfaceTools<double,3,false>(m);
    gbs_bind_surfaceTools<double,1,true>(m);
    gbs_bind_surfaceTools<double,2,true>(m);
    gbs_bind_surfaceTools<double,3,true>(m);
}