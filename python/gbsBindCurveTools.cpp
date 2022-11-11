#include "gbsBindCurveTools.h"

void gbs_bind_curveTools(py::module &m)
{

    gbs_bind_curveTools<double,1,false>(m);
    gbs_bind_curveTools<double,2,false>(m);
    gbs_bind_curveTools<double,3,false>(m);
    gbs_bind_curveTools<double,1,true>(m);
    gbs_bind_curveTools<double,2,true>(m);
    gbs_bind_curveTools<double,3,true>(m);
}