#include "gbsBindCurveTools.h"

void gbs_bind_curveTools(py::module &m)
{
    gbs_bind_blossom<double, 1>(m);
    gbs_bind_blossom<double, 2>(m);
    gbs_bind_blossom<double, 3>(m);
    gbs_bind_blossom<double, 4>(m);

    gbs_bind_curveTools<double,1,false>(m);
    gbs_bind_curveTools<double,2,false>(m);
    gbs_bind_curveTools<double,3,false>(m);
    gbs_bind_curveTools<double,1,true>(m);
    gbs_bind_curveTools<double,2,true>(m);
    gbs_bind_curveTools<double,3,true>(m);
}