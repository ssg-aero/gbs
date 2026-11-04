#include "gbsBindSurfaceTools.h"
#include <gbs/bsstools.h>
#include <gbs/bssjoin.h>

#include <string>

void gbs_bind_surfaceTools(py::module &m)
{

    py::enum_<gbs::SurfaceBound>(m, "SurfaceBound", py::arithmetic())
        .value("U_START", gbs::SurfaceBound::U_START)
        .value("V_START", gbs::SurfaceBound::V_START)
        .value("V_END", gbs::SurfaceBound::V_END)
        .value("U_END", gbs::SurfaceBound::U_END)
    ;

    // join_surfaces, join_curves (gbs/bssjoin.h): non-rational B-splines
    {
        using T = double;
        py::class_<SurfaceJoin<T, 3>>(m, "SurfaceJoin3d", "What join_surfaces gives")
            .def_readonly("surface", &SurfaceJoin<T, 3>::surface, "On the parameters of a, b beyond the side they share")
            .def_readonly("map", &SurfaceJoin<T, 3>::map, "(u, v) of b -> (m0 u + m1 v + m2, m3 u + m4 v + m5) on surface")
            .def_readonly("knots_removed", &SurfaceJoin<T, 3>::knots_removed, "Times the knot of the junction was removed: 0 for C0, 1 for C1, ...")
            .def_readonly("deviation", &SurfaceJoin<T, 3>::deviation, "Bound of the distance between surface and a or b at the same parameters");
        m.def("join_surfaces", &join_surfaces<T, 3>,
              "a and b joined into one surface along a side they share within tol, a keeping its parameters, the knot of the "
              "junction removed as far as tol allows; None if they share no side, are not clamped there, or face opposite ways",
              py::arg("a"), py::arg("b"), py::arg("tol"));
        auto curve_join = [&]<std::size_t dim>() {
            const auto name = "CurveJoin" + std::to_string(dim) + "d";
            py::class_<CurveJoin<T, dim>>(m, name.c_str(), "What join_curves gives")
                .def_readonly("curve", &CurveJoin<T, dim>::curve, "a, then b shifted to start where a ends")
                .def_readonly("knots_removed", &CurveJoin<T, dim>::knots_removed, "Times the knot of the junction was removed")
                .def_readonly("deviation", &CurveJoin<T, dim>::deviation, "Bound of the distance between curve and a or b");
            m.def("join_curves", &join_curves<T, dim>,
                  "a then b (its start within tol of a's end) as one curve, each at its own speed, the knot of the junction "
                  "removed as far as tol allows; None if they do not meet, or are not clamped there",
                  py::arg("a"), py::arg("b"), py::arg("tol"));
        };
        curve_join.template operator()<2>();
        curve_join.template operator()<3>();
    }

    gbs_bind_surfaceTools<double,1,false>(m);
    gbs_bind_surfaceTools<double,2,false>(m);
    gbs_bind_surfaceTools<double,3,false>(m);
    gbs_bind_surfaceTools<double,1,true>(m);
    gbs_bind_surfaceTools<double,2,true>(m);
    gbs_bind_surfaceTools<double,3,true>(m);
}