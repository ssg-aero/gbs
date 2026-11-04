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
        py::enum_<SurfaceSide>(m, "SurfaceSide", "A side of a surface, at one of its bounds")
            .value("u_min", SurfaceSide::u_min)
            .value("u_max", SurfaceSide::u_max)
            .value("v_min", SurfaceSide::v_min)
            .value("v_max", SurfaceSide::v_max);
        py::class_<SurfaceJoin<T, 3>>(m, "SurfaceJoin3d", "What join_surfaces gives")
            .def_readonly("surface", &SurfaceJoin<T, 3>::surface, "On the parameters of a, b beyond the side they share")
            .def_readonly("map", &SurfaceJoin<T, 3>::map, "(u, v) of b -> (m0 u + m1 v + m2, m3 u + m4 v + m5) on surface")
            .def_readonly("knots_removed", &SurfaceJoin<T, 3>::knots_removed, "Times the knot of the junction was removed: 0 for C0, 1 for C1, ...")
            .def_readonly("deviation", &SurfaceJoin<T, 3>::deviation, "Bound of the distance between surface and a or b at the same parameters (<= tol)")
            .def_readonly("side_a", &SurfaceJoin<T, 3>::side_a, "The side of a joined")
            .def_readonly("side_b", &SurfaceJoin<T, 3>::side_b, "The side of b joined");
        using S = BSSurface<T, 3>;
        m.def("join_surfaces", py::overload_cast<const S &, const S &, T>(&join_surfaces<T, 3>),
              "a and b joined into one surface along a side they share within tol (the smoothest join of all such sides; a "
              "side collapsed to a point is none), a keeping its parameters, the knot of the junction removed as far as tol "
              "allows; None if they share no side, are not clamped, or face opposite ways",
              py::arg("a"), py::arg("b"), py::arg("tol"));
        m.def("join_surfaces", py::overload_cast<const S &, const S &, SurfaceSide, SurfaceSide, T>(&join_surfaces<T, 3>),
              "a and b joined into one surface along the side side_a of a and side_b of b; None if those sides are not one "
              "curve within tol",
              py::arg("a"), py::arg("b"), py::arg("side_a"), py::arg("side_b"), py::arg("tol"));
        auto curve_join = [&]<std::size_t dim>() {
            const auto name = "CurveJoin" + std::to_string(dim) + "d";
            py::class_<CurveJoin<T, dim>>(m, name.c_str(), "What join_curves gives")
                .def_readonly("curve", &CurveJoin<T, dim>::curve, "a, then b shifted to start where a ends")
                .def_readonly("knots_removed", &CurveJoin<T, dim>::knots_removed, "Times the knot of the junction was removed")
                .def_readonly("deviation", &CurveJoin<T, dim>::deviation, "Bound of the distance between curve and a or b");
            m.def("join_curves", &join_curves<T, dim>,
                  "a then b (its start within tol of a's end) as one curve, each at its own speed, the knot of the junction "
                  "removed as far as tol allows; None if they do not meet, or are not clamped",
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