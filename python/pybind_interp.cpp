#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gbs/bscinterp.h>
PYBIND11_MODULE(pygbs, m) {

    using T = double;
        
        m.def( "interpolate_function",
                py::overload_cast<const std::vector<T> &,
                const std::vector<T> &,
                size_t>(&gbs::interpolate<T>),
                py::arg("Q"),py::arg("u"),py::arg("p")
        );

         py::enum_<gbs::KnotsCalcMode>(m, "KnotsCalcMode", py::arithmetic())
        .value("EQUALLY_SPACED", gbs::KnotsCalcMode::EQUALLY_SPACED)
        .value("CHORD_LENGTH", gbs::KnotsCalcMode::CHORD_LENGTH)
        .value("CENTRIPETAL", gbs::KnotsCalcMode::CENTRIPETAL)
        ;

        m.def(  "interpolate_cn_3d", 
                py::overload_cast<const std::vector< gbs::constrType<T, 3, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<T,3>), 
                "Cn interpolation");
        m.def(  "interpolate_cn_2d", 
                py::overload_cast<const std::vector< gbs::constrType<T, 2, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<T,2>), 
                "Cn interpolation");
        m.def(  "interpolate_cn_1d", 
                py::overload_cast<const std::vector< gbs::constrType<T, 1, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<T,1>), 
                "Cn interpolation");
}