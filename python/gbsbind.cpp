#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gbslib/bscurve.h>

namespace py = pybind11;

PYBIND11_MODULE(gbs, m) {
    py::class_<gbs::BSCurve3d_d>(m, "BSCurve3d_d")
        .def(py::init<
                const gbs::points_vector<double,3> &,
                const std::vector<double> &,
                size_t 
                >())
        .def("value", &gbs::BSCurve3d_d::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("begin", &gbs::BSCurve3d_d::begin,"Curve evaluation at begin",py::arg("d") = 0)
        .def("end", &gbs::BSCurve3d_d::end,"Curve evaluation at end",py::arg("d") = 0)
        .def("degree", &gbs::BSCurve3d_d::degree,"Curve's degree")
        .def("knotsFlats", &gbs::BSCurve3d_d::knotsFlats,"Curve's knots")
        .def("insertKnot", &gbs::BSCurve3d_d::insertKnot,"Insert knot with the given multiplicity",py::arg("u"),py::arg("m") = 1)
        .def("removeKnot", &gbs::BSCurve3d_d::removeKnot,"Try to remove m times the given knot",py::arg("u"),py::arg("tol"),py::arg("m") = 1)
        .def("poles", &gbs::BSCurve3d_d::poles,"Curve's poles")
        .def("reverse", &gbs::BSCurve3d_d::reverse,"reverse curve orientation")
        .def("trim", &gbs::BSCurve3d_d::trim,"Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones")
        .def("changeBounds", &gbs::BSCurve3d_d::changeBounds,"Change parametrization to fit between k1 and k2")
        .def("bounds", &gbs::BSCurve3d_d::bounds,"Returns curves's start stop values")
        ;

    py::class_<gbs::BSCurve2d_d>(m, "BSCurve2d_d")
        .def(py::init<
                const gbs::points_vector<double,2> &,
                const std::vector<double> &,
                size_t 
                >())
        .def("value", &gbs::BSCurve2d_d::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("begin", &gbs::BSCurve2d_d::begin,"Curve evaluation at begin",py::arg("d") = 0)
        .def("end", &gbs::BSCurve2d_d::end,"Curve evaluation at end",py::arg("d") = 0)
        .def("degree", &gbs::BSCurve2d_d::degree,"Curve's degree")
        .def("knotsFlats", &gbs::BSCurve2d_d::knotsFlats,"Curve's knots")
        .def("insertKnot", &gbs::BSCurve2d_d::insertKnot,"Insert knot with the given multiplicity",py::arg("u"),py::arg("m") = 1)
        .def("removeKnot", &gbs::BSCurve2d_d::removeKnot,"Try to remove m times the given knot",py::arg("u"),py::arg("tol"),py::arg("m") = 1)
        .def("poles", &gbs::BSCurve2d_d::poles,"Curve's poles")
        .def("reverse", &gbs::BSCurve2d_d::reverse,"reverse curve orientation")
        .def("trim", &gbs::BSCurve2d_d::trim,"Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones")
        .def("changeBounds", &gbs::BSCurve2d_d::changeBounds,"Change parametrization to fit between k1 and k2")
        .def("bounds", &gbs::BSCurve2d_d::bounds,"Returns curves's start stop values")
        ;

        py::class_<gbs::BSCurve<double,1>>(m, "BSCurve1d_d")
        .def(py::init<
                const gbs::points_vector<double,1> &,
                const std::vector<double> &,
                size_t 
                >())
        .def("value", &gbs::BSCurve<double,1>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("begin", &gbs::BSCurve<double,1>::begin,"Curve evaluation at begin",py::arg("d") = 0)
        .def("end", &gbs::BSCurve<double,1>::end,"Curve evaluation at end",py::arg("d") = 0)
        .def("degree", &gbs::BSCurve<double,1>::degree,"Curve's degree")
        .def("knotsFlats", &gbs::BSCurve<double,1>::knotsFlats,"Curve's knots")
        .def("insertKnot", &gbs::BSCurve<double,1>::insertKnot,"Insert knot with the given multiplicity",py::arg("u"),py::arg("m") = 1)
        .def("removeKnot", &gbs::BSCurve<double,1>::removeKnot,"Try to remove m times the given knot",py::arg("u"),py::arg("tol"),py::arg("m") = 1)
        .def("poles", &gbs::BSCurve<double,1>::poles,"Curve's poles")
        .def("reverse", &gbs::BSCurve<double,1>::reverse,"reverse curve orientation")
        .def("trim", &gbs::BSCurve<double,1>::trim,"Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones")
        .def("changeBounds", &gbs::BSCurve<double,1>::changeBounds,"Change parametrization to fit between k1 and k2")
        .def("bounds", &gbs::BSCurve<double,1>::bounds,"Returns curves's start stop values")
        ;
}