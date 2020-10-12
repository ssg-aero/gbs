#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gbslib/bscurve.h>
#include <gbslib/bscinterp.h>

namespace py = pybind11;

// https://stackoverflow.com/questions/47487888/pybind11-template-class-of-many-types
template <typename T, size_t dim, bool rational>
void declare_curve(py::module &m)
{
       
        using Class = typename std::conditional<rational, gbs::BSCurveRational<T, dim>, gbs::BSCurve<T, dim>>::type;

        std::string typestr;
        std::string rationalstr;
        if(rational) rationalstr="Rational";
        T x;
        auto typeid_name = typeid(x).name();

        if     ( strcmp (typeid_name,"float")  == 0 )  typestr = "f";
        else if( strcmp (typeid_name,"double") == 0 ) typestr = "d";
        else                                  typestr = typeid_name ;

        std::string pyclass_name = std::string("BSCurve")+ rationalstr + std::to_string(dim) + "d_" + typestr;

        py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<
                const gbs::points_vector<T,dim+rational> &,
                const std::vector<T> &,
                size_t 
                >())
        .def("value", &Class::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("begin", &Class::begin,"Curve evaluation at begin",py::arg("d") = 0)
        .def("end", &Class::end,"Curve evaluation at end",py::arg("d") = 0)
        .def("degree", &Class::degree,"Curve's degree")
        .def("knotsFlats", &Class::knotsFlats,"Curve's knots")
        .def("insertKnot", &Class::insertKnot,"Insert knot with the given multiplicity",py::arg("u"),py::arg("m") = 1)
        .def("removeKnot", &Class::removeKnot,"Try to remove m times the given knot",py::arg("u"),py::arg("tol"),py::arg("m") = 1)
        .def("poles", &Class::poles,"Curve's poles")
        .def("reverse", &Class::reverse,"reverse curve orientation")
        .def("trim", &Class::trim,"Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones")
        .def("changeBounds", &Class::changeBounds,"Change parametrization to fit between k1 and k2")
        .def("bounds", &Class::bounds,"Returns curves's start stop values")
        ;
}

PYBIND11_MODULE(gbs, m) {

        declare_curve<double,3,false>(m);
        declare_curve<double,2,false>(m);
        declare_curve<double,1,false>(m);

        declare_curve<float,3,false>(m);
        declare_curve<float,2,false>(m);
        declare_curve<float,1,false>(m);

        declare_curve<double,3,true>(m);
        declare_curve<double,2,true>(m);
        declare_curve<double,1,true>(m);

         py::enum_<gbs::KnotsCalcMode>(m, "KnotsCalcMode", py::arithmetic())
        .value("EQUALY_SPACED", gbs::KnotsCalcMode::EQUALY_SPACED)
        .value("CHORD_LENGTH", gbs::KnotsCalcMode::CHORD_LENGTH)
        .value("CENTRIPETAL", gbs::KnotsCalcMode::CENTRIPETAL)
        ;

        m.def("interpolate_cn_3d_d", &gbs::interpolate<double,3>, "Cn interpolation");
        m.def("interpolate_cn_2d_d", &gbs::interpolate<double,2>, "Cn interpolation");
        m.def("interpolate_cn_1d_d", &gbs::interpolate<double,1>, "Cn interpolation");
        m.def("interpolate_cn_3d_f", &gbs::interpolate<float,3>, "Cn interpolation");
        m.def("interpolate_cn_2d_f", &gbs::interpolate<float,2>, "Cn interpolation");
        m.def("interpolate_cn_1d_f", &gbs::interpolate<float,1>, "Cn interpolation");
}