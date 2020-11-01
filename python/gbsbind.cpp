#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/extrema.h>

#include <tuple>

namespace py = pybind11;

// template <typename T, size_t dim, bool rational>
// inline void makePythonName(py::module &m)
// {

// }

// https://stackoverflow.com/questions/47487888/pybind11-template-class-of-many-types
template <typename T, size_t dim, bool rational>
inline void declare_bscurve(py::module &m)
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
        .def(py::init<const Class &>())
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
        .def("changeBounds",py::overload_cast<T,T>(&Class::changeBounds),"Change parametrization to fit between k1 and k2")
        .def("changeBounds",py::overload_cast<const std::array<T,2>&>(&Class::changeBounds),"Change parametrization to fit between k1 and k2")
        .def("bounds", &Class::bounds,"Returns curves's start stop values")
        ;
}

// template <typename T, size_t dim, bool rational>
// inline void declare(py::module &m,ptfc)
// {

// }
template <typename T, size_t i>
auto extrema_PC_(gbs::extrema_PC_result<double> &res,bool &ok,py::args args)
{
        try
        {
                auto p_crv = py::cast<const gbs::BSCurve<T, i> &>(args[0]);
                auto p_pnt = py::cast<const std::array<T, i>>(args[1]);
                auto tol = py::cast<T>(args[2]);
                auto r = gbs::extrema_PC<T, i>(p_crv, p_pnt, tol);
                res = {r.u,r.d};
                ok = true;
        }
        catch (const std::exception &e)
        {
        }
};

auto extrema_PC_(py::args args) -> gbs::extrema_PC_result<double>
{
        gbs::extrema_PC_result<double> res;
        bool ok = false;

        extrema_PC_<double, 1>(res, ok, args);
        extrema_PC_<double, 2>(res, ok, args);
        extrema_PC_<double, 3>(res, ok, args);

        extrema_PC_<float, 1>(res, ok, args);
        extrema_PC_<float, 2>(res, ok, args);
        extrema_PC_<float, 3>(res, ok, args);

        if (!ok)
        {
                throw std::runtime_error("wrong argument type");
        }
        return res;
}

PYBIND11_MODULE(pygbs, m) {

        declare_bscurve<double,3,false>(m);
        declare_bscurve<double,2,false>(m);
        declare_bscurve<double,1,false>(m);

        declare_bscurve<float,3,false>(m);
        declare_bscurve<float,2,false>(m);
        declare_bscurve<float,1,false>(m);

        declare_bscurve<double,3,true>(m);
        declare_bscurve<double,2,true>(m);
        declare_bscurve<double,1,true>(m);

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

        py::class_<gbs::extrema_PC_result<double> >(m,"extrema_PC_result")
        .def_readwrite("d", &gbs::extrema_PC_result<double>::d)
        .def_readwrite("u", &gbs::extrema_PC_result<double>::u)
        ;

        m.def("extrema_PC", &extrema_PC_);
}