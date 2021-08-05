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
//https://pybind11.readthedocs.io/en/stable/advanced/classes.html#overriding-virtual-functions-in-python
// template < size_t dim>
// class PyCurve : public gbs::Curve<double,dim> {
// public:
//     /* Inherit the constructors */
//     using gbs::Curve<double,dim>::Curve;

//     /* Trampoline (need one for each virtual function) */
//     std::array<double, dim> value(double u, size_t d = 0) override {
//         PYBIND11_OVERRIDE_PURE(
//             std::array<double, dim>, /* Return type */
//             gbs::Curve<double,dim>,      /* Parent class */
//             value,          /* Name of function in C++ (must match Python name) */
//             u,d      /* Argument(s) */
//         );
//     }
// };
// https://stackoverflow.com/questions/47487888/pybind11-template-class-of-many-types
template <typename T, size_t dim, bool rational>
inline void declare_bscurve(py::module_ &m)
{
       
        using Class = typename std::conditional<rational, gbs::BSCurveRational<T, dim>, gbs::BSCurve<T, dim>>::type;
        using ClassBase =  gbs::Curve<T,dim>;

        std::string typestr;
        std::string rationalstr;
        if(rational) rationalstr="Rational";
        T x;
        auto typeid_name = typeid(x).name();

        if     ( strcmp (typeid_name,"float")  == 0 )  typestr = "_f";
        else if( strcmp (typeid_name,"double") == 0 ) typestr = "";
        else                                  typestr = typeid_name ;

        std::string pyclass_name = std::string("BSCurve")+ rationalstr + std::to_string(dim) + "d" + typestr;



        py::class_<Class, ClassBase>(m, pyclass_name.c_str())
        // py::class_<Class>(m, pyclass_name.c_str())
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
        .def("increaseDegree", &Class::increaseDegree,"Increment curve's degree of 1")
        ;
}




// template <typename T, size_t i>
// auto extrema_PC_(std::array<double,2> &res, bool &ok, py::args args)
// {
//         try
//         {
//                 auto p_crv = py::cast<const gbs::BSCurve<T, i> &>(args[0]);
//                 auto p_pnt = py::cast<const std::array<T, i>>(args[1]);
//                 auto tol = py::cast<T>(args[2]);
//                 // auto r = gbs::extrema_curve_point<T, i>(p_crv, p_pnt, tol);
//                 // res = {r.u, r.d};
//                 res = gbs::extrema_curve_point<T, i>(p_crv, p_pnt, tol);
//                 ok = true;
//         }
//         catch (const std::exception &e)
//         {
//         }
// };
// auto extrema_PC_(py::args args) -> std::array<double,2>
// {
//         std::array<double,2> res;
//         bool ok = false;

//         extrema_PC_<double, 1>(res, ok, args);
//         extrema_PC_<double, 2>(res, ok, args);
//         extrema_PC_<double, 3>(res, ok, args);

//         if (!ok)
//         {
//                 throw std::runtime_error("wrong argument type");
//         }
//         return res;
// }


auto extrema_curve_point_3d(py::args args) -> std::array<double,2>
{
        return gbs::extrema_curve_point<double,3>(
                py::cast<const gbs::Curve<double, 3> &>(args[0]),
                py::cast<const std::array<double, 3> >(args[1]),
                py::cast<double>(args[2])
                // ,
                // py::cast<nlopt::algorithm>(args[3]),
                // py::cast<size_t>(args[4])
        );
}

template <typename T, size_t d>
auto extrema_curve_point_(py::args args)
{
        return gbs::extrema_curve_point<T, d>(
            py::cast<const gbs::Curve<T, d> &>(args[0]),
            py::cast<const std::array<T, d>>(args[1]),
            py::cast<T>(args[2])
            // ,
            // py::cast<nlopt::algorithm>(args[3]),
            // py::cast<size_t>(args[4])
        );
}

template<typename T,size_t d>
auto conv_arr(const std::array<T,d> &arr_) -> std::array<double,d>
{
        std::array<double,d> res;
        std::transform(
                arr_.begin(),
                arr_.end(),
                res.begin(),
                [](const auto &v)
                {
                        return static_cast<double>(v);
                }

        );
        return res;
}

auto extrema_curve_point(py::args args) -> std::array<double,2>
{
        try
        {
                return extrema_curve_point_<double,3>(args);
        }
        catch (const std::exception &e)
        {
        }
        try
        {
                return extrema_curve_point_<double,2>(args);
        }
        catch (const std::exception &e)
        {
        }
        try
        {
                return extrema_curve_point_<double,1>(args);
        }
        catch (const std::exception &e)
        {
        }
        // try
        // {
        //         return conv_arr(extrema_curve_point_<float,3>(args));
        // }
        // catch (const std::exception &e)
        // {
        // }
        // try
        // {
        //         return conv_arr(extrema_curve_point_<float,2>(args));
        // }
        // catch (const std::exception &e)
        // {
        // }
        // try
        // {
        //         return conv_arr(extrema_curve_point_<float,1>(args));
        // }
        // catch (const std::exception &e)
        // {
        // }
        throw std::runtime_error("wrong argument type");
}

#include <gbs-render/vtkcurvesrender.h>
// inline void f_plot_curves_2d(const std::vector<std::shared_ptr<gbs::Curve<double,2>>> &crv_lst){gbs::plot(crv_lst);};
// inline void f_plot_curves(const std::vector<std::shared_ptr<gbs::Curve<double,3>>> &crv_lst){gbs::plot(crv_lst);};
inline void f_plot_curves_2d(const std::vector<gbs::BSCurve<double,2>> &crv_lst){gbs::plot(crv_lst);};
inline void f_plot_curves(const std::vector<gbs::BSCurve<double,3>> &crv_lst){gbs::plot(crv_lst);};


PYBIND11_MODULE(pygbs, m) {

        // py::class_<gbs::Curve<double,3> >(m, "Curve3d");
        py::class_<gbs::Curve<double,3> >(m, "Curve3d")
                .def("value", &gbs::Curve<double,3>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);
        py::class_<gbs::Curve<double,2> >(m, "Curve2d")
                .def("value", &gbs::Curve<double,2>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);
        py::class_<gbs::Curve<double,1> >(m, "Curve1d")
                .def("value", &gbs::Curve<double,1>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);

        declare_bscurve<double,3,false>(m);
        declare_bscurve<double,2,false>(m);
        declare_bscurve<double,1,false>(m);

        declare_bscurve<double,3,true>(m);
        declare_bscurve<double,2,true>(m);
        declare_bscurve<double,1,true>(m);

        // py::class_<gbs::Curve<float,3> >(m, "Curve3d_f")
        //         .def("value", &gbs::Curve<float,3>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);
        // py::class_<gbs::Curve<float,2> >(m, "Curve2d_f")
        //         .def("value", &gbs::Curve<float,2>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);
        // py::class_<gbs::Curve<float,1> >(m, "Curve1d_f")
        //         .def("value", &gbs::Curve<float,1>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0);

        // declare_bscurve<float,3,false>(m);
        // declare_bscurve<float,2,false>(m);
        // declare_bscurve<float,1,false>(m);

        // declare_bscurve<float,3,true>(m);
        // declare_bscurve<float,2,true>(m);
        // declare_bscurve<float,1,true>(m);

        py::class_<gbs::BSCfunction<double>>(m,"BSCfunction")
                .def(py::init<const gbs::BSCfunction<double> &>())
                .def("value",&gbs::BSCfunction<double>::value,"Function evaluation at givent parameter",py::arg("u"),py::arg("d") = 0)
                .def("basisCurve",&gbs::BSCfunction<double>::basisCurve )
                .def("bounds",&gbs::BSCfunction<double>::bounds )
                .def("__call__",&gbs::BSCfunction<double>::operator(),"Function evaluation at givent parameter",py::arg("u"),py::arg("d") = 0)
                // .def("__reduce__", [](gbs::BSCfunction<double> const &self) { // for pickle https://github.com/pybind/pybind11/issues/1261
                //         return py::make_tuple(py::cpp_function([](){return gbs::BSCfunction<double>();}), py::make_tuple());
                // })
                // .def(py::pickle(
                //         [](const gbs::BSCfunction<double> &f) {return py::make_tuple(f.basisCurve());},
                //         [](py::tuple t){
                //                 if (t.size() != 1)
                //                         throw std::runtime_error("Invalid state!");
                //                 gbs::BSCfunction<double> f{ t[0].cast<gbs::BSCurve<double,1>>};
                //                 return f;
                //         }
                // ))
                .def("__copy__",  [](const  gbs::BSCfunction<double> &self) {
                        return  gbs::BSCfunction<double>(self);
                })
                ;

         py::enum_<gbs::KnotsCalcMode>(m, "KnotsCalcMode", py::arithmetic())
        .value("EQUALY_SPACED", gbs::KnotsCalcMode::EQUALY_SPACED)
        .value("CHORD_LENGTH", gbs::KnotsCalcMode::CHORD_LENGTH)
        .value("CENTRIPETAL", gbs::KnotsCalcMode::CENTRIPETAL)
        ;

        m.def(  "interpolate_cn_3d", 
                py::overload_cast<const std::vector< gbs::constrType<double, 3, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<double,3>), 
                "Cn interpolation");
        m.def(  "interpolate_cn_2d", 
                py::overload_cast<const std::vector< gbs::constrType<double, 2, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<double,2>), 
                "Cn interpolation");
        m.def(  "interpolate_cn_1d", 
                py::overload_cast<const std::vector< gbs::constrType<double, 1, 1> > &, 
                size_t , gbs::KnotsCalcMode>(&gbs::interpolate<double,1>), 
                "Cn interpolation");

        m.def( "interpolate_cn_function",
                py::overload_cast<const std::vector<double> &,
                const std::vector<double> &,
                size_t>(&gbs::interpolate<double>),
                py::arg("Q"),py::arg("u"),py::arg("p")
        );

        // py::class_<gbs::extrema_PC_result<double> >(m,"extrema_PC_result")
        // .def_readwrite("d", &gbs::extrema_PC_result<double>::d)
        // .def_readwrite("u", &gbs::extrema_PC_result<double>::u)
        // ;

        // m.def("extrema_curve_point_3d", &extrema_curve_point<double,3>);
        // m.def("extrema_curve_point_2d", &extrema_curve_point<double,2>);
        // m.def("extrema_curve_point_1d", &extrema_curve_point<double,1>);
        m.def("extrema_curve_point", &extrema_curve_point);
        
        m.def( "plot_curves_2d",
        &f_plot_curves_2d);
        m.def( "plot_curves",
        &f_plot_curves);
}

// #include <gbs-render/vtkcurvesrender.h>
// inline void f_plot_curves_2d(const std::vector<std::shared_ptr<gbs::Curve<double,2>>> &crv_lst){gbs::plot(crv_lst);};
// inline void f_plot_curves(const std::vector<std::shared_ptr<gbs::Curve<double,3>>> &crv_lst){gbs::plot(crv_lst);};

// PYBIND11_MODULE(pygbs_render, m) {
//         // auto f_plot_curves_2d = [](const std::vector<gbs::Curve<double,2>> &crv_lst){gbs::plot(crv_lst);};
//         m.def( "plot_curves_2d",
//         &f_plot_curves_2d);
//         m.def( "plot_curves",
//         &f_plot_curves);
// }