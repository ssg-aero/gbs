#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/bscinterp.h>
#include <gbs/bscbuild.h>
#include <gbs/bssanalysis.h>
#include <gbs/bscanalysis.h>
#include <gbs/bssbuild.h>
#include <gbs/bsctools.h>
#include <gbs/bsstools.h>
#include <gbs/bscapprox.h>
#include <gbs/extrema.h>
#include <gbs/transform.h>
#include <gbs-render/vtkcurvesrender.h>
#include <gbs-render/vtkgridrender.h>
#include <gbs-mesh/tfi.h>

#include "gbsbindapprox.h"
#include "gbsbindextrema.h"

#include <vtk_bind.h>
#include <tuple>
#include <functional>

namespace py = pybind11;

using gbs::operator-;
using gbs::operator+;
using gbs::operator*;
using gbs::operator/;

void gbs_bind_mesh(py::module &m);
void gbs_bind_render(py::module &m);
void gbs_bind_interp(py::module &m);

// static const std::array<size_t, 3> dims{1, 2, 3};
// using T = double;

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
inline auto declare_bscurve(py::module_ &m)
{
       
        using Class = typename std::conditional<rational, gbs::BSCurveRational<T, dim>, gbs::BSCurve<T, dim>>::type;
        using ClassBase =  gbs::Curve<T,dim>;

        std::string typestr;
        std::string rationalstr;
        if(rational) rationalstr="Rational";
        T x;
        // auto typeid_name = typeid(x).name();

        // if     ( strcmp (typeid_name,"float")  == 0 )  typestr = "_f";
        // else if( strcmp (typeid_name,"double") == 0 ) typestr = "";
        // else                                  typestr = typeid_name ;

        std::string pyclass_name = std::string("BSCurve")+ rationalstr + std::to_string(dim) + "d";// + typestr;

        auto cls = py::class_<Class, std::shared_ptr<Class>, ClassBase>(m, pyclass_name.c_str())
            // py::class_<Class>(m, pyclass_name.c_str())
            .def(py::init<const gbs::points_vector<T, dim + rational> &, const std::vector<T> &, size_t>())
            .def(py::init<const gbs::points_vector<T, dim + rational> &, const std::vector<T> &, const std::vector<size_t> &, size_t>())
            .def(py::init<const Class &>())
        //     .def("value", &Class::value, "Curve evaluation at given parameter", py::arg("u"), py::arg("d") = 0)
            .def("begin", &Class::begin, "Curve evaluation at begin", py::arg("d") = 0)
            .def("end", &Class::end, "Curve evaluation at end", py::arg("d") = 0)
            .def("degree", &Class::degree, "Curve's degree")
            .def("knotsFlats", &Class::knotsFlats, "Curve's flat knots")
            .def("knots", &Class::knots, "Curve's knots")
            .def("mults", &Class::mults, "Curve's knots multiplicities")
            .def("insertKnot", &Class::insertKnot, "Insert knot with the given multiplicity", py::arg("u"), py::arg("m") = 1)
            .def("removeKnot", &Class::removeKnot, "Try to remove m times the given knot", py::arg("u"), py::arg("tol"), py::arg("m") = 1)
            .def("poles", &Class::poles, "Curve's poles")
            .def("pole",py::overload_cast<size_t>(&Class::pole),"Edit specified pole")
            .def("pole",py::overload_cast<size_t>(&Class::pole,py::const_),"Access specified pole")
            .def("copyKnots",&Class::copyKnots,"Replace curve's knots")
            .def("reverse", &Class::reverse, "reverse curve orientation")
            .def("trim", &Class::trim, "Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones",py::arg("u1"),py::arg("u2"),py::arg("permanently")=true)
            .def("changeBounds", py::overload_cast<T, T>(&Class::changeBounds), "Change parametrization to fit between k1 and k2")
            .def("changeBounds", py::overload_cast<const std::array<T, 2> &>(&Class::changeBounds), "Change parametrization to fit between k1 and k2")
        //     .def("bounds", &Class::bounds, "Returns curves's start stop values")
            .def("increaseDegree", &Class::increaseDegree, "Increment curve's degree of 1")
            .def("__copy__", [](const Class &self)
                 { return Class(self); })
        //     .def("__call__", &Class::operator(), "Curve evaluation at given parameter", py::arg("u"), py::arg("d") = 0)
        ;
        // if( rational )
        // {
        //         cls.def(
        //                 py::init<const gbs::points_vector<T, dim> &, const std::vector<T> &, const std::vector<T> &, size_t>()
        //                 , py::arg("poles"), py::arg("knots_flats"), py::arg("weights"), py::arg("degree")
        //         )
        //         .def("polesProjected",&Class::polesProjected)
        //         .def("weights",&Class::weights)
        //         ;

        // }
        return cls;
}

template <typename T, size_t dim, bool rational>
inline void declare_bssurface(py::module_ &m)
{
       
        using Class = typename std::conditional<rational, gbs::BSSurfaceRational<T, dim>, gbs::BSSurface<T, dim>>::type;
        using ClassBase =  gbs::Surface<T,dim>;

        std::string typestr;
        std::string rationalstr;
        if(rational) rationalstr="Rational";
        T x;
        // auto typeid_name = typeid(x).name();

        // if     ( strcmp (typeid_name,"float")  == 0 )  typestr = "_f";
        // else if( strcmp (typeid_name,"double") == 0 ) typestr = "";
        // else                                  typestr = typeid_name ;

        std::string pyclass_name = std::string("BSSurface")+ rationalstr + std::to_string(dim) + "d";// + typestr;

        py::class_<Class,std::shared_ptr<Class>, ClassBase>(m, pyclass_name.c_str())

            .def(py::init<
                 const gbs::points_vector<T, dim + rational> &,
                 const std::vector<T> &,
                 const std::vector<T> &,
                 size_t,
                 size_t>())
        //     .def(py::init<const Class &>())
        //     .def("value", &Class::value, "Surface evaluation at given parameter", py::arg("u"), py::arg("v"), py::arg("du") = 0, py::arg("dv") = 0)
            // .def("begin", &Class::begin,"Curve evaluation at begin",py::arg("d") = 0)
            // .def("end", &Class::end,"Curve evaluation at end",py::arg("d") = 0)
            .def("degreeU", &Class::degreeU,"Surface's U degree")
            .def("degreeV", &Class::degreeV,"Surface's V degree")
            .def("knotsFlatsU", &Class::knotsFlatsU,"Surface's U knots")
            .def("knotsFlatsV", &Class::knotsFlatsV,"Surface's V knots")
            // .def("insertKnot", &Class::insertKnot,"Insert knot with the given multiplicity",py::arg("u"),py::arg("m") = 1)
            // .def("removeKnot", &Class::removeKnot,"Try to remove m times the given knot",py::arg("u"),py::arg("tol"),py::arg("m") = 1)
            .def("poles", &Class::poles,"Surface's poles")
            // .def("reverse", &Class::reverse,"reverse curve orientation")
            // .def("trim", &Class::trim,"Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones")
            // .def("changeBounds",py::overload_cast<T,T>(&Class::changeBounds),"Change parametrization to fit between k1 and k2")
            // .def("changeBounds",py::overload_cast<const std::array<T,2>&>(&Class::changeBounds),"Change parametrization to fit between k1 and k2")
            .def("bounds", &Class::bounds,"Returns surface's start stop values")
            // .def("increaseDegree", &Class::increaseDegree,"Increment curve's degree of 1")
            .def("__copy__", [](const Class &self)
                 { return Class(self); })
            .def("__call__", &Class::operator(), "Surface evaluation at given parameter", py::arg("u"), py::arg("v"), py::arg("du") = 0, py::arg("dv") = 0);
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



// PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(gbs, m) {

        const size_t dim1 = 2;
        const size_t dim2 = 3;
        
        // m.def("add",
        //       py::overload_cast<const gbs::point<double, 3> &, const gbs::point<double, 3> &>(&gbs::operator+<double, 3>),
        //       "Add points/vector",py::arg("v1"), py::arg("v2"));
        // m.def("add",
        //       py::overload_cast<const gbs::point<double, 2> &, const gbs::point<double, 2> &>(&gbs::operator+<double, 2>),
        //       "Add points/vector",py::arg("v1"), py::arg("v2"));
        // m.def("__add__", py::overload_cast<const gbs::point<double, 3> &, const gbs::point<double, 3> &>(&gbs::operator+<double, 3>),
        m.def("add", [](const gbs::point<double, 2> &v1,const gbs::point<double, 2> &v2){return v1+v2;},
                "Add points/vector",py::arg("v1"), py::arg("v2"));
        m.def("add", [](const gbs::point<double, 3> &v1,const gbs::point<double, 3> &v2){return v1+v2;},
                "Add points/vector",py::arg("v1"), py::arg("v2"));
        m.def("sub", [](const gbs::point<double, 2> &v1,const gbs::point<double, 2> &v2){return v1-v2;},
                "Substract points/vector",py::arg("v1"), py::arg("v2"));
        m.def("sub", [](const gbs::point<double, 3> &v1,const gbs::point<double, 3> &v2){return v1-v2;},
                "Substract points/vector",py::arg("v1"), py::arg("v2"));
        m.def("mult", [](double s, const gbs::point<double, 3> &v){return s*v;},
                "Add points/vector",py::arg("s"), py::arg("v"));
        m.def("mult", [](double s, const gbs::point<double, 2> &v){return s*v;},
                "Add points/vector",py::arg("s"), py::arg("v"));
        m.def("neg",[](const gbs::point<double, 2> &v){return -1.*v;});
        m.def("neg",[](const gbs::point<double, 3> &v){return -1.*v;});
        m.def("adim",[](const gbs::point<double, 2> &v){return v/gbs::norm(v);});
        m.def("adim",[](const gbs::point<double, 3> &v){return v/gbs::norm(v);});
        m.def("norm",[](const gbs::point<double, 2> &v){return gbs::norm(v);});
        m.def("norm",[](const gbs::point<double, 3> &v){return gbs::norm(v);});
        m.def("dist",[](const gbs::point<double, 3> &a, const gbs::point<double, 3> &b){return gbs::distance(a,b);});
        m.def("dist",[](const gbs::point<double, 2> &a, const gbs::point<double, 2> &b){return gbs::distance(a,b);});
        // m.def("__mult__", py::overload_cast<double, const gbs::point<double, 3> &>(&gbs::operator*<double, double, 3>),
        //         "Multiply points/vector",py::arg("v1"), py::arg("v2"));

        // py::class_<gbs::Curve<double,3> >(m, "Curve3d");
        py::class_<gbs::Curve<double,3>, std::shared_ptr<gbs::Curve<double,3>> >(m, "Curve3d")
        .def("value", &gbs::Curve<double,3>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("bounds", &gbs::Curve<double,3>::bounds, "Returns curves's start stop values")
        .def("midPosition", &gbs::Curve<double,3>::midPosition, "Returns curves's start stop average value")
        .def("begin", &gbs::Curve<double,3>::begin, "Curve evaluation at begin", py::arg("d") = 0)
        .def("end", &gbs::Curve<double,3>::end, "Curve evaluation at end", py::arg("d") = 0)
        .def("__call__",&gbs::Curve<double,3>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        ;
        py::class_<gbs::Curve<double,2>, std::shared_ptr<gbs::Curve<double,2>> >(m, "Curve2d")
        .def("value", &gbs::Curve<double,2>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("bounds", &gbs::Curve<double,2>::bounds, "Returns curves's start stop values")
        .def("midPosition", &gbs::Curve<double,2>::midPosition, "Returns curves's start stop average value")
        .def("begin", &gbs::Curve<double,2>::begin, "Curve evaluation at begin", py::arg("d") = 0)
        .def("end", &gbs::Curve<double,2>::end, "Curve evaluation at end", py::arg("d") = 0)
        .def("__call__",&gbs::Curve<double,2>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        ;
        py::class_<gbs::Curve<double,1>, std::shared_ptr<gbs::Curve<double,1>> >(m, "Curve1d")
        .def("value", &gbs::Curve<double,1>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("bounds", &gbs::Curve<double,1>::bounds, "Returns curves's start stop values")
        .def("midPosition", &gbs::Curve<double,1>::midPosition, "Returns curves's start stop average value")
        .def("begin", &gbs::Curve<double,1>::begin, "Curve evaluation at begin", py::arg("d") = 0)
        .def("end", &gbs::Curve<double,1>::end, "Curve evaluation at end", py::arg("d") = 0)
        .def("__call__",&gbs::Curve<double,1>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        ;

        py::class_<gbs::Surface<double,3>, std::shared_ptr<gbs::Surface<double,3>> >(m, "Surface3d")
        .def("value", &gbs::Surface<double,3>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        .def("bounds", &gbs::Surface<double,3>::bounds, "Returns surface's start stop values")
        .def("__call__",&gbs::Surface<double,3>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        ;
        py::class_<gbs::Surface<double,2>, std::shared_ptr<gbs::Surface<double,2>> >(m, "Surface2d")
        .def("value", &gbs::Surface<double,2>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        .def("bounds", &gbs::Surface<double,2>::bounds, "Returns surface's start stop values")
        .def("__call__",&gbs::Surface<double,2>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        ;
        py::class_<gbs::Surface<double,1>, std::shared_ptr<gbs::Surface<double,1>> >(m, "Surface1d")
        .def("value", &gbs::Surface<double,1>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        .def("bounds", &gbs::Surface<double,1>::bounds, "Returns surface's start stop values")
        .def("__call__",&gbs::Surface<double,1>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
        ;
        
        py::class_<gbs::Line<double,2>,  std::shared_ptr<gbs::Line<double,2>> >(m, "Line2d")
        .def(py::init<const gbs::point<double, 2> &, const gbs::point<double, 2>&>())
        .def(py::init<const gbs::ax1<double, 2> &>())
        .def("value", &gbs::Curve<double,2>::value,"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        .def("__call__",&gbs::Curve<double,2>::operator(),"Curve evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
        ;

        py::class_<gbs::CurveOnSurface<double,3>, std::shared_ptr<gbs::CurveOnSurface<double,3>>, gbs::Curve<double,3>>(m, "CurveOnSurface3d")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const std::shared_ptr<gbs::Surface<double, 3>>&>())
        // .def(py::init<const gbs::BSCurve<double, 2> &, const gbs::BSSurface<double, 3>&>())
        // .def(py::init<const gbs::BSCurve<double, 2> &, const gbs::SurfaceOfRevolution<double>&>())
        // .def(py::init<const gbs::CurveOffset<double, 2,gbs::BSCfunction<double>> &, const gbs::SurfaceOfRevolution<double>&>())
        .def(py::init<const gbs::CurveOnSurface<double,3> &>())
        ;
        py::class_<gbs::CurveOnSurface<double,2>, std::shared_ptr<gbs::CurveOnSurface<double,2>>, gbs::Curve<double,2>>(m, "CurveOnSurface2d")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const std::shared_ptr<gbs::Surface<double, 2>>&>())
        // .def(py::init<const gbs::BSCurve<double, 2> &, const gbs::BSSurface<double, 2>&>())
        .def(py::init<const gbs::CurveOnSurface<double,2> &>())
        ;
        py::class_<gbs::CurveOffset<double,2,gbs::BSCfunction<double>>,std::shared_ptr<gbs::CurveOffset<double,2,gbs::BSCfunction<double>>>, gbs::Curve<double,2> >(m, "CurveOffset2d_bs")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const gbs::BSCfunction<double> &>())
        // .def(py::init<const gbs::BSCurve<double, 2> &, const gbs::BSCfunction<double> &>())
        // .def(py::init<const gbs::CurveOffset<double,2,gbs::BSCfunction<double>> &>())
        .def( "basisCurve", &gbs::CurveOffset< double, 2, gbs::BSCfunction<double> >::basisCurve )
        .def( "offset", &gbs::CurveOffset< double, 2, gbs::BSCfunction<double> >::offset )
        ;
        py::class_<gbs::CurveOffset<double,2,std::function<double(double,size_t)>>, std::shared_ptr<gbs::CurveOffset<double,2,std::function<double(double,size_t)>>>, gbs::Curve<double,2> >(m, "CurveOffset2d_func")
        .def(py::init<const std::shared_ptr<gbs::Curve<double, 2>> &, const std::function<double(double,size_t)> &>())
        // .def(py::init<const gbs::BSCurve<double, 2> &, const std::function<double(double,size_t)> &>())
        // .def(py::init<const gbs::CurveOffset<double,2,std::function<double(double,size_t)>> &>())
        .def( "basisCurve", &gbs::CurveOffset< double, 2, std::function<double(double,size_t)> >::basisCurve )
        .def( "offset", &gbs::CurveOffset< double, 2, std::function<double(double,size_t)> >::offset )
        ;
        // py::class_<gbs::CurveOffset<double,2,const py::object &>, gbs::Curve<double,2> >(m, "CurveOffset2d_func")
        // .def(py::init<const gbs::BSCurve<double, 2> &, const  py::object & >())
        // // .def(py::init<const gbs::CurveOffset<double,2,const py::object &>())
        // ;

        gbs::ax2<double,3> ax2_z {
                gbs::point<double,3>{0., 0., 0.}, 
                gbs::point<double,3>{0., 0., 1.}, 
                gbs::point<double,3>{1., 0., 0.}
        };
        py::class_<gbs::SurfaceOfRevolution<double>, std::shared_ptr<gbs::SurfaceOfRevolution<double>>, gbs::Surface<double, 3>>(m, "SurfaceOfRevolution")
        .def(py::init< std::shared_ptr<gbs::Curve<double, 2>> &, const gbs::ax2<double, 3>, double, double>(),
                py::arg("crv"), py::arg("ax") = ax2_z, py::arg("a1") = 0., py::arg("a2") = 2 * std::numbers::pi)
        // .def(py::init<gbs::BSCurve<double, 2> &, const gbs::ax2<double, 3>, double, double>(),
        //         py::arg("crv"), py::arg("ax") = ax2_z, py::arg("a1") = 0., py::arg("a2") = std::numbers::pi)
        // .def(py::init<gbs::CurveOnSurface<double, 2> &, const gbs::ax2<double, 3>, double, double>(),
        //         py::arg("crv"), py::arg("ax") = ax2_z, py::arg("a1") = 0., py::arg("a2") = std::numbers::pi)
        ;

        declare_bscurve<double,3,false>(m);
        declare_bscurve<double,2,false>(m);
        declare_bscurve<double,1,false>(m);

        declare_bscurve<double,3,true>(m).def(
                py::init<const gbs::points_vector<double,3> &, const std::vector<double> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots_flats"), py::arg("weights"), py::arg("degree")
        )
        .def(
                py::init<const gbs::points_vector<double,3> &, const std::vector<double> &, const std::vector<size_t> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots"), py::arg("mults"), py::arg("weights"), py::arg("degree")
        )
        .def("polesProjected",&gbs::BSCurveRational<double,3>::polesProjected)
        .def("weights",&gbs::BSCurveRational<double,3>::weights)
        ;
        declare_bscurve<double,2,true>(m).def(
                py::init<const gbs::points_vector<double,2> &, const std::vector<double> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots_flats"), py::arg("weights"), py::arg("degree")
        )
        .def(
                py::init<const gbs::points_vector<double,2> &, const std::vector<double> &, const std::vector<size_t> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots"), py::arg("mults"), py::arg("weights"), py::arg("degree")
        )
        .def("polesProjected",&gbs::BSCurveRational<double,2>::polesProjected)
        .def("weights",&gbs::BSCurveRational<double,2>::weights)
        ;
        declare_bscurve<double,1,true>(m).def(
                py::init<const gbs::points_vector<double,1> &, const std::vector<double> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots_flats"), py::arg("weights"), py::arg("degree")
        )
        .def(
                py::init<const gbs::points_vector<double,1> &, const std::vector<double> &, const std::vector<size_t> &, const std::vector<double> &, size_t>()
                , py::arg("poles"), py::arg("knots"), py::arg("mults"), py::arg("weights"), py::arg("degree")
        )
        .def("polesProjected",&gbs::BSCurveRational<double,1>::polesProjected)
        .def("weights",&gbs::BSCurveRational<double,1>::weights)
        ;

        declare_bssurface<double,3,false>(m);
        declare_bssurface<double,2,false>(m);
        declare_bssurface<double,1,false>(m);

        declare_bssurface<double,3,true>(m);
        declare_bssurface<double,2,true>(m);
        declare_bssurface<double,1,true>(m);

        py::class_<gbs::BSCfunction<double>>(m,"BSCfunction")
                .def(py::init<const gbs::BSCfunction<double> &>())
                .def(py::init<const std::vector<std::array<double, 1>> &, const std::vector<double> &, size_t >(),
                        py::arg("poles"),py::arg("knots_flats"),py::arg("deg"))
                .def(py::init<const std::vector<double> &, const std::vector<double> &, size_t >(),
                        py::arg("poles"),py::arg("knots_flats"),py::arg("deg"))
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<size_t> &, size_t >(),
                        py::arg("poles"), py::arg("knots"), py::arg("mults"), py::arg("deg"))
                .def("value",&gbs::BSCfunction<double>::value,"Function evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
                .def("basisCurve",&gbs::BSCfunction<double>::basisCurve )
                .def("bounds",&gbs::BSCfunction<double>::bounds )
                .def("__call__",&gbs::BSCfunction<double>::operator(),"Function evaluation at given parameter",py::arg("u"),py::arg("d") = 0)
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
        py::class_<gbs::BSSfunction<double>>(m,"BSSfunction")
                .def(py::init<const gbs::BSSfunction<double> &>())
                .def(py::init<const gbs::BSSurface<double, 1> &>())
                .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &, size_t, size_t>())
                .def("basisSurface",&gbs::BSSfunction<double>::basisSurface)
                .def("bounds",&gbs::BSSfunction<double>::bounds )
                .def("__call__",&gbs::BSSfunction<double>::operator(),"Function evaluation at given parameter",py::arg("u"),py::arg("v"),py::arg("du") = 0,py::arg("dv") = 0)
                .def("__copy__",  [](const  gbs::BSSfunction<double> &self) {
                        return  gbs::BSSfunction<double>(self);
                })
                ;

        m.def("eval",
                py::overload_cast<
                        double,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 3> &,
                        size_t,
                        size_t,
                        bool
                >(&gbs::eval_value_decasteljau<double,3>),
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );
        m.def("eval",
                py::overload_cast<
                        double,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 2> &,
                        size_t,
                        size_t,
                        bool
                >(&gbs::eval_value_decasteljau<double,2>),
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );
        m.def("eval",
                py::overload_cast<
                        double,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 1> &,
                        size_t,
                        size_t,
                        bool
                >(&gbs::eval_value_decasteljau<double,1>),
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );

        // for( const auto dim : dims)
        // {
        //         m.def("eval",
        //                 [] ( double u, const std::vector<T> &k, const std::vector<std::array<T, dims[1]>> &poles, const std::vector<T> &weights, size_t p, size_t d = 0,bool use_span_reduction =true )
        //                 {
        //                         return gbs::eval_rational_value_simple(u, k, gbs::add_weights_coord(poles, weights), p , d, use_span_reduction);
        //                 }
        //         ,
        //                 "Non rational BSplineCurve evaluation",
        //                 py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("weights"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        //         );
        // }

        m.def("eval",
                [] ( double u, const std::vector<double> &k, const std::vector<std::array<double, 3>> &poles, const std::vector<double> &weights, size_t p, size_t d = 0,bool use_span_reduction =true )
                {
                        return gbs::eval_rational_value_simple<double,3>(u, k, gbs::add_weights_coord(poles, weights), p , d, use_span_reduction);
                }
        ,
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("weights"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );
        m.def("eval",
                [] ( double u, const std::vector<double> &k, const std::vector<std::array<double, 2>> &poles, const std::vector<double> &weights, size_t p, size_t d = 0,bool use_span_reduction =true )
                {
                        return gbs::eval_rational_value_simple<double,2>(u, k, gbs::add_weights_coord(poles, weights), p , d, use_span_reduction);
                }
                ,
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("weights"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );
        m.def("eval",
                [] ( double u, const std::vector<double> &k, const std::vector<std::array<double, 1>> &poles, const std::vector<double> &weights, size_t p, size_t d = 0,bool use_span_reduction =true )
                {
                        return gbs::eval_rational_value_simple<double,1>(u, k, gbs::add_weights_coord(poles, weights), p , d, use_span_reduction);
                }
                ,
                "Non rational BSplineCurve evaluation",
                py::arg("u"), py::arg("knots_flat"), py::arg("poles"), py::arg("weights"), py::arg("degree"), py::arg("derivative") = 0, py::arg("use_span_reduction") = true
        );

        m.def("eval",
                py::overload_cast<
                        double,
                        double,
                        const std::vector<double> &,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 3> &,
                        size_t,
                        size_t,
                        size_t,
                        size_t
                >(&gbs::eval_value_decasteljau<double,3>),
                "Non rational BSplineSurface evaluation",
                py::arg("u"), py::arg("v"), py::arg("ku"), py::arg("kv"), py::arg("poles"), py::arg("degreeU"), py::arg("degreeV"), py::arg("du") = 0, py::arg("dv") = 0
        );
        m.def("eval",
                py::overload_cast<
                        double,
                        double,
                        const std::vector<double> &,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 2> &,
                        size_t,
                        size_t,
                        size_t,
                        size_t
                >(&gbs::eval_value_decasteljau<double,2>),
                "Non rational BSplineSurface evaluation",
                py::arg("u"), py::arg("v"), py::arg("ku"), py::arg("kv"), py::arg("poles"), py::arg("degreeU"), py::arg("degreeV"), py::arg("du") = 0, py::arg("dv") = 0
        );
        m.def("eval",
                py::overload_cast<
                        double,
                        double,
                        const std::vector<double> &,
                        const std::vector<double> &,
                        const gbs::points_vector<double, 1> &,
                        size_t,
                        size_t,
                        size_t,
                        size_t
                >(&gbs::eval_value_decasteljau<double,1>),
                "Non rational BSplineSurface evaluation",
                py::arg("u"), py::arg("v"), py::arg("ku"), py::arg("kv"), py::arg("poles"), py::arg("degreeU"), py::arg("degreeV"), py::arg("du") = 0, py::arg("dv") = 0
        );

        m.def("unflat_knots",
                py::overload_cast<const std::vector<double> &>(&gbs::unflat_knots<double>),
                "Get unique knots and theire multiplicities",
                py::arg("knots_flat")
        );

        m.def(
                "join",
                [](const gbs::BSCurve<double,3> &c1, gbs::BSCurve<double,3> &c2){return gbs::join<double, 3, false, false>(c1,c2);},
                "Join 2 curves into one",
                py::arg("crv1"), py::arg("crv2")
        );
        m.def(
                "join",
                [](const gbs::BSCurve<double,2> &c1, gbs::BSCurve<double,2> &c2){return gbs::join<double, 2, false, false>(c1,c2);},
                "Join 2 curves into one",
                py::arg("crv1"), py::arg("crv2")
        );
        m.def(
                "join",
                [](const gbs::BSCurve<double,1> &c1, gbs::BSCurve<double,1> &c2){return gbs::join<double, 1, false, false>(c1,c2);},
                "Join 2 curves into one",
                py::arg("crv1"), py::arg("crv2")
        );

        m.def("make_shared",
                // &std::make_shared<gbs::SurfaceOfRevolution<double>>
                [](const gbs::SurfaceOfRevolution<double> & s){return std::make_shared<gbs::SurfaceOfRevolution<double>>(s);}
        );
        // py::class_<gbs::BSCurve<double,2>, std::shared_ptr<gbs::BSCurve<double,2>>>(m, "shr_BSCurve2d")
        // py::class_<std::shared_ptr<gbs::BSCurve<double,2>>>(m, "shr_BSCurve2d");
        // .def(py::init<const gbs::BSCurve<double,2> &>());
        m.def("make_shared",
                [](const gbs::BSCurve<double,2> & s){return std::make_shared<gbs::BSCurve<double,2>>(s);}
        );
        m.def("make_shared",
                [](const gbs::BSCurve<double,3> & s){return std::make_shared<gbs::BSCurve<double,3>>(s);}
        );
        m.def("make_shared",
                [](const gbs::CurveOffset<double,2,gbs::BSCfunction<double>> & s){return std::make_shared<gbs::CurveOffset<double,2,gbs::BSCfunction<double>>>(s);}
        );
        m.def("make_shared",
                [](const gbs::CurveOffset<double,2,std::function<double(double,size_t)>> & s){return std::make_shared<gbs::CurveOffset<double,2,std::function<double(double,size_t)>>>(s);}
        );
        

        py::enum_<gbs::KnotsCalcMode>(m, "KnotsCalcMode", py::arithmetic())
            .value("EQUALY_SPACED", gbs::KnotsCalcMode::EQUALY_SPACED)
            .value("CHORD_LENGTH", gbs::KnotsCalcMode::CHORD_LENGTH)
            .value("CENTRIPETAL", gbs::KnotsCalcMode::CENTRIPETAL);

        m.def("curve_parametrization",
                py::overload_cast<
                const std::vector<std::array<double,3> > &,
                gbs::KnotsCalcMode,
                bool>(&gbs::curve_parametrization<double,3>),
                "Builds curve's parametrization from passing points, the result cal be set to range from 0. to 1.",
                py::arg("pts"), py::arg("mode"), py::arg("adimensional") = false
        );
        m.def("curve_parametrization",
                py::overload_cast<
                const std::vector<std::array<double,2> > &,
                gbs::KnotsCalcMode,
                bool>(&gbs::curve_parametrization<double,2>),
                "Builds curve's parametrization from passing points, the result cal be set to range from 0. to 1.",
                py::arg("pts"), py::arg("mode"), py::arg("adimensional") = false
        );
        m.def("curve_parametrization",
                py::overload_cast<
                const std::vector<std::array<double,1> > &,
                gbs::KnotsCalcMode,
                bool>(&gbs::curve_parametrization<double,1>),
                "Builds curve's parametrization from passing points, the result cal be set to range from 0. to 1.",
                py::arg("pts"), py::arg("mode"), py::arg("adimensional") = false
        );

        m.def("build_segment", py::overload_cast<const gbs::point<double,3> &, const gbs::point<double,3> &, bool>(&gbs::build_segment<double,3>),
                py::arg("p1"),py::arg("p2"),py::arg("normalized_param") = false
        );
        m.def("build_segment", py::overload_cast<const gbs::point<double,2> &, const gbs::point<double,2> &, bool>(&gbs::build_segment<double,2>),
                py::arg("p1"),py::arg("p2"),py::arg("normalized_param") = false
        );

        gbs_bind_interp(m);

        m.def("to_bscurve_3d",
                [](const gbs::BSCurve<double, 2> &crv,double z){return gbs::add_dimension(crv,z);},
                "Convert 2d curve to 3d curve",
                py::arg("crv"), py::arg("z") = 0.
        );
        m.def("to_bssurface_3d",
                [](const gbs::BSSurface<double, 2> &crv,double z){return gbs::add_dimension(crv,z);},
                "Convert 2d curve to 3d curve",
                py::arg("crv"), py::arg("z") = 0.
        );
        // m.def("loft",
        //         py::overload_cast<const std::list<gbs::BSCurve<double, 2>> &, size_t>(&gbs::loft<double,2>),
        //         py::arg("bs_lst"), py::arg("v_degree_max") = 3
        // );
        // m.def("loft",
        //         py::overload_cast<const std::list<gbs::BSCurve<double, 3>> &, size_t>(&gbs::loft<double,3>),
        //         py::arg("bs_lst"), py::arg("v_degree_max") = 3
        // );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 1>>> &, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,1>),
                py::arg("bs_lst"), py::arg("v_degree_max") = 3, py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,2>),
                py::arg("bs_lst"), py::arg("v_degree_max") = 3, py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,3>),
                py::arg("bs_lst"), py::arg("v_degree_max") = 3, py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 1>>> &, 
                        const std::vector<double>&, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,1>),
                py::arg("bs_lst"), py::arg("v"), py::arg("q"), py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &, 
                        const std::vector<double>&, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,2>),
                py::arg("bs_lst"), py::arg("v"), py::arg("q"), py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        m.def("loft",
                py::overload_cast<
                        const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &, 
                        const std::vector<double>&, 
                        size_t, 
                        double, 
                        size_t, 
                        size_t
                >(&gbs::loft<double,3>),
                py::arg("bs_lst"), py::arg("v"), py::arg("q"), py::arg("dev")=0.01, py::arg("np")=100, py::arg("deg_approx")=5
        );
        // m.def("loft",
        //         py::overload_cast<const std::list<gbs::BSCurve<double, 3>> &, const gbs::BSCurve<double, 3> &, size_t>(&gbs::loft<double,3>),
        //         py::arg("bs_lst"), py::arg("spine"), py::arg("v_degree_max") = 3
        // );
        // APPROX
        gbs_bind_approx<double,1>(m);
        gbs_bind_approx<double,2>(m);
        gbs_bind_approx<double,3>(m);

        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,3> &, size_t>(&gbs::abs_curv<double,3,100>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"),py::arg("n")=30
        );
        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,2> &, size_t>(&gbs::abs_curv<double,2,100>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"),py::arg("n")=30
        );
        // m.def("len_curv",
        //         py::overload_cast<const gbs::Curve<double,3> &, size_t>(&gbs::length<double,3,250>),
        //         "Precise curve length using 250 gauss integration points",
        //         py::arg("crv"),py::arg("d")=0
        // );
        // m.def("len_curv",
        //         py::overload_cast<const gbs::Curve<double,2> &, size_t>(&gbs::length<double,2,250>),
        //         "Precise curve length using 250 gauss integration points",
        //         py::arg("crv"),py::arg("d")=0
        // );
        // m.def("len_curv_fast",
        //         py::overload_cast<const gbs::Curve<double,3> &, size_t>(&gbs::length<double,3,10>),
        //         "Precise curve length using 10 gauss integration points",
        //         py::arg("crv"),py::arg("d")=0
        // );
        // m.def("len_curv_fast",
        //         py::overload_cast<const gbs::Curve<double,2> &, size_t>(&gbs::length<double,2,10>),
        //         "Precise curve length using 10 gauss integration points",
        //         py::arg("crv"),py::arg("d")=0
        // );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,3> &, size_t>(&gbs::length<double,3,250>),
                "Precise curve length using 250 gauss integration points",
                py::arg("crv"),py::arg("d")=0
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,2> &, size_t>(&gbs::length<double,2,250>),
                "Precise curve length using 250 gauss integration points",
                py::arg("crv"),py::arg("d")=0
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,3> &, double, double, size_t>(&gbs::length<double,3,250>),
                "Precise curve length using 250 gauss integration points between u1 and u2",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,2> &, double, double, size_t>(&gbs::length<double,2,250>),
                "Precise curve length using 250 gauss integration points between u1 and u2",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0
        );
        m.def("length",
              [](const gbs::point<double, 3> &p1, const gbs::point<double, 3> &p2){return gbs::norm<double,3>(p2-p1);},
              "Distance between 2 points",py::arg("p1"), py::arg("p2"));
        m.def("length",
              [](const gbs::point<double, 2> &p1, const gbs::point<double, 2> &p2){return gbs::norm<double,2>(p2-p1);},
              "Distance between 2 points",py::arg("p1"), py::arg("p2"));
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,3> &, size_t>(&gbs::length<double,3,10>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("d")=0
        );
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,2> &, size_t>(&gbs::length<double,2,10>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("d")=0
        );
        m.def("normal_direction",
                py::overload_cast<const gbs::Curve<double,3> &, double>(&gbs::normal_direction<double>),
                "Curve normal direction at given parameter",
                py::arg("crv"),py::arg("u")
        );
        m.def("normal_direction",
                py::overload_cast<const gbs::Curve<double,2> &, double>(&gbs::normal_direction<double>),
                "Curve normal direction at given parameter",
                py::arg("crv"),py::arg("u")
        );
        m.def("tangential_direction",
                py::overload_cast<const gbs::Curve<double,3> &, double>(&gbs::tangential_direction<double,3>),
                "Curve tangential direction at given parameter",
                py::arg("crv"),py::arg("u")
        );
        m.def("tangential_direction",
                py::overload_cast<const gbs::Curve<double,2> &, double>(&gbs::tangential_direction<double,2>),
                "Curve tangential direction at given parameter",
                py::arg("crv"),py::arg("u")
        );
        m.def("translate",
              py::overload_cast<gbs::BSCurve<double, 2> &, const std::array<double, 2> &>(&gbs::translate<double, 2>),
              "Translate geom",
              py::arg("crv"), py::arg("vec"));
        m.def("translate",
              py::overload_cast<gbs::BSCurve<double, 3> &, const std::array<double, 3> &>(&gbs::translate<double, 3>),
              "Translate geom",
              py::arg("crv"), py::arg("vec"));
        m.def("rotate",
              py::overload_cast<gbs::point<double, 2> &, double,const gbs::point<double, 2>&>(&gbs::rotate<double,2>),
              "Rotate geom",
              py::arg("pt"), py::arg("angle"), py::arg("center")= gbs::point<double, 2>{0.,0.});
        m.def("rotated",
              py::overload_cast<const gbs::point<double, 2> &, double,const gbs::point<double, 2>&>(&gbs::rotated<double,2>),
              "Rotate geom",
              py::arg("pt"), py::arg("angle"), py::arg("center")= gbs::point<double, 2>{0.,0.});
        m.def("rotate",
              py::overload_cast<gbs::BSCurve<double, 2> &, double>(&gbs::rotate<gbs::BSCurve<double, 2>, double>),
              "Rotate geom",
              py::arg("crv"), py::arg("angle"));
        m.def("rotate",
              py::overload_cast<gbs::BSCurve<double, 3> &, double, const std::array<double, 3> &>(&gbs::rotate<gbs::BSCurve<double, 3>, double>),
              "Rotate geom",
              py::arg("crv"), py::arg("angle"), py::arg("axis"));
        m.def("rotate",
              py::overload_cast<gbs::BSSurface<double, 3> &, double, const std::array<double, 3> &>(&gbs::rotate<gbs::BSSurface<double, 3>, double>),
              "Rotate geom",
              py::arg("srf"), py::arg("angle"), py::arg("axis"));
        m.def("scaled",
                py::overload_cast<const gbs::point<double,3> &, double, gbs::point<double,3> >(&gbs::scaled<double,3>),
                py::arg("pt"), py::arg("s"), py::arg("center")=gbs::point<double,3>{});
        m.def("scaled",
                py::overload_cast<const gbs::point<double,2> &, double, gbs::point<double,2> >(&gbs::scaled<double,2>),
                py::arg("pt"), py::arg("s"), py::arg("center")=gbs::point<double,2>{});
        m.def("scale",
                py::overload_cast<gbs::point<double,3> &, double, gbs::point<double,3> >(&gbs::scale<double,3>),
                py::arg("pt"), py::arg("s"), py::arg("center")=gbs::point<double,3>{});
        m.def("scale",
                py::overload_cast<gbs::point<double,2> &, double, gbs::point<double,2> >(&gbs::scale<double,2>),
                py::arg("pt"), py::arg("s"), py::arg("center")=gbs::point<double,2>{});
        // m.def("to_bscurverational_3d",
        //         [](const gbs::BSCurveRational<double, 2> &crv,double z){return gbs::add_dimension(crv,z);},
        //         "Convert 2d curve to 3d curve",
        //         py::arg("crv"), py::arg("z") = 0.
        // );
         // py::class_<gbs::extrema_PC_result<double> >(m,"extrema_PC_result")
         // .def_readwrite("d", &gbs::extrema_PC_result<double>::d)
         // .def_readwrite("u", &gbs::extrema_PC_result<double>::u)
         // ;

        ///// EXTREMA
        py::enum_<nlopt::algorithm>(m, "nlopt_algorithm", py::arithmetic())
            .value("LN_PRAXIS", nlopt::algorithm::LN_PRAXIS)
            .value("LN_COBYLA", nlopt::algorithm::LN_COBYLA);
        gbs_bind_extrema<double,1>(m);
        gbs_bind_extrema<double,2>(m);
        gbs_bind_extrema<double,3>(m);
        m.def("extrema_curve_curve",
                py::overload_cast<const gbs::Line<double, 2>&, const gbs::Line<double, 2>&>(&gbs::extrema_curve_curve<double>),
                py::arg("crv1"),  py::arg("crv2")
        );


         // m.def("discretize_curve",&f_discretize_curve);
        m.def("discretize_curve_unif",
               py::overload_cast<const gbs::Curve<double, 3> &, size_t>(&gbs::discretize<double, 3>),
               "Uniformly spaced discretization",
               py::arg("crv"), py::arg("np"));
        m.def("discretize_curve_unif",
               py::overload_cast<const gbs::Curve<double, 2> &, size_t>(&gbs::discretize<double, 2>),
               "Uniformly spaced discretization",
               py::arg("crv"), py::arg("np"));
         m.def("discretize_curve",
               py::overload_cast<const gbs::Curve<double, 3> &, size_t, double, size_t>(&gbs::discretize<double, 3>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
         m.def("discretize_curve",
               py::overload_cast<const gbs::Curve<double, 2> &, size_t, double, size_t>(&gbs::discretize<double, 2>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
         m.def("deviation_based_params",
               py::overload_cast<const gbs::Curve<double, 3> &, size_t, double, size_t>(&gbs::deviation_based_params<double, 3>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
         m.def("deviation_based_params",
               py::overload_cast<const gbs::Curve<double, 2> &, size_t, double, size_t>(&gbs::deviation_based_params<double, 2>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
        m.def("deviation_based_params",
               py::overload_cast<const gbs::Curve<double, 3> &, double, double, size_t, double, size_t>(&gbs::deviation_based_params<double, 3>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
         m.def("deviation_based_params",
               py::overload_cast<const gbs::Curve<double, 2> &, double, double, size_t, double, size_t>(&gbs::deviation_based_params<double, 2>),
               "Curve discretization based on deviation",
               py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);
         // m.def("discretize_surface",
         //       py::overload_cast<const gbs::Surface<double, 3> &, size_t, size_t>(&gbs::discretize<double, 3>),
         //       " ",
         //       py::arg("srf"), py::arg("nu"), py::arg("nv"));
        m.def("normal_direction",
                py::overload_cast<const gbs::Curve<double,2>&, double>(&gbs::normal_direction<double>),"Compute normalized normal direction of the curve using curvature, thus is the later is null result is +/-infinity",
                py::arg("crv"), py::arg("u") );
        m.def("normal_direction",
                py::overload_cast<const gbs::Curve<double,3>&, double>(&gbs::normal_direction<double>),"Compute normalized normal direction of the curve using curvature, thus is the later is null result is +/-infinity",
                py::arg("crv"), py::arg("u") );
        m.def("tangential_direction",
                py::overload_cast<const gbs::Curve<double,2>&, double>(&gbs::normal_direction<double>),"Compute normalized tangential direction of the curve",
                py::arg("crv"), py::arg("u") );
        m.def("tangential_direction",
                py::overload_cast<const gbs::Curve<double,3>&, double>(&gbs::normal_direction<double>),"Compute normalized tangential direction of the curve",
                py::arg("crv"), py::arg("u") );
        //////// RENDER
        gbs_bind_render(m);
        //////////// MESH
        gbs_bind_mesh(m);
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