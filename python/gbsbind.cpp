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
#include <gbs-render/vtkGbsRender.h>
#include <gbs-render/vtkgridrender.h>
#include <gbs-mesh/tfi.h>


#include "gbsbindextrema.h" // TODO remove

#include <vtk_bind.h>
#include <tuple>
#include <functional>

namespace py = pybind11;

using gbs::operator-;
using gbs::operator+;
using gbs::operator*;
using gbs::operator/;
using gbs::operator^;

void gbs_bind_curves(py::module &m);
void gbs_bind_surfaces(py::module &m);
void gbs_bind_mesh(py::module &m);
void gbs_bind_render(py::module &m);
void gbs_bind_interp(py::module &m);
void gbs_bind_approx(py::module &m);
void gbs_bind_build_curve(py::module &m);
void gbs_bind_curveTools(py::module &m);
void gbs_bind_surfaceTools(py::module &m);
void gbs_bind_shaping(py::module &m);

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

        // const size_t dim1 = 2;
        // const size_t dim2 = 3;
        using T = double;
        
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
        m.def("cross", [](const gbs::point<double, 3> &v1,const gbs::point<double, 3> &v2){return v1^v2;},
                "Cross product points/vector",py::arg("v1"), py::arg("v2"));
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

        gbs_bind_curves(m);

        gbs_bind_surfaces(m);


        m.def("basis_function",
                py::overload_cast<
                double,
                size_t,
                size_t,
                size_t,
                const std::vector<double> &
                >(&gbs::basis_function<double>),
                "BSpline basis function",
                py::arg("u"),
                py::arg("i"),
                py::arg("p"),
                py::arg("d"),
                py::arg("k")
        );

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

        gbs_bind_curveTools(m);

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
                [](const gbs::CurveOffset2D<double, gbs::BSCfunction<double>> & s){return std::make_shared<gbs::CurveOffset2D<double, gbs::BSCfunction<double>>>(s);}
        );
        m.def("make_shared",
                [](const gbs::CurveOffset2D<double, std::function<double(double,size_t)>> & s){return std::make_shared<gbs::CurveOffset2D<double, std::function<double(double,size_t)>>>(s);}
        );
        m.def("build_simple_mult_flat_knots",
                py::overload_cast<double, double, size_t, size_t>(&gbs::build_simple_mult_flat_knots<double>)
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

        gbs_bind_build_curve(m);

        gbs_bind_surfaceTools(m);

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

        gbs_bind_approx(m);

        // Curve analysys
        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,3> &, size_t, size_t>(&gbs::abs_curv<double,3,10>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"),py::arg("n")=30,py::arg("p")=3
        );
        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,2> &, size_t, size_t>(&gbs::abs_curv<double,2,10>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"),py::arg("n")=30,py::arg("p")=3
        );
        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,3> &, double, double, size_t, size_t>(&gbs::abs_curv<double,3,10>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"), py::arg("u1"), py::arg("u2"),py::arg("n")=30,py::arg("p")=3
        );
        m.def("abs_curv",
                py::overload_cast<const gbs::Curve<double,2> &, double, double, size_t, size_t>(&gbs::abs_curv<double,2,10>),
                "Builds a function returning curve's parameter corresponding to the curvilinear abscissa",
                py::arg("crv"), py::arg("u1"), py::arg("u2"),py::arg("n")=30,py::arg("p")=3
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
        const size_t n_int_pt{31};
        const size_t n_int_fast_pt{15};
        m.def("length",
                py::overload_cast<const gbs::Curve<double,3> &, size_t, bool>(&gbs::length<double,3,n_int_pt>),
                "Precise curve length using 30 gauss_kronrods integration points",
                py::arg("crv"),py::arg("d")=0, py::arg("adaptive")=true
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,2> &, size_t, bool>(&gbs::length<double,2,n_int_pt>),
                "Precise curve length using 30 gauss_kronrods integration points",
                py::arg("crv"),py::arg("d")=0, py::arg("adaptive")=true
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,3> &, double, double, size_t, bool>(&gbs::length<double,3,n_int_pt>),
                "Precise curve length using 30 gauss_kronrods integration points between u1 and u2",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0, py::arg("adaptive")=true
        );
        m.def("length",
                py::overload_cast<const gbs::Curve<double,2> &, double, double, size_t, bool>(&gbs::length<double,2,n_int_pt>),
                "Precise curve length using 30 gauss_kronrod integration points between u1 and u2",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0, py::arg("adaptive")=true
        );
        m.def("length",
              [](const gbs::point<double, 3> &p1, const gbs::point<double, 3> &p2){return gbs::norm<double,3>(p2-p1);},
              "Distance between 2 points",py::arg("p1"), py::arg("p2"));
        m.def("length",
              [](const gbs::point<double, 2> &p1, const gbs::point<double, 2> &p2){return gbs::norm<double,2>(p2-p1);},
              "Distance between 2 points",py::arg("p1"), py::arg("p2"));
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,3> &, size_t, bool>(&gbs::length<double,3,n_int_fast_pt>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("d")=0, py::arg("adaptive")=false
        );
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,2> &, size_t, bool>(&gbs::length<double,2,n_int_fast_pt>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("d")=0, py::arg("adaptive")=false
        );
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,3> &, double, double, size_t, bool>(&gbs::length<double,3,n_int_fast_pt>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0, py::arg("adaptive")=false
        );
        m.def("length_fast",
                py::overload_cast<const gbs::Curve<double,2> &, double, double, size_t, bool>(&gbs::length<double,2,n_int_fast_pt>),
                "Precise curve length using 10 gauss integration points",
                py::arg("crv"),py::arg("u1"),py::arg("u2"),py::arg("d")=0, py::arg("adaptive")=false
        );
        ///////////////////
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
        m.def("translated",
              py::overload_cast<const std::array<double, 3> &, const std::array<double, 3> &>(&gbs::translated<double, 3>),
              "Translate geom",
              py::arg("pt"), py::arg("vec")); 
        m.def("rotate",
              py::overload_cast<gbs::point<double, 2> &, double,const gbs::point<double, 2>&>(&gbs::rotate<double,2>),
              "Rotate geom",
              py::arg("pt"), py::arg("angle"), py::arg("center")= gbs::point<double, 2>{0.,0.});
        m.def("rotated",
              py::overload_cast<const gbs::point<double, 2> &, double,const gbs::point<double, 2>&>(&gbs::rotated<double,2>),
              "Rotate geom",
              py::arg("pt"), py::arg("angle"), py::arg("center")= gbs::point<double, 2>{0.,0.});
        m.def("rotated",
              py::overload_cast<const gbs::point<double, 3> &, double,const std::array<double, 3>&>(&gbs::rotated<double>),
              "Rotate geom",
              py::arg("pt"), py::arg("angle"), py::arg("ax")= gbs::point<double, 3>{0.,0.,1.});
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

         m.def("uniform_distrib_params",
               py::overload_cast<const gbs::Curve<double, 3> &, double, double, size_t, size_t>(&gbs::uniform_distrib_params<double, 3>),
               "Uniform Curve discretization",
               py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np"), py::arg("n_law") = 30);
         m.def("uniform_distrib_params",
               py::overload_cast<const gbs::Curve<double, 2> &, double, double, size_t, size_t>(&gbs::uniform_distrib_params<double, 2>),
               "Uniform Curve discretization",
               py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np"), py::arg("n_law") = 30);

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
        // SHAPING
        gbs_bind_shaping(m);
        //////// RENDER
        gbs_bind_render(m);
        //////////// MESH
        gbs_bind_mesh(m);
}

// #include <gbs-render/vtkGbsRender.h>
// inline void f_plot_curves_2d(const std::vector<std::shared_ptr<gbs::Curve<double,2>>> &crv_lst){gbs::plot(crv_lst);};
// inline void f_plot_curves(const std::vector<std::shared_ptr<gbs::Curve<double,3>>> &crv_lst){gbs::plot(crv_lst);};

// PYBIND11_MODULE(pygbs_render, m) {
//         // auto f_plot_curves_2d = [](const std::vector<gbs::Curve<double,2>> &crv_lst){gbs::plot(crv_lst);};
//         m.def( "plot_curves_2d",
//         &f_plot_curves_2d);
//         m.def( "plot_curves",
//         &f_plot_curves);
// }