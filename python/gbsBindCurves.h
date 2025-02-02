#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "repr.h"
#include <gbs/curves>
#include <gbs-io/tojson.h>
using namespace gbs;

template <size_t dim>
inline auto add_ext(const char *func_name)
{
    return std::string( func_name ) + std::to_string( dim) + std::string("d");
}

template <typename T, size_t dim, bool rational>
inline auto declare_bscurve(py::module_ &m)
{
    using Class = typename std::conditional<rational, gbs::BSCurveRational<T, dim>, gbs::BSCurve<T, dim>>::type;
    using ClassBase =  gbs::Curve<T,dim>;

    std::string typestr;
    std::string rationalstr;
    if(rational) rationalstr="Rational";

    std::string pyclass_name = std::string("BSCurve")+ rationalstr + std::to_string(dim) + "d";// + typestr;
    std::string pyclass_base_name = std::string("BSCurveBase")+ rationalstr + std::to_string(dim) + "d";// + typestr;

    py::class_<BSCurveGeneral<T, dim, rational>, std::shared_ptr<BSCurveGeneral<T, dim, rational>>, ClassBase>(m, pyclass_base_name.c_str());

    auto cls = py::class_<Class, std::shared_ptr<Class>, BSCurveGeneral<T, dim, rational> >(m, pyclass_name.c_str())
        .def(py::init<const gbs::points_vector<T, dim + rational> &, const std::vector<T> &, size_t>(),
            py::arg("poles"), 
            py::arg("knots_flats"), 
            py::arg("degree")
        )
        .def(py::init<const gbs::points_vector<T, dim + rational> &, const std::vector<T> &, const std::vector<size_t> &, size_t>(),
            py::arg("poles"), py::arg("knots"), py::arg("mults"), py::arg("degree")
        )
        .def(py::init<const Class &>(),py::arg("crv"))
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
        .def("copyPoles",&Class::copyPoles,"Replace curve's poles")
        .def("reverse", &Class::reverse, "Reverse curve orientation")
        .def("reversed", &Class::reversed, "Return curve with reversed orientation")
        .def("trim", &Class::trim, "Permanently trim curve between u1 and u2 by inserting knots and dropping useless ones",py::arg("u1"),py::arg("u2"),py::arg("permanently")=true)
        .def("changeBounds", py::overload_cast<T, T>(&Class::changeBounds), "Change parametrization to fit between k1 and k2")
        .def("changeBounds", py::overload_cast<const std::array<T, 2> &>(&Class::changeBounds), "Change parametrization to fit between k1 and k2")
        .def("increaseDegree", &Class::increaseDegree, "Increment curve's degree of step", py::arg("step")=1)
        .def("__copy__", [](const Class &self)
                { return Class(self); })
        .def(
            py::pickle(
                    [](const Class &crv) {
                            return py::make_tuple(// __getstate__
                                    crv.poles(),
                                    crv.knots(),
                                    crv.mults(),
                                    crv.degree()
                            );
                    },
                    [](py::tuple t) { // __setstate__
                            if (t.size() != 4)
                                    throw std::runtime_error("Invalid state!");
                            return Class{
                                    t[0].cast<gbs::points_vector<T,dim+rational>>(),
                                    t[1].cast<std::vector<T>>(),
                                    t[2].cast<std::vector<size_t>>(),
                                    t[3].cast<size_t>(),
                            };
                    }
            )
        )
        .def("__repr__", [](const Class &self) { return build_rep( self ); } )
    ;
    return cls;
}

template <typename T, size_t dim>
inline void gbs_bind_curves(py::module &m)
{
    //////////////////

    py::class_<Curve<T, dim>, std::shared_ptr<Curve<T, dim>>, gbs::Geom<T, dim>>(m, add_ext<dim>("Curve").c_str())
        .def(
            "value",
            &Curve<T, dim>::value,
            "Curve evaluation at given parameter at derivative order d",
            py::arg("u"),
            py::arg("d") = 0)
        .def(
            "value",
            [](const Curve<T, dim> &self, const std::vector<T> &u_lst, size_t d)
            {
                py::array points = py::cast(self.values(u_lst, d));
                return points;
            },
            py::arg("u_lst"),
            py::arg("d") = 0)
        .def(
            "bounds",
            &Curve<T, dim>::bounds,
            "Returns curves's start stop positions")
        .def(
            "midPosition",
            &Curve<T, dim>::midPosition,
            "Returns curves's start stop average position")
        .def(
            "begin",
            &Curve<T, dim>::begin,
            "Curve evaluation at begin at derivative order d",
            py::arg("d") = 0)
        .def(
            "end",
            &Curve<T, dim>::end,
            "Curve evaluation at end at derivative order d",
            py::arg("d") = 0)
        .def(
            "d_dm",
            &Curve<T, dim>::d_dm,
            "Curve first derivative respectively to the curvilinear abscissa",
            py::arg("m"))
        .def(
            "d_dm2",
            &Curve<T, dim>::d_dm,
            "Curve second derivative respectively to the curvilinear abscissa",
            py::arg("m"))
        .def(
            "__call__",
            &Curve<T, dim>::operator(),
            "Curve evaluation at given parameter at derivative order d",
            py::arg("u"),
            py::arg("d") = 0)
            .def(
                "__call__",
                [](const Curve<T, dim>& self, const std::vector<T>& u_lst, size_t d) {
                        py::array points = py::cast(self.values(u_lst, d));
                        return  points;
                },
                py::arg("u_lst"),
                py::arg("d") = 0
            )
        ;

    /////////////////

    py::class_<Line<T, dim>,  std::shared_ptr<Line<T, dim>>, Curve<T, dim> >(m,  add_ext<dim>("Line").c_str())
    .def(py::init<const point<T, dim> &, const point<T, dim>&>())
    .def(py::init<const ax1<T, dim> &>())
    .def(py::init<const Line<T, dim> &>())
    .def("__repr__", [](const Line<T, dim> &self) { return build_rep( self ); } )
    .def("__copy__", [](const Line<T, dim> &self) { return Line<T, dim>(self); })
    .def(py::pickle(
        [](const Line<T, dim> &crv) {
                return py::make_tuple(// __getstate__
                        crv.getAx()
                );
        },
        [](py::tuple t) { // __setstate__
                if (t.size() != 1)
                        throw std::runtime_error("Invalid state!");
                return Line<T, dim>{
                        t[0].cast<gbs::ax1<T,dim>>()
                };
        }
    ) )
    ;

    /////////////////

    py::class_<gbs::CurveOnSurface<T, dim>, std::shared_ptr<gbs::CurveOnSurface<T, dim>>, gbs::Curve<T, dim>>(m, add_ext<dim>("CurveOnSurface").c_str())
    .def(py::init<const std::shared_ptr<gbs::Curve<T, 2>> &, const std::shared_ptr<gbs::Surface<T, dim>>&>())
    .def(py::init<const gbs::CurveOnSurface<T, dim> &>())
    .def("__copy__", [](const gbs::CurveOnSurface<T, dim> &self){ return gbs::CurveOnSurface<T, dim>(self); })
    .def("__repr__", [](const gbs::CurveOnSurface<T, dim> &self) { return build_rep( self ); } )
    .def("basisCurve",&gbs::CurveOnSurface<T, dim>::basisCurve)
    .def("basisSurface",&gbs::CurveOnSurface<T, dim>::basisSurface)
    .def(py::pickle(
            [](CurveOnSurface<T, dim> &crv) {
                    // auto crv_ = crv.p_basisCurve();
                    // auto srf_ = crv.p_basisSurface();
                    // return py::make_tuple(// __getstate__
                    //         crv_,
                    //         srf_
                    // );
                    return py::make_tuple(// __getstate__
                            crv.basisCurve(),
                            crv.basisSurface()
                    );
            },
            [](py::tuple t) { // __setstate__
                    if (t.size() != 2)
                            throw std::runtime_error("Invalid state!");
                    return CurveOnSurface<T, dim>{
                            t[0].cast<std::shared_ptr<gbs::Curve<T, 2>>>(),
                            t[1].cast<std::shared_ptr<gbs::Surface<T, dim>>>()
                    };
            }
    ))
    ;


    //////////////////

    declare_bscurve<T, dim, false>(m);
    declare_bscurve<T, dim, true>(m).def(
            py::init<const gbs::points_vector<T, dim> &, const std::vector<T> &, const std::vector<T> &, size_t>(), 
            py::arg("poles"), 
            py::arg("knots_flats"), 
            py::arg("weights"), 
            py::arg("degree")
    )
    .def(
        py::init<const gbs::points_vector<T, dim> &, const std::vector<T> &, const std::vector<size_t> &, const std::vector<T> &, size_t>(), 
        py::arg("poles"), 
        py::arg("knots"), 
        py::arg("mults"), 
        py::arg("weights"), 
        py::arg("degree")
    )
    .def("polesProjected",&gbs::BSCurveRational<T, dim>::polesProjected)
    .def("weights",&gbs::BSCurveRational<T, dim>::weights)
    ;

    py::class_<gbs::CurveTrimmed<T, dim>, std::shared_ptr<gbs::CurveTrimmed<T, dim>>, gbs::Curve<T, dim>>(m, add_ext<dim>("CurveTrimmed").c_str())
    .def(py::init<const std::shared_ptr<gbs::Curve<T, dim>> &, T, T>(), py::arg("crv"), py::arg("u1"), py::arg("u2"))
    .def(py::init<const gbs::CurveTrimmed<T, dim> &>())
    // .def("__copy__", [](const gbs::CurveOnSurface<T, dim> &self){ return gbs::CurveOnSurface<T, dim>(self); })
    // .def("__repr__", [](const gbs::CurveOnSurface<T, dim> &self) { return build_rep( self ); } )
    ;


    
    py::class_<
            gbs::CurveOffset<T, dim, gbs::BSCfunction<T>>, 
            std::shared_ptr<gbs::CurveOffset<T, dim, gbs::BSCfunction<T>>>, 
            gbs::Curve<T, dim>
    >(m, add_ext<dim>("CurveOffset_bs_base").c_str() );


    py::class_<
            gbs::CurveExtended<T, dim>, 
            std::shared_ptr<gbs::CurveExtended<T, dim>>, 
            gbs::Curve<T, dim>
    >(m, add_ext<dim>("CurveExtended").c_str() )
        .def(py::init<std::shared_ptr<gbs::Curve<T, dim>>, bool>()
        , py::arg("crv")
        , py::arg("clamped") = false
        )
    ;


    py::class_<
            gbs::CurveOffset<T, dim, std::function<T(T, size_t)>>, 
            std::shared_ptr<gbs::CurveOffset<T, dim, std::function<T(T, size_t)>>>, 
            gbs::Curve<T, dim>
    >(m, add_ext<dim>("CurveOffset_func_base").c_str() );

}