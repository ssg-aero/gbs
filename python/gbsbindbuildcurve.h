#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


#include <gbs/bscbuild.h>
using namespace gbs;

template <size_t dim>
inline auto add_ext(const char *func_name)
{
    return std::string( func_name ) + std::to_string( dim) + std::string("d");
}

template <typename T, size_t dim>
inline void gbs_bind_build_curve(py::module &m)
{
    m.def( add_ext<dim>("build_segment").c_str(),
        py::overload_cast<
            const point<T, dim> &,
            const point<T, dim> &,
            bool
        >(&build_segment<T,dim>),
        "Build BSCurve build segment between point p1 and point p2",
        py::arg("p1"),
        py::arg("p2"),
        py::arg("normalized_param")=false
    );

    m.def(add_ext<dim>("build_ellipse").c_str(),
        &build_ellipse<T,dim>,
        "Build NURBS definition of an ellipse, if dim > 2 the ellipse is in the xy plane",
        py::arg("radius1"),
        py::arg("radius2"),
        py::arg("center") = std::array<T, dim>{} 
    );

    m.def(add_ext<dim>("build_circle").c_str(),
        &build_circle<T,dim>,
        "Build NURBS definition of a circle, if dim > 2 the circle is in the xy plane",
        py::arg("radius"),
        py::arg("center") = std::array<T, dim>{} 
    );

    m.def("build_derivate",
        &build_derivate<T,dim>,
        "Builds the derivate curve, aka \n  d C(u) \n  ------  \n   d u",
        py::arg("crv")
    );

    m.def("build_integrate",
        &build_integrate<T,dim>,
        "Builds the integral curve from point P0 aka \n  /  \n  |   C(u) d u \n  /",
        py::arg("crv"),
        py::arg("P0") = std::array<T,dim>{}
    );

}