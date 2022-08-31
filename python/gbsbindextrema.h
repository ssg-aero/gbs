#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/extrema.h>
using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_extrema(py::module &m)
{

    m.def("extrema_curve_point",
          py::overload_cast<
            const Curve<T, dim> &,
            const std::array<T, dim> &,
            T,
            nlopt::algorithm,
            size_t
        >(&extrema_curve_point<T, dim>),
        "Curve-Point minimum finding",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("tol_u"),
        py::arg("solver") = default_nlopt_algo,
        py::arg("n_bracket") = 30
    );

    m.def("extrema_curve_point",
          py::overload_cast<
            const Curve<T, dim> &,
            const std::array<T, dim> &,
            T,
            T,
            nlopt::algorithm
        >(&extrema_curve_point<T, dim>),
        "Curve-Point minimum finding",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("u0"),
        py::arg("tol_u"),
        py::arg("solver") = default_nlopt_algo
    );

    /* if nlopt::algorithm is still problematic */
    // m.def("extrema_curve_point",
    //     [](const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T tol_u)
    //     {
    //         return extrema_curve_point<T, dim>(crv, pnt, tol_u);
    //     },
    //     "Curve-Point minimum finding",
    //     py::arg("crv"),
    //     py::arg("pnt"),
    //     py::arg("tol_u")
    // );

    m.def(
        "extrema_surf_pnt", 
        py::overload_cast<
            const Surface<T, dim> &,
            const std::array<T, dim> &,
            T, 
            T,
            T,
            nlopt::algorithm
        >(&extrema_surf_pnt<T,dim>),
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("u0")=0.,
        py::arg("v0")=0.,
        py::arg("tol_x")=1e-6,
        py::arg("solver") = default_nlopt_algo
    );

    m.def(
        "extrema_surf_pnt", 
        py::overload_cast<
            const Surface<T, dim> &,
            const std::array<T, dim> &,
            T,
            nlopt::algorithm,
            size_t,
            size_t
        >(&extrema_surf_pnt<T,dim>),
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = default_nlopt_algo,
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30
    );

    m.def(
        "extrema_curve_curve",
        py::overload_cast<
            const Curve<T, dim> &,
            const Curve<T, dim> &,
            T,
            nlopt::algorithm,
            size_t,
            size_t
        >(&extrema_curve_curve<T,dim>),
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = default_nlopt_algo,
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30
    );

    m.def(
        "extrema_curve_curve",
        py::overload_cast<
            const Curve<T, dim> &,
            const Curve<T, dim> &,
            T,
            T,
            T,
            nlopt::algorithm
        >(&extrema_curve_curve<T,dim>),
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("u10"),
        py::arg("u20"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = default_nlopt_algo
    );

}