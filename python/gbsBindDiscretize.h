#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include <gbs/bscanalysis.h>
#include <gbs/bssanalysis.h>
using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_discretize(py::module &m)
{
    m.def(
        "discretize_curve_unif",
        [](const Curve<T, dim> &crv, size_t np)
        {
            py::array points = py::cast(discretize(crv, np));
            return points;
        },
        "Uniformly spaced discretization",
        py::arg("crv"), py::arg("np"));

    m.def(
        "discretize_curve",
        [](const Curve<T, dim> &crv, size_t np, T dev, size_t n_max_pts)
        {
            py::array points = py::cast(discretize(crv, np, dev));
            return points;
        },
        "Curve discretization based on deviation",
        py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);

    m.def(
        "deviation_based_params",
        [](const gbs::Curve<T, dim> &crv, size_t np, T dev, size_t n_max_pts)
        {
            py::array u_arr = py::cast(deviation_based_params(crv, np, dev));
            return u_arr;
        },
        "Curve discretization parameters based on deviation",
        py::arg("crv"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);

    m.def(
        "deviation_based_params",
        [](const gbs::Curve<T, dim> &crv, T u1, T u2, size_t np, T dev, size_t n_max_pts)
        {
            py::array u_arr = py::cast(deviation_based_params(crv, u1, u2, np, dev));
            return u_arr;
        },
        "Curve discretization based on deviation",
        py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np") = 30, py::arg("dev") = 0.01, py::arg("n_max_pts") = 5000);

    m.def(
        "uniform_distrib_params",
        [](const gbs::Curve<T, dim> &crv, T u1, T u2, size_t np, size_t n_law)
        {
            py::array u_arr = py::cast(uniform_distrib_params(crv, u1, u2, np, n_law));
            return u_arr;
        },
        "Uniform Curve discretization",
        py::arg("crv"), py::arg("u1"), py::arg("u2"), py::arg("np"), py::arg("n_law") = 30);

    m.def(
        "discretize_surface_unif",
        //    py::overload_cast<const gbs::Surface<double, 3> &, size_t, size_t>(&gbs::discretize<double, 3>),
        [](const gbs::Surface<T, dim> &srf, size_t nu, size_t nv)
        {
            py::array points = py::cast(discretize(srf, nu, nv));
            return points;
        },
        " ",
        py::arg("srf"), py::arg("nu"), py::arg("nv"));
}