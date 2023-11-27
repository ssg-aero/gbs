#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bssinterp.h>
using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_interp_surface(py::module &m)
{
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<std::array<T, dim>> &, size_t, size_t, size_t, KnotsCalcMode>(&interpolate<T, dim>),
          "Cn Surface interpolation",
          py::arg("Q"), py::arg("n_poles_v"), py::arg("p"), py::arg("q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
}