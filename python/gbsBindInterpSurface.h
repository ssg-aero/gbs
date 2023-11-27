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

    /////////////
    // Utility //
    /////////////
    m.def(
        "build_poles",
        py::overload_cast<
            const std::vector<std::array<T, dim>> &, 
            const std::vector<T> &, 
            const std::vector<T> &, 
            const std::vector<T> &, 
            const std::vector<T> &, 
            size_t,
            size_t
        >(&build_poles<T, dim>),
        py::arg("Q"), py::arg("k_flat_u"), py::arg("k_flat_v"), py::arg("u"), py::arg("v"), py::arg("p"), py::arg("q")
    );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<bss_constraint<T, dim>> &, const std::vector<T> &, const std::vector<T> &, size_t, size_t>(&build_poles<T, dim>),
        py::arg("Q"), py::arg("k_flat_u"), py::arg("k_flat_v"), py::arg("p"), py::arg("q")
    );

}