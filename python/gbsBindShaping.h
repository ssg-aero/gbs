#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bscshaping.h>
#include <gbs/bssshaping.h>

using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_shaping(py::module &m)
{
    m.def( "move_to_point_delta",
        py::overload_cast<
            points_vector<T, dim> &,
            const point<T, dim> &,
            const std::vector<T> &,
            size_t,
            T,
            size_t,
            size_t
        >(&move_to_point_delta<T,dim>),
        "Edit poles to pass through pt at parameter u",
        py::arg("poles"),
        py::arg("D"),
        py::arg("k"),
        py::arg("p"),
        py::arg("u"),
        py::arg("i1"),
        py::arg("i2")
    );

    m.def( "moved_to_point",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const point<T, dim> &,
            T
        >(&moved_to_point<T,dim>),
        "Builds crv copy and edit poles to pass through pt at parameter u",
        py::arg("poles"),
        py::arg("pt"),
        py::arg("u")
    );

    m.def( "move_to_point",
        py::overload_cast<
            BSCurve<T,dim> &,
            const point<T, dim> &,
            T
        >(&move_to_point<T,dim>),
        "Edit curve's poles to pass through pt at parameter u",
        py::arg("poles"),
        py::arg("pt"),
        py::arg("u")
    );

    m.def( "move_to_constraints_delta",
        py::overload_cast<
            points_vector<T, dim> &,
            const std::vector<bsc_constraint<T,dim>> &,
            const std::vector<T> &,
            size_t,
            size_t,
            size_t
        >(&move_to_constraints_delta<T,dim>),
        "Edit poles to pass through constraints",
        py::arg("poles"),
        py::arg("constraints_delta"),
        py::arg("k"),
        py::arg("p"),
        py::arg("i1"),
        py::arg("i2")
    );

    m.def( "moved_to_constraints",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &
        >(&moved_to_constraints<T,dim>),
        "Builds crv copy and edit poles to match constraints",
        py::arg("crv"),
        py::arg("constraints")
    );

    m.def( "move_to_constraints",
        py::overload_cast<
            BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &
        >(&move_to_constraints<T,dim>),
        "Edit poles to match constraints",
        py::arg("crv"),
        py::arg("constraints")
    );
}