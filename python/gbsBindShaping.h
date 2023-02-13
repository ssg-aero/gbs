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
    // Curves
    m.def( "find_span_range",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &
        >(&find_span_range<T,dim>),
        "Find pole range affected by constrains",
        py::arg("crv"),
        py::arg("constraints")
    );

    m.def( "convert_constraints_to_delta",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &
        >(&convert_constraints_to_delta<T,dim>),
        "Convert points constrains into delta",
        py::arg("crv"),
        py::arg("constraints")
    );

    m.def( "moved_to_point",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const point<T, dim> &,
            T,
            bool
        >(&moved_to_point<T,dim>),
        "Builds crv copy and edit poles to pass through pt at parameter u",
        py::arg("poles"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("fix_bounds")=false
    );

    m.def( "move_to_point",
        py::overload_cast<
            BSCurve<T,dim> &,
            const point<T, dim> &,
            T,
            bool
        >(&move_to_point<T,dim>),
        "Edit curve's poles to pass through pt at parameter u",
        py::arg("poles"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("fix_bounds")=false
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

    m.def( "moved_to_constraints_delta",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &,
            size_t,
            size_t
        >(&moved_to_constraints_delta<T,dim>),
        "Edit poles to pass through constraints",
        py::arg("crv"),
        py::arg("constraints_delta"),
        py::arg("i1"),
        py::arg("i2")
    );

    m.def( "moved_to_constraints",
        py::overload_cast<
            const BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &,
            bool
        >(&moved_to_constraints<T,dim>),
        "Builds crv copy and edit poles to match constraints",
        py::arg("crv"),
        py::arg("constraints"),
        py::arg("fix_bounds")=false
    );

    m.def( "move_to_constraints",
        py::overload_cast<
            BSCurve<T,dim> &,
            const std::vector<bsc_constraint<T,dim>> &,
            bool
        >(&move_to_constraints<T,dim>),
        "Edit poles to match constraints",
        py::arg("crv"),
        py::arg("constraints"),
        py::arg("fix_bounds")=false
    );
    // Surfaces
    // m.def( "move_to_constraints_delta",
    //     py::overload_cast<
    //         points_vector<T, dim> &,
    //         const std::vector<bss_constraint<T,dim>> &,
    //         const std::vector<T> &,
    //         const std::vector<T> &,
    //         size_t,
    //         size_t,
    //         size_t,
    //         size_t,
    //         size_t,
    //         size_t
    //     >(&move_to_constraints_delta<T,dim>),
    //     "Edit poles to match constraints",
    //     py::arg("poles"),
    //     py::arg("constraints_delta"),
    //     py::arg("ku"),
    //     py::arg("kv"),
    //     py::arg("p"),
    //     py::arg("q"),
    //     py::arg("i1"),
    //     py::arg("i2"),
    //     py::arg("j1"),
    //     py::arg("j2")
    // );

    m.def( "moved_to_constraints_delta",
        py::overload_cast<
            const BSSurface<T,dim> &,
            const std::vector<bss_constraint<T,dim>> &,
            size_t,
            size_t,
            size_t,
            size_t
        >(&moved_to_constraints_delta<T,dim>),
        "Edit poles to match constraints",
        py::arg("poles"),
        py::arg("constraints_delta"),
        py::arg("i1"),
        py::arg("i2"),
        py::arg("j1"),
        py::arg("j2")
    );

    m.def( "moved_to_constraints",
        py::overload_cast<
            const BSSurface<T,dim> &,
            const std::vector<bss_constraint<T,dim>> &,
            bool
        >(&moved_to_constraints<T,dim>),
        "Builds srf copy and edit poles to match constraints",
        py::arg("srf"),
        py::arg("constraints"),
        py::arg("fix_bounds")=false
    );

    m.def( "moved_to_point",
        py::overload_cast<
            const BSSurface<T,dim> &,
            const point<T, dim> &,
            T,
            T,
            bool
        >(&moved_to_point<T,dim>),
        "Builds srf copy and edit poles to match constraints",
        py::arg("srf"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("v"),
        py::arg("fix_bounds")=false
    );

    m.def( "move_to_point",
        py::overload_cast<
            BSSurface<T,dim> &,
            const point<T, dim> &,
            T,
            T,
            bool
        >(&move_to_point<T,dim>),
        "Builds srf copy and edit poles to match constraints",
        py::arg("srf"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("v"),
        py::arg("fix_bounds")=false
    );

}