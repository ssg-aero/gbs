#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bsctools.h>
using namespace gbs;

template <typename T, size_t dim, bool rational>
inline void gbs_bind_curveTools(py::module &m)
{
    
    m.def(
        "join",
        py::overload_cast<
            const BSCurveGeneral<T, dim, rational>&,
            const BSCurveGeneral<T, dim, rational>&
        >(&join<T, dim, rational, rational>),
        "Join 2 curves into one",
        py::arg("crv1"), py::arg("crv2")
    );

    m.def(
        "extention_to_point",
        py::overload_cast<
            const BSCurveGeneral<T, dim, rational>&,
            const point<T,dim> &,
            T,
            bool,
            std::optional<size_t>
        >(&extention_to_point<T, dim, rational>),
        "Build a curve as an extention of the curve crv to the point pt from position u following curve's direction",
        py::arg("crv"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );

    m.def(
        "extention_to_point",
        py::overload_cast<
            const BSCurveGeneral<T, dim, rational>&,
            const point<T,dim> &,
            T,
            bool,
            std::optional<size_t>
        >(&extention_to_point<T, dim, rational>),
        "Build a curve as an extention of the curve crv to the point pt from position u following curve's direction",
        py::arg("crv"),
        py::arg("pt"),
        py::arg("u"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );

    m.def(
        "extended_to_point",
        py::overload_cast<
            const BSCurveGeneral<T, dim, rational>&,
            const point<T,dim> &,
            bool,
            bool,
            std::optional<size_t>
        >(&extended_to_point<T, dim, rational>),
        "Extend curve's end/start to point, the curve definition is kept, i.e. between original curve's bound definition are identical",
        py::arg("crv"),
        py::arg("pt"),
        py::arg("at_end"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );

    m.def(
        "extended",
        py::overload_cast<
            const BSCurveGeneral<T, dim, rational>&,
            T,
            bool,
            bool,
            bool,
            std::optional<size_t>
        >(&extended<T, dim, rational>),
        "Extend curve's end by length, the curve definition is kept, i.e. between original curve's bound definition are identical",
        py::arg("crv"),
        py::arg("l"),
        py::arg("at_end"),
        py::arg("relative"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );
}