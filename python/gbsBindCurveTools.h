#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bsctools.h>
using namespace gbs;

/// `blossom` of the poles of d coordinates: d = 1 to 3 for the curves, 2 to 4 for the
/// homogeneous poles of the rational ones.
template <typename T, size_t d>
inline void gbs_bind_blossom(py::module &m)
{
    m.def(
        "blossom",
        [](const std::vector<std::array<T, d>> &poles, const std::vector<T> &knots, size_t degree, size_t span, const std::vector<T> &t) {
            return blossom(poles, knots, degree, span, t);
        },
        "The polar form (blossom) of the polynomial that the B-spline of poles, flat knots and degree follows on its span "
        "[knots[span], knots[span + 1]], at the degree values t (de Boor's algorithm, a value per level): symmetric, affine "
        "in each value, the curve where they are equal; pole i is the blossom at knots[i + 1 : i + degree + 1]. Of a "
        "rational curve, its homogeneous poles (the weight last): the blossom is homogeneous. ValueError for knots not "
        "of len(poles) + degree + 1 values, or decreasing, a span outside [degree, len(poles) - 1] or empty, t not of "
        "degree values",
        py::arg("poles"),
        py::arg("knots"),
        py::arg("degree"),
        py::arg("span"),
        py::arg("t")
    );
}

template <typename T, size_t dim, bool rational>
inline void gbs_bind_curveTools(py::module &m)
{
    m.def(
        "blossom",
        [](const BSCurveGeneral<T, dim, rational> &crv, size_t span, const std::vector<T> &t) { return blossom(crv, span, t); },
        "The blossom of the polynomial of the span [k[span], k[span + 1]] of the curve (k its flat knots) at the degree "
        "values t; homogeneous (the weight last) for a rational curve",
        py::arg("crv"),
        py::arg("span"),
        py::arg("t")
    );

    
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