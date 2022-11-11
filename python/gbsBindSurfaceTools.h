#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bsstools.h>
using namespace gbs;

template <typename T, size_t dim, bool rational>
inline void gbs_bind_surfaceTools(py::module &m)
{
    m.def(
        "add_dimension",
        py::overload_cast<
            const BSSurfaceGeneral<T, dim, rational> &,
            T
        >( &add_dimension<T, dim, rational> ),
        "Create a surface's copy with and additional dimension. The value of the dimension can be specified, the default is 0.",
        py::arg("srf"),
        py::arg("val")=T(0)
    );

    m.def(
        "extention_to_curve", &extention_to_curve<T, dim, rational>,
        py::arg("srf"), py::arg("crv"), py::arg("natural_end"), py::arg("max_cont")=std::nullopt
    );

    m.def(
        "join",
        py::overload_cast<
            const BSSurfaceGeneral<T,dim,rational> &,
            const BSSurfaceGeneral<T,dim,rational> &
        >(
            &join<T, dim, rational, rational>
        ),
        "Join surfaces in v direction",
        py::arg("srf1"),
        py::arg("srf2")
    );

    m.def(
        "split_u",
        py::overload_cast<
            const BSSurfaceGeneral<T,dim,rational> &,
            T
        >(
            &split_u<T, dim, rational>
        ),
        "Split surface at u parameter, aka in v direction",
        py::arg("srf"),
        py::arg("u")
    );

    m.def(
        "split_v",
        py::overload_cast<
            const BSSurfaceGeneral<T,dim,rational> &,
            T
        >(
            &split_v<T, dim, rational>
        ),
        "Split surface at v parameter, aka in u direction",
        py::arg("srf"),
        py::arg("v")
    );

    m.def(
        "extended_to_curve",
        py::overload_cast<
            const BSSurfaceGeneral<T,dim,rational> &,
            const BSCurveGeneral<T,dim,rational> &,
            bool,
            std::optional<size_t>
        >(extended_to_curve<T, dim, rational>),
        "Surface extention to curve",
        py::arg("srf"),
        py::arg("crv"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );

    m.def(
        "extended",
        py::overload_cast<
            const BSSurfaceGeneral<T,dim,rational> &,
            T,
            SurfaceBound,
            bool,
            std::optional<size_t>
        >(&extention<T, dim, rational>),
        "Surface extension in the given direction",
        py::arg("srf"),
        py::arg("l_ext"),
        py::arg("location"),
        py::arg("natural_end"),
        py::arg("max_cont") = std::nullopt
    );


}