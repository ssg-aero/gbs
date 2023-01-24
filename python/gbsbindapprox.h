#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bscapprox.h>
#include <gbs/bssapprox.h>
using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_approx(py::module &m)
{
    m.def(
        "approx",
        py::overload_cast<
            const points_vector<T, dim> &, 
            size_t,
            KnotsCalcMode,
            bool,
            T,
            T,
            size_t
        >(&approx<T, dim>),
        "Points curves fitting",
        py::arg("pts"), 
        py::arg("deg"), 
        py::arg("mode") = KnotsCalcMode::CHORD_LENGTH,
        py::arg("fix_bound")= true, 
        py::arg("d_max")= 1e-3, 
        py::arg("d_avg")= 1e-4, 
        py::arg("n_max")= 200
    );

    m.def(
        "approx",
        py::overload_cast<
            const points_vector<T, dim> &, 
            size_t,
            size_t, 
            KnotsCalcMode
        >(&approx<T, dim>),
        "Points curves fitting",
        py::arg("pts"), 
        py::arg("deg"), 
        py::arg("n_poles"), 
        py::arg("mode") = KnotsCalcMode::CHORD_LENGTH
    );

    // m.def(
    //     "refine_approx",
    //     py::overload_cast<
    //         const points_vector<T, dim> &,
    //         const std::vector<T> &,
    //         const gbs::BSCurve<T, dim> &,
    //         bool,
    //         T,
    //         T,
    //         size_t
    //     >(&approx<T, dim>),
    //     py::arg("pts"),
    //     py::arg("pts"),
    //     py::arg("pts"),
    //     py::arg("pts"),
    //     py::arg("pts"),
    //     py::arg("pts"),
    //     py::arg("pts"),
    // );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::Curve<T,dim> &,
            T,
            T,
            size_t,
            size_t
        >(&gbs::approx<T,dim>),
        "Approximate curve respecting original curve's parametrization",
        py::arg("crv"), 
        py::arg("deviation"), 
        py::arg("tol"),
        py::arg("p"),
        py::arg("bp") = 30
    );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::Curve<T,dim> &,
            const std::function<std::array<T,dim>(const std::array<T,dim>&)> &,
            T,
            T,
            size_t,
            size_t
        >(&gbs::approx<T,dim>),
        "Approximate curve's transformation respecting original curve's parametrization",
        py::arg("crv"), 
        py::arg("trf_func"), 
        py::arg("deviation"), 
        py::arg("tol"),
        py::arg("p"),
        py::arg("bp") = 30
    );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::Curve<T,dim> &,
            T,
            T,
            size_t,
            size_t,
            size_t
        >(&gbs::approx<T,dim>),
        "Approximate curve respecting original curve's parametrization between 2 parameters",
        py::arg("crv"),
        py::arg("deviation"),
        py::arg("tol"),
        py::arg("p"),
        py::arg("n_poles"),
        py::arg("bp") = 30
    );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::Curve<T,dim> &,
            T,
            T,
            T,
            T,
            size_t,
            size_t
        >(&gbs::approx<T,dim>),
        "Approximate curve respecting original curve's parametrization",
        py::arg("crv"),
        py::arg("u1"),
        py::arg("u2"),
        py::arg("deviation"),
        py::arg("tol"),
        py::arg("p"),
        py::arg("bp") = 30
    );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::Curve<T,dim> &,
            size_t,
            size_t
        >(&gbs::approx<T,dim>),
        "Approximate curve with uniform point prepartition",
        py::arg("crv"),
        py::arg("p"),
        py::arg("np")
    );

    m.def(
        "approx",
        py::overload_cast<
            const gbs::points_vector<T,dim> &,
            size_t,
            size_t,
            const std::vector<T> &,
            bool
        >(&gbs::approx<T,dim>),
        "Approximate curve passing thru points at givent parameter",
        py::arg("pts"),
        py::arg("p"),
        py::arg("n_poles"),
        py::arg("u"),
        py::arg("fix_bound") = true
    );

    m.def(
        "approx_bound_fixed",
        py::overload_cast<
            const gbs::points_vector<T,dim> &,
            size_t,
            size_t,
            const std::vector<T> &,
            const std::vector<T> &
        >(&gbs::approx_bound_fixed<T,dim>),
        "Approximate curve passing thru points at givent parameter",
        py::arg("pts"),
        py::arg("p"),
        py::arg("n_poles"),
        py::arg("u"),
        py::arg("k_flat")
    );

    m.def(
        "approx",
        py::overload_cast<
            const std::vector<bss_constrain<T, dim>> &, 
            const std::vector<T> &,
            const std::vector<T> &,
            size_t,
            size_t>(&approx<T, dim>),
        "Constrains surface fitting",
        py::arg("constrains"), 
        py::arg("k_flaut_u"), 
        py::arg("k_flaut_v"), 
        py::arg("p"), 
        py::arg("q")
    );
}