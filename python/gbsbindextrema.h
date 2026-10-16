#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
namespace py = pybind11;

#include <gbs/extrema.h>
using namespace gbs;

// The extrema are computed by Newton's method: `solver` is ignored and deprecated, it is
// removed with NLopt (docs/sources/design/nlopt_removal.md).
inline void warn_solver_deprecated(const std::optional<nlopt::algorithm> &solver)
{
    if (solver)
        PyErr_WarnEx(PyExc_DeprecationWarning, "extrema: 'solver' is ignored (Newton's method) and will be removed", 1);
}

template <typename T, size_t dim>
inline void gbs_bind_extrema(py::module &m)
{
    using Solver = std::optional<nlopt::algorithm>;

    m.def("extrema_curve_point",
        [](const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T tol_u, Solver solver, size_t n_bracket) {
            warn_solver_deprecated(solver);
            return extrema_curve_point<T, dim>(crv, pnt, tol_u, default_nlopt_algo, n_bracket);
        },
        "Curve-Point minimum finding: {u, distance}",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("tol_u"),
        py::arg("solver") = py::none(),
        py::arg("n_bracket") = 30
    );

    m.def("extrema_curve_point",
        [](const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T u0, T tol_u, Solver solver) {
            warn_solver_deprecated(solver);
            return extrema_curve_point<T, dim>(crv, pnt, u0, tol_u);
        },
        "Curve-Point minimum finding from u0: {u, distance}",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("u0"),
        py::arg("tol_u"),
        py::arg("solver") = py::none()
    );

    m.def(
        "extrema_surf_pnt",
        [](const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T u0, T v0, T tol_x, Solver solver) {
            warn_solver_deprecated(solver);
            return extrema_surf_pnt<T, dim>(srf, pnt, u0, v0, tol_x);
        },
        "Surface-Point minimum finding from (u0, v0): {u, v, distance}",
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("u0"),
        py::arg("v0"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none()
    );

    m.def(
        "extrema_surf_pnt",
        [](const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, Solver solver, size_t n_bracket_u, size_t n_bracket_v) {
            warn_solver_deprecated(solver);
            return extrema_surf_pnt<T, dim>(srf, pnt, tol_x, default_nlopt_algo, n_bracket_u, n_bracket_v);
        },
        "Surface-Point minimum finding: {u, v, distance}",
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none(),
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30
    );

    m.def(
        "extrema_curve_curve",
        [](const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T tol_x, Solver solver, size_t n_bracket_u, size_t n_bracket_v) {
            warn_solver_deprecated(solver);
            return extrema_curve_curve<T, dim>(crv1, crv2, tol_x, default_nlopt_algo, n_bracket_u, n_bracket_v);
        },
        "Curve-Curve minimum finding: {u1, u2, distance}",
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none(),
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30
    );

    m.def(
        "extrema_curve_curve",
        [](const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20, T tol_x, Solver solver) {
            warn_solver_deprecated(solver);
            return extrema_curve_curve<T, dim>(crv1, crv2, u10, u20, tol_x);
        },
        "Curve-Curve minimum finding from (u10, u20): {u1, u2, distance}",
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("u10"),
        py::arg("u20"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none()
    );

    m.def(
        "extrema_surf_curve",
        [](const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, Solver solver, size_t n_bracket_u, size_t n_bracket_v, size_t n_bracket_w) {
            warn_solver_deprecated(solver);
            return extrema_surf_curve<T, dim>(srf, crv, tol_x, default_nlopt_algo, n_bracket_u, n_bracket_v, n_bracket_w);
        },
        "Surface-Curve minimum finding: {u, v, u_c, distance}",
        py::arg("srf"),
        py::arg("crv"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none(),
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30,
        py::arg("n_bracket_w") = 30
    );

    m.def(
        "extrema_surf_curve",
        [](const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x, Solver solver) {
            warn_solver_deprecated(solver);
            return extrema_surf_curve<T, dim>(srf, crv, u_c0, u_s0, v_s0, tol_x);
        },
        "Surface-Curve minimum finding from (u_c0, u_s0, v_s0): {u, v, u_c, distance}",
        py::arg("srf"),
        py::arg("crv"),
        py::arg("u_c0"),
        py::arg("u_s0"),
        py::arg("v_s0"),
        py::arg("tol_x")=1e-6,
        py::arg("solver") = py::none()
    );

    m.def(
        "solve_coordinate",
        &solve_coordinate<T, dim>,
        py::arg("crv"),
        py::arg("iCoord"),
        py::arg("Xi"),
        py::arg("u0")= std::nullopt,
        py::arg("digits")=20,
        py::arg("maxit")=20
    );

    m.def(
        "max_coordinate",
        &max_coordinate<T, dim>,
        py::arg("crv"),
        py::arg("iCoord"),
        py::arg("Xi"),
        // py::arg("u0")= std::nullopt,
        py::arg("digits")=20,
        py::arg("maxit")=20
    );

    m.def(
        "min_coordinate",
        &min_coordinate<T, dim>,
        py::arg("crv"),
        py::arg("iCoord"),
        py::arg("Xi"),
        // py::arg("u0")= std::nullopt,
        py::arg("digits")=20,
        py::arg("maxit")=20
    );
}