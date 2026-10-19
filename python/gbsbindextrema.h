#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/extrema.h>
using namespace gbs;

// The former forms with numbers of samples, extrema_curve_point(crv, pnt, tol, 30) or
// extrema_surf_pnt(srf, pnt, tol, 30, 30), have the arity of the forms with a start and a
// tolerance: they would run Newton from u0 = tol with a tolerance of 30. pybind11 3 accepts an
// int for a float even without conversion, so these overloads, registered first, catch strict
// ints there, positionally only (the grid of the surfaces is passed by keyword), and raise an
// explicit TypeError.
[[noreturn]] inline void removed_samples(const char *what)
{
    throw py::type_error(what);
}

template <typename T, size_t dim>
inline void gbs_bind_extrema_removed(py::module &m)
{
    m.def("extrema_curve_point",
        [](const Curve<T, dim> &, const std::array<T, dim> &, T, long) -> std::array<T, 2> {
            removed_samples("extrema_curve_point: n_bracket was removed (global minimum by subdivision); with a start, pass u0 and tol_u as floats");
        },
        py::arg("crv"), py::arg("pnt"), py::arg("tol_u"), py::arg("n_bracket").noconvert(), py::pos_only());
    m.def("extrema_curve_curve",
        [](const Curve<T, dim> &, const Curve<T, dim> &, T, long, long) -> std::array<T, 3> {
            removed_samples("extrema_curve_curve: n_bracket_u, n_bracket_v were removed (global minimum by subdivision); with a start, pass u10, u20 and tol_x as floats");
        },
        py::arg("crv1"), py::arg("crv2"), py::arg("tol_x"), py::arg("n_bracket_u").noconvert(), py::arg("n_bracket_v").noconvert(), py::pos_only());
    m.def("extrema_surf_pnt",
        [](const Surface<T, dim> &, const std::array<T, dim> &, T, long, long) -> std::array<T, 3> {
            removed_samples("extrema_surf_pnt: pass the grid by keyword (n_bracket_u=..., n_bracket_v=...); with a start, pass u0, v0 and tol_x as floats");
        },
        py::arg("srf"), py::arg("pnt"), py::arg("tol_x"), py::arg("n_bracket_u").noconvert(), py::arg("n_bracket_v").noconvert(), py::pos_only());
    m.def("extrema_surf_curve",
        [](const Surface<T, dim> &, const Curve<T, dim> &, T, long, long, long) -> std::array<T, 4> {
            removed_samples("extrema_surf_curve: pass the grid by keyword (n_bracket_u=..., n_bracket_v=..., n_bracket_w=...); with a start, pass u_c0, u_s0, v_s0 and tol_x as floats");
        },
        py::arg("srf"), py::arg("crv"), py::arg("tol_x"), py::arg("n_bracket_u").noconvert(), py::arg("n_bracket_v").noconvert(),
        py::arg("n_bracket_w").noconvert(), py::pos_only());
}

template <typename T, size_t dim>
inline void gbs_bind_extrema(py::module &m)
{
    gbs_bind_extrema_removed<T, dim>(m);

    m.def("extrema_curve_point",
        [](const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T tol_u) {
            return extrema_curve_point<T, dim>(crv, pnt, tol_u);
        },
        "Curve-Point minimum finding: {u, distance}",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("tol_u")
    );

    m.def("extrema_curve_point",
        [](const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T u0, T tol_u) {
            return extrema_curve_point<T, dim>(crv, pnt, u0, tol_u);
        },
        "Curve-Point minimum finding from u0: {u, distance}",
        py::arg("crv"),
        py::arg("pnt"),
        py::arg("u0"),
        py::arg("tol_u")
    );

    m.def(
        "extrema_surf_pnt",
        [](const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T u0, T v0, T tol_x) {
            return extrema_surf_pnt<T, dim>(srf, pnt, u0, v0, tol_x);
        },
        "Surface-Point minimum finding from (u0, v0): {u, v, distance}",
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("u0"),
        py::arg("v0"),
        py::arg("tol_x")=1e-6
    );

    m.def(
        "extrema_surf_pnt",
        [](const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, size_t n_bracket_u, size_t n_bracket_v) {
            return extrema_surf_pnt<T, dim>(srf, pnt, tol_x, {n_bracket_u, n_bracket_v});
        },
        "Surface-Point minimum finding: {u, v, distance}; the grid (surfaces outside 3D) is keyword-only",
        py::arg("srf"),
        py::arg("pnt"),
        py::arg("tol_x")=1e-6,
        py::kw_only(),
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30
    );

    m.def(
        "extrema_curve_curve",
        [](const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T tol_x) {
            return extrema_curve_curve<T, dim>(crv1, crv2, tol_x);
        },
        "Curve-Curve minimum finding: {u1, u2, distance}",
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("tol_x")=1e-6
    );

    m.def(
        "extrema_curve_curve",
        [](const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20, T tol_x) {
            return extrema_curve_curve<T, dim>(crv1, crv2, u10, u20, tol_x);
        },
        "Curve-Curve minimum finding from (u10, u20): {u1, u2, distance}",
        py::arg("crv1"),
        py::arg("crv2"),
        py::arg("u10"),
        py::arg("u20"),
        py::arg("tol_x")=1e-6
    );

    m.def(
        "extrema_surf_curve",
        [](const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, size_t n_bracket_u, size_t n_bracket_v, size_t n_bracket_w) {
            return extrema_surf_curve<T, dim>(srf, crv, tol_x, {n_bracket_u, n_bracket_v, n_bracket_w});
        },
        "Surface-Curve minimum finding: {u, v, u_c, distance}; the grid (surfaces outside 3D) is keyword-only",
        py::arg("srf"),
        py::arg("crv"),
        py::arg("tol_x")=1e-6,
        py::kw_only(),
        py::arg("n_bracket_u") = 30,
        py::arg("n_bracket_v") = 30,
        py::arg("n_bracket_w") = 30
    );

    m.def(
        "extrema_surf_curve",
        [](const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x) {
            return extrema_surf_curve<T, dim>(srf, crv, u_c0, u_s0, v_s0, tol_x);
        },
        "Surface-Curve minimum finding from (u_c0, u_s0, v_s0): {u, v, u_c, distance}",
        py::arg("srf"),
        py::arg("crv"),
        py::arg("u_c0"),
        py::arg("u_s0"),
        py::arg("v_s0"),
        py::arg("tol_x")=1e-6
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