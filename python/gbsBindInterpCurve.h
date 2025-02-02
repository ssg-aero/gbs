#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bscinterp.h>
using namespace gbs;

template <typename T, size_t dim>
inline void gbs_bind_interp_curve(py::module &m)
{
    ///////////////////
    // Mx continuity //
    ///////////////////
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<T, dim> &, size_t, KnotsCalcMode>(&interpolate<T, dim>),
          "Cn Curve interpolation, parametrization is computed according mode specification",
          py::arg("pts"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def("interpolate_cn",
          py::overload_cast<const points_vector<T, dim> &, const std::vector<T> &, size_t>(&interpolate<T, dim>),
          "Cn Curve interpolation, with specified parametrization",
          py::arg("pts"), py::arg("u"), py::arg("p"));

    // TODO: check if the folowing biding is used
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<constrType<T, dim, 1>> &, size_t, KnotsCalcMode>(&interpolate<T, dim>),
          "Cn Curve interpolation, parametrization is computed according mode specification",
          py::arg("Q"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    ////////////////
    // Bezier seg //
    ////////////////
    // C2
    m.def(
        "interpolate_c2",
        [](const std::vector<gbs::constrType<T, dim, 2>> &Q, gbs::KnotsCalcMode mode)
        {
            return interpolate<T, dim, 2>(Q, mode, true);
        },
        "C2 Curve interpolation with natural curvature, parametrization is computed according mode specification",
        py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def(
        "interpolate_c2",
        [](const std::vector<gbs::constrType<T, dim, 3>> &Q, gbs::KnotsCalcMode mode)
        {
            return interpolate<T, dim, 3>(Q, mode, false);
        },
        "C2 Curve interpolation with specified curvature, parametrization is computed according mode specification",
        py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def(
        "interpolate_c2",
        [](const std::vector<gbs::constrType<T, dim, 2>> &Q, const std::vector<double> &u)
        {
            return interpolate<T, dim, 2>(Q, u, true);
        },
        "C2 Curve interpolation with natural curvature, with specified parametrization",
        py::arg("Q"), py::arg("u") );

    m.def(
        "interpolate_c2",
        [](const std::vector<gbs::constrType<T, dim, 3>> &Q, const std::vector<double> &u)
        {
            return interpolate<T, dim, 3>(Q, u, false);
        },
        "C2 Curve interpolation with specified curvature, with specified parametrization",
        py::arg("Q"), py::arg("u") );

    // C1
    m.def(
        "interpolate_c1",
        [](const std::vector<gbs::constrType<T, dim, 1>> &Q, gbs::KnotsCalcMode mode)
        {
            return interpolate<T, dim, 1>(Q, mode, true);
        },
        "C1 Curve interpolation with natural curvature, parametrization is computed according mode specification",
        py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def(
        "interpolate_c1",
        [](const std::vector<std::array<T, dim>> &pts, gbs::KnotsCalcMode mode)
        {
            std::vector<constrType<T, dim, 1>> Q(pts.size());
            std::ranges::transform(pts, Q.begin(),[](const auto &pt){return constrType<T, dim, 1>{pt};});
            return interpolate<T, dim, 1>(Q, mode, true);
        },
        "C1 Curve interpolation with natural curvature, parametrization is computed according mode specification",
        py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def(
        "interpolate_c1",
        [](const std::vector<gbs::constrType<T, dim, 2>> &Q, gbs::KnotsCalcMode mode)
        {
            return interpolate<T, dim, 2>(Q, mode, false);
        },
        "C1 Curve interpolation with specified curvature, parametrization is computed according mode specification",
        py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def(
        "interpolate_c1",
        [](const std::vector<gbs::constrType<T, dim, 1>> &Q, const std::vector<double> &u)
        {
            return interpolate<T, dim, 1>(Q, u, true);
        },
        "C1 Curve interpolation with natural curvature, with specified parametrization",
        py::arg("Q"), py::arg("u"));

    m.def(
        "interpolate_c1",
        [](const std::vector<std::array<T, dim>> &pts, const std::vector<double> &u)
        {
            std::vector<constrType<T, dim, 1>> Q(pts.size());
            std::ranges::transform(pts, Q.begin(),[](const auto &pt){return constrType<T, dim, 1>{pt};});
            return interpolate<T, dim, 1>(Q, u, true);
        },
        "C1 Curve interpolation with natural curvature, with specified parametrization",
        py::arg("Q"), py::arg("u"));

    m.def(
        "interpolate_c1",
        [](const std::vector<gbs::constrType<T, dim, 2>> &Q, const std::vector<double> &u)
        {
            return interpolate<T, dim, 2>(Q, u, false);
        },
        "C1 Curve interpolation with specified curvature, with specified parametrization",
        py::arg("Q"), py::arg("u"));

    ///////////////////////////
    // General interpolation //
    ///////////////////////////

    m.def("interpolate",
          py::overload_cast<const bsc_bound<T, dim> &, const bsc_bound<T, dim> &, const std::vector<bsc_constraint<T, dim>> &, size_t>(&interpolate<T, dim>),
          "Arbitrary constrained interpolation, constrins are specified (u,pt,derivate_order), except at bounds the later are specifed (u,pt), derivatives then ca be added",
          py::arg("pt_begin"), py::arg("pt_end"), py::arg("cstr_lst"), py::arg("p"));


    //////////////
    // Utility //
    //////////////

    m.def("build_3pt_tg_vec",
            py::overload_cast<
                  const std::vector<std::array<T, dim>> &,
                  const std::vector<T> &
                  >(&build_3pt_tg_vec<T, dim>),"Bessel's method",
            py::arg("pts"), py::arg("u")
      );

    py::class_<constrPoint<T, dim>>(m, ("constrPoint"+ std::to_string( dim) + std::string("d")).c_str())
    .def(py::init<const std::array<T, dim> &, T, size_t>(), py::arg("pt"), py::arg("u"), py::arg("d") );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<constrPoint<T, dim>> &, const std::vector<T> &, size_t>(&build_poles<T, dim>),
        py::arg("Q"), py::arg("k_flat"), py::arg("d")
    );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<constrType<T,dim,1> > &, const std::vector<T> &, const std::vector<T> &, size_t>(&build_poles<T, dim, 1>),
        py::arg("Q"), py::arg("k_flat"), py::arg("u"), py::arg("d")
    );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<constrType<T,dim,2> > &, const std::vector<T> &, const std::vector<T> &, size_t>(&build_poles<T, dim, 2>),
        py::arg("Q"), py::arg("k_flat"), py::arg("u"), py::arg("d")
    );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<constrType<T,dim,3> > &, const std::vector<T> &, const std::vector<T> &, size_t>(&build_poles<T, dim, 3>),
        py::arg("Q"), py::arg("k_flat"), py::arg("u"), py::arg("d")
    );

    m.def(
        "build_poles",
        py::overload_cast<const std::vector<std::array<T, dim>> &, const std::vector<T> &, const std::vector<T> &, size_t>(&build_poles<T, dim>),
        py::arg("pts"), py::arg("k_flat"), py::arg("u"), py::arg("d")
    );

}