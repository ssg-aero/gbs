#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
using namespace gbs;

void gbs_bind_interp(py::module &m)
{

    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 3> &, size_t, KnotsCalcMode>(&interpolate<double, 3>),
          "Cn interpolation",
          py::arg("pts"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 2> &, size_t, KnotsCalcMode>(&interpolate<double, 2>),
          "Cn interpolation",
          py::arg("pts"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 1> &, size_t, KnotsCalcMode>(&interpolate<double, 1>),
          "Cn interpolation",
          py::arg("pts"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<double> &, const std::vector<double> &, size_t>(&interpolate<double>),
          "CN interpolation",
          py::arg("Q"), py::arg("u"), py::arg("p"));
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 3> &, const std::vector<double> &, size_t>(&interpolate<double, 3>),
          "Cn interpolation",
          py::arg("pts"), py::arg("u"), py::arg("p"));
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 2> &, const std::vector<double> &, size_t>(&interpolate<double, 2>),
          "Cn interpolation",
          py::arg("pts"), py::arg("u"), py::arg("p"));
    m.def("interpolate_cn",
          py::overload_cast<const points_vector<double, 1> &, const std::vector<double> &, size_t>(&interpolate<double, 1>),
          "Cn interpolation",
          py::arg("pts"), py::arg("u"), py::arg("p"));

    m.def("interpolate_cn",
          py::overload_cast<const std::vector<constrType<double, 3, 1>> &, size_t, KnotsCalcMode>(&interpolate<double, 3>),
          "Cn interpolation",
          py::arg("Q"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<constrType<double, 2, 1>> &, size_t, KnotsCalcMode>(&interpolate<double, 2>),
          "Cn interpolation",
          py::arg("Q"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<constrType<double, 1, 1>> &, size_t, KnotsCalcMode>(&interpolate<double, 1>),
          "Cn interpolation",
          py::arg("Q"), py::arg("p"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def("interpolate_cn",
          py::overload_cast<const std::vector<std::array<double, 3>> &, size_t, size_t, size_t, KnotsCalcMode>(&interpolate<double, 3>),
          "Cn interpolation",
          py::arg("Q"), py::arg("n_poles_v"), py::arg("p"), py::arg("q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<std::array<double, 2>> &, size_t, size_t, size_t, KnotsCalcMode>(&interpolate<double, 2>),
          "Cn interpolation",
          py::arg("Q"), py::arg("n_poles_v"), py::arg("p"), py::arg("q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_cn",
          py::overload_cast<const std::vector<std::array<double, 1>> &, size_t, size_t, size_t, KnotsCalcMode>(&interpolate<double, 1>),
          "Cn interpolation",
          py::arg("Q"), py::arg("n_poles_v"), py::arg("p"), py::arg("q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);

    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 3, 2>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 3, 2>( Q, mode, true);
                }
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 2, 2>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 2, 2>( Q, mode, true);
                }  
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 3, 3>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 3, 3>( Q, mode, false);
                }
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 2, 3>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 2, 3>( Q, mode, false);
                }  
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
/////////////
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 3, 2>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 3, 2>( Q, mode, false);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 2, 2>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 2, 2>( Q, mode, false);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 3, 1>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 3, 1>( Q, mode, true);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 2, 1>> &Q, gbs::KnotsCalcMode mode)
                {
                      return interpolate<double, 2, 1>( Q, mode, true);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("mode") = KnotsCalcMode::CHORD_LENGTH);
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 3, 2>> &Q, const std::vector<double> &u)
                {
                      return interpolate<double, 3, 2>( Q, u, false);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("u"));
    m.def("interpolate_c1",
                [](const std::vector<gbs::constrType<double, 2, 2>> &Q, const std::vector<double> &u)
                {
                      return interpolate<double, 2, 2>( Q, u, false);
                }  
          ,
          "C1 interpolation",
          py::arg("Q"), py::arg("u"));
//////////////
    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 3, 2>> &Q, const std::vector<double> &u)
                {
                      return interpolate<double, 3, 2>( Q, u, true);
                }  
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("u"));
    m.def("interpolate_c2",
                [](const std::vector<gbs::constrType<double, 2, 2>> &Q, const std::vector<double> &u)
                {
                      return interpolate<double, 2, 2>( Q, u, true);
                }  
          ,
          "C2 interpolation",
          py::arg("Q"), py::arg("u"));
    m.def("interpolate",
          py::overload_cast<const bsc_bound<double, 3> &, const bsc_bound<double, 3> &, const std::vector<bsc_constraint<double, 3>> &, size_t>(&interpolate<double, 3>),
          "Arbitrary constrained interpolation, constrins are specified (u,pt,derivate_order), except at bounds the later are specifed (u,pt), derivatives then ca be added",
          py::arg("pt_begin"), py::arg("pt_end"), py::arg("cstr_lst"), py::arg("p"));
    m.def("interpolate",
          py::overload_cast<const bsc_bound<double, 2> &, const bsc_bound<double, 2> &, const std::vector<bsc_constraint<double, 2>> &, size_t>(&interpolate<double, 2>),
          "Arbitrary constrained interpolation, constrins are specified (u,pt,derivate_order), except at bounds the later are specifed (u,pt), derivatives then ca be added",
          py::arg("pt_begin"), py::arg("pt_end"), py::arg("cstr_lst"), py::arg("p"));
    m.def("interpolate",
          py::overload_cast<const bsc_bound<double, 1> &, const bsc_bound<double, 1> &, const std::vector<bsc_constraint<double, 1>> &, size_t>(&interpolate<double, 1>),
          "Arbitrary constrained interpolation, constrins are specified (u,pt,derivate_order), except at bounds the later are specifed (u,pt), derivatives then ca be added",
          py::arg("pt_begin"), py::arg("pt_end"), py::arg("cstr_lst"), py::arg("p"));

      m.def("build_3pt_tg_vec",
            py::overload_cast<
                  const std::vector<std::array<double, 3>> &,
                  const std::vector<double> &
                  >(&build_3pt_tg_vec<double,3>),"Bessel's method",
            py::arg("pts"), py::arg("u")
      );
      m.def("build_3pt_tg_vec",
            py::overload_cast<
                  const std::vector<std::array<double, 2>> &,
                  const std::vector<double> &
                  >(&build_3pt_tg_vec<double,2>),"Bessel's method",
            py::arg("pts"), py::arg("u")
      );
      m.def("build_3pt_tg_vec",
            py::overload_cast<
                  const std::vector<std::array<double, 1>> &,
                  const std::vector<double> &
                  >(&build_3pt_tg_vec<double,1>),"Bessel's method",
            py::arg("pts"), py::arg("u")
      );
}