#include <gbs/render/vtkcurvesrender.h>
#include <gbs/render/vtkgridrender.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include "vtk_bind.h"

namespace py = pybind11;
// inline void f_plot_curves_2d(const
// std::vector<std::shared_ptr<gbs::Curve<double,2>>>
// &crv_lst){gbs::plot(crv_lst);}; inline void f_plot_curves(const
// std::vector<std::shared_ptr<gbs::Curve<double,3>>>
// &crv_lst){gbs::plot(crv_lst);};
inline void
f_plot_curves_2d(const std::vector<gbs::BSCurve<double, 2>> &crv_lst) {
  gbs::plot(crv_lst);
};
inline void f_plot_curves(const std::vector<gbs::BSCurve<double, 3>> &crv_lst) {
  gbs::plot(crv_lst);
};
// inline auto f_make_curve3d_actor(const gbs::Curve<double, 3>& crv,
// std::array<double,3>  col, size_t np){return
// py::cast(gbs::make_actor(crv,col,np));}
inline vtkSmartPointer<vtkActor>
f_make_curve3d_actor(const gbs::Curve<double, 3> &crv,
                     std::array<double, 3> col, size_t np, double dev) {
  return gbs::make_actor(crv, col, np, dev);
}
inline vtkSmartPointer<vtkActor>
f_make_surf3d_actor(const gbs::Surface<double, 3> &srf,
                    std::array<double, 3> col, size_t n1, size_t n2) {
  return gbs::make_actor(srf, col, n1, n2);
}
inline vtkSmartPointer<vtkActor>
f_make_curve2d_actor(const gbs::Curve<double, 2> &crv,
                     std::array<double, 3> col, size_t np, double dev) {
  return gbs::make_actor(crv, col, np, dev);
}
inline vtkSmartPointer<vtkActor>
f_make_surf2d_actor(const gbs::Surface<double, 2> &srf,
                    std::array<double, 3> col, size_t n1, size_t n2) {
  return gbs::make_actor(srf, col, n1, n2);
}
// inline auto f_discretize_curve(const gbs::Curve<double,3> &crv, size_t n,
// double dev_max, size_t n_max_pts){return
// discretize(crv,n,dev_max,n_max_pts);}

void gbs_bind_render(py::module &m) {
  m.def("plot_curves_2d", &f_plot_curves_2d);
  m.def("plot_curves", &f_plot_curves);
  m.def("make_curve3d_actor", &f_make_curve3d_actor, py::arg("crv"),
        py::arg("col") =
            std::array<double, 3>{51. / 255., 161. / 255., 201. / 255.},
        py::arg("np"), py::arg("dev") = 0.01);
  m.def("make_surf3d_actor", &f_make_surf3d_actor, py::arg("srf"),
        py::arg("col") =
            std::array<double, 3>{51. / 255., 161. / 255., 201. / 255.},
        py::arg("nu") = 200, py::arg("nv") = 200);

  m.def("make_actor",
        py::overload_cast<const gbs::Curve<double, 3> &, std::array<double, 3>,
                          size_t, double>(&gbs::make_actor<double, 3>),
        "Discretize curve and make vkt actor", py::arg("bsc"),
        py::arg("col") =
            std::array<double, 3>{255. / 255., 99. / 255., 71. / 255},
        py::arg("np") = 100, py::arg("dev") = 0.01);
  m.def("make_actor",
        py::overload_cast<const gbs::Curve<double, 2> &, std::array<double, 3>,
                          size_t, double>(&gbs::make_actor<double, 2>),
        "Discretize curve and make vkt actor", py::arg("bsc"),
        py::arg("col") =
            std::array<double, 3>{255. / 255., 99. / 255., 71. / 255},
        py::arg("np") = 100, py::arg("dev") = 0.01);
  m.def(
      "make_actor",
      py::overload_cast<const gbs::Surface<double, 3> &, std::array<double, 3>,
                        size_t, size_t>(&gbs::make_actor<double, 3>),
      "Discretize surface and make actor", py::arg("srf"),
      py::arg("col") =
          std::array<double, 3>{51. / 255., 161. / 255., 201. / 255.},
      py::arg("n1") = 200, py::arg("n2") = 200);
  m.def(
      "make_actor",
      py::overload_cast<const gbs::Surface<double, 2> &, std::array<double, 3>,
                        size_t, size_t>(&gbs::make_actor<double, 2>),
      "Discretize surface and make actor", py::arg("srf"),
      py::arg("col") =
          std::array<double, 3>{51. / 255., 161. / 255., 201. / 255.},
      py::arg("n1") = 200, py::arg("n2") = 200);
  m.def("make_actor",
        py::overload_cast<const gbs::points_vector<double, 3> &, double, bool,
                          std::array<double, 3> &>(&gbs::make_actor<double, 3>),
        "Actor from points", py::arg("pts"), py::arg("pt_size") = 5.,
        py::arg("render_as_sphere") = true,
        py::arg("col") = std::array<double, 3>{0.3, 0.3, 0.3});
  m.def("make_actor",
        py::overload_cast<const gbs::points_vector<double, 2> &, double, bool,
                          std::array<double, 3> &>(&gbs::make_actor<double, 2>),
        "Actor from points", py::arg("pts"), py::arg("pt_size") = 5.,
        py::arg("render_as_sphere") = true,
        py::arg("col") = std::array<double, 3>{0.3, 0.3, 0.3});
}