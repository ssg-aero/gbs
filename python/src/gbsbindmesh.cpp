#include <gbs/mesh/tfi.h>
#include <gbs/render/vtkgridrender.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vtkProperty.h>
#include "vtk_bind.h"

namespace py = pybind11;

void gbs_bind_mesh(py::module &m) {
  m.def("msh_curves_lattice",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const std::vector<size_t> &, const std::vector<size_t> &,
            const std::shared_ptr<gbs::Surface<double, 3>>>(
            &gbs::msh_curves_lattice<double, 3, 1, 1>),
        "Compute bi-directional curve 3d lattice", py::arg("iso_ksi"),
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("eth_j"),
        py::arg("n_iso_ksi"), py::arg("n_iso_eth"), py::arg("p_srf") = nullptr);
  m.def("msh_curves_lattice",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const std::vector<size_t> &, const std::vector<size_t> &,
            const std::shared_ptr<gbs::Surface<double, 2>>>(
            &gbs::msh_curves_lattice<double, 2, 1, 1>),
        "Compute bi-directional curve 2d lattice", py::arg("iso_ksi"),
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("eth_j"),
        py::arg("n_iso_ksi"), py::arg("n_iso_eth"), py::arg("p_srf") = nullptr);
  m.def("msh_curves_lattice",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 1>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 1>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const std::vector<size_t> &, const std::vector<size_t> &,
            const std::shared_ptr<gbs::Surface<double, 1>>>(
            &gbs::msh_curves_lattice<double, 1, 1, 1>),
        "Compute bi-directional curve 1d lattice", py::arg("iso_ksi"),
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("eth_j"),
        py::arg("n_iso_ksi"), py::arg("n_iso_eth"), py::arg("p_srf") = nullptr);

  m.def("msh_curves_lattice",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &,
            const std::vector<std::vector<double>> &,
            const std::vector<std::vector<double>> &,
            const std::vector<size_t> &, const std::vector<size_t> &>(
            &gbs::msh_curves_lattice<double, 3, 1, 1>),
        "Compute bi-directional curve 3d lattice", py::arg("iso_ksi"),
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("eth_j"),
        py::arg("n_iso_ksi"), py::arg("n_iso_eth"));
  m.def("msh_curves_lattice",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &,
            const std::vector<std::vector<double>> &,
            const std::vector<std::vector<double>> &,
            const std::vector<size_t> &, const std::vector<size_t> &>(
            &gbs::msh_curves_lattice<double, 2, 1, 1>),
        "Compute bi-directional curve 2d lattice", py::arg("iso_ksi"),
        py::arg("iso_eth"), py::arg("ksi_i"), py::arg("eth_j"),
        py::arg("n_iso_ksi"), py::arg("n_iso_eth"));

  m.def("tfi_mesh",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 3>>> &, size_t,
            size_t, double>(&gbs::tfi_mesh_2d<double, 3, 1, 1, true>),
        py::arg("iso_ksi"), py::arg("iso_eth"), py::arg("n_ksi"),
        py::arg("n_eth"), py::arg("tol"));

  m.def("tfi_mesh",
        py::overload_cast<
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &,
            const std::vector<std::shared_ptr<gbs::Curve<double, 2>>> &, size_t,
            size_t, double>(&gbs::tfi_mesh_2d<double, 2, 1, 1, true>),
        py::arg("iso_ksi"), py::arg("iso_eth"), py::arg("n_ksi"),
        py::arg("n_eth"), py::arg("tol"));

  m.def("tfi_mesh",
        py::overload_cast<
            const std::vector<std::vector<std::array<gbs::point<double, 3>, 1>>>
                &,
            const std::vector<std::vector<std::array<gbs::point<double, 3>, 1>>>
                &,
            const std::vector<std::vector<
                std::array<std::array<gbs::point<double, 3>, 1>, 1>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const gbs::BSCfunction<double> &, const gbs::BSCfunction<double> &>(
            &gbs::tfi_mesh_2d<double, 3, 1, 1, true>),
        "Transfinite interpolation of point set", py::arg("X_ksi"),
        py::arg("X_eth"), py::arg("X_ksi_eth"), py::arg("ksi_i"),
        py::arg("eth_j"), py::arg("ksi"), py::arg("eth"));
  m.def("tfi_mesh",
        py::overload_cast<
            const std::vector<std::vector<std::array<gbs::point<double, 2>, 1>>>
                &,
            const std::vector<std::vector<std::array<gbs::point<double, 2>, 1>>>
                &,
            const std::vector<std::vector<
                std::array<std::array<gbs::point<double, 2>, 1>, 1>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const gbs::BSCfunction<double> &, const gbs::BSCfunction<double> &>(
            &gbs::tfi_mesh_2d<double, 2, 1, 1, true>),
        "Transfinite interpolation of point set", py::arg("X_ksi"),
        py::arg("X_eth"), py::arg("X_ksi_eth"), py::arg("ksi_i"),
        py::arg("eth_j"), py::arg("ksi"), py::arg("eth"));
  m.def("tfi_mesh",
        py::overload_cast<
            const std::vector<std::vector<std::array<gbs::point<double, 1>, 1>>>
                &,
            const std::vector<std::vector<std::array<gbs::point<double, 1>, 1>>>
                &,
            const std::vector<std::vector<
                std::array<std::array<gbs::point<double, 1>, 1>, 1>>> &,
            const std::vector<double> &, const std::vector<double> &,
            const gbs::BSCfunction<double> &, const gbs::BSCfunction<double> &>(
            &gbs::tfi_mesh_2d<double, 1, 1, 1, true>),
        "Transfinite interpolation of point set", py::arg("X_ksi"),
        py::arg("X_eth"), py::arg("X_ksi_eth"), py::arg("ksi_i"),
        py::arg("eth_j"), py::arg("ksi"), py::arg("eth"));
  m.def("tfi_mesh",
        py::overload_cast<const std::shared_ptr<gbs::Surface<double, 3>> &,
                          const std::vector<double> &,
                          const std::vector<double> &, size_t, size_t>(
            &gbs::tfi_mesh_2d<double, 3, 1, 1, true>),
        "Transfinite mesh on surface");
  m.def("tfi_mesh",
        py::overload_cast<const std::shared_ptr<gbs::Surface<double, 2>> &,
                          const std::vector<double> &,
                          const std::vector<double> &, size_t, size_t>(
            &gbs::tfi_mesh_2d<double, 2, 1, 1, true>),
        "Transfinite mesh on surface");
  m.def(
      "make_structuredgrid",
      py::overload_cast<const gbs::points_vector<double, 3> &, size_t, size_t>(
          &gbs::make_structuredgrid<double, 3>),
      "Make structured grid", py::arg("pts"), py::arg("ni"), py::arg("nj"));
  m.def(
      "make_structuredgrid",
      py::overload_cast<const gbs::points_vector<double, 2> &, size_t, size_t>(
          &gbs::make_structuredgrid<double, 2>),
      "Make structured grid", py::arg("pts"), py::arg("ni"), py::arg("nj"));
  m.def("project_points",
        py::overload_cast<const gbs::Surface<double, 3> &,
                          gbs::points_vector<double, 3> &, double>(
            &gbs::project_points_on_surface<double, 3>),
        "Project points on surface", py::arg("srf"), py::arg("pts"),
        py::arg("tol_projection") = 1e-6);
}