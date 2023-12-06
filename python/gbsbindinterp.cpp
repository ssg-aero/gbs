#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>

#include "gbsBindInterpCurve.h"
#include "gbsBindInterpSurface.h"

using namespace gbs;

void gbs_bind_interp(py::module &m)
{

    m.def("interpolate_cn",
          py::overload_cast<const std::vector<double> &, const std::vector<double> &, size_t>(&interpolate<double>),
          "CN interpolation",
          py::arg("Q"), py::arg("u"), py::arg("p"));

      gbs_bind_interp_curve<double,1>(m);
      gbs_bind_interp_curve<double,2>(m);
      gbs_bind_interp_curve<double,3>(m);

      gbs_bind_interp_surface<double,1>(m);
      gbs_bind_interp_surface<double,2>(m);
      gbs_bind_interp_surface<double,3>(m);

}