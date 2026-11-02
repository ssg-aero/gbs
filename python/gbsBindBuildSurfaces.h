#pragma
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gbs/bssbuild.h>

namespace py = pybind11;

template <typename T, size_t dim>
inline void gbs_bind_build_surface(py::module &m)
{

    using namespace gbs;

    m.def(
        "loftbs",
        [](const std::list<gbs::BSCurve<T, dim>> &bs_lst, size_t v_degree_max)
        {
            return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), v_degree_max);
        },
        "Loft Surface using BSCurve.",
        py::arg("bs_lst"), py::arg("v_degree_max"));

    m.def(
        "loftbs",
        [](const std::list<gbs::BSCurve<T, dim>> &bs_lst, const std::vector<T> &v, size_t v_degree)
        {
            return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), v, v_degree);
        },
        "Loft Surface using BSCurve as v isoparametrics curves. Each isoV(v[j]) is coincident with curve bs_lst[j].",
        py::arg("bs_lst"), py::arg("v"), py::arg("v_degree"));

    m.def(
        "loftbs",
        [](const std::list<gbs::BSCurve<T, dim>> &bs_lst, const std::vector<T> &v, const std::vector<double> &flat_v, size_t v_degree)
        {
            return loft_generic<T, dim>(bs_lst.begin(), bs_lst.end(), v, flat_v, v_degree);
        },
        "Loft Surface using BSCurve with imposed v direction parametrization. Each isoV(v[j]) is coincident with curve bs_lst[j].",
        py::arg("bs_lst"), py::arg("v"), py::arg("flat_v"), py::arg("v_degree"));

    m.def("loft",
          py::overload_cast<
              const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &,
              size_t,
              double,
              size_t,
              size_t>(&gbs::loft<T, dim>),
          "Loft Surface using BSCurve approximations of curves (if base type is not BSCurve).",
          py::arg("crv_lst"), py::arg("v_degree_max") = 3, py::arg("dev") = 0.01, py::arg("np") = 100, py::arg("deg_approx") = 5);

    m.def("loft",
          py::overload_cast<
              const std::vector<std::shared_ptr<gbs::Curve<T, dim>>> &,
              const std::vector<double> &,
              size_t,
              double,
              size_t,
              size_t>(&gbs::loft<T, dim>),
          "Loft Surface using BSCurve approximations of curves (if base type is not BSCurve). Each isoV(v[j]) is coincident with curve crv_lst[j].",
          py::arg("crv_lst"), py::arg("v"), py::arg("q"), py::arg("dev") = 0.01, py::arg("np") = 100, py::arg("deg_approx") = 5);

    m.def(
        "loft_pieces",
        [](const std::vector<std::vector<gbs::BSCurve<T, dim>>> &pieces, const std::vector<T> &v, size_t q, T tol)
        {
            return gbs::loft_pieces(pieces, v, q, tol);
        },
        "Skinned surfaces through contours cut in pieces: pieces[k][i] the piece k of the contour i, each lofted on its own "
        "at the parameters v. Consecutive pieces sharing their end poles (within tol) give surfaces sharing those boundaries exactly.",
        py::arg("pieces"), py::arg("v"), py::arg("q") = 3, py::arg("tol") = gbs::brep_default_tolerance<T>);

    m.def(
        "loft_pieces",
        [](const std::vector<std::vector<gbs::BSCurveRational<T, dim>>> &pieces, const std::vector<T> &v, size_t q, T tol)
        {
            return gbs::loft_pieces(pieces, v, q, tol);
        },
        "Skinned rational surfaces through contours cut in rational pieces (end poles and weights shared).",
        py::arg("pieces"), py::arg("v"), py::arg("q") = 3, py::arg("tol") = gbs::brep_default_tolerance<T>);

    m.def(
        "gordonbs",
        [](const std::list<gbs::BSCurve<T, dim>> &u_crv_lst, const std::list<gbs::BSCurve<T, dim>> &v_crv_lst, T tol)
        {
            return gbs::gordon<T, dim>(u_crv_lst.begin(), u_crv_lst.end(), v_crv_lst.begin(), v_crv_lst.end(), tol);
        },
        "Bidirectional interpolation using Gordon's algorythm.\nCurves must be compatibles, aka curve lattice has to intersect at the same parameters values",
        py::arg("u_crv_lst"), py::arg("v_crv_lst"), py::arg("tol") = 1e-6);
}
