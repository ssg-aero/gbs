#pragma once
#include <array>
#include <execution>
#include <fstream>
#include <gbs/bscinterp.h>
#include <gbs/curves>
#include <gbs/io/tojson.h>
#include <gbs/surfaces>
#include <iostream>
#include <rapidjson/document.h>
#include <rapidjson/rapidjson.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <gbs/third_parties/magic_enum.hpp>

namespace gbs {
static const char *kTypeNames[] = {"Null",  "False",  "True",  "Object",
                                   "Array", "String", "Number"};

auto print_js_obj_content(const auto &js_obj) {
  for (auto &m : js_obj)
    printf("Type of member %s is %s\n", m.name.GetString(),
           kTypeNames[m.value.GetType()]);
}
template <typename T> auto get_val(const rapidjson::Value &val) -> T {
  throw std::invalid_argument(
      "auto get_val(const rapidjson::Value &val) -> T unsupported type");
  return T{};
}

template <> inline auto get_val<double>(const rapidjson::Value &val) -> double {
  if (!val.IsDouble())
    throw std::invalid_argument(
        "auto get_val(const rapidjson::Value &val) wrong type" +
        std::to_string(val.GetType()));
  return val.GetDouble();
}

template <> inline auto get_val<float>(const rapidjson::Value &val) -> float {
  if (!val.IsFloat())
    throw std::invalid_argument(
        "auto get_val(const rapidjson::Value &val) wrong type");
  return val.GetFloat();
}

template <> inline auto get_val<int>(const rapidjson::Value &val) -> int {
  if (!val.IsInt())
    throw std::invalid_argument(
        "auto get_val(const rapidjson::Value &val) wrong type");
  return val.GetInt();
}

template <> inline auto get_val<size_t>(const rapidjson::Value &val) -> size_t {
  if (!val.IsUint64())
    throw std::invalid_argument(
        "auto get_val(const rapidjson::Value &val) wrong type");
  return val.GetUint64();
}

template <>
inline auto get_val<std::string>(const rapidjson::Value &val) -> std::string {
  if (!val.IsString())
    throw std::invalid_argument(
        "auto get_val(const rapidjson::Value &val) wrong type");
  return std::string{val.GetString()};
}

template <typename T>
auto make_vec(const rapidjson::Value &a) -> std::vector<T> {
  if (!a.IsArray()) {
    throw std::invalid_argument("auto make_vec(const rapidjson::Value &a) -> "
                                "std::vector<T> not and array");
  }
  std::vector<T> v_(a.Size());
  std::transform(
      // std::execution::par,
      a.Begin(), a.End(), v_.begin(),
      [](const auto &val) { return get_val<T>(val); });
  // std::vector<T> v_;
  // for(const auto &val : a)
  // {
  //     v_.push_back(get_val<T>(val));
  // }
  return v_;
}

template <typename T, size_t dim>
auto make_array(const rapidjson::Value &a) -> std::array<T, dim> {
  if (!a.IsArray()) {
    throw std::invalid_argument("auto make_array(const rapidjson::Value &a) -> "
                                "std::array<T,dim> not and array");
  }
  if (a.Size() != dim) {
    throw std::invalid_argument("auto make_array(const rapidjson::Value &a) -> "
                                "std::array<T,dim> wrong size");
  }
  std::array<T, dim> v_;
  std::transform(
      // std::execution::par,
      a.Begin(), a.End(), v_.begin(),
      [](const auto &val) { return get_val<T>(val); });
  return v_;
}

template <typename T, size_t dim>
auto make_ax2(const rapidjson::Value &a) -> ax2<T, dim> {
  if (!a.IsArray()) {
    throw std::invalid_argument("auto make_array(const rapidjson::Value &a) -> "
                                "std::array<T,dim> not and array");
  }
  if (a.Size() != 3) {
    throw std::invalid_argument("auto make_array(const rapidjson::Value &a) -> "
                                "std::array<T,dim> wrong size");
  }
  auto arr = a.GetArray();
  return {
      make_array<T, dim>(arr[0]),
      make_array<T, dim>(arr[1]),
      make_array<T, dim>(arr[2]),
  };
}

template <typename T, size_t dim>
auto make_point_vec(const rapidjson::Value &a)
    -> std::vector<std::array<T, dim>> {
  if (!a.IsArray()) {
    throw std::invalid_argument(
        "auto make_point_vec(const rapidjson::Value &a) -> std::vector< "
        "std::array<T,dim> > not and array");
  }
  std::vector<std::array<T, dim>> v_(a.Size());
  std::transform(
      // std::execution::par,
      a.Begin(), a.End(), v_.begin(),
      [](const auto &val) { return make_array<T, dim>(val); });
  return v_;
}

template <typename T, size_t dim, size_t nc>
auto make_constrains_vec(const rapidjson::Value &a)
    -> std::vector<gbs::constrType<T, dim, nc>> {
  if (!a.IsArray()) {
    throw std::invalid_argument(
        "auto make_point_vec(const rapidjson::Value &a) -> std::vector< "
        "std::array<T,dim> > not and array");
  }
  auto n = a.GetArray()[0].GetArray().Size();

  std::vector<gbs::constrType<T, dim, nc>> v_(n);
  for (auto j = 0; j < n; j++) {
    for (auto i = 0; i < nc; i++) {
      v_[j][i] = make_array<T, dim>(a.GetArray()[i].GetArray()[j]);
      // std::transform(
      //     std::execution::par,
      //     a.GetArray()[i].GetArray().Begin(),
      //     a.GetArray()[i].GetArray().End(),
      //     v_.begin(),
      //     v_.begin(),
      //     [](const auto &val, const auto &constain) { return make_array<T,
      //     dim>(val) + constain; });
    }
  }
  return v_;
}

template <typename T, size_t dim>
auto bscurve_direct(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> {
  // std::cerr << "Curve name: " << a["name"].GetString() << std::endl;
  assert(std::strcmp(a["type"].GetString(), "bscurve") == 0);
  auto knots = make_vec<T>(a["knots"]);
  auto deg = get_val<size_t>(a["deg"]);
  auto poles = make_point_vec<T, dim>(a["poles"]);

  return gbs::BSCurve<T, dim>(poles, knots, deg);
}

template <typename T, size_t dim>
auto bscurve_interp_cn(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> {
  assert(std::strcmp(a["type"].GetString(), "bscurve_interp_cn") == 0);
  if (a.HasMember("params")) {
    auto u = make_vec<T>(a["params"]);
    auto deg = get_val<size_t>(a["deg"]);
    auto points = make_point_vec<T, dim>(a["points"]);
    return gbs::interpolate(points, u, deg);
  } else if (a.HasMember("param_calc_mode")) {
    auto deg = get_val<size_t>(a["deg"]);
    auto points = make_point_vec<T, dim>(a["points"]);
    auto enum_str = std::string(a["param_calc_mode"].GetString());
    auto mode = magic_enum::enum_cast<gbs::KnotsCalcMode>(enum_str);
    return gbs::interpolate(points, deg, mode.value());
  } else {
    auto deg = get_val<size_t>(a["deg"]);
    auto points = make_point_vec<T, dim>(a["points"]);
    return gbs::interpolate(points, deg, gbs::KnotsCalcMode::CHORD_LENGTH);
  }
}

template <typename T, size_t dim>
auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> {
  auto mode = gbs::KnotsCalcMode::CHORD_LENGTH;
  std::vector<T> u;
  if (a.HasMember("param_calc_mode") && a["param_calc_mode"].IsString()) {
    if (strcmp(a["param_calc_mode"].GetString(), "SPECIFIED") == 0) {
      u = gbs::make_vec<T>(a["params"]);
    } else {
      auto enum_str = std::string(a["param_calc_mode"].GetString());
      mode = magic_enum::enum_cast<gbs::KnotsCalcMode>(enum_str).value();
    }
  }
  if (!a.HasMember("constrains"))
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (!a["constrains"].IsArray())
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  auto nc = a["constrains"].GetArray().Size();
  if (nc <= 0)
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (!a["constrains"].GetArray()[0].IsArray())
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (nc == 1) {
    auto Q = make_constrains_vec<T, dim, 1>(a["constrains"]);
    if (u.size())
      return gbs::interpolate<T, dim>(Q, u);
    else
      return gbs::interpolate<T, dim>(Q, mode);
  } else if (nc == 2) {
    auto Q = make_constrains_vec<T, dim, 2>(a["constrains"]);
    if (u.size())
      return gbs::interpolate<T, dim>(Q, u);
    else
      return gbs::interpolate<T, dim>(Q, mode);
  } else {
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> number of constrains not implemented");
    return gbs::BSCurve<T, dim>{}; // tur off warning
  }
}

template <typename T, size_t dim>
auto make_bscurve(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> {
  if (std::strcmp(a["type"].GetString(), "bscurve_interp_cn") == 0) {
    return bscurve_interp_cn<T, dim>(a);
  } else if (std::strcmp(a["type"].GetString(), "bscurve") == 0) {
    return bscurve_direct<T, dim>(a);
  } else if (std::strcmp(a["type"].GetString(), "bscurve_interp") == 0) {

    return bscurve_interp<T, dim>(a);
  } else {
    throw std::invalid_argument(
        "auto make_curve(const rapidjson::Value &a) unsupported type");
    return gbs::BSCurve<T, dim>{}; // tur off warning
  }
}

template <typename T, size_t dim>
auto bscurve_interp(const rapidjson::Value &a, const std::vector<T> &u)
    -> gbs::BSCurve<T, dim> {
  if (!a.HasMember("constrains"))
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (!a["constrains"].IsArray())
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  auto nc = a["constrains"].GetArray().Size();
  if (nc <= 0)
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (!a["constrains"].GetArray()[0].IsArray())
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> constrains not defined");
  if (nc == 1) {
    auto Q = make_constrains_vec<T, dim, 1>(a["constrains"]);
    return gbs::interpolate<T, dim>(Q, u);
  } else if (nc == 2) {
    auto Q = make_constrains_vec<T, dim, 2>(a["constrains"]);
    return gbs::interpolate<T, dim>(Q, u);
  } else {
    throw std::invalid_argument(
        "auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, "
        "dim> number of constrains not implemented");
    return gbs::BSCurve<T, dim>{}; // tur off warning
  }
}

template <typename T, size_t dim>
auto make_bscurve(const rapidjson::Value &a, const std::vector<T> &u)
    -> gbs::BSCurve<T, dim> {
  if (std::strcmp(a["type"].GetString(), "bscurve_interp") == 0) {
    return bscurve_interp<T, dim>(a, u);
  } else {
    throw std::invalid_argument(
        "auto make_curve(const rapidjson::Value &a) unsupported type");
    return gbs::BSCurve<T, dim>{}; // tur off warning
  }
}

template <typename T, size_t dim>
auto make_surface(const rapidjson::Value &val)
    -> std::shared_ptr<Surface<T, dim>>;

template <typename T, size_t dim>
auto make_curve(const rapidjson::Value &val) -> std::shared_ptr<Curve<T, dim>> {
  switch (val["type"].GetInt()) {
  case static_cast<int>(entity_type::BSCurve):
    return std::make_shared<BSCurve<T, dim>>(BSCurve<T, dim>{
        make_point_vec<T, dim>(val["poles"]),
        flat_knots(make_vec<T>(val["knots"]), make_vec<size_t>(val["mults"])),
        val["degree"].GetUint64()});
    break;
  case static_cast<int>(entity_type::BSCurveRational):
    return std::make_shared<BSCurveRational<T, dim>>(BSCurveRational<T, dim>{
        make_point_vec<T, dim>(val["poles"]),
        flat_knots(make_vec<T>(val["knots"]), make_vec<size_t>(val["mults"])),
        make_vec<T>(val["weights"]), val["degree"].GetUint64()});
    break;
  case static_cast<int>(entity_type::CurveOnSurface):
    return std::make_shared<CurveOnSurface<T, dim>>(
        CurveOnSurface<T, dim>{make_curve<T, 2>(val["curve2d"]),
                               make_surface<T, dim>(val["surface"])});
    break;
  default:
    throw std::invalid_argument("Unsupported curve type");
    break;
  }
}

template <typename T, size_t dim>
auto make_surface(const rapidjson::Value &val)
    -> std::shared_ptr<Surface<T, dim>> {
  switch (val["type"].GetInt()) {
  case static_cast<int>(entity_type::BSSurface):
    return std::make_shared<BSSurface<T, dim>>(BSSurface<T, dim>{
        make_point_vec<T, dim>(val["poles"]),
        flat_knots(make_vec<T>(val["knotsU"]), make_vec<size_t>(val["multsU"])),
        flat_knots(make_vec<T>(val["knotsV"]), make_vec<size_t>(val["multsV"])),
        val["degreeU"].GetUint64(), val["degreeV"].GetUint64()});
    break;
  case static_cast<int>(entity_type::BSSurfaceRational):
    return std::make_shared<BSSurfaceRational<T, dim>>(
        BSSurfaceRational<T, dim>{
            add_weights_coord(make_point_vec<T, dim>(val["poles"]),
                              make_vec<T>(val["weights"])),
            flat_knots(make_vec<T>(val["knotsU"]),
                       make_vec<size_t>(val["multsU"])),
            flat_knots(make_vec<T>(val["knotsV"]),
                       make_vec<size_t>(val["multsV"])),
            val["degreeU"].GetUint64(), val["degreeV"].GetUint64()});
    break;
  default:
    throw std::invalid_argument("Unsupported surface type");
    break;
  }
}

template <typename T>
auto make_surface(const rapidjson::Value &val)
    -> std::shared_ptr<Surface<T, 3>> {
  switch (val["type"].GetInt()) {
  case static_cast<int>(entity_type::SurfaceOfRevolution): {
    return std::make_shared<SurfaceOfRevolution<T>>(SurfaceOfRevolution<T>{
        make_curve<T, 2>(val["curve"]), make_ax2<T, 3>(val["axis2"]),
        get_val<T>(val["theta1"]), get_val<T>(val["theta2"])});
    break;
  }
  default:
    return make_surface<T, 3>(val);
    break;
  }
}

inline auto parse_file(const char *fname, rapidjson::Document &document) {
  std::ifstream f(fname);
  std::string str;
  if (f) {
    std::ostringstream ss;
    ss << f.rdbuf();
    str = ss.str();
  }

  document.Parse(str.data());
}
} // namespace gbs