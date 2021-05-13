#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <rapidjson/rapidjson.h>
#include <rapidjson/document.h>
#include <execution>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <array>
#include <tools/magic_enum.hpp>
#include <stdexcept>

namespace gbs
{
    static const char* kTypeNames[] = 
        { "Null", "False", "True", "Object", "Array", "String", "Number" };

    auto print_js_obj_content(const auto &js_obj)
    {
        for (auto &m : js_obj)
            printf("Type of member %s is %s\n",
                   m.name.GetString(), kTypeNames[m.value.GetType()]);
    }
    template <typename T>
    auto get_val(const rapidjson::Value &val) -> T
    {
        throw std::exception("auto get_val(const rapidjson::Value &val) -> T unsupported type");
        return T{};
    }

    template <>
    inline auto get_val<double>(const rapidjson::Value &val) -> double
    {
        if (!val.IsDouble())
            throw std::exception("auto get_val(const rapidjson::Value &val) wrong type");
        return val.GetDouble();
    }

    template <>
    inline auto get_val<float>(const rapidjson::Value &val) -> float
    {
        if (!val.IsFloat())
            throw std::exception("auto get_val(const rapidjson::Value &val) wrong type");
        return val.GetFloat();
    }

    template <>
    inline auto get_val<int>(const rapidjson::Value &val) -> int
    {
        if (!val.IsInt())
            throw std::exception("auto get_val(const rapidjson::Value &val) wrong type");
        return val.GetInt();
    }

    template <>
    inline auto get_val<size_t>(const rapidjson::Value &val) -> size_t
    {
        if (!val.IsUint64())
            throw std::exception("auto get_val(const rapidjson::Value &val) wrong type");
        return val.GetUint64();
    }

    template <>
    inline auto get_val<std::string>(const rapidjson::Value &val) -> std::string
    {
        if (!val.IsString())
            throw std::exception("auto get_val(const rapidjson::Value &val) wrong type");
        return std::string{val.GetString()};
    }

    template <typename T>
    auto make_vec(const rapidjson::Value &a) -> std::vector<T>
    {
        if (!a.IsArray())
        {
            throw std::exception("auto make_vec(const rapidjson::Value &a) -> std::vector<T> not and array");
        }
        std::vector<T> v_(a.Size());
        std::transform(
            std::execution::par,
            a.Begin(),
            a.End(),
            v_.begin(),
            [](const auto &val) { return get_val<T>(val); });
        return v_;
    }

    template <typename T, size_t dim>
    auto make_array(const rapidjson::Value &a) -> std::array<T, dim>
    {
        if (!a.IsArray())
        {
            throw std::exception("auto make_array(const rapidjson::Value &a) -> std::array<T,dim> not and array");
        }
        if (a.Size() != dim)
        {
            throw std::exception("auto make_array(const rapidjson::Value &a) -> std::array<T,dim> wrong size");
        }
        std::array<T, dim> v_;
        std::transform(
            std::execution::par,
            a.Begin(),
            a.End(),
            v_.begin(),
            [](const auto &val) { return get_val<T>(val); });
        return v_;
    }

    template <typename T, size_t dim>
    auto make_point_vec(const rapidjson::Value &a) -> std::vector<std::array<T, dim>>
    {
        if (!a.IsArray())
        {
            throw std::exception("auto make_point_vec(const rapidjson::Value &a) -> std::vector< std::array<T,dim> > not and array");
        }
        std::vector<std::array<T, dim>> v_(a.Size());
        std::transform(
            std::execution::par,
            a.Begin(),
            a.End(),
            v_.begin(),
            [](const auto &val) { return make_array<T, dim>(val); });
        return v_;
    }

    template <typename T, size_t dim, size_t nc>
    auto make_constrains_vec(const rapidjson::Value &a) -> std::vector<gbs::constrType<T, dim, nc>>
    {
        if (!a.IsArray())
        {
            throw std::exception("auto make_point_vec(const rapidjson::Value &a) -> std::vector< std::array<T,dim> > not and array");
        }
        auto n = a.GetArray()[0].GetArray().Size();

        std::vector<gbs::constrType<T, dim, nc>> v_(n);
        for (auto j = 0; j < n; j++)
        {
            for (auto i = 0; i < nc; i++)
            {
                v_[j][i] =  make_array<T, dim>( a.GetArray()[i].GetArray()[j] );
                // std::transform(
                //     std::execution::par,
                //     a.GetArray()[i].GetArray().Begin(),
                //     a.GetArray()[i].GetArray().End(),
                //     v_.begin(),
                //     v_.begin(),
                //     [](const auto &val, const auto &constain) { return make_array<T, dim>(val) + constain; });
            }
        }
        return v_;
    }

    template <typename T, size_t dim>
    auto bscurve_direct(const rapidjson::Value &a) -> gbs::BSCurve<T, dim>
    {
        // std::cerr << "Curve name: " << a["name"].GetString() << std::endl;
        assert(std::strcmp(a["type"].GetString(), "bscurve") == 0);
        auto knots = make_vec<T>(a["knots"]);
        auto deg = get_val<size_t>(a["deg"]);
        auto poles = make_point_vec<T, dim>(a["poles"]);

        return gbs::BSCurve<T, dim>(poles, knots, deg);
    }

    template <typename T, size_t dim>
    auto bscurve_interp_cn(const rapidjson::Value &a) -> gbs::BSCurve<T, dim>
    {
        assert(std::strcmp(a["type"].GetString(), "bscurve_interp_cn") == 0);
        if (a.HasMember("params"))
        {
            auto u = make_vec<T>(a["params"]);
            auto deg = get_val<size_t>(a["deg"]);
            auto points = make_point_vec<T, dim>(a["points"]);
            return gbs::interpolate(points, u, deg);
        }
        else if (a.HasMember("param_calc_mode"))
        {
            auto deg = get_val<size_t>(a["deg"]);
            auto points = make_point_vec<T, dim>(a["points"]);
            auto enum_str = std::string(a["param_calc_mode"].GetString());
            auto mode = magic_enum::enum_cast<gbs::KnotsCalcMode>(enum_str);
            return gbs::interpolate(points, deg, mode.value());
        }
        else
        {
            auto deg = get_val<size_t>(a["deg"]);
            auto points = make_point_vec<T, dim>(a["points"]);
            return gbs::interpolate(points, deg, gbs::KnotsCalcMode::CHORD_LENGTH);
        }
    }

    template <typename T, size_t dim>
    auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim>
    {
        auto mode = gbs::KnotsCalcMode::CHORD_LENGTH;
        if (a.HasMember("param_calc_mode") && a["param_calc_mode"].IsString())
        {
            auto enum_str = std::string(a["param_calc_mode"].GetString());
            mode = magic_enum::enum_cast<gbs::KnotsCalcMode>(enum_str).value();
        }
        if (!a.HasMember("constrains"))
            throw std::exception("auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> constrains not defined");
        if (!a["constrains"].IsArray())
            throw std::exception("auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> constrains not defined");
        auto nc = a["constrains"].GetArray().Size();
        if (nc <= 0)
            throw std::exception("auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> constrains not defined");
        if (!a["constrains"].GetArray()[0].IsArray())
            throw std::exception("auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> constrains not defined");
        if (nc == 1)
        {
            auto Q = make_constrains_vec<T, dim, 1>(a["constrains"]);
            return gbs::interpolate<T, dim>(Q, mode);
        }
        else if (nc == 2)
        {
            auto Q = make_constrains_vec<T, dim, 2>(a["constrains"]);
            return gbs::interpolate<T, dim>(Q, mode);
        }
        else
        {
            throw std::exception("auto bscurve_interp(const rapidjson::Value &a) -> gbs::BSCurve<T, dim> number of constrains not implemented");
            return gbs::BSCurve<T, dim>{}; // tur off warning
        }
    }

    template <typename T, size_t dim>
    auto make_bscurve(const rapidjson::Value &a) -> gbs::BSCurve<T, dim>
    {
        if( std::strcmp(a["type"].GetString(),"bscurve_interp_cn") == 0)
        {
            return bscurve_interp_cn<T,dim>(a);
        }
        else if( std::strcmp(a["type"].GetString(),"bscurve")  == 0)
        {
            return bscurve_direct<T,dim>(a);
        }
        else if( std::strcmp(a["type"].GetString(),"bscurve_interp")  == 0)
        {

            return bscurve_interp<T,dim>(a);
        }
        else
        {
            throw std::exception("auto make_curve(const rapidjson::Value &a) unsupported type");
            return gbs::BSCurve<T, dim>{}; // tur off warning
        }
    }


    inline auto parse_file(const char *fname, rapidjson::Document &document)
    {
        std::ifstream f(fname);
        std::string str;
        if (f)
        {
            std::ostringstream ss;
            ss << f.rdbuf();
            str = ss.str();
        }

        document.Parse(str.data());
    }
}