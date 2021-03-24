#pragma once

#include <gbs/gbslib.h>
#include <gbs/bscurve.h>
#include <version>
#include <iostream>
#include <algorithm>
#ifndef __cpp_lib_format
#include <fmt/core.h>
using fmt::format;
#else
#include <format>
using std::format;
#endif
namespace gbs
{

    const char * format_Float = "{:*< 8.6}";
    const char * format_uInt = " {:>7}";
    const char * format_Int = " {:> 7}";

    auto print(double x)
    {
        std::cout << format(format_Float, x);
    }
    auto print(float x)
    {
        std::cout << format(format_Float, x);
    }
    auto print(size_t i)
    {
        std::cout << format(format_uInt, i);
    }
    template <typename T, size_t dim>
    auto print(const point<T, dim> &pt) -> void
    {
        std::cout << "[";
        std::for_each(
            pt.begin(),
            std::next(pt.end(), -1),
            [](const auto x_) {
                std::cout << " ";
                print(x_);
                std::cout  << " |";
            });
        std::cout << " ";
        print(pt.back());
        std::cout << " ]";
    }

    template <typename T, typename L>
    auto print(const std::pair<T,L> p) ->void
    {
        std::cout << "{ ";
        print(p.first);
        std::cout  << " | ";
        print(p.second);
        std::cout << " }";
    }


    template <typename T, size_t dim>
    auto print(const std::vector<T> &v) -> void
    {
        std::cout << "[" << std::endl;
        for (const auto &v_ : v)
        {
            std::cout << "\t";
            print(v_);
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    template <typename T, typename L>
    auto print(const std::vector<std::pair<T,L>> &v) -> void
    {
        std::cout << "[" << std::endl;
        for (const auto &v_ : v)
        {
            std::cout << "\t";
            print(v_);
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    template <typename T, size_t dim>
    auto print(const points_vector<T, dim> &pts) -> void
    {
        std::cout << "[" << std::endl;
        for (const auto &pt : pts)
        {
            std::cout << "\t";
            print(pt);
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    template <typename T, size_t dim>
    auto print(const BSCurve<T, dim> &crv) -> void
    {
        std::cout << "Degree: ";
        print( crv.degree() );
        std::cout << std::endl;
        std::cout << "Poles:" << std::endl;
        print(crv.poles());
        std::cout << "Knots anf mult:" << std::endl;
        auto k_m = unflat_knots(crv.knotsFlats());
        print(k_m); 
    }

} // namespace gbs