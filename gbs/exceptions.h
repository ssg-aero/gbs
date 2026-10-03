#pragma once
#include <stdexcept>
#include <array>
#include <string>

namespace gbs{
    template<typename T>
    class OutOfBoundsEval : public std::domain_error
    {
    public:
        explicit OutOfBoundsEval(T v, const std::array<T, 2> &bounds, const char *eval_msg ="Eval ") :
            std::domain_error{ std::string{eval_msg} +  std::to_string(v) + " out of bounds [ " +  std::to_string(bounds[0]) + " , "  + std::to_string(bounds[1]) + " ] error."}
            { }

    };

    /**
     * @brief Error raised by the native BREP core (gbs-brep): invalid handle,
     * violated precondition of a builder, inconsistent topology.
     */
    class BRepError : public std::runtime_error
    {
    public:
        explicit BRepError(const std::string &msg) : std::runtime_error{"brep: " + msg} {}
    };
}
