#pragma once
#include <stdexcept>
#include <array>
namespace gbs{
    template<typename T>
    class OutOfBoundsEval : public std::domain_error
    {
    public:
        explicit OutOfBoundsEval(T v, const std::array<T, 2> &bounds, const char *eval_msg ="Eval ") : 
            std::domain_error{ "Eval "+  std::to_string(v) + " out of bounds [ " +  std::to_string(bounds[0]) + " , "  + std::to_string(bounds[1]) + " ] error."}
            { }

    };
}