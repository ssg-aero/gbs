#pragma once
#include <stdexcept>
#include <array>
namespace gbs{
    template<typename T>
    class OutOfBoundsEval : public std::out_of_range
    {
        T m_val;
        std::array<T,2> m_bounds;

    public:
        explicit OutOfBoundsEval(T v, const std::array<T, 2> &bounds, const char *eval_msg ="Eval ") : 
            m_val{v}, 
            m_bounds{bounds}, 
            std::out_of_range{ eval_msg + msg()}
        { }
        auto val() const noexcept {return m_val;}
        const auto& bounds() const noexcept {return m_bounds;}
        auto msg()
        {
            return std::to_string(val()) + " out of bounds [ " +  std::to_string(bounds()[0]) + " , "  + std::to_string(bounds()[1]) + " ] error.";
        }
    };
}