#pragma once
#include "gbslib.h"

namespace gbs
{
/**
 * @brief Check if curves insice a container have the same bounds
 * 
 * @param p_crv_begin 
 * @param p_crv_end 
 * @return true 
 * @return false 
 */
    template<typename T, size_t dim>
    auto check_curves_bounds(const auto &crv_begin, const auto &crv_end) -> bool
    {
        auto [u1, u2] = crv_begin->bounds();
        return std::all_of(
            // std::execution::par,
            std::next(crv_begin,1),crv_end,
            [u1, u2](const auto &crv)
            {
                auto [u1_, u2_] = crv.bounds();
                return ( (abs(u1_-u1) < knot_eps<T>) && (abs(u2_-u2) < knot_eps<T>) );
            }
        );
    }
/**
 * @brief Check if curves inside a pointers container have the same bounds
 * 
 * @param p_crv_begin 
 * @param p_crv_end 
 * @return true 
 * @return false 
 */
    template<typename T, size_t dim>
    auto check_p_curves_bounds(const auto &p_crv_begin, const auto &p_crv_end) -> bool
    {
        auto [u1, u2] = (*p_crv_begin)->bounds();
        return std::all_of(
            // std::execution::par,
            std::next(p_crv_begin,1),p_crv_end,
            [u1, u2](const auto &p_crv)
            {
                auto [u1_, u2_] = p_crv->bounds();
                return ( (abs(u1_-u1) < knot_eps<T>) && (abs(u2_-u2) < knot_eps<T>) );
            }
        );
    }
}