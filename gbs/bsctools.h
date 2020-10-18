#pragma once
#include <gbs/bscurve.h>

namespace gbs

{
    /**
     * @brief Trim curve's head
     * 
     * @tparam T 
     * @tparam dim 
     * @param p 
     * @param knots_flats 
     * @param poles 
     * @param u 
     */
    template <typename T, size_t dim>
    auto trim_begin(size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles, T u) -> void
    {
        if(u-knots_flats.front()<knot_eps) return;
        
        for (auto i = 0; i < p; i++)
        {
            insert_knot(u, p, knots_flats, poles);
        }

        auto it_l = std::lower_bound(knots_flats.begin(), knots_flats.end(), u);
        auto i_l = it_l - knots_flats.begin();
        knots_flats.erase(knots_flats.begin(), it_l);
        knots_flats.insert(knots_flats.begin(), knots_flats.front());
        poles.erase(poles.begin(), std::next(poles.begin(), i_l - 1));

    }
    /**
     * @brief Trim curve's tail
     * 
     * @tparam T 
     * @tparam dim 
     * @param p 
     * @param knots_flats 
     * @param poles 
     * @param u 
     */
    template <typename T, size_t dim>
    auto trim_end(size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles, T u) -> void
    {
        if(knots_flats.back()-u<knot_eps) return;
        
        for (auto i = 0; i < p; i++)
        {
            insert_knot(u, p, knots_flats, poles);
        }

        auto it_h = std::lower_bound(knots_flats.begin(), knots_flats.end(), u);
        auto i_h = it_h - knots_flats.begin();
        knots_flats.erase(std::next(it_h, p), knots_flats.end());
        knots_flats.push_back(knots_flats.back());
        poles.erase(std::next(poles.begin(), i_h), poles.end());

    }
    /**
     * @brief Trim between u1 and u2
     * 
     * @tparam T 
     * @tparam dim 
     * @param p 
     * @param knots_flats 
     * @param poles 
     * @param u1 
     * @param u2 
     */
    template <typename T, size_t dim>
    auto trim(size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles, T u1, T u2) -> void
    {
        trim_begin(p,knots_flats,poles,fmin(u1,u2));
        trim_end(p,knots_flats,poles,fmax(u1,u2));
    }
} // namespace occt_utils