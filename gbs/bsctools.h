#pragma once
#include <gbs/bscurve.h>

namespace gbs

{
    /**
     * @brief Set all bspline curves to the same degree
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_degree(Container &bs_lst) -> void
    {
        auto &C = bs_lst.front();
        // Set first curve at maximun degree
        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                while (C.degree()< C_.degree())
                {
                    C.increaseDegree();
                } 
            });
        // Set the others at first curve's degree
        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                while (C_.degree()< C.degree())
                {
                    C_.increaseDegree();
                } 
            });
    }
    /**
     * @brief unify the knots of a bspline curve set
     * 
     * @tparam Container 
     * @param bs_lst 
     */
    template <typename Container>
    auto unify_knots(Container &bs_lst) -> void
    {
        auto p = bs_lst.front().degree();
        std::for_each(bs_lst.begin(),bs_lst.end(),[&p](const auto &C_)
        {
            if (C_.degree() != p)
                throw std::exception("unify_knots: need curves with same degree");
        });

        auto &C = bs_lst.front();

        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                C_.changeBounds(C.bounds());
                auto km_ = unflat_knots(C_.knotsFlats());
                C.insertKnots( km_ );
            });

        auto km = unflat_knots(C.knotsFlats());

        std::for_each(
            std::next(bs_lst.begin()),
            bs_lst.end(),
            [&](auto &C_) {
                C_.insertKnots( km );
            });
    }
} // namespace gbs