#pragma once
#include <gbs/bscurve.h>
#include <gbs/vecop.h>

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

    // template <typename T, size_t dim, bool rational1, bool rational2>
    // auto join(  const gbs::BSCurveGeneral<T, dim, rational1> &crv1, 
    //             const gbs::BSCurveGeneral<T, dim, rational2> &crv2) 
    //             -> gbs::BSCurveGeneral<T, dim, rational1 || rational2>

    template <typename T,size_t dim>
     auto join(const BSCurve<T, dim> &crv1, const BSCurve<T, dim> &crv2)
    {
        auto poles1(crv1.poles());
        auto poles2(crv2.poles());
        // TODO transform to rational if required
        auto k1(crv1.knotsFlats());
        auto k2(crv2.knotsFlats());

        auto p1 = crv1.degree();
        auto p2 = crv2.degree();
        while (p1<p2)
        {
            increase_degree(k1,poles1,p1);
            p1++;
        }
        while (p2<p1)
        {
            increase_degree(k2,poles2,p2);
            p2++;
        }

        poles1.back() = 0.5 * (poles1.back()+poles2.front());
        poles1.insert(poles1.end(), std::next(poles2.begin()), poles2.end());

        auto k1_end = k1.back();
        auto k2_start = k2.front();
        std::transform(
            k2.begin(),
            k2.end(),
            k2.begin(),
            [&k2_start,&k1_end](const auto &k_)
            {
                return k_ + k1_end - k2_start;
            }
        );
        // k1.pop_back();
        k1.erase( k1.end() - p1, k1.end() );
        k1.insert(k1.end() , std::next(k2.begin(),p2), k2.end());

        // if(rational1 || rational2)
        //     return gbs::BSCurveRational<T,dim>(poles1,k1,p1);
        // else
            return 
            BSCurve<T,dim>(poles1,k1,p1);
    }

} // namespace gbs