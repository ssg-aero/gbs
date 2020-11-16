#pragma once
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
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

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto join(const gbs::BSCurveGeneral<T, dim, rational1> *crv1,
              const gbs::BSCurveGeneral<T, dim, rational2> *crv2)// -> std::unique_ptr<gbs::BSCurveGeneral<T, dim, rational1 || rational2>>
    {
        // make curves copies uniform in definition
        typedef std::conditional<rational1 || rational2,gbs::BSCurveRational<T,dim>,gbs::BSCurve<T,dim>>::type crvType;
        std::unique_ptr< gbs::BSCurveGeneral<T, dim, rational1 || rational2> > crv1_cp,crv2_cp;
        crv1_cp = std::make_unique<crvType>(*crv1); 
        crv2_cp = std::make_unique<crvType>(*crv2);
        // recover data of same dimension
        auto poles1 = crv1_cp->poles();
        auto poles2 = crv2_cp->poles();
        auto k1(crv1_cp->knotsFlats());
        auto k2(crv2_cp->knotsFlats());
        auto p1 = crv1_cp->degree();
        auto p2 = crv2_cp->degree();      
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
        // join poles
        if(rational1 || rational2)
        {
            //set tail/head with a weight of 1
            scale_poles(poles1,1. / poles1.back().back());
            scale_poles(poles2,1. / poles2.front().back());
            auto pt = T(0.5) * (weight_projection(poles1.back())+weight_projection(poles2.front()));
            // poles1.back() = add_weight(pt,poles1.back().back());
        }
        else
        {
            poles1.back() = T(0.5) * (poles1.back() + poles2.front());
        }
        poles1.insert(poles1.end(), std::next(poles2.begin()), poles2.end());
        // join knots
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
        // create result
        return std::make_unique<crvType>(poles1,k1,p1);
    }

    template <typename T, size_t dim, bool rational1, bool rational2>
    auto c2_connect(const BSCurveGeneral<T, dim, rational1> &crv1,
                    const BSCurveGeneral<T, dim, rational2> &crv2) -> BSCurve<T,dim>
    {
    std::vector<gbs::constrType<T,dim,3>> Q =
    {
        {crv1.value(1.),crv1.value(1.,1),crv1.value(1.,2)},
        {crv2.value(0.),crv2.value(0.,1),crv2.value(0.,2)}
    };

    return interpolate(Q,gbs::KnotsCalcMode::CHORD_LENGTH);
    }

} // namespace gbs