#pragma once
#include <gbslib/gbslib.h>
#include <gbslib/vecop.h>

namespace gbs
{
    /**
     * @brief Builds and array of knots' single values and store multiplicity
     * 
     * @tparam T 
     * @tparam L 
     * @param u 
     * @param m 
     * @param k 
     */
    template <class T, class L>
    inline void unflat_knots(const std::vector<T> &u, std::vector<L> &m, std::vector<T> &k)
    {
        m.clear();
        k.clear();
        m.push_back(1);
        k.push_back(u.front());

        // std::foreach(u.begin())

        auto n = u.size();
        for (size_t i = 1; i < n; i++)
        {
            if (fabs(u[i] - k.back()) < knot_eps)
            {
                m.back()++;
            }
            else
            {
                m.push_back(1);
                k.push_back(u[i]);
            }
        }
    }
    /**
     * @brief Build flat knots vector from knots values and multiplicity
     * 
     * @tparam T 
     * @tparam L 
     * @param k 
     * @param m 
     * @return std::vector<T> 
     */
    template <class T, class L>
    inline std::vector<T> flat_knots(const std::vector<T> &k, const std::vector<L> &m)
    {
        if(k.size()!=m.size())
        {
            throw std::length_error( "multiplicities must have same length as knots" );
        }
        // test_size(k, m);

        std::vector<T> flat_k;
        auto nk = k.size();
        for (size_t i = 0; i < nk; i++)
        {
            L n = m[i];
            T kv = k[i];
            for (L j = 0; j < n; j++)
            {
                flat_k.push_back(kv);
            }
        }
        flat_k.shrink_to_fit();
        return flat_k;
    }
    /**
     * @brief Eval flat knots multiplicities
     * 
     * @tparam T 
     * @param u 
     * @param k_flat 
     * @return auto 
     */
    template< typename T>
    auto multiplicity(const T & u,const std::vector<T> &k_flat)
    {
        return std::upper_bound(k_flat.begin(), k_flat.end(),u) - std::lower_bound(k_flat.begin(), k_flat.end(),u);
    }
    /**
     * @brief set the values ranging from 0. to 1.
     * 
     * @tparam T 
     * @param k 
     */
    template< typename T>
    auto adimension(std::vector<T> &k) -> void
    {
        auto d = k.back() - k.front();
        if(fabs(d)>knot_eps)
        {
            std::transform(k.begin(),k.end(),k.begin(),[&d](auto const &k_){return k_/d;});
        }
    }

    enum class KnotsCalcMode { EQUALY_SPACED , CHORD_LENGTH, CENTRIPETAL};
    /**
     * @brief Builds curve's parametrization from passing points, the result cal be set to range from 0. to 1.
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts 
     * @param mode 
     * @param adimensionnal 
     * @return std::vector<T> 
     */
    template< typename T, size_t dim> 
    auto curve_parametrization( const std::vector<std::array<T,dim> > &pts, KnotsCalcMode mode, bool adimensionnal =false ) -> std::vector<T>
    {

        std::vector<T> k(pts.size());
        k[0] = 0.;
        switch (mode)
        {
        case KnotsCalcMode::CHORD_LENGTH:
            std::transform(++pts.begin(), pts.end(),pts.begin(), ++k.begin(),
                           [&, k_ = 0.](const auto &pt1,const auto &pt2) mutable {
                               k_ += gbs::distance(pt1, pt2);
                               return k_;
                           });
            break;
        case KnotsCalcMode::CENTRIPETAL:
            std::transform(++pts.begin(), pts.end(),pts.begin(), ++k.begin(),
                           [&, k_ = 0.](const auto &pt1,const auto &pt2) mutable {
                               k_ += sqrt(gbs::distance(pt1, pt2));
                               return k_;
                           });
            break;
        default: // KnotsCalcMode::EQUALY_SPACED:
            auto step = 1. / (pts.end() - pts.begin() - 1); //TODO, make/use range func
            std::generate(++k.begin(),k.end(),
                           [&step,k_=0.]() mutable { return k_+=step; });
                break;
        }

        if(adimensionnal )
        {
            adimension(k);
        }

        return k;
    }
    /**
     * @brief Insert knot value
     * 
     * @tparam T 
     * @tparam dim 
     * @param u 
     * @param p 
     * @param knots_flats 
     * @param poles 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto insert_knot(T u, size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles)
    {

        //Ckeck if knot can be inserted
        std::vector<int> mult;
        std::vector<double> knots;
        gbs::unflat_knots(knots_flats, mult, knots);
        auto iu = std::find_if(knots.begin(),knots.end(),[&](const auto u_){return fabs(u_-u)<knot_eps;}) - knots.begin();

        if(iu<knots.size() && mult[iu]>=p) return;

        // Start inserting knot

        auto knots_flats_ = knots_flats; //copy/move for failproof
        std::vector<std::array<T, dim>> Q(poles.size() + 1);

        auto k = std::lower_bound(knots_flats_.begin(), knots_flats_.end(), u);
        k = knots_flats_.insert(k, 1, u);
        auto ik = (k - knots_flats_.begin()) - 1;

        Q.front() = poles.front();
        Q.back() = poles.back();
        T alpha;
        auto count = Q.size() - 1;
        for (auto i = 1; i < count; i++)
        {
            if (i <= ik - p)
                alpha = 1.;
            else if (i > ik)
                alpha = 0.;
            else
                alpha = (u - knots_flats[i]) / (knots_flats[i + p] - knots_flats[i]);
            Q[i] = alpha * poles[i] + (1 - alpha) * poles[i - 1];
        }

        //move
        knots_flats = std::move(knots_flats_);
        poles = std::move(Q);
    }
    /**
     * @brief Remove knot if possible
     * 
     * @tparam T 
     * @tparam dim 
     * @param u 
     * @param p 
     * @param knots_flats 
     * @param P 
     * @param tol 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto remove_knot(T u, size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &P, T tol)
    {
        // auto knots_flats_ = knots_flats; //copy/move for fail saife
        auto Pi_ = P;
        auto Pj_ = P;

        // std::vector<T> knots;
        // std::vector<size_t> m;
        // unflat_knots(knots_flats_,m,knots);

        // auto k = std::next(std::lower_bound(knots.begin(), knots.end(), u));
        // auto r = std::next(std::lower_bound(knots_flats_.begin(), knots_flats_.end(), *k)-1);

        // auto r = std::next(std::upper_bound(knots_flats_.begin(), knots_flats_.end(),u),-1)-knots_flats_.begin();
        // auto s = multiplicity(u,knots_flats_);

        auto k = std::next(std::upper_bound(knots_flats.begin(), knots_flats.end(),u),-1);
        auto r = k - knots_flats.begin();
        auto s = multiplicity(u,knots_flats);

        Pi_.erase(std::next(Pi_.begin(),(2*r-s-p)/2-1));// less mem and less copy possible
        Pj_.erase(std::next(Pj_.begin(),(2*r-s-p)/2-1));

        size_t i = r-p;
        size_t j = r-s;
        int t = 0;
        int d = j-i;
        while (d>t)
        {
            auto ai = (u-knots_flats[i]   )/(knots_flats[i+p+1+t]-knots_flats[i]);
            auto aj = (u-knots_flats[j-t])/(knots_flats[i+p+1]-knots_flats[j-t]);
            Pi_[i] = (P[i] - (1 - ai) * P[i - 1]) / ai;
            Pj_[j] = (P[j] - (1 - aj) * P[j + 1]) / aj;
            i++;
            j--;
            d=j-i;
        }

        auto D = delta(Pi_, Pj_);
        // if (std::reduce(
        //         std::execution::par,
        //         D.begin(), D.end(),
        //         T(0.),
        //         [](const auto &t_, const auto &v_) { return t_ + sq_norm(v_); }) < tol)
        // {
            knots_flats.erase(k);
            P = Pi_;
        // }
    }

} // namespace gbs