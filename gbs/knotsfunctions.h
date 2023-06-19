#pragma once
#include <gbs/gbslib.h>
#include <gbs/vecop.h>
#include <gbs/basisfunctions.h>
#include <gbs/maths.h>

#include <Eigen/Dense>

namespace gbs
{
    using namespace Eigen;

    template <typename T>
    using VectorX = Matrix<T, Dynamic, 1>;
    template <typename T>
    using MatrixX = Matrix<T, Dynamic, Dynamic>;

    template <typename T>
    auto build_simple_mult_flat_knots(const std::vector<T> &u, size_t p) -> std::vector<T>
    {

        auto n = u.size();
        auto nk = n + p + 1;
        auto u1 = u.front();
        auto u2 = u.back();

        std::vector<T> k_flat(nk);
        std::fill(k_flat.begin(), std::next(k_flat.begin(), p+1), u1);
        std::fill(std::next(k_flat.begin(), nk - p - 1), k_flat.end(), u2);

        auto delta_ = u2 - u1;
        if(u.size()>2)
        {
            for (int j = 1; j < n - p; j++) // TODO use std algo
            {
                k_flat[j + p] = 0;
                for (int i = j; i <= j + p - 1 ; i++)
                {
                    k_flat[j + p] += u[i] / p;
                }

                // k_flat[j + p] = j / T(n - p) * delta_ + u1;
            }
        }

        return k_flat;
    } 

    /**
     * @brief 
     * 
     * @tparam T 
     * @param u1 : Knot start value
     * @param u2 : Knot send value
     * @param n  : Poles number
     * @param p  : Degree
     * @return std::vector<T> 
     */

    template <typename T>
    auto build_simple_mult_flat_knots(T u1, T u2, size_t n, size_t p) -> std::vector<T>
    {

        auto nk = n + p + 1;
        
        std::vector<T> k_flat(nk);
        std::fill(k_flat.begin(), std::next(k_flat.begin(), p+1), u1);
        std::fill(std::next(k_flat.begin(), nk  - p - 1 ), k_flat.end(), u2);
        auto delta_ = u2 - u1;

        for (int j = 1; j < n - p; j++)
        {
            k_flat[j + p] = u1 + delta_ * j / T(n - p + 1);
        }

        return k_flat;
    } 
    /**
     * @brief Builds and array of knots' single values and store multiplicity
     * 
     * @tparam T 
     * @tparam L 
     * @param knots_flat 
     * @param mult 
     * @param knots 
     */
    template <class T, class L>
    inline void unflat_knots(const std::vector<T> &knots_flat, std::vector<L> &mult, std::vector<T> &knots)
    {
        mult.clear();
        knots.clear();
        mult.push_back(1);
        knots.push_back(knots_flat.front());

        // std::foreach(knots_flat.begin())

        auto n = knots_flat.size();
        for (size_t i = 1; i < n; i++)
        {
            if (fabs(knots_flat[i] - knots.back()) <= knot_eps<T>)
            {
                mult.back()++;
            }
            else
            {
                mult.push_back(1);
                knots.push_back(knots_flat[i]);
            }
        }
    }
    /**
     * @brief Builds and array of knots' single values
     * 
     * @tparam T 
     * @tparam L 
     * @param knots_flat 
     * @param knots 
     */
    template <class T>
    inline void unflat_knots(const std::vector<T> &knots_flat, std::vector<T> &knots)
    {
        knots.clear();
        knots.push_back(knots_flat.front());

        // std::foreach(knots_flat.begin())

        auto n = knots_flat.size();
        for (size_t i = 1; i < n; i++)
        {
            if (fabs(knots_flat[i] - knots.back()) > knot_eps<T>)
            {
                knots.push_back(knots_flat[i]);
            }
        }
    }
    /**
     * @brief Extract knots and multiplicity from flat knots
     * 
     * @tparam T 
     * @param knots_flat 
     * @return std::pair< std::vector<T> , std::vector<size_t> > 
     */
    template <class T>
    auto knots_and_mults(const std::vector<T> &knots_flat) -> std::pair< std::vector<T> , std::vector<size_t> >
    {
        std::vector<size_t> mult;
        std::vector<T> knots;
        unflat_knots(knots_flat,mult,knots);
        return std::make_pair(knots,mult);
    }
    /**
     * @brief Retunr unique knot with it's multiplicity
     * 
     * @tparam T 
     * @param knots_flat 
     * @return std::vector<std::pair<T, size_t>> 
     */
    template <class T>
    auto unflat_knots(const std::vector<T> &knots_flat) -> std::vector<std::pair<T, size_t>>
    {
        std::vector<std::pair<T, size_t>> result;
        result.push_back(std::make_pair(knots_flat.front(), 1));

        auto n = knots_flat.size();
        for (size_t i = 1; i < n; i++)
        {
            if (fabs(knots_flat[i] - result.back().first) < knot_eps<T>)
            {
                result.back().second++;
            }
            else
            {
                result.push_back(std::make_pair(knots_flat[i], 1));
            }
        }
        return result;
    }

    template <class T>
    auto multiplicity(const std::vector<T> &knots_flat, T knot) -> size_t
    {
        size_t m = 0;
        std::for_each(knots_flat.begin(),knots_flat.end(),[&m,&knot](const auto &k_){if(fabs(k_-knot)<knot_eps<T>) m++;});
        return m;
    }

    /**
     * @brief Build flat knots vector from knots values and multiplicity
     * 
     * @tparam T 
     * @tparam L 
     * @param knots
     * @param mult 
     * @return std::vector<T> 
     */
    template <class T, class L>
    inline std::vector<T> flat_knots(const std::vector<T> &knots, const std::vector<L> &mult)
    {
        if(knots.size()!=mult.size())
        {
            throw std::length_error( "multiplicities must have same length as knots" );
        }
        // test_size(knots, mult);

        std::vector<T> knots_flat;
        auto nk = knots.size();
        for (size_t i = 0; i < nk; i++)
        {
            L n = mult[i];
            T kv = knots[i];
            for (L j = 0; j < n; j++)
            {
                knots_flat.push_back(kv);
            }
        }
        knots_flat.shrink_to_fit();
        return knots_flat;
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
     * @brief Build flat knots with uniform multiplicity, k is supposed to be ordered
     * 
     * @tparam T 
     * @param k : knots
     * @param p : degree
     * @param mult : multiplicity
     * @return std::vector<T> 
     */
    template <typename T> 
    auto build_mult_flat_knots(const std::vector<T> &k, size_t p, size_t mult) -> std::vector<T>
    {
        std::vector<size_t> m(k.size());
        m.front()=p+1;
        std::fill(++m.begin(),--m.end(),mult);
        m.back()=p+1;
        return gbs::flat_knots(k, m);
    }
    template <typename T> 
    auto build_mult_flat_knots(T u1, T u2, size_t n, size_t p, size_t mult) -> std::vector<T>
    {
        auto k = make_range<T>(u1, u2, n);
        std::vector<size_t> m(n);
        m.front()=p+1;
        std::fill(++m.begin(),--m.end(),mult);
        m.back()=p+1;
        return gbs::flat_knots(k, m);
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
        if(fabs(d)>knot_eps<T>)
        {
            std::transform(k.begin(),k.end(),k.begin(),[&d](auto const &k_){return k_/d;});
        }
    }

    enum class KnotsCalcMode { EQUALLY_SPACED , CHORD_LENGTH, CENTRIPETAL};
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
                           [&, k_ = T(0.)](const auto &pt1,const auto &pt2) mutable {
                               k_ += gbs::distance(pt1, pt2);
                               return k_;
                           });
            break;
        case KnotsCalcMode::CENTRIPETAL:
            std::transform(++pts.begin(), pts.end(),pts.begin(), ++k.begin(),
                           [&, k_ = T(0.)](const auto &pt1,const auto &pt2) mutable {
                               k_ += sqrt(gbs::distance(pt1, pt2));
                               return k_;
                           });
            break;
        default: // KnotsCalcMode::EQUALLY_SPACED:
            T step = 1. / (pts.end() - pts.begin() - 1); //TODO, make/use range func
            std::generate(++k.begin(),k.end(),
                           [&step,k_=T(0.)]() mutable { return k_+=step; });
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
        std::vector<T> knots;
        gbs::unflat_knots(knots_flats, mult, knots);
        auto iu = std::find_if(
            knots.begin(),
            knots.end(),
            [&](const auto u_){return fabs(u_-u)<knot_eps<T>;}
            ) - knots.begin();



        auto knots_flats_ = knots_flats; //copy/move for failproof
        std::vector<std::array<T, dim>> Q(poles.size() + 1);

        auto k = std::lower_bound(knots_flats_.begin(), knots_flats_.end(), u);
        k = knots_flats_.insert(k, 1, u);
        auto ik = (k - knots_flats_.begin()) - 1;
        if(iu<knots.size() && mult[iu]>=p) return iu == 0 ? 0 : ik;

        // Start inserting knot
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
        return ik;
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
    /**
     * @brief Change parametrization to fit between k1 and k2
     * 
     * @param k1 
     * @param k2 
     */
    template <typename T>
    auto change_bounds(T k1, T k2, std::vector<T> &knots) -> void
    {
        auto k1_ = knots.front();
        auto k2_ = knots.back();
        std::transform(knots.begin(), knots.end(), knots.begin(),
                       [&](const auto &k_) {
                           return k1 + (k2 - k1) * (k_ - k1_) / (k2_ - k1_);
                       });
    }

    template <typename T, size_t nc>
    auto build_poles_matix(const std::vector<T> &k_flat, const std::vector<T> &u, size_t deg, size_t n_poles, MatrixX<T> &N) -> void
    {
        // auto n_pt = Q.size();
        // auto n_poles = int(Q.size() * nc);
        // auto n_params = int(n_poles / nc);
        auto n_params = int(u.size());
        // Eigen::MatrixX<T> N(n_poles, n_poles);

        for (int i = 0; i < n_params; i++)
        {
            for (int j = 0; j < n_poles; j++)
            {
                for (int deriv = 0; deriv < nc; deriv++)
                {
                    N(nc * i + deriv, j) = gbs::basis_function(u[i], j, deg, deriv, k_flat);
                }
            }
        }
        // return N;
    }

    template <typename T, size_t dim>
    auto increase_degree(std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles,size_t p)
    {
        // TODO use Pieg94more complex but faster method, this  needs nurbs' Bezier segment extraction
        // TODO: ckeck curve is not degenerated i.e. first and last mult = deg +1
        std::vector<T> knots;
        std::vector<size_t> mult;
        unflat_knots(knots_flats,mult,knots);
        std::transform( 
            mult.begin() , 
            mult.end(), 
            mult.begin(), 
            [](const auto &m_){return m_+1;} 
            );
        auto new_knots_flat = flat_knots(knots,mult);
        auto new_nb_poles = new_knots_flat.size() - p - 1 - 1; // TODO process periodic curve

        auto u = make_range(knots_flats.front(),knots_flats.back(),new_nb_poles);
        MatrixX<T> N(new_nb_poles, new_nb_poles);
        build_poles_matix<T,1>(new_knots_flat,u,p+1,new_nb_poles,N);
        auto N_inv = N.partialPivLu();

        points_vector<T,dim> Q{new_nb_poles};
        std::transform(
            u.begin(),
            u.end(),
            Q.begin(),
            [&](const auto &u_){return eval_value_decasteljau(u_,knots_flats,poles,p);}
        );

        std::vector<std::array<T, dim>> new_poles(new_nb_poles);
        VectorX<T> b(new_nb_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < new_nb_poles; i++)
            {
                    b(i) = Q[i][d];//Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);
            
            for (int i = 0; i < new_nb_poles; i++)
            {
                new_poles[i][d] = x(i);
            }
        }

        knots_flats = new_knots_flat;
        poles       = new_poles;
        
    }
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
        if(u-knots_flats.front()<knot_eps<T>) return;
        
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
        if(knots_flats.back()-u<knot_eps<T>) return;
        
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
        trim_begin<T,dim>(p,knots_flats,poles,fmin(u1,u2));
        trim_end<T,dim>(p,knots_flats,poles,fmax(u1,u2));
    }
/**
 * Inserts a value into a sorted vector in an ordered and unique manner, with an optional tolerance for equality checking.
 *
 * The function inserts the specified value into the sorted vector while maintaining the order of the elements.
 * If the value is already present in the vector or is within the specified tolerance, it will not be inserted again.
 *
 * @tparam T The type of the elements in the vector.
 * @param vec The vector to insert the value into. It must already be sorted.
 * @param value The value to be inserted into the vector.
 * @param tolerance The tolerance value for equality checking. Defaults to gbs::knot_eps<T>.
 */
    template <typename T>
    void insert_unique_ordered(std::vector<T>& vec, T value, T tolerance = knot_eps<T>) {
        auto iter = std::lower_bound(vec.begin(), vec.end(), value);
        if (iter == vec.end() || std::abs(value - *iter) > tolerance) {
            vec.insert(iter, value);
        } // else the value is already in the vector or close enough
    }
/**
 * Inserts a value into a sorted vector in an ordered manner.
 *
 * The function inserts the specified value into the sorted vector while maintaining the order of the elements.
 * If the value is already present in the vector, it will still be inserted at the appropriate position.
 *
 * @tparam T The type of the elements in the vector.
 * @param vec The vector to insert the value into. It must already be sorted.
 * @param value The value to be inserted into the vector.
 */
    template <typename T>
    void insert_ordered(std::vector<T>& vec, T value) {
        auto iter = std::lower_bound(vec.begin(), vec.end(), value);
        vec.insert(iter, value);
    }

} // namespace gbs