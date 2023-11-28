#pragma once
#include <gbs/gbslib.h>
#include <gbs/vecop.h>
#include <gbs/basisfunctions.h>
#include <gbs/maths.h>
#include <list>
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
        return flat_knots(k, m);
    }
    template <typename T> 
    auto build_mult_flat_knots(T u1, T u2, size_t n, size_t p, size_t mult) -> std::vector<T>
    {
        auto k = make_range<T>(u1, u2, n);
        std::vector<size_t> m(n);
        m.front()=p+1;
        std::fill(++m.begin(),--m.end(),mult);
        m.back()=p+1;
        return flat_knots(k, m);
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
                               k_ += distance(pt1, pt2);
                               return k_;
                           });
            break;
        case KnotsCalcMode::CENTRIPETAL:
            std::transform(++pts.begin(), pts.end(),pts.begin(), ++k.begin(),
                           [&, k_ = T(0.)](const auto &pt1,const auto &pt2) mutable {
                               k_ += sqrt(distance(pt1, pt2));
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
     * Inserts new knots, if possible, into a B-spline or NURBS knot vector and updates the control points.
     *
     * Template Parameters:
     *   T: The type of the knot values and control point coordinates (typically float or double).
     *   dim: The dimension of the control points (e.g., 2 for 2D points, 3 for 3D points).
     *
     * Parameters:
     *   u: The knot value to be inserted.
     *   p: The degree of the B-spline or NURBS curve.
     *   r: The number of times the knot u is to be inserted, adjusted based on existing multiplicity.
     *   UP: The original knot vector, which will be updated with the new knots.
     *   P: The original control points, which will be updated to the new control points.
     *
     * Returns:
     *   The index in the updated knot vector where the last insertion of the knot u occurred.
     *   This value is adjusted based on the degree of the curve and the current multiplicity of u.
     *
     * Notes:
     *   - The function updates both UP and P in place.
     *   - The function checks if there is enough room to insert the new knots given their current multiplicity.
     *   - If there is no room for the additional knots (s >= p), the function returns an adjusted index.
     *   - The value of r is minimized to not exceed p - s, ensuring valid knot insertion.
     */
    template <typename T, size_t dim>
    auto insert_knots(T u, size_t p, size_t r, std::vector<T> &UP, std::vector<std::array<T, dim>> &P)
    {

        auto np   = P.size()-1;
        auto mp   = np + p + 1;
        auto nq   = np + r;
        // Knots are added after the potential existing value
        auto k_it = std::upper_bound(UP.begin(), UP.end(), u); 
        auto k    = std::distance(UP.begin(), k_it) - 1;
        // Checking the multiplicity
        auto s    = std::distance(std::lower_bound(UP.begin(), UP.end(), u), k_it);
        // Check if there is room for insertion
        if(s>=p) // No room
        {
            auto at_end = k_it == UP.end() ? 1 : 0;
            return k-p-at_end;
        }
        r = std::min(r, p-s);

        // Allocate the new knots and vector vectors
        std::vector<std::array<T, dim>> Q(P.size() + r);
        std::vector<T> UQ(UP.size() + r);

        // Create the new knots vector
        for( size_t i{}; i<=k; i++) UQ[i] = UP[i];
        for( size_t i{1}; i<=r; i++) UQ[i+k] = u;
        for( size_t i= k+1; i<=mp; i++) UQ[i+r] = UP[i];
        // Prepare new poles vector
        for(size_t i{}; i <=k-p; i++) Q[i] = P[i];
        for(size_t i =k-s; i <=np; i++) Q[i+r] = P[i];
        // Temporary storage
        std::vector<std::array<T, dim>> R(p+1);
        for(size_t i{}; i<=p-s; i++) R[i] = P[k-p+i];

        // Insert knots
        size_t L{};
        for(size_t j{1}; j <= r ; j++)
        {
            L = k-p+j;
            for(size_t i{}; i<=p-j-s; i++)
            {
                auto alpha = (u-UP[L+i])/(UP[i+k+1]-UP[L+i]);
                R[i] = alpha*R[i+1]+(1.-alpha)*R[i];
            }
            Q[L]=R[0];
            Q[k+r-j-s] = R[p-j-s];
        }
        for(auto i{L+1}; i < k-s; i++)
        {
            Q[i] = R[i-L];
        }
        //move
        UP = std::move(UQ);
        P  = std::move(Q);

        return k+1+r-(p+1);
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
        return insert_knots(u, p , 1, knots_flats, poles);
    }

    // /**
    //  * @brief Remove knot if possible
    //  * 
    //  * @tparam T 
    //  * @tparam dim 
    //  * @param u 
    //  * @param p 
    //  * @param knots_flats 
    //  * @param P 
    //  * @param tol 
    //  * @return auto 
    //  */
    // template <typename T, size_t dim>
    // auto remove_knot(T u, size_t p, std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &P, T tol)
    // {
    //     // auto knots_flats_ = knots_flats; //copy/move for fail saife
    //     auto Pi_ = P;
    //     auto Pj_ = P;

    //     // std::vector<T> knots;
    //     // std::vector<size_t> m;
    //     // unflat_knots(knots_flats_,m,knots);

    //     // auto k = std::next(std::lower_bound(knots.begin(), knots.end(), u));
    //     // auto r = std::next(std::lower_bound(knots_flats_.begin(), knots_flats_.end(), *k)-1);

    //     // auto r = std::next(std::upper_bound(knots_flats_.begin(), knots_flats_.end(),u),-1)-knots_flats_.begin();
    //     // auto s = multiplicity(u,knots_flats_);

    //     auto k = std::next(std::upper_bound(knots_flats.begin(), knots_flats.end(),u),-1);
    //     auto r = k - knots_flats.begin();
    //     auto s = multiplicity(u,knots_flats);

    //     Pi_.erase(std::next(Pi_.begin(),(2*r-s-p)/2-1));// less mem and less copy possible
    //     Pj_.erase(std::next(Pj_.begin(),(2*r-s-p)/2-1));

    //     size_t i = r-p;
    //     size_t j = r-s;
    //     int t = 0;
    //     int d = j-i;
    //     while (d>t)
    //     {
    //         auto ai = (u-knots_flats[i]   )/(knots_flats[i+p+1+t]-knots_flats[i]);
    //         auto aj = (u-knots_flats[j-t])/(knots_flats[j+p+1]-knots_flats[j-t]);
    //         Pi_[i] = (P[i] - (1 - ai) * P[i - 1]) / ai;
    //         Pj_[j] = (P[j] - aj * P[j + 1]) / (1 - aj);
    //         i++;
    //         j--;
    //         d=j-i;
    //     }

    //     auto D = delta(Pi_, Pj_);
    //     // if (std::reduce(
    //     //         std::execution::par,
    //     //         D.begin(), D.end(),
    //     //         T(0.),
    //     //         [](const auto &t_, const auto &v_) { return t_ + sq_norm(v_); }) < tol)
    //     // {
    //         knots_flats.erase(k);
    //         P = Pi_;
    //     // }
    // }

/**
  \brief Removes an internal knot from a curve.
  This is A5.8 on p185 from the NURB book modified to not check for 
  tolerance before removing the knot.

  \param r  the knot to remove
  \param s  the multiplicity of the knot
  \param num  the number of times to try to remove the knot

  \warning r \b must be an internal knot.

  \author Philippe Lavoie 
  \date 24 January 1997
**/

/**
 * @brief Removes an internal knot from a curve.
  This is A5.8 on p185 from the NURB book modified to not check for 
  tolerance before removing the knot.
 * 
 * @tparam T 
 * @tparam dim 
 * @param u The knot value to remove
 * @param num the number of times to try to remove the knot
 * @param U Flat knots vector
 * @param P Poles
 * @param p Degree
 * @return void 
 */
    template <typename T, size_t dim>
    void remove_knot(T u, size_t num, std::vector<T> &U, std::vector<std::array<T, dim>> &P, size_t p)
    {
        // Determine the knot interval that contains u
        auto r = std::distance(U.begin(), std::ranges::upper_bound(U, u)) - 1;
        // if u is not an internal knot
        if(std::abs(U[r]-u)>knot_eps<T> || r ==0 || r == U.size())
        {
            return ;
        }
        // Find the multiplicity of the knot u
        auto s = multiplicity(u, U);

        int m = U.size();
        int ord = p + 1;
        int fout = (2 * r - s - p) / 2;
        int last = r - s;
        int first = r - p;
        T alfi, alfj;
        int i, j, k, ii, jj, off;

        std::vector<std::array<T, dim>> temp(2 * p + 1);

        if (num < 1)
        {
            return;
        }

        int t;
        for (t = 0; t < num; ++t)
        {
            off = first - 1;
            temp[0] = P[off];
            temp[last + 1 - off] = P[last + 1];
            i = first;
            j = last;
            ii = 1;
            jj = last - off;
            while (j - i > t)
            {
                alfi = (u - U[i]) / (U[i + ord + t] - U[i]);
                alfj = (u - U[j - t]) / (U[j + ord] - U[j - t]);
                temp[ii] = (P[i] - (1.0 - alfi) * temp[ii - 1]) / alfi;
                temp[jj] = (P[j] - alfj * temp[jj + 1]) / (1.0 - alfj);
                ++i;
                ++ii;
                --j;
                --jj;
            }
            i = first;
            j = last;
            while (j - i > t)
            {
                P[i] = temp[i - off];
                P[j] = temp[j - off];
                ++i;
                --j;
            }
            --first;
            ++last;
        }
        if (t == 0)
        {
            return;
        }

        for (k = r + 1; k < m; ++k)
            U[k - t] = U[k];
        j = fout;
        i = j; // Pj thru Pi will be overwritten
        for (k = 1; k < t; k++)
            if ((k % 2) == 1)
                ++i;
            else
                --j;
        for (k = i + 1; k < P.size(); k++)
        { // Shift
            P[j++] = P[k];
        }

        P.resize(P.size()-t);
        U.resize(U.size()-t);

    }

    /**
     * @brief This function removes a knot from a B-spline.
     *
     * @param u The value of the knot to remove.
     * @param p The order of the B-spline.
     * @param num The number of times the knot u should be removed.
     * @param U A vector of knots.
     * @param Pw A vector of control points.
     * @param tol The tolerance for knot removal.
     *
     * @return The actual number of times the knot was removed.
     */
    template <typename T, size_t dim>
    auto remove_knot(T u, size_t p, size_t num, std::vector<T> &U, std::vector<std::array<T, dim>> &Pw, T tol)
    {
        // Determine the knot interval that contains u
        auto r = std::distance(U.begin(), std::ranges::upper_bound(U, u)) - 1;
        // if u is not a knot
        if(std::abs(U[r]-u)>knot_eps<T>)
        {
            // return false;
        }
        // Find the multiplicity of the knot u
        auto s = multiplicity(u, U);
        auto m = U.size();
        auto n = Pw.size();
        
        // Determine the first and last affected control points
        auto first = r - p;
        auto last = r - s;
        
        size_t t{0};
        std::list<std::array<T, dim>> Pi;
        std::list<std::array<T, dim>> Pj;
        
        // Try to remove the knot u num times
        for (; t < num; t++)
        {
            auto i{first};
            auto j{last};
            
            // Initialize lists of control points Pi and Pj
            Pi = {Pw[first - 1]};
            Pj = {Pw[last + 1]};
            
            while (j > i + t)
            {
                // Compute blending factors
                auto ai = (u - U[i]) / (U[i + p + 1 + t] - U[i]);
                auto aj = (u - U[j - t]) / (U[j + p + 1] - U[j - t]);

                // Compute new control points
                Pi.push_back((Pw[i] - (1 - ai) * Pi.back()) / ai);
                Pj.push_front((Pw[j] - aj * Pj.front()) / (1 - aj));
                i++;
                j--;
            }
            
            bool remove_flag = false;
            if (j < i + t)
            {
                // If the knot can be removed, update Pi and Pj
                auto d = distance(Pi.back(), Pj.front());
                if (d < tol)
                {
                    remove_flag = true;
                    Pj.pop_front();
                }
            }
            else
            {
                // Check if knot removal is possible
                auto ai = (u - U[i]) / (U[i + p + 1 + t] - U[i]);
                auto d = distance(Pw[i], ai * Pi.back() + (1 - ai) * Pj.front());
                if (d < tol)
                {
                    remove_flag = true;
                    Pj.pop_front();
                    Pi.pop_back();
                }
            }
            
            first--;
            last++;

            if (!remove_flag)
            {
                break;
            }

        }
        if (t > 0)
        {
            // If at least one knot has been removed, update the control points Pw
            auto Pw_beg = Pw.begin();
            std::vector<std::array<T, dim>> head{Pw_beg, std::next(Pw_beg, first)};
            std::vector<std::array<T, dim>> tail{std::next(Pw_beg, last + 1), Pw.end()};
            Pw = std::move(head);
            Pw.insert(Pw.end(), Pi.begin(), Pi.end());
            Pw.insert(Pw.end(), Pj.begin(), Pj.end());
            Pw.insert(Pw.end(), tail.begin(), tail.end());

            U.erase(std::next(U.begin(),r-t+1), std::next(U.begin(),r+1));

            assert(U.size()-p-1==Pw.size());
        }

        // Return the number of times the knot was actually removed
        return t;
    }
    /**
     * @brief This function removes a knot from a B-spline.
     *
     * @param u The value of the knot to remove.
     * @param p The order of the B-spline.
     * @param U A vector of knots.
     * @param Pw A vector of control points.
     * @param tol The tolerance for knot removal.
     *
     * @return If knot was removed.
     */
    template <typename T, size_t dim>
    auto remove_knot(T u, size_t p, std::vector<T> &U, std::vector<std::array<T, dim>> &Pw, T tol)
    {
        // Determine the knot interval that contains u
        auto r = std::distance(U.begin(), std::ranges::upper_bound(U, u)) - 1;

        // if u is not a knot
        if(std::abs(U[r]-u)>knot_eps<T>)
        {
            return false;
        }
        
        // Find the multiplicity of the knot u
        auto s = multiplicity(u, U);
        auto m = U.size();
        auto n = Pw.size();
        
        // Determine the first and last affected control points
        auto first = r - p;
        auto last = r - s;
        
        std::list<std::array<T, dim>> Pi;
        std::list<std::array<T, dim>> Pj;
        

        auto i{first};
        auto j{last};
        
        // Initialize lists of control points Pi and Pj
        Pi = {Pw[first - 1]};
        Pj = {Pw[last + 1]};
        
        while (j > i)
        {
            // Compute blending factors
            auto ai = (u - U[i]) / (U[i + p + 1 ] - U[i]);
            auto aj = (u - U[j]) / (U[j + p + 1] - U[j]);

            // Compute new control points
            Pi.push_back((Pw[i] - (1 - ai) * Pi.back()) / ai);
            Pj.push_front((Pw[j] - aj * Pj.front()) / (1 - aj));
            i++;
            j--;
        }
        
        bool remove_flag = false;
        if (j < i)
        {
            // If the knot can be removed, update Pi and Pj
            auto d = distance(Pi.back(), Pj.front());
            if (d < tol)
            {
                remove_flag = true;
                Pj.pop_front();
            }
        }
        else
        {
            // Check if knot removal is possible
            auto ai = (u - U[i]) / (U[i + p + 1] - U[i]);
            auto d = distance(Pw[i], ai * Pi.back() + (1 - ai) * Pj.front());
            if (d < tol)
            {
                remove_flag = true;
                Pj.pop_front();
                Pi.pop_back();
            }
        }

        if (remove_flag)
        {
            auto Pw_beg = Pw.begin();
            std::vector<std::array<T, dim>> head{Pw_beg, std::next(Pw_beg, first)};
            std::vector<std::array<T, dim>> tail{std::next(Pw_beg, last + 1), Pw.end()};
            Pw = std::move(head);
            Pw.insert(Pw.end(), Pi.begin(), Pi.end());
            Pw.insert(Pw.end(), Pj.begin(), Pj.end());
            Pw.insert(Pw.end(), tail.begin(), tail.end());

            U.erase(std::next(U.begin(),r), std::next(U.begin(),r+1));
        }


        return remove_flag;
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

/**
 * Increases the degree of a Bezier curve represented by control points P
 * using degree elevation technique.
 *
 * @tparam T         The data type of the control points.
 * @tparam dim       The dimension of each control point.
 * @param P          The original control points of the Bezier curve.
 * @param p          The original degree of the Bezier curve.
 * @param t          The desired degree increase (new_degree = p + t).
 * @return           The control points of the Bezier curve with the increased degree.
 */
    template<typename T, size_t dim>
    auto increase_bezier_degree(const std::vector<std::array<T,dim>> &P, size_t p, size_t t)
    {
        std::vector<std::array<T,dim>> Pt(p+1+t);
        for(int i{}; i <= p+t; i++)
        {
            Pt[i] = std::array<T,dim>{};
            // Calculate the coefficients using binomial law for Bezier curve degree elevation.
            for(auto j{std::max<int>(0,i-t)}; j <= std::min<int>(p,i); j++)
            {
                auto C = binomial_law<T,int>(p,j)*binomial_law<T,int>(t, i-j) / binomial_law<T,int>(p+t,i);
                Pt[i] = Pt[i] + C * P[j];
            }
        }
        return Pt;
    }

/**
 * @brief Split a NURBS curve into its corresponding Bezier segments
 * @tparam T Numerical type of the curve's parameters and coordinates
 * @tparam dim Dimensionality of the curve
 * @param k Knot vector of the NURBS curve
 * @param poles Control points of the NURBS curve
 * @param p Degree of the NURBS curve
 * @return Vector of pairs, where each pair contains a Bezier segment's control points and its parameter range
 */
    template<typename T, size_t dim>
    auto bezier_segments(std::vector<T> k,  std::vector<std::array<T,dim>> poles, size_t p) 
        -> std::vector< std::pair< std::vector<std::array<T,dim>>, std::pair<T,T> > >
    {
        // Insert knots to obtain segments of degree p with multiplicity p
        {
            auto [U, M] = knots_and_mults(k);
            auto m = U.size() - 1;
            for (size_t i{1}; i < m; i++)
            {
                insert_knots(U[i], p, p - M[i], k, poles);
            }
        }

        // Convert NURBS segments to Bezier form
        std::vector< std::pair< std::vector<std::array<T,dim>>, std::pair<T,T> > > bezier_seg;
        {
            auto [U, M] = knots_and_mults(k);
            
            assert(poles.size() == p * (U.size() - 1)+1); // sanity check: the number of poles should match the number of segments

            auto step = p+1;
            auto itp{poles.begin()};
            auto itu{U.begin()};
            auto itp_ {itp};
            auto itu_ {itu};
            
            // iterate through each segment
            while (itp_ != poles.end())
            {
                itp_ = next(itp, step); // advance step positions
                itu_ = next(itu); // advance one position

                std::vector<std::array<T, dim>> poles_bezier{itp, itp_}; // extract the control points for the current Bezier segment

                // append the current segment (control points and parameter range) to the output list
                bezier_seg.push_back(std::make_pair(poles_bezier, std::make_pair(*itu, *itu_)));

                itp=next(itp_,-1); // back up one position
                itu = itu_;
            }
        }

        return bezier_seg; // return the list of Bezier segments
    }

    template <typename T, size_t nc>
    auto build_poles_matrix(const std::vector<T> &k_flat, const std::vector<T> &u, size_t deg, size_t n_poles, MatrixX<T> &N) -> void
    {
        auto n_params = int(u.size());

        for (int i = 0; i < n_params; i++)
        {
            for (int j = 0; j < n_poles; j++)
            {
                for (int deriv = 0; deriv < nc; deriv++)
                {
                    N(nc * i + deriv, j) = basis_function(u[i], j, deg, deriv, k_flat);
                }
            }
        }
    }
    template <typename T, size_t dim>
    auto increase_degree(std::vector<T> &k, std::vector<std::array<T, dim>> &poles,size_t p, size_t t)
    {
        auto [U, M] = knots_and_mults(k);

        // Split bezier
        auto bezier_seg = bezier_segments(k, poles, p);

        // Increase Bezier degree
        transform(
            bezier_seg.begin(), 
            bezier_seg.end(), 
            bezier_seg.begin(), 
            [p,t](const auto &pa)
            {
                return make_pair(
                    increase_bezier_degree(pa.first, p, t),
                    pa.second
                );
            }
        );
        size_t p_new = p+t;
        // Join poles and knots
        std::vector<std::array<T,dim>> poles_new{bezier_seg.front().first.begin(), bezier_seg.front().first.end()};
        std::vector<size_t> M_new{p_new+1};
        std::vector<T> U_new{bezier_seg.front().second.first,bezier_seg.front().second.second};
        for(size_t s{1}; s<bezier_seg.size(); s++)
        {
            const auto & seg = bezier_seg[s];
            for(size_t i{1}; i <= p_new; i++)
            {
                poles_new.push_back(seg.first[i]);
            }
            M_new.push_back(p_new);
            U_new.push_back(seg.second.second);
        }
        M_new.push_back(p_new+1);
        // Remove knots
        auto k_new = flat_knots(U_new, M_new);
        auto m = U.size()-1;
        for(size_t i{1}; i < m; i++){
            remove_knot(U[i], p_new-(M[i]+t), k_new, poles_new, p_new);
        }
        k=std::move(k_new);
        poles=std::move(poles_new);
    }

    // template <typename T, size_t dim>
    // auto increase_degree(std::vector<T> &knots_flats, std::vector<std::array<T, dim>> &poles,size_t p)
    // {
    //     // TODO use Pieg94more complex but faster method, this  needs nurbs' Bezier segment extraction
    //     // TODO: ckeck curve is not degenerated i.e. first and last mult = deg +1
    //     std::vector<T> knots;
    //     std::vector<size_t> mult;
    //     unflat_knots(knots_flats,mult,knots);
    //     std::transform( 
    //         mult.begin() , 
    //         mult.end(), 
    //         mult.begin(), 
    //         [](const auto &m_){return m_+1;} 
    //         );
    //     auto new_knots_flat = flat_knots(knots,mult);
    //     auto new_nb_poles = new_knots_flat.size() - p - 1 - 1; // TODO process periodic curve

    //     auto u = make_range(knots_flats.front(),knots_flats.back(),new_nb_poles);
    //     MatrixX<T> N(new_nb_poles, new_nb_poles);
    //     build_poles_matrix<T,1>(new_knots_flat,u,p+1,new_nb_poles,N);
    //     auto N_inv = N.partialPivLu();

    //     points_vector<T,dim> Q{new_nb_poles};
    //     std::transform(
    //         u.begin(),
    //         u.end(),
    //         Q.begin(),
    //         [&](const auto &u_){return eval_value_decasteljau(u_,knots_flats,poles,p);}
    //     );

    //     std::vector<std::array<T, dim>> new_poles(new_nb_poles);
    //     VectorX<T> b(new_nb_poles);
    //     for (int d = 0; d < dim; d++)
    //     {
    //         for (int i = 0; i < new_nb_poles; i++)
    //         {
    //                 b(i) = Q[i][d];//Pas top au niveau de la localisation mémoire
    //         }

    //         auto x = N_inv.solve(b);
            
    //         for (int i = 0; i < new_nb_poles; i++)
    //         {
    //             new_poles[i][d] = x(i);
    //         }
    //     }

    //     knots_flats = new_knots_flat;
    //     poles       = new_poles;
        
    // }
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
        
        insert_knots(u, p, p, knots_flats, poles);

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
        
        insert_knots(u, p, p, knots_flats, poles);

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
 * @param tolerance The tolerance value for equality checking. Defaults to knot_eps<T>.
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