#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <Eigen/Dense>
#include <gbs/bssinterp.h>
#include <gbs/bscanalysis.h>

namespace gbs
{
    /**
     * @brief Reject ill-posed approximation sizes, consistently for every point-set
     * entry point (curve overloads previously checked this in only one place).
     *
     * A degree-p B-spline needs at least @c p+1 control points (same requirement as
     * @c check_curve), and an approximation cannot have more poles than data points
     * (it would no longer be a least-squares problem). Surfaces enforce the analogous
     * checks in bssapprox.h.
     */
    inline void check_approx_sizes(size_t n_poles, size_t p, size_t n_pts)
    {
        if (n_poles < p + 1)
        {
            throw std::length_error(
                "approx: a degree-" + std::to_string(p) + " B-spline needs at least " +
                std::to_string(p + 1) + " poles (got " + std::to_string(n_poles) + ")");
        }
        if (n_poles > n_pts)
        {
            throw std::length_error(
                "approx: n_poles (" + std::to_string(n_poles) +
                ") must not exceed the number of points (" + std::to_string(n_pts) + ")");
        }
    }
    template <typename T>
    void removeRow(MatrixX<T> &matrix, unsigned int rowToRemove)
    {
        unsigned int numRows = matrix.rows() - 1;
        unsigned int numCols = matrix.cols();

        if (rowToRemove < numRows)
            matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

        matrix.conservativeResize(numRows, numCols);
    }

    template <typename T>
    void removeColumn(MatrixX<T> &matrix, unsigned int colToRemove)
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols() - 1;

        if (colToRemove < numCols)
            matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

        matrix.conservativeResize(numRows, numCols);
    }
    /**
     * @brief Approximate a set of points with exact matchin at bounds
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts     : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param k_flat  : curve's parametrization
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx_bound_fixed(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, const std::vector<T> &k_flat) -> BSCurve<T, dim>
    {
        check_approx_sizes(n_poles, p, pts.size());
        auto n_params = int(u.size());

        // Banded one-pass assembly (A2.3): each interior collocation row has only
        // p+1 non-zero basis values (columns [span-p, span]); fill_basis_row writes
        // exactly those via a single ders_basis_funs pass instead of n_poles
        // recursive basis_function calls. We fill the full-width rows, then split off
        // the two fixed boundary columns (poles P0 and P_{n-1}) into the RHS.
        MatrixX<T> N_full(n_params - 2, n_poles);
        N_full.setZero();
        for (int i = 0; i < n_params - 2; i++)
            fill_basis_row(N_full, Eigen::Index(i), u[i + 1], k_flat, p, size_t{0});

        MatrixX<T> N = N_full.middleCols(1, n_poles - 2);
        auto Nbegin = N_full.col(0);
        auto Nend = N_full.col(n_poles - 1);

        auto p_begin = pts.front();
        auto p_end = pts.back();

        // auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto N_inv = N.colPivHouseholderQr();

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        VectorX<T> b(n_pt-2);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt-2; i++)
            {
                b(i) = pts[i+1][d] - p_begin[d] * Nbegin(i) - p_end[d] * Nend(i); //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 1; i < n_poles - 1; i++)
            {
                poles[i][d] = x(i - 1);
            }
        }
        poles.front() = pts.front();
        poles.back() = pts.back();

        return BSCurve<T,dim>(poles, k_flat, p);
    }
    /**
     * @brief Approximate a point set
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param k_flat  : curve's parametrization
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, const std::vector<T> &k_flat) -> BSCurve<T, dim>
    {
        check_approx_sizes(n_poles, p, pts.size());
        MatrixX<T> N(u.size(), n_poles);
        build_poles_matrix<T, 1>(k_flat, u, p, n_poles, N);

        // auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto N_inv = N.colPivHouseholderQr();

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        VectorX<T> b(n_pt);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt; i++)
            {
                b(i) = pts[i][d]; //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return BSCurve<T,dim>(poles, k_flat, p);
    }
    /**
     * @brief Approximate a point set, the curve's parametrization is automaticaly computed
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts     : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param fix_bound : force match on bounds
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, bool fix_bound) -> BSCurve<T, dim>
    {
        // Reject ill-posed sizes before building the knot vector: with n_poles < p+1
        // the (n-p) term in build_simple_mult_flat_knots underflows (size_t).
        check_approx_sizes(n_poles, p, pts.size());

        auto k_flat = build_simple_mult_flat_knots(u.front(),u.back(),n_poles,p);

        if (fix_bound)
        {
            return approx_bound_fixed(pts, p, n_poles, u, k_flat);
        }
        else
        {
            return approx(pts, p, n_poles, u, k_flat);
        }
    }

    template <typename T, size_t dim>
    auto refine_approx(const points_vector<T, dim> &pts, const std::vector<T> &u,const BSCurve<T, dim> &crv,bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200)  -> BSCurve<T, dim>
    {
        auto n_pts = pts.size();
        auto n_poles = crv.poles().size();
        auto p = crv.degree();
        BSCurve<T, dim> crv_refined{crv};
        // A candidate knot is only inserted when the worst point sits "well inside"
        // a span, i.e. at least this fraction of the span away from both its knots.
        // This avoids stacking a new knot right on top of an existing one; it is a
        // placement filter only and must NOT influence the stop criterion.
        constexpr T min_span_fraction = T(0.33);
        for (size_t i {} ; i < n_max; i++)
        {
            // Stop-criterion statistics: the TRUE max/average deviation over ALL
            // points (C1: d_avg_ used to accumulate only record-breaking distances;
            // C2: d_max_ used to ignore points near a knot, masking unmet tolerance).
            T d_avg_{}, d_max_{};
            // Knot-insertion candidate: the worst point that is well inside a span.
            T d_insert{}, u_max{-1};
            std::vector<T> knots{crv_refined.knotsFlats()};

            for(size_t j {} ; j < n_pts ; j++)
            {
                auto u_ = u[j];
                auto d = norm(crv_refined(u_)-pts[j]);

                // (a) true deviation over every point -> drives the stop criterion
                d_avg_ += d;
                if (d > d_max_)
                    d_max_ = d;

                // (b) where to insert next: worst point not too close to a knot
                if (d > d_insert)
                {
                    auto it = std::lower_bound(knots.begin(), knots.end(), u_);
                    if (it == knots.end() || it == knots.begin())
                        continue;
                    auto uh = *it;
                    auto ul = *std::prev(it);
                    if (uh <= ul)
                        continue;
                    auto dul = (u_ - ul) / (uh - ul);
                    auto duh = (uh - u_) / (uh - ul);
                    if (dul > min_span_fraction && duh > min_span_fraction)
                    {
                        d_insert = d;
                        u_max = u_;
                    }
                }
            }

            d_avg_ /= pts.size();

            // No span is "open enough" for a new knot: we cannot refine further,
            // regardless of whether the tolerance is met -> do not claim convergence,
            // just stop (the returned curve reflects the best we can place).
            if (u_max < 0)
                break;

            auto it = std::lower_bound(knots.begin(), knots.end(), u_max);
            knots.insert(it, u_max);
            n_poles++;
            if (fix_bound)
            {
                crv_refined = approx_bound_fixed(pts, p, n_poles, u, knots);
            }
            else
            {
                crv_refined = approx(pts, p, n_poles, u, knots);
            }
            // std::cout << "d_avg: " << d_avg_ << ", d_max: " << d_max_ << ", u_max:" << u_max << std::endl;
            if (
                (d_max_ < d_max && d_avg_ < d_avg)
                || n_poles>= pts.size() - 5 // The idea of approximation id to get less pole tha  interpolation
                )
                break;
            // std::cout << "n poles: " << crv_refined.poles().size() << ", n_flat: " << crv_refined.knotsFlats().size() << std::endl;
        }
        return crv_refined;
    }
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u, size_t p, bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200) -> BSCurve<T, dim>
    {
        auto n_poles = p * 2;
        // auto n_poles = pts.size() / 5;
        auto crv = approx(pts, p, n_poles, u, fix_bound);
        return refine_approx(pts,u,crv,fix_bound,d_max,d_avg,n_max);
    }
    /**
     * @brief Approximate a point set, the curve's poles' number and parametrization are automaticaly computed  
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts : point set to approximate
     * @param p   : curve's degree
     * @param mode : curve parametrization mode
     * @param fix_bound : force match on bounds
     * @param adimensionnal : makes bound going from 0. to 1.
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, KnotsCalcMode mode, bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200, bool adimensionnal =false) -> BSCurve<T, dim>
    {
        auto u = curve_parametrization(pts, mode, adimensionnal);
        return approx(pts,u,p,fix_bound,d_max,d_avg,n_max);
    }

    /**
     * @brief  Approximate a point set, the curve's parametrization is automaticaly computed and the bounds are matched
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts  : point set to approximate
     * @param p   : curve's degree
     * @param n_poles : desired poles' number
     * @param mode   : curve parametrization mode
     * @param adimensionnal : makes bound going from 0. to 1.
     * @return auto 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, KnotsCalcMode mode, bool adimensionnal =false) // -> BSCurve<T,dim>
    {
        auto u = curve_parametrization(pts, mode, adimensionnal);
        return approx(pts, p, n_poles, u, true);
    }
    /**
     * @brief Build a BSCurve approximation of crv, bounds approximation match crv's bounds.
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param deviation 
     * @param n_poles 
     * @param p 
     * @param mode 
     * @param np : Number of point for preliminary curve's discretization 
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation, size_t p, KnotsCalcMode mode, size_t np = 30) -> BSCurve<T, dim>
    {
        auto pts = discretize(crv,np,deviation);
        return approx(pts,p,mode,true);
    }

    /**
     * @brief  Build a BSCurve approximation of crv, bounds approximation match crv's bounds.
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param deviation 
     * @param n_poles 
     * @param p 
     * @param mode 
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation,  size_t n_poles, size_t p, KnotsCalcMode mode) -> BSCurve<T, dim>
    {
        if(n_poles < p + 1)
        {
            throw std::length_error("More poles needed for approximation");
        }
        auto np = n_poles*10; // gives a reasonable number of points for starting discretization
        auto pts = discretize(crv,np,deviation);
        return approx(pts,p,n_poles,mode);
    }

    /**
     * @brief Approximate curve with knots imposed
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param p 
     * @param n_poles 
     * @param k_flat 
     * @param deviation 
     * @param np 
     * @param n_max_pts 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T, dim> &crv, size_t p, size_t n_poles, std::vector<T> k_flat, T deviation, size_t np, size_t n_max_pts=5000)
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation,n_max_pts);
        std::vector<T> u(u_lst.size());
        std::transform(
            GBS_PAR_EXEC
            u_lst.begin(),u_lst.end(),u.begin(),
            [](const auto &u_){return u_;}
        );
        auto pts = make_points(crv,u_lst);
        return approx_bound_fixed<T,dim>(pts,p,n_poles,u,k_flat);
    }
    /**
     * @brief Approximate curve respecting original curve's parametrization
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param deviation 
     * @param tol : global tolerance, max tol = 10 * tol
     * @param p   : approximation curve's degree
     * @param np  : number of point for first approximation
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation,T tol, size_t p, size_t np = 30) -> BSCurve<T, dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation);
        auto pts = make_points(crv,u_lst);
        size_t n_poles = std::max<size_t>(np / 3, p+1 );
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        auto crv_approx = approx(pts, p, n_poles, u, true);
        return refine_approx<T,dim>(pts,u,crv_approx,true,10.*tol,tol);
    }
    /**
     * @brief Approximate curve respecting original curve's parametrization
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param deviation 
     * @param tol : global tolerance, max tol = 10 * tol
     * @param p   : approximation curve's degree
     * @param np  : number of point for first approximation
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation,T tol, size_t p, size_t n_poles, size_t np) -> BSCurve<T, dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation);
        auto pts = make_points(crv,u_lst);
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        auto crv_approx = approx(pts, p, n_poles, u, true);
        return refine_approx<T,dim>(pts,u,crv_approx,true,10.*tol,tol);
    }
    /**
     * @brief Approximate curve respecting original curve's parametrization between params u1 and u2
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param u1 
     * @param u2 
     * @param deviation 
     * @param tol 
     * @param p 
     * @param np 
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T u1, T u2, T deviation,T tol, size_t p, size_t np = 30) -> BSCurve<T, dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation);
        auto pts = make_points(crv,u_lst);
        auto n_poles = p * 2;
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        auto crv_approx = approx(pts, p, n_poles, u, true);
        return refine_approx<T,dim>(pts,u,crv_approx,true,10.*tol,tol);
    }

    template <typename T, size_t dim>
    auto approx(
        const Curve<T,dim> &crv,
        const std::function<std::array<T,dim>(const std::array<T,dim>&)> &trf_func, 
        T deviation,
        T tol, 
        size_t p, 
        size_t np = 30
        ) -> BSCurve<T, dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation);
        auto pts = make_points(crv,u_lst);
        std::transform(
            GBS_PAR_EXEC
            pts.begin(),
            pts.end(),
            pts.begin(),
            trf_func
        );
        auto n_poles = p * 2;
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        auto crv_approx = approx(pts, p, n_poles, u, true);
        return refine_approx<T,dim>(pts,u,crv_approx,true,10.*tol,tol);
    }

    /**
     * @brief Approximate curve with uniformly spaced points
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv: curve 
     * @param p  : degree
     * @param np : nomber of points
     * @return BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, size_t p, size_t np ) -> BSCurve<T, dim>
    {
        auto u_lst = uniform_distrib_params<T, dim>(crv, np);
        auto pts = make_points(crv,u_lst);
        auto n_poles = p * 2;
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        return approx(pts, p, n_poles, u, true);
    }
} // namespace gbs