#pragma once
#include <vector>
#include <Eigen/Dense>
#include <gbs/bssinterp.h>
#include <gbs/bscanalysis.h>

namespace gbs
{
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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx_bound_fixed(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, std::vector<T> k_flat) -> gbs::BSCurve<T, dim>
    {
        auto n_params = int(u.size());
        MatrixX<T> N(n_params-2, n_poles-2);


        for (int i = 0; i < n_params-2; i++) // TODO check if matrix is row or col dominant
        {
            for (int j = 0; j < n_poles-2; j++)
            {
                    N(i , j) = gbs::basis_function(u[i+1], j+1, p, 0, k_flat);
            }
        }

        auto p_begin = pts.front();
        auto p_end = pts.back();
        std::vector<T> Nbegin(n_params-2),Nend(n_params-2);
        for (int i = 0; i < n_params-2; i++)
        {
            Nbegin[i] = gbs::basis_function(u[i+1], 0, p, 0, k_flat);
            Nend[i] = gbs::basis_function(u[i+1], n_poles-1, p, 0, k_flat);
        }


        // auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto N_inv = N.colPivHouseholderQr();

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        VectorX<T> b(n_pt-2);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt-2; i++)
            {
                b(i) = pts[i+1][d] - p_begin[d] * Nbegin[i] - p_end[d] * Nend[i]; //Pas top au niveau de la localisation mémoire
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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, std::vector<T> k_flat) -> gbs::BSCurve<T, dim>
    {
        MatrixX<T> N(u.size(), n_poles);
        build_poles_matix<T, 1>(k_flat, u, p, n_poles, N);

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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, bool fix_bound) -> gbs::BSCurve<T, dim>
    {

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
    auto refine_approx(const points_vector<T, dim> &pts, const std::vector<T> &u,const gbs::BSCurve<T, dim> &crv,bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200)  -> gbs::BSCurve<T, dim>
    {
        auto n_pts = pts.size();
        auto n_poles = crv.poles().size();
        auto p = crv.degree();
        gbs::BSCurve<T, dim> crv_refined{crv};
        auto j_ = make_range<size_t>(size_t{},n_pts);
        for (size_t i {} ; i < n_max; i++)
        {
            auto d_avg_ = 0., d_max_ = 0., u_max = -1.;
            std::vector<T> knots{crv_refined.knotsFlats()};

            for(size_t j {} ; j < n_pts ; j++)
            // std::for_each(std::execution::par, j_.begin(),j_.end(),
            // [&](auto j)
            {
                auto u_ = u[j];
                auto d = norm(crv_refined(u_)-pts[j]);
                if(d>d_max_)
                {
                    auto it = std::lower_bound(knots.begin(), knots.end(), u_);
                    auto uh = *it;
                    if (it == knots.begin())
                        it = std::next(it);
                    auto ul = *(std::next(it, -1));
                    auto dul = (u_ - ul) / (uh - ul);
                    auto duh = (uh - u_) / (uh - ul);
                    if (dul > 0.33 && duh > 0.33)
                    {
                        d_max_ =d;
                        u_max = u_;
                    }
                    d_avg_ += d;
                }
            }
            // );


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
            d_avg_ /= pts.size();
            // std::cout << "d_avg: " << d_avg_ << ", d_max: " << d_max_ << ", u_max:" << u_max << std::endl;
            if (
                u_max < 0 // Nothing needed
                || (d_max_ < d_max && d_avg_ < d_avg)
                || n_poles>= pts.size() - 5 // The idea of approximation id to get less pole tha  interpolation
                )
                break;
            // std::cout << "n poles: " << crv_refined.poles().size() << ", n_flat: " << crv_refined.knotsFlats().size() << std::endl;
        }
        return crv_refined;
    }
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u, size_t p, bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200) -> gbs::BSCurve<T, dim>
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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, gbs::KnotsCalcMode mode, bool fix_bound, T d_max = 1e-3, T d_avg= 1e-4, size_t n_max = 200) -> gbs::BSCurve<T, dim>
    {
        auto u = gbs::curve_parametrization(pts, mode, true);
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
     * @return auto 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, gbs::KnotsCalcMode mode) // -> gbs::BSCurve<T,dim>
    {
        auto u = gbs::curve_parametrization(pts, mode, true);
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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation, size_t p, gbs::KnotsCalcMode mode, size_t np = 30) -> gbs::BSCurve<T, dim>
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
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const Curve<T,dim> &crv, T deviation,  size_t n_poles, size_t p, gbs::KnotsCalcMode mode) -> gbs::BSCurve<T, dim>
    {
        if(n_poles < p + 1)
        {
            throw std::exception("More poles needed for approximation");
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
    auto approx(const gbs::Curve<T, dim> &crv, size_t p, size_t n_poles, std::vector<T> k_flat, T deviation, size_t np, size_t n_max_pts=5000)
    {
        auto u_lst = deviation_based_params<T, dim>(crv, np,deviation,n_max_pts);
        std::vector<T> u(u_lst.size());
        std::transform(
            std::execution::par,
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
        auto n_poles = p * 2;
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
    /**
     * @brief Approximate curve with uniforly spaced points
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param p 
     * @param np 
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