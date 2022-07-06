#pragma once
#include <nlopt.hpp>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/extrema.h>
#include <boost/math/quadrature/gauss.hpp>
#include <list>

namespace{
    const size_t N_gauss_pt = 250;
}
namespace gbs
{
/**
 * @brief Computes the deviation from a set of points
 * 
 * @tparam T 
 * @tparam dim 
 * @param points 
 * @param crv 
 * @return (pos of max, max distance, average distance) 
 */
    template <typename T, size_t dim>
    auto dev_from_points(const std::vector<std::array<T, dim>> &points, const BSCurve<T, dim> &crv)
    {
        auto d_avg = 0., d_max = 0., u_max = 0.;
        auto u0 = crv.knotsFlats().front();
        std::for_each(
            points.begin(),
            points.end(),
            [&](const auto &pnt) {
                auto [res_u,res_d] = extrema_curve_point(crv, pnt, u0, 1e-6);
                u0 = res_u;
                if (res_d > d_max)
                {
                    d_max = res_d;
                    u_max = res_u;
                }
                d_avg += res_d;
            });
        d_avg /= points.size();
        return std::make_tuple(u_max, d_max, d_avg);
    }
/**
 * @brief Compute Curve length betwee bounds
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam rational 
 * @param crv : The curve
 * @param u1  : Starting point
 * @param u2  : End point
 * @param d  : derivate order max 2 respective to u2
 * @return T 
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto length(const Curve<T,dim> &crv,T u1 , T u2, size_t d = 0) -> T
    {
        using namespace boost::math::quadrature;

        auto f = [&](const auto& u) { return norm(crv.value(u,1)); };
        switch (d)
        {
        case 0:
            return  gauss<T, N>::integrate(f, u1, u2); // TODO: check if integration points ok
            break;
        case 1:
            return f(u2);
            break;
        case 2:
            return (crv.value(u2,1) * crv.value(u2,2)) / f(u2);
            break;
        default:
            {
                throw std::runtime_error("Not implemented.");
                return  0.;
            }
            break;
        }

    }
/**
 * @brief Compute full curve length
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param d  : derivate order max 2 respective to u2
 * @return T 
 */
    template <typename T, size_t dim, size_t N = N_gauss_pt>
    auto length(const Curve<T,dim> &crv, size_t d = 0) -> T
    {
        auto [u1,u2] = crv.bounds();
        return length<T,dim,N>(crv,u1,u2,d);
    }
/**
 * @brief For internal use, be cautious out of [u1,u2] function is not valid
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam N 
 * @param crv 
 * @param u1 
 * @param u2 
 * @param n 
 * @return BSCfunction<T> 
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv(const Curve<T, dim> &crv, T u1, T u2, size_t n = 30) -> BSCfunction<T>
    {
        points_vector<T, 1> u = make_range(point<T, 1>{u1}, point<T, 1>{u2}, n);

        std::vector<T> dm(n - 1);
        std::transform(
            std::execution::par,
            u.begin(), std::next(u.end(), -1), std::next(u.begin(), 1), dm.begin(),
            [&crv](const auto &u1_,const auto &u2_) {
                return length<T, dim, N>(crv, u1_[0], u2_[0]);
            });

        std::vector<T> m(n);
        m.front() = 0.;
        std::transform(dm.begin(), dm.end(), m.begin(), std::next(m.begin()),
                       [](T dm_, T sum_) {
                           return sum_ + dm_;
                       });


        return  BSCfunction<T>{ interpolate<T, 1>(u, m, fmin(3, n)) };
    }
/**
 * @brief Construct an inverse function returning the paramenter on curve corresponding to the curvilinear abscissa
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam N    Number of points used for Gauss integration of length
 * @param crv   Curve
 * @param n     Number of points to create function interpolation
 * @return BSCurve<T,1> 
 */
    template <typename T, size_t dim, size_t N = 10>
    auto abs_curv(const Curve<T, dim> &crv, size_t n = 30) -> BSCfunction<T>
    {
        auto [u1, u2] = crv.bounds();
        return abs_curv<T,dim,N>(crv,u1,u2,n);
    }
/**
 * @brief Create a list of parameters uniformly spaced on curve between 2 parameters, can raise if n_law the point number to build the curvilinear law is too small
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam N    Number of points used for Gauss integration of length
 * @param crv   Curve
 * @param u1    Start parameter
 * @param u2    End parameter
 * @param n     Number of points
 * @param n_law Points uses for abs curve law
 * @return std::list<T> 
 */
    template <typename T, size_t dim, size_t N = 10>
    auto uniform_distrib_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, size_t n_law = 30) -> std::list<T>
    {
        std::list<T> u_lst(n);
        u_lst.front() = u1;
        u_lst.back() = u2;

        auto f_u = abs_curv<T,dim,N>(crv,u1,u2,n_law);
        auto dm = f_u.bounds()[1] / ( n - T(1) );

        std::generate(
            std::next(u_lst.begin(), 1),
            std::next(u_lst.end()  ,-1),
            [&,m_=T(0.)]() mutable {
                m_ += dm;
                return f_u(m_);
                });

        if(!std::is_sorted(u_lst.begin(),u_lst.end()))
            throw std::length_error("Building abs curve fails, please refine n_law");
        return u_lst;
    }
/**
 * @brief Create a list of parameters uniformly spaced on curve, can raise if n_law the point number to build the curvilinear law is too small
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam N    Number of points used for Gauss integration of length
 * @param crv Curve
 * @param n   Number points for distribution
 * @return std::list<T> 
 */
    template <typename T, size_t dim, size_t N = 10>
    auto uniform_distrib_params(const Curve<T, dim> &crv, size_t n, size_t n_law = 30) -> std::list<T>
    {
        auto [u1, u2] = crv.bounds();
        return uniform_distrib_params(crv,u1,u2,n,n_law);
    }

    template <typename T, size_t dim>
    auto deviation_based_params(const Curve<T, dim> &crv, T u1, T u2, size_t n, T dev_max, size_t n_max_pts=5000) -> std::list<T>
    {
        std::list<T> u_lst;
        for(size_t i{}; i < n ; i++)
            u_lst.push_back( u1+(u2-u1)*i/(n-1.) );


        // refine
        bool inserted = true;
        while (inserted && u_lst.size() < n_max_pts)
        {
            auto it_u1 = u_lst.begin();
            auto it_u3 = std::next(it_u1);
            inserted = false;
            while (it_u3 != u_lst.end())
            {
                auto u_ = 0.5 * (*it_u1 + *it_u3);
                auto v1 = crv(u_) - crv(*it_u1);
                auto v2 = crv(*it_u3) - crv(*it_u1);
                auto dev_ =  norm(v1 ^ v2) / norm(v1) / norm(v2) ;

                if (dev_ > dev_max)
                {
                    inserted = true;
                    u_lst.insert(it_u3, u_);
                }

                it_u1 = it_u3;
                it_u3 = std::next(it_u1);

            }
        }
        return u_lst;
    }
/**
 * @brief 
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @param dev_max 
 * @param n_max_pts 
 * @return std::list<T> 
 */
    template <typename T, size_t dim>
    auto deviation_based_params(const Curve<T, dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000) -> std::list<T>
    {

        auto [u1, u2] = crv.bounds();
        return deviation_based_params(crv, u1, u2, n, dev_max, n_max_pts );
    }
/**
 * @brief Create a vector of points at curve's parameters positions
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam _Container 
 * @param crv Curve
 * @param u_lst Poisitions on curve
 * @return gbs::points_vector<T,dim> 
 */
    template <typename T, size_t dim, typename _Container>
    auto make_points(const Curve<T,dim> &crv,const _Container &u_lst, size_t d = 0) -> gbs::points_vector<T,dim>
    {
        points_vector<T,dim> points(u_lst.size());

        std::transform(
            std::execution::par,
            u_lst.begin(),u_lst.end(),points.begin(),
            [d,&crv](T u_){return crv(u_,d);}
        );

        return points;
    }
/**
 * @brief Uniformly spaced discretization
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @return PointArray<T,dim> 
 */
    template <typename T, size_t dim>
    auto discretize(const Curve<T,dim> &crv, size_t n) -> gbs::points_vector<T,dim>
    {
        return make_points(crv,uniform_distrib_params<T,dim>(crv,n));
    }
    /**
     * @brief Build curve's discretization with a nim/max number of points and refined to respect the maximal deviation
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv curve
     * @param n minimum number of points
     * @param dev_max maximal deviation
     * @param n_max_pts maximal number of points
     * @return gbs::points_vector<T,dim> 
     */
    template <typename T, size_t dim>
    auto discretize(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000) -> gbs::points_vector<T,dim>
    {
        auto u_lst = deviation_based_params<T, dim>(crv, n,dev_max,n_max_pts);
        // build points
        return make_points(crv,u_lst);

    }
/**
 * @brief Build curve's discretization with a nim/max number of points and refined to respect the maximal deviation, both points and corresponding paramets are returned
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param n 
 * @param dev_max 
 * @param n_max_pts 
 * @return (gbs::points_vector<T,dim>, std::vector<T>) 
 */
    template <typename T, size_t dim>
    auto discretize_with_params(const Curve<T,dim> &crv, size_t n, T dev_max, size_t n_max_pts=5000)
    {
        auto u_lst = deviation_based_params<T, dim>(crv, n,dev_max,n_max_pts);
        std::vector<T> u{ std::begin(u_lst), std::end(u_lst) };
        // build points
        return std::make_tuple( make_points(crv,u_lst) , u );
    }
/**
 * @brief Compute normalized tangential direction of the curve
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return point<T,dim> 
 */
    template <typename T, size_t dim>
    auto tangential_direction(const Curve<T,dim> &crv,T u) -> point<T,dim>
    {
        auto tg = crv(u, 1);
        return tg / norm(tg);
    }
/**
 * @brief Compute curve direction information at given parameter
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return ax1<T,dim> 
 */
    template <typename T,size_t dim>
    auto tangential_line(const Curve<T,dim> &crv,T u) -> ax1<T,dim>
    {
        return {crv(u),tangential_direction(crv,u)};
    }
/**
 * @brief Compute normalized normal direction of the curve
 * 
 * @tparam T 
 * @param crv 
 * @param u 
 * @return point<T,2> 
 */
    template <typename T>
    auto normal_direction(const Curve<T,2> &crv,T u) -> point<T,2>
    {
        auto tg = crv(u, 1);
        point<T, 2> n;
        n[0] = -tg[1];
        n[1] =  tg[0];
        return n / norm(n);
    }
/**
 * @brief Compute normalized normal direction of the curve using curvature, thus is the later is null result is +/-infinity
 * 
 * @tparam T 
 * @param crv 
 * @param u 
 * @return point<T,3> 
 */
    template <typename T>
    auto normal_direction(const Curve<T,3> &crv,T u) -> point<T,3>
    {
        auto tg = crv(u, 1);
        auto cu = crv(u, 2);
        auto n_pln = tg^cu; // normal osculating plane
        n_pln = n_pln / norm(n_pln);
        return tg^n_pln / norm(tg);
    }
/**
 * @brief Compute curve's normal direction information at given parameter
 * 
 * @tparam T 
 * @tparam dim 
 * @param crv 
 * @param u 
 * @return ax1<T,dim> 
 */
    template <typename T,size_t dim>
    auto normal_line(const Curve<T,dim> &crv,T u) -> ax1<T,dim>
    {
        return {crv(u),normal_direction(crv,u)};
    }
/**
 * @brief find the closest curve to point
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt : the point
 * @param curve_begin 
 * @param curve_end 
 * @param tol 
 * @return auto : curve's iterator
 */
    template <typename T, size_t dim>
    auto closest_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = 1e-6)
    {
        return std::min_element(
            std::execution::par,
            curve_begin, curve_end,
            [&pt, tol](const auto &c1, const auto &c2)
            {
                return extrema_curve_point(c1, pt, c1.midPoint(), tol)[1] < extrema_curve_point(c2, pt, c2.midPoint(), tol)[1];
            });
    }
/**
 * @brief  find the closest curve's pointer to point
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt : the point
 * @param curve_begin 
 * @param curve_end 
 * @param tol 
 * @return auto : curve pointer's iterator
 */
    template <typename T, size_t dim>
    auto closest_p_curve(const point<T, dim> pt, const auto &curve_begin, const auto &curve_end, T tol = 1e-6)
    {
        return std::min_element(
            std::execution::par,
            curve_begin, curve_end,
            [&pt, tol](const auto &c1, const auto &c2)
            {
                return extrema_curve_point(*c1, pt, c1->midPoint(), tol)[1] < extrema_curve_point(*c2, pt, c2->midPoint(), tol)[1];
            });
    }

    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T u1, T u2, T tol_x, nlopt::algorithm solver=nlopt::LN_COBYLA)
    {
        auto f = [&crv](const std::vector<double> &x)
        {
            return std::vector<double>{1./gbs::sq_norm(crv(x[0],2))};
        };
        auto gf = [&crv](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto f2 = crv(x[0],2);
            auto f3 = crv(x[0],3);
            return std::vector<double>{ -2 * ( f2 * f3 ) / ( ( f2 * f2 ) * ( f2 * f2 ) )};
        };
        std::vector<T> x{ 0.5 * (u1 + u2)};
        std::vector<T> lb{u1};
        std::vector<T> hb{u2};
        auto minf = gbs::solve_D_nlop(
            f,gf,
            x, lb, hb,
            tol_x*tol_x,
            solver
        );

        return std::array<T,2> { x[0], 1 / sqrt(minf)};
    }

    template <typename T, size_t dim>
    auto max_curvature_pos(const Curve<T,dim> &crv, T tol_x, nlopt::algorithm solver=nlopt::LN_COBYLA)
    {
        auto [ u1, u2 ] = crv.bounds();
        return max_curvature_pos(crv, u1, u2, tol_x, solver);

    }
} // namespace gbs