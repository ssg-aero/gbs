#pragma once
#include <gbs/basisfunctions.h>
#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>

namespace gbs
{
template <typename T,size_t dim,size_t nc>
    using constrType = std::array<std::array<T,dim>,nc >;

template <typename T, size_t dim>
struct constrPoint
{
    std::array<T, dim> v = std::array<T, dim>{0};
    T u;
    size_t d;
};

template <typename T,size_t dim,size_t nc>
    auto build_poles(const std::vector<constrType<T,dim,nc> > &Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> std::vector<std::array<T,dim> >
{
    auto n_pt = Q.size();
    auto n_poles = int(Q.size() * nc);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);


    build_poles_matix<T,nc>(k_flat,u,deg,n_poles,N);
    auto N_inv = N.partialPivLu(); //TODO solve block system
    // auto N_inv = N.colPivHouseholderQr(); //TODO solve block system

    // std::cout << N << std::endl;

    std::vector<std::array<T, dim>> poles(n_poles);

    VectorX<T> b(n_poles);
    for (int d = 0; d < dim; d++)
    {
        for (int i = 0; i < n_pt; i++)
        {
            for (int deriv = 0; deriv < nc; deriv++)
            {
                b(nc * i + deriv) = Q[i][deriv][d];//Pas top au niveau de la localisation mémoire
            }
        }

        auto x = N_inv.solve(b);

        for (int i = 0; i < n_poles; i++)
        {
            poles[i][d] = x(i);
        }
    }

    return poles;
}

template <typename T, size_t dim>
auto build_poles(const std::vector<constrPoint<T, dim>> &Q, const std::vector<T> &k_flat, size_t deg) -> std::vector<std::array<T, dim>>
{
    //TODO sort Q by contrain order to get a block systen
    auto n_poles = Q.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);
    for (int i = 0; i < n_poles; i++) // Eigen::ColMajor is default
    {
        for (int j = 0; j < n_poles; j++)
        {
            N(i, j) = gbs::basis_function(Q[i].u, j, deg, Q[i].d, k_flat);
        }
    }
    auto N_inv = N.partialPivLu();
    VectorX<T> b(n_poles);
    std::vector<std::array<T, dim>> poles(n_poles);
    for (int d = 0; d < dim; d++)
    {
        for (int i = 0; i < n_poles; i++)
        {
            b(i) = Q[i].v[d];
        }

        auto x = N_inv.solve(b);

        for (int i = 0; i < n_poles; i++)
        {
            poles[i][d] = x(i);
        }
    }

    return poles;
}

template <typename T, size_t dim,typename _FwdIt>
auto get_constrain(const _FwdIt &begin, const _FwdIt &end, size_t order) -> std::vector<std::array<T, dim>> 
{
    std::vector<std::array<T, dim>> pts(end-begin);
    std::transform(begin,end, pts.begin(), [&order](const auto &q_) { return q_[order]; });
    return pts;
}

template <typename T, size_t dim, size_t nc>
auto get_constrain(const std::vector<gbs::constrType<T, dim, nc>> &Q, size_t order) -> std::vector<std::array<T, dim>> 
{
    return get_constrain<T,dim>(Q.begin(),Q.end(),order);
}

template <typename T, size_t dim, size_t nc>
auto interpolate(const std::vector<gbs::constrType<T, dim, nc>> &Q, const std::vector<T> &u)->gbs::BSCurve<T,dim>
{
    size_t p = 2 * nc - 1;
    
    std::vector<size_t> m(Q.size());
    m.front()=p+1;
    std::fill(++m.begin(),--m.end(),nc);
    m.back()=p+1;

    auto k_flat = gbs::flat_knots(u, m);

    auto poles = gbs::build_poles(Q, k_flat, u, p);

    return gbs::BSCurve<T,dim>(poles, k_flat, p);
}

template <typename T, size_t dim, size_t nc>
auto interpolate(const std::vector<gbs::constrType<T, dim, nc>> &Q, gbs::KnotsCalcMode mode)->gbs::BSCurve<T,dim>
{
    
    auto pts = get_constrain(Q,0);
    auto u = gbs::curve_parametrization(pts, mode);

    return interpolate<T,dim,nc>(Q,u);
}

template <typename T>
auto build_simple_mult_flat_knots(const std::vector<T> &u, size_t p) -> std::vector<T>
{

    auto n = u.size();
    auto nk = n + p + 1;
    auto u1 = u.front();
    auto u2 = u.back();

    std::vector<T> k_flat(nk);
    std::fill(k_flat.begin(), std::next(k_flat.begin(), p), u1);
    std::fill(std::next(k_flat.begin(), nk - 1 - p), k_flat.end(), u2);

    auto delta_ = u2 - u1;

    for (int j = 1; j < n - p; j++) // TODO use std algo
    {
        k_flat[j + p] = 0;
        for (int i = j; i <= j + p - 1; i++)
        {
            k_flat[j + p] += u[i] / p;
        }

        // k_flat[j + p] = j / T(n - p) * delta_ + u1;
    }

    return k_flat;
} 

template <typename T>
auto build_simple_mult_flat_knots(T u1, T u2, size_t n, size_t p) -> std::vector<T>
{

    auto nk = n + p + 1;
    
    std::vector<T> k_flat(nk);
    std::fill(k_flat.begin(), std::next(k_flat.begin(), p), u1);
    std::fill(std::next(k_flat.begin(), nk - 1 - p), k_flat.end(), u2);
    auto delta_ = u2 - u1;

    for (int j = 1; j < n - p; j++)
    {
        k_flat[j + p] = u1 + delta_ * j / T(n - p);
    }

    return k_flat;
} 

// interp cn
//Nurbs book p365
template <typename T, size_t dim>
auto interpolate(const std::vector<gbs::constrType<T, dim, 1>> &Q,const std::vector<T> &u, size_t p) -> gbs::BSCurve<T, dim>
{ 
    auto k_flat = build_simple_mult_flat_knots<T>(u,p);
    auto poles = gbs::build_poles<T,dim,1>(Q, k_flat, u, p);
    return gbs::BSCurve<T,dim>(poles, k_flat, p);
}
template <typename T, size_t dim>
auto interpolate(const std::vector<gbs::constrType<T, dim, 1>> &Q, size_t p, gbs::KnotsCalcMode mode ) -> gbs::BSCurve<T, dim>
{

    if(p>=Q.size()) throw std::domain_error("Degree must be strictly inferior to points number");
    auto pts = get_constrain(Q, 0);
    auto u = gbs::curve_parametrization(pts, mode, true);
    return interpolate(Q,u,p);
}

template <typename T, size_t dim>
auto interpolate(const points_vector<T,dim> &Q, size_t p, gbs::KnotsCalcMode mode ) -> gbs::BSCurve<T, dim>
{
    std::vector<gbs::constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return gbs::constrType<T, dim, 1>{pt_};});
    return interpolate<T,dim>(Q_,p,mode);
}

template <typename T, size_t dim>
auto interpolate(const points_vector<T,dim> &Q,const std::vector<T> &u, size_t p ) -> gbs::BSCurve<T, dim>
{ 
    std::vector<gbs::constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return gbs::constrType<T, dim, 1>{pt_};});
    return interpolate<T,dim>(Q_,u,p);
}

template < typename T, size_t dim>
using bsc_constrain = std::tuple<T,point<T,dim>,size_t>;
template < typename T, size_t dim>
using bsc_bound = std::pair<T,point<T,dim>>;

/**
 * @brief Most general interpolation
 * 
 * @tparam T 
 * @tparam dim 
 * @param pt_begin 
 * @param pt_end 
 * @param cstr_lst 
 * @param p 
 * @return auto 
 */
template < typename T, size_t dim>
auto interpolate(const bsc_bound<T,dim> &pt_begin,const bsc_bound<T,dim> &pt_end,const std::vector<bsc_constrain<T,dim>> &cstr_lst,size_t p)
{
    auto nc= cstr_lst.size();
    auto n = nc + 2;

    assert(n >= p + 1);

    auto k_flat = build_simple_mult_flat_knots(std::get<0>(pt_begin),std::get<0>(pt_end),n,p);

    MatrixX<T> N(n,n);
    for (auto j = 0; j < n; j++)
    {
        N(0, j) = basis_function(std::get<0>(pt_begin), j, p,0, k_flat);
    }

    for (int i = 0; i < nc; i++) // Eigen::ColMajor is default
    {
        for (auto j = 0; j < n; j++)
        {
            N(i + 1, j) = basis_function(std::get<0>(cstr_lst[i]), j, p, std::get<2>(cstr_lst[i]), k_flat);
        }
    }


    for (auto j = 0; j < n; j++)
    {
        N(n-1, j) = basis_function(std::get<0>(pt_end), j, p,0, k_flat);
    }

    auto N_inv = N.partialPivLu();
    VectorX<T> b(n);
    std::vector<std::array<T, dim>> poles(n);
    for (int d = 0; d < dim; d++)
    {
        b(0) = std::get<1>(pt_begin)[d];
        for (int i = 0; i < nc; i++)
        {
            b(i+1) = std::get<1>(cstr_lst[i])[d];
        }
        b(n-1) = std::get<1>(pt_end)[d];

        auto x = N_inv.solve(b);

        for (int i = 0; i < n; i++)
        {
            poles[i][d] = x(i);
        }
    }

    return BSCurve<T,dim>(poles, k_flat, p);

}


template <typename T, size_t dim>
auto build_3pt_tg_vec(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u)
{
    // Nurbs book p 386, Bessel's method
    auto d_u = delta(u);
    auto q   = delta(pts);
    auto d   = q/d_u;

    std::vector<T> alpha(u.size() - 2);
    std::transform(
        std::execution::par,
        d_u.begin(), std::next(d_u.end(),-1),
        std::next(d_u.begin()),
        alpha.begin(),
        [](const auto &d_u1, const auto &d_u2) { return d_u1 / (d_u2+d_u1); });

    std::vector<std::array<T, dim>> D(u.size());

    blend(d,alpha,D,
        d.begin(),std::next(d.end(),-1),
        alpha.begin(),
        std::next(D.begin())
        );

    D.front() = 2.*q.front() - *std::next(D.begin());
    D.back()  = 2.*q.back()  - *std::next(D.end(),-1);

    return D;
}

template <typename T, size_t dim>
auto build_3pt_tg_dir(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u)
{
    // Nurbs book p 386, Bessel's method
    auto d_u = delta(u);
    auto q   = delta(pts);

    std::vector<T> alpha(u.size() - 2);
    std::transform(
        std::execution::par,
        d_u.begin(), std::next(d_u.end(),-1),
        std::next(d_u.begin()),
        alpha.begin(),
        [](const auto &d_u1, const auto &d_u2) { return d_u1 / (d_u2+d_u1); });

    std::vector<std::array<T, dim>> D(u.size());

    blend(q,alpha,D,
        q.begin(),std::next(q.end(),-1),
        alpha.begin(),
        std::next(D.begin())
        );

    D.front() = 2.*q.front() - *std::next(D.begin());
    D.back()  = 2.*q.back()  - *std::next(D.end(),-1);

    adim(D);

    return D;
}
} // namespace gbs