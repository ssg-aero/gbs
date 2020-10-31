#pragma once
#include <gbs/basisfunctions.h>
#include <gbs/knotsfunctions.h>
#include <gbs/bscurve.h>

namespace gbs
{
template <typename T,size_t dim,size_t nc>
    using constrType = std::array<std::array<T,dim>,nc >;

template <typename T,size_t dim,size_t nc>
    auto build_poles(const std::vector<constrType<T,dim,nc> > Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> std::vector<std::array<T,dim> >
{
    auto n_pt = Q.size();
    auto n_poles = int(Q.size() * nc);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);


    build_poles_matix<T,nc>(k_flat,u,deg,n_poles,N);
    auto N_inv = N.partialPivLu(); //TODO solve banded system

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
auto interpolate(const std::vector<gbs::constrType<T, dim, nc>> &Q, gbs::KnotsCalcMode mode)->gbs::BSCurve<T,dim>
{
    
    auto pts = get_constrain(Q,0);

    size_t p = 2 * nc - 1;
    
    std::vector<size_t> m(Q.size());
    m.front()=p+1;
    std::fill(++m.begin(),--m.end(),nc);
    m.back()=p+1;

    auto u = gbs::curve_parametrization(pts, mode);
    auto k_flat = gbs::flat_knots(u, m);

    auto poles = gbs::build_poles(Q, k_flat, u, p);

    return gbs::BSCurve<T,dim>(poles, k_flat, p);
}

template <typename T>
auto build_simple_mult_flat_knots(const std::vector<T> &u, size_t n, size_t p) -> std::vector<T>
{

    auto nk = n + p + 1;
    
    std::vector<T> k_flat(nk);
    std::fill(k_flat.begin(), std::next(k_flat.begin(), p), T(0.));
    std::fill(std::next(k_flat.begin(), nk - 1 - p), k_flat.end(), u.back()-u.front());


    for (int j = 1; j < n - p; j++) // TODO use std algo
    {
        // k_flat[j + p] = 0;
        // for (int i = j; i <= j + p - 1; i++)
        // {
        //     k_flat[j + p] += u[i] / p;
        // }
        k_flat[j + p] = j / T(n-p);
    }

    return k_flat;
} 

// interp cn
//Nurbs book p365
template <typename T, size_t dim>
auto interpolate(const std::vector<gbs::constrType<T, dim, 1>> &Q, size_t p, gbs::KnotsCalcMode mode ) -> gbs::BSCurve<T, dim>
{

    if(p>=Q.size()) throw std::domain_error("Degree must be strictly inferior to points number");
    auto pts = get_constrain(Q, 0);
    auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    auto k_flat = build_simple_mult_flat_knots<T>(u,pts.size(),p);
    auto poles = gbs::build_poles<T,dim,1>(Q, k_flat, u, p);
    return gbs::BSCurve<T,dim>(poles, k_flat, p);
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