#pragma once

#include <gbs/execution.h>
#include <Eigen/Dense>
#include <gbs/bscurve.h>

#ifdef GBS_USE_MODULES
    import knots_functions;
    import basis_functions;
    import vecop;
#else
    #include "vecop.ixx"
    #include "knotsfunctions.ixx"
    #include "basisfunctions.ixx"
#endif
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

template <typename T, size_t dim>
auto build_3pt_tg_vec(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u)
{
    // Nurbs book p 386, Bessel's method
    auto d_u = delta(u);
    auto q   = delta(pts);
    auto d   = q/d_u;

    std::vector<T> alpha(u.size() - 2);
    std::transform(
        GBS_PAR_EXEC
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

    D.front() = T(2.)*q.front() - *std::next(D.begin());
    D.back()  = T(2.)*q.back()  - *std::next(D.end(),-1);

    return D;
}

template <typename T, size_t dim>
auto build_3pt_tg_dir(const std::vector<std::array<T, dim>> &pts,const std::vector<T> &u)
{
    if(pts.size()<3)
    {
        throw std::invalid_argument("build_3pt_tg_dir: at least 3 pts needs to be provided");
    }
    // Nurbs book p 386, Bessel's method
    auto d_u = delta(u);
    auto q   = delta(pts);

    std::vector<T> alpha(u.size() - 2);
    std::transform(
        GBS_PAR_EXEC
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

template <typename T,size_t dim,size_t nc>
    auto build_poles(const std::vector<constrType<T,dim,nc> > &Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> std::vector<std::array<T,dim> >
{
    auto n_pt = Q.size();
    auto n_poles = int(Q.size() * nc);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);


    build_poles_matrix<T,nc>(k_flat,u,deg,n_poles,N);
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
    // Banded one-pass assembly: each constraint row has only deg+1 non-zero
    // entries (the basis support), filled in a single basis pass.
    N.setZero();
    for (size_t i = 0; i < n_poles; ++i)
        fill_basis_row(N, Eigen::Index(i), Q[i].u, k_flat, deg, Q[i].d);
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

template <typename T,size_t dim>
auto build_poles(const points_vector<T,dim> &Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> points_vector<T,dim>
{
    std::vector<constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return constrType<T, dim, 1>{pt_};});
    return build_poles(Q_,k_flat,u,deg);
}

template <typename T, size_t dim,typename _FwdIt>
auto get_constrain(const _FwdIt &begin, const _FwdIt &end, size_t order) -> std::vector<std::array<T, dim>> 
{
    std::vector<std::array<T, dim>> pts(end-begin);
    std::transform(begin,end, pts.begin(), [&order](const auto &q_) { return q_[order]; });
    return pts;
}

template <typename T, size_t dim, size_t nc>
auto get_constrain(const std::vector<constrType<T, dim, nc>> &Q, size_t order) -> std::vector<std::array<T, dim>> 
{
    return get_constrain<T,dim>(Q.begin(),Q.end(),order);
}

template <typename T, size_t dim, size_t nc>
auto interpolate(const std::vector<constrType<T, dim, nc>> &Q, const std::vector<T> &u, bool add_computed=false)->BSCurve<T,dim>
{

    auto build_crv = [&u](size_t nc_, const auto &Q)
    {
        size_t p = 2 * nc_ - 1;

        std::vector<size_t> m(Q.size());
        m.front() = p + 1;
        std::fill(++m.begin(), --m.end(), nc_);
        m.back() = p + 1;

        auto k_flat = flat_knots(u, m);

        auto poles = build_poles(Q, k_flat, u, p);

        return BSCurve<T, dim>(poles, k_flat, p);
    };

    if(add_computed)
    {
        auto n = Q.size();
        points_vector<T,dim> v(n);
        std::transform(
            Q.begin(), Q.end(),
            v.begin(),
            [](const constrType<T, dim, nc> &q){return q[nc-1];}
        );
        auto w = build_3pt_tg_vec(v,u);
        std::vector<constrType<T, dim, nc+1>> Q_(Q.size());
        for( size_t i{} ; i < nc ; i++)
        {
            for(size_t j{}; j < n ; j++)
            {
                Q_[j][i] = Q[j][i];
            }
        }
        for(size_t j{}; j < n ; j++)
        {
            Q_[j][nc] = w[j];
        }
        // Natural BSCurve if nc == 2
        Q_[0][nc]  = point<T,dim>{};
        Q_[n-1][nc] = point<T,dim>{};
        return build_crv(nc+1, Q_);
    }
    else
    {
        return build_crv(nc, Q);
    }
}

template <typename T, size_t dim, size_t nc>
auto interpolate(const std::vector<constrType<T, dim, nc>> &Q, KnotsCalcMode mode, bool add_computed=false)->BSCurve<T,dim>
{
    
    auto pts = get_constrain(Q,0);
    auto u = curve_parametrization(pts, mode);

    return interpolate<T,dim,nc>(Q,u, add_computed);
}

// interp cn
//Nurbs book p365
template <typename T, size_t dim>
auto interpolate(const std::vector<constrType<T, dim, 1>> &Q,const std::vector<T> &u, size_t p) -> BSCurve<T, dim>
{ 
    auto k_flat = build_simple_mult_flat_knots<T>(u,p);
    auto poles = build_poles<T,dim,1>(Q, k_flat, u, p);
    return BSCurve<T,dim>(poles, k_flat, p);
}
template <typename T, size_t dim>
auto interpolate(const std::vector<constrType<T, dim, 1>> &Q, size_t p, KnotsCalcMode mode ) -> BSCurve<T, dim>
{

    if(p>=Q.size()) throw std::domain_error("Degree must be strictly inferior to points number");
    auto pts = get_constrain(Q, 0);
    auto u = curve_parametrization(pts, mode, true);
    return interpolate(Q,u,p);
}

template <typename T, size_t dim>
auto interpolate(const points_vector<T,dim> &Q, size_t p, KnotsCalcMode mode ) -> BSCurve<T, dim>
{
    std::vector<constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return constrType<T, dim, 1>{pt_};});
    return interpolate<T,dim>(Q_,p,mode);
}

template <typename T, size_t dim>
auto interpolate(const points_vector<T,dim> &Q,const std::vector<T> &u, size_t p ) -> BSCurve<T, dim>
{ 
    std::vector<constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return constrType<T, dim, 1>{pt_};});
    return interpolate<T,dim>(Q_,u,p);
}

template <typename T>
auto interpolate(const std::vector<T> &Q,const std::vector<T> &u, size_t p) -> BSCfunction<T>
{ 
    auto k_flat = build_simple_mult_flat_knots<T>(u,p);
    std::vector<constrType<T, 1, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &v_){return constrType<T, 1, 1>{{v_}};});
    return BSCfunction<T>{ interpolate<T,1>(Q_,u,p) };
}

template < typename T, size_t dim>
using bsc_constraint = std::tuple<T,point<T,dim>,size_t>;
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
auto interpolate(const bsc_bound<T,dim> &pt_begin,const bsc_bound<T,dim> &pt_end,const std::vector<bsc_constraint<T,dim>> &cstr_lst,size_t p)
{
    auto nc= cstr_lst.size();
    auto n = nc + 2;

    if(n < p + 1)
    {
        throw std::invalid_argument{"constraints number shall be greater than degree - 1"};
    }

    std::vector<T> u{std::get<0>(pt_begin),std::get<0>(pt_end)};
    for(const auto &cstr : cstr_lst)
    {
        insert_ordered(u, std::get<0>(cstr));
    }

    auto k_flat = build_simple_mult_flat_knots(u,p);

    // Banded one-pass assembly: each row (begin point, the derivative
    // constraints, end point) has only p+1 non-zero entries.
    MatrixX<T> N(n,n);
    N.setZero();
    fill_basis_row(N, Eigen::Index(0), std::get<0>(pt_begin), k_flat, p, size_t{0});
    for (size_t i = 0; i < nc; ++i)
        fill_basis_row(N, Eigen::Index(i + 1), std::get<0>(cstr_lst[i]), k_flat, p, std::get<2>(cstr_lst[i]));
    fill_basis_row(N, Eigen::Index(n - 1), std::get<0>(pt_end), k_flat, p, size_t{0});

    auto N_inv = N.partialPivLu();
    // auto N_inv = N.colPivHouseholderQr();
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


} // namespace gbs