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

template <typename T,size_t dim>
auto build_poles(const points_vector<T,dim> &Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> points_vector<T,dim>
{
    std::vector<gbs::constrType<T, dim, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &pt_){return gbs::constrType<T, dim, 1>{pt_};});
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
auto get_constrain(const std::vector<gbs::constrType<T, dim, nc>> &Q, size_t order) -> std::vector<std::array<T, dim>> 
{
    return get_constrain<T,dim>(Q.begin(),Q.end(),order);
}

template <typename T, size_t dim, size_t nc>
auto interpolate(const std::vector<gbs::constrType<T, dim, nc>> &Q, const std::vector<T> &u, bool add_computed=false)->gbs::BSCurve<T,dim>
{

    auto build_crv = [&u](size_t nc_, const auto &Q)
    {
        size_t p = 2 * nc_ - 1;

        std::vector<size_t> m(Q.size());
        m.front() = p + 1;
        std::fill(++m.begin(), --m.end(), nc_);
        m.back() = p + 1;

        auto k_flat = gbs::flat_knots(u, m);

        auto poles = gbs::build_poles(Q, k_flat, u, p);

        return gbs::BSCurve<T, dim>(poles, k_flat, p);
    };

    if(add_computed)
    {
        auto n = Q.size();
        points_vector<T,dim> v(n);
        std::transform(
            Q.begin(), Q.end(),
            v.begin(),
            [](const gbs::constrType<T, dim, nc> &q){return q[nc-1];}
        );
        auto w = build_3pt_tg_vec(v,u);
        std::vector<gbs::constrType<T, dim, nc+1>> Q_(Q.size());
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
auto interpolate(const std::vector<gbs::constrType<T, dim, nc>> &Q, gbs::KnotsCalcMode mode, bool add_computed=false)->gbs::BSCurve<T,dim>
{
    
    auto pts = get_constrain(Q,0);
    auto u = gbs::curve_parametrization(pts, mode);

    return interpolate<T,dim,nc>(Q,u, add_computed);
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

template <typename T>
auto interpolate(const std::vector<T> &Q,const std::vector<T> &u, size_t p) -> gbs::BSCfunction<T>
{ 
    auto k_flat = build_simple_mult_flat_knots<T>(u,p);
    std::vector<gbs::constrType<T, 1, 1>> Q_(Q.size());
    std::transform(Q.begin(),Q.end(),Q_.begin(),[](const auto &v_){return gbs::constrType<T, 1, 1>{{v_}};});
    return gbs::BSCfunction<T>{ interpolate<T,1>(Q_,u,p) };
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

    MatrixX<T> N(n,n);
    for (auto j = 0; j < n; j++)
    {
        N(0, j) = basis_function(std::get<0>(pt_begin), j, p,0, k_flat);
    }

    for (int i = 0; i < nc; i++) // Eigen::ColMajor is default
    {
        auto u = std::get<0>(cstr_lst[i]);
        auto d = std::get<2>(cstr_lst[i]);
        for (auto j = 0; j < n; j++)
        {
            N(i + 1, j) = basis_function(u, j, p, d, k_flat);
        }
    }


    for (auto j = 0; j < n; j++)
    {
        N(n-1, j) = basis_function(std::get<0>(pt_end), j, p,0, k_flat);
    }

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