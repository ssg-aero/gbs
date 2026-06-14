#pragma once
#include <gbs/bscinterp.h>
#include <gbs/bssurf.h>
#include <exception>
namespace gbs
{

    // ( u, v, x, du, dv)
    template < typename T, size_t dim>
    using bss_constraint = std::tuple<T,T,point<T,dim>,size_t,size_t>; 

    template < typename T, size_t dim>
    auto get_constraints_bounds(const std::vector<bss_constraint<T, dim>> &Q)
    {
        auto u_min = std::get<0>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<0>(Q1) < std::get<0>(Q2);}));
        auto u_max = std::get<0>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<0>(Q1) > std::get<0>(Q2);}));
        auto v_min = std::get<1>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<1>(Q1) < std::get<1>(Q2);}));
        auto v_max = std::get<1>(*std::min_element(Q.begin(),Q.end(),[](auto Q1, auto Q2){return std::get<1>(Q1) > std::get<1>(Q2);}));
        return std::array<T,4>{u_min, u_max, v_min, v_max};
    }

    template <typename Container>
    auto extract_U(size_t i_V,const Container &Q, size_t n_col) -> Container
    {
        if (Q.size() % n_col)
        {
            throw std::length_error("size error");
        }
        Container row( n_col );
        std::copy(std::next(Q.begin(),i_V*n_col),std::next(Q.begin(),(i_V+1)*n_col),row.begin());
        return row;
    }

    template <typename Container>
    auto extract_V(size_t i_U,const Container &Q, size_t n_col) -> Container
    {
        if (Q.size() % n_col)
        {
            throw std::length_error("size error");
        }
        auto n_row = Q.size() / n_col;
        Container col( n_row );
        for(auto j = 0 ; j < n_row; j++)
        {
            col[j] = Q[i_U+n_col*j];
        }
        return col;
    }

    template <typename T>
    auto fill_poles_matrix(MatrixX<T> &N,const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q)
    {
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = k_flat_u.size() - (p+1); 
        auto n_poles_v = k_flat_v.size() - (q+1); 
        size_t iRow, iCol;
        for (size_t iv = 0; iv < n_params_v; iv++)
        {
            for (size_t iu = 0; iu < n_params_u; iu++)
            {
                iRow = iu + n_params_u * iv;
                for (size_t j = 0; j < n_poles_v; j++)
                {
                    auto Nv = basis_function(v[iv], std::next(k_flat_v.begin(), j), q, k_flat_v.end());
                    for (size_t i = 0; i < n_poles_u; i++)
                    {
                        iCol = i + n_poles_u * j;
                        N(iRow, iCol) = basis_function(u[iu], std::next(k_flat_u.begin(), i), p, k_flat_u.end()) * Nv;
                    }
                }
            }
        }
    }

    template <typename T, size_t dim>
    auto build_poles(const auto &N_inv, const std::vector<std::array<T, dim>> &Q)
    {
        auto n_poles = N_inv.cols();
        auto n_pt = Q.size();
        VectorX<T> b(n_pt);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt; i++)
            {
                b(i) = Q[i][d]; //Pas top au niveau de la localisation mémoire
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
    auto build_poles(const std::vector<std::array<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        auto n_poles = Q.size();
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = n_params_u;
        auto n_poles_v = n_params_v;
        if (n_poles != n_params_u * n_params_v)
        {
            throw std::length_error("size error");
        }

        // Tensor-product (grid) interpolation: the (nu*nv)x(nu*nv) collocation
        // matrix is the Kronecker product Mv (x) Mu, so the system separates as
        //   Mu * P * Mv^T = Q   =>   P = Mu^-1 * Q * Mv^-T,
        // with Mu (nu x nu) and Mv (nv x nv) the 1-D collocation matrices. We
        // factor the two small matrices once and solve as two sequences of 1-D
        // solves instead of one dense O((nu*nv)^3) LU (issue #34). Poles and Q
        // are stored U-first: index i + n_poles_u * j.
        MatrixX<T> Mu(n_poles_u, n_poles_u);
        Mu.setZero();
        for (size_t iu = 0; iu < n_params_u; ++iu)
            fill_basis_row(Mu, Eigen::Index(iu), u[iu], k_flat_u, p, size_t{0});

        MatrixX<T> Mv(n_poles_v, n_poles_v);
        Mv.setZero();
        for (size_t iv = 0; iv < n_params_v; ++iv)
            fill_basis_row(Mv, Eigen::Index(iv), v[iv], k_flat_v, q, size_t{0});

        auto lu_u = Mu.partialPivLu();
        auto lu_v = Mv.partialPivLu();

        std::vector<std::array<T, dim>> poles(n_poles);
        MatrixX<T> Qd(n_poles_u, n_poles_v);
        for (int d = 0; d < dim; ++d)
        {
            for (size_t iv = 0; iv < n_params_v; ++iv)
                for (size_t iu = 0; iu < n_params_u; ++iu)
                    Qd(iu, iv) = Q[iu + n_params_u * iv][d];

            MatrixX<T> A  = lu_u.solve(Qd);              // Mu^-1 Q          (nu x nv)
            MatrixX<T> Pt = lu_v.solve(A.transpose());   // P^T = Mv^-1 A^T  (nv x nu)

            for (size_t j = 0; j < n_poles_v; ++j)
                for (size_t i = 0; i < n_poles_u; ++i)
                    poles[i + n_poles_u * j][d] = Pt(j, i);
        }

        return poles;
    }

    template <typename T, size_t dim>
    auto interpolate(const std::vector<std::array<T, dim>> &Q, size_t n_poles_v, size_t p, size_t q, KnotsCalcMode mode) -> BSSurface<T, dim>
    {

        if (Q.size() % n_poles_v)
        {
            throw std::length_error("size error");
        }

        size_t n_poles_u = Q.size() / n_poles_v;

        auto avg_p = [&](auto extract_f, size_t ni, size_t nj) mutable {
            
            std::vector<T> params(nj , T(0));
            for (auto i = 0; i < ni; i++)
            {
                auto pts = extract_f(i);
                auto params_i = curve_parametrization(pts, mode, true);
                std::transform(
                    GBS_PAR_EXEC
                    params_i.begin(),
                    params_i.end(),
                    params.begin(),
                    params.begin(),
                    [&](const auto p_new, const auto p_) {
                        return p_ + p_new / ni;
                    });
            }
            return params;
        };

        auto u = avg_p(
            [&](size_t i) {
                return extract_U(i, Q, n_poles_u);
            },
            n_poles_v, n_poles_u);
        adimension(u);
        auto ku = build_simple_mult_flat_knots<T>(u, p);
        auto v = avg_p(
            [&](size_t i) {
                return extract_V(i, Q, n_poles_u);
            },
            n_poles_u, n_poles_v);
        adimension(v);
        auto kv = build_simple_mult_flat_knots<T>(v, q);
        
        return BSSurface<T, dim>(build_poles(Q,ku,kv,u,v,p,q),ku,kv,p,q);

    }
    /**
     * @brief Internal use for template <typename T, size_t dim>
     * auto interpolate(const std::vector<bss_constraint<T,dim>> &Q, size_t n_polesv, size_t p, size_t q, const std::array<T,4> &bounds = {0.,1.,0.,1.}, const std::array<size_t,2> &mults = {1,1} ) -> BSSurface<T, dim>
     * 
     * Use it at your own risks!
     * 
     * @tparam T 
     * @tparam dim 
     * @param Q 
     * @param k_flat_u 
     * @param k_flat_v 
     * @param p 
     * @param q 
     * @return std::vector<std::array<T, dim>> 
     */
    template <typename T, size_t dim>
    auto build_poles(const std::vector<bss_constraint<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        // TODO sort Q by contrain order to get a block systen
        auto n_poles = Q.size();
        auto n = k_flat_u.size() - p - 1;
        auto m = k_flat_v.size() - q - 1;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);
        // Banded one-pass assembly: each constraint row is the outer product of
        // the (p+1) non-zero u-basis derivatives and the (q+1) non-zero v-basis
        // derivatives, evaluated with one basis pass per direction instead of
        // n*m recursive basis_function pairs.
        N.setZero();
        for (size_t k = 0; k < n_poles; ++k)
        {
            const auto &[u, v, x, du, dv] = Q[k];
            if (du > p || dv > q)
                continue; // identically-zero row
            const size_t span_u = std::min<size_t>(find_span(n, p, u, k_flat_u) - k_flat_u.begin(), n - 1);
            const size_t span_v = std::min<size_t>(find_span(m, q, v, k_flat_v) - k_flat_v.begin(), m - 1);
            const size_t i_min = (span_u >= p) ? span_u - p : 0;
            const size_t j_min = (span_v >= q) ? span_v - q : 0;
            const size_t ii_max = (i_min + p < n) ? p : n - 1 - i_min; // clamp bands to poles
            const size_t jj_max = (j_min + q < m) ? q : m - 1 - j_min;
            if (p <= bspline_stack_max_degree && q <= bspline_stack_max_degree)
            {
                T Nu[bspline_stack_max_degree + 1];
                T Nv[bspline_stack_max_degree + 1];
                basis_ders(span_u, p, du, u, k_flat_u, Nu);
                basis_ders(span_v, q, dv, v, k_flat_v, Nv);
                for (size_t jj = 0; jj <= jj_max; ++jj)
                    for (size_t ii = 0; ii <= ii_max; ++ii)
                        N(Eigen::Index(k), Eigen::Index((i_min + ii) + n * (j_min + jj))) = Nu[ii] * Nv[jj];
            }
            else
            {
                // Pathological degree: recursive fallback (still banded).
                for (size_t jj = 0; jj <= jj_max; ++jj)
                    for (size_t ii = 0; ii <= ii_max; ++ii)
                        N(Eigen::Index(k), Eigen::Index((i_min + ii) + n * (j_min + jj))) =
                            basis_function(u, i_min + ii, p, du, k_flat_u) *
                            basis_function(v, j_min + jj, q, dv, k_flat_v);
            }
        }

        auto N_inv = N.partialPivLu();

        VectorX<T> b(n_poles);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_poles; i++)
            {
                b(i) = std::get<2>(Q[i])[d];
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return poles;
    }
    /**
     * @brief Build BSSurface<T,dim> satisfying the constraints vector provided: 
     *          * constraints parameters u,v are within bounds
     *          * Q.size() % n_polesv == 0 ( constraints have to be properly shaped aka Q.size() = n_polesv* n_polesv)
     *          * (n_polesu-p-1) % mu == 0 and (n_polesv-q-1) % mv == 0
     * 
     * @tparam T 
     * @tparam dim 
     * @param Q 
     * @param n_polesv 
     * @param p : degree U
     * @param q : degree V
     * @param bounds 
     * @param mults 
     * @return BSSurface<T, dim> 
     */
    template <typename T, size_t dim>
    auto interpolate(const std::vector<bss_constraint<T,dim>> &Q, size_t n_polesv, size_t p, size_t q, const std::array<T,4> &bounds = {0.,1.,0.,1.}, const std::array<size_t,2> &mults = {1,1} ) -> BSSurface<T, dim>
    {
        auto [u1, u2, v1, v2] = bounds;
        auto [mu, mv        ] = mults;

        auto [u_min, u_max, v_min, v_max] = get_constraints_bounds(Q);
        if(u_min < u1 || u_max > u2 || v_min < v1 || v_max > v2)
        {
            throw std::invalid_argument("constains must be specified within bounds");
        }

        if( Q.size() % n_polesv !=0)
        {
            throw std::length_error("interpolate, Q.size() = n_polesv * n_polesu required");
        }
        auto n_polesu = Q.size() / n_polesv;
        auto ku = build_mult_flat_knots<T>(u1, u2, 2 + (n_polesu-p-1)/mu, p, mu);
        auto kv = build_mult_flat_knots<T>(v1, v2, 2 + (n_polesv-q-1)/mv, q, mv);
        if(ku.size() - p - 1 != n_polesu)
        {
            throw std::invalid_argument("(n_polesu-p-1) % mu == 0 required");
        }
        if(kv.size() - q - 1 != n_polesv)
        {
            throw std::invalid_argument("(n_polesv-q-1) % mv == 0 required");
        }
        return BSSurface<T, dim>(build_poles(Q,ku,kv,p,q),ku,kv,p,q);

    }
} // namespace gbs