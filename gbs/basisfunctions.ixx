#ifdef GBS_USE_MODULES
    module;
#else
    #pragma once
#endif
#include <gbs/gbslib.h>
#include  <vector>
#include  <utility>
#include <gbs/execution.h>
#include  <algorithm>
#include  <cmath>

#ifdef GBS_USE_MODULES
    export module basis_functions;

    import vecop;
    import math;
#else
    #include "vecop.ixx"
    #include "math.ixx"
#endif

#ifdef GBS_USE_MODULES
    export namespace gbs
#else
    namespace gbs
#endif
{
    /**
     * @brief Basis function used to compute BSpline
     *
     * @tparam InIt
     * @tparam T
     * @param u  Parameter on BSpline object
     * @param it Flat knots iterator
     * @param p  BSpline object degree
     * @param _Last Container's end
     * @return T
     */
    template <std::input_iterator InIt, std::floating_point T>
    T basis_function(T u, const InIt &it, size_t p, const InIt &last)
    {
        auto u_last = *std::next(last, -1);

        T ui = *it;
        T ui1 = *std::next(it, 1);
        if (p == 0)
        {
            return ((ui <= u) && (u < ui1)) || (std::abs(ui1 - u_last) < knot_eps<T> && std::abs(u - u_last) < knot_eps<T>)
                       ? T(1.)
                       : T(0.);
        }
        else
        {
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps<T>)
            {
                C1 = (u - ui) / C1;
                C1 *= basis_function(u, it, p - 1, last);
            }
            if (C2 > knot_eps<T>)
            {
                C2 = (ui1p - u) / C2;
                C2 *= basis_function(u, std::next(it, 1), p - 1, last);
            }
            return C1 + C2;
        }
    }

    /**
     * @brief Basis function used to compute BSpline's derivatives
     *
     * @tparam InIt
     * @tparam T
     * @param u   Parameter on BSpline object
     * @param it  Flat knots iterator
     * @param p   BSpline object degree
     * @param d   Derivative order
     * @param _Last
     * @return T
     */
    template <std::input_iterator InIt, std::floating_point T>
    T basis_function(T u, const InIt &it, size_t p, size_t d, const InIt &last)
    {
        if (d == 0)
        {
            return basis_function(u, it, p, last);
        }
        else if (d > p)
        {
            return T(0.);
        }
        else
        {
            T ui = *it;
            T ui1 = *std::next(it, 1);
            T uip = *std::next(it, p);
            T ui1p = *std::next(it, p + 1);
            T C1 = (uip - ui);
            T C2 = (ui1p - ui1);
            if (C1 > knot_eps<T>)
            {
                C1 = basis_function(u, it, p - 1, d - 1, last) / C1;
            }
            if (C2 > knot_eps<T>)
            {
                C2 = basis_function(u, std::next(it, 1), p - 1, d - 1, last) / C2;
            }
            return p * (C1 - C2);
        }
    }

    template <typename T>
    auto basis_funcs(size_t i, size_t p, T u, const std::vector<T> &k, std::vector<T>& N)
    {

        std::vector<T> left(p+1), right(p+1);
        T saved{}, temp{};
        N[0] = static_cast<T>(1.);
        for(size_t j{1}; j <= p; j++)
        {
            left[j]  = u - k[i+1-j];
            right[j] = k[i+j] - u;
            saved = 0.0;
            for( size_t r{0}; r <j; r++)
            {
                temp = N[r] / (right[r+1]+left[j-r]);
                N[r] = saved + right[r+1]*temp;
                saved = left[j-r]*temp;
            }
            N[j] = saved;
        }
    }

    // Largest B-spline degree the allocation-free evaluators handle on the
    // stack. Real CAD/CAM degrees are well under this; beyond it the evaluators
    // fall back to the (correct but slow) recursive basis_function so behaviour
    // never silently breaks. Stack cost is O((max+1)^2) scalars.
    inline constexpr size_t bspline_stack_max_degree = 24;

    /**
     * @brief Allocation-free B-spline basis functions (Piegl & Tiller A2.2).
     *
     * Writes the p+1 non-zero basis functions of knot span @p i into @p N, i.e.
     * N[r] = N_{i-p+r,p}(u) for r in [0, p]. Working storage lives on the stack;
     * requires p <= bspline_stack_max_degree.
     */
    template <typename T>
    void basis_funcs(size_t i, size_t p, T u, const std::vector<T> &k, T *N)
    {
        T left[bspline_stack_max_degree + 1];
        T right[bspline_stack_max_degree + 1];
        N[0] = T(1);
        for (size_t j = 1; j <= p; ++j)
        {
            left[j]  = u - k[i + 1 - j];
            right[j] = k[i + j] - u;
            T saved = T(0);
            for (size_t r = 0; r < j; ++r)
            {
                T temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
    }

    /**
     * @brief Allocation-free basis functions and derivatives (Piegl & Tiller
     *        A2.3, "The NURBS Book").
     *
     * Computes every derivative order 0..n of the p+1 non-zero basis functions
     * of knot span @p i in a single O(p*(p+n)) pass. Output @p ders is row-major
     * of shape (n+1) x (p+1): ders[r*(p+1)+j] = d^r N_{i-p+j,p}/du^r. All working
     * storage is on the stack; requires p <= bspline_stack_max_degree.
     */
    template <typename T>
    void ders_basis_funs(size_t i, size_t p, size_t n, T u, const std::vector<T> &k, T *ders)
    {
        const size_t s = p + 1; // row stride of the triangular work arrays
        T ndu[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
        T a[2 * (bspline_stack_max_degree + 1)];
        T left[bspline_stack_max_degree + 1];
        T right[bspline_stack_max_degree + 1];
        auto NDU = [&](size_t r, size_t c) -> T & { return ndu[r * s + c]; };
        auto A   = [&](size_t r, size_t c) -> T & { return a[r * s + c]; };

        NDU(0, 0) = T(1);
        for (size_t j = 1; j <= p; ++j)
        {
            left[j]  = u - k[i + 1 - j];
            right[j] = k[i + j] - u;
            T saved = T(0);
            for (size_t r = 0; r < j; ++r)
            {
                NDU(j, r) = right[r + 1] + left[j - r];   // lower triangle
                T temp = NDU(r, j - 1) / NDU(j, r);
                NDU(r, j) = saved + right[r + 1] * temp;  // upper triangle
                saved = left[j - r] * temp;
            }
            NDU(j, j) = saved;
        }

        // Order-0 derivatives are the basis functions themselves.
        for (size_t j = 0; j <= p; ++j)
            ders[j] = NDU(j, p);

        // Higher orders via Eq. (2.9): for each basis function index r, build the
        // coefficient table a[] up to order n.
        for (size_t r = 0; r <= p; ++r)
        {
            size_t s1 = 0, s2 = 1; // alternating rows of a[]
            A(0, 0) = T(1);
            for (size_t kk = 1; kk <= n; ++kk)
            {
                T d = T(0);
                const int rk = int(r) - int(kk);
                const int pk = int(p) - int(kk);
                if (r >= kk)
                {
                    A(s2, 0) = A(s1, 0) / NDU(size_t(pk + 1), size_t(rk));
                    d = A(s2, 0) * NDU(size_t(rk), size_t(pk));
                }
                const int j1 = (rk >= -1) ? 1 : -rk;
                const int j2 = (int(r) - 1 <= pk) ? int(kk) - 1 : int(p) - int(r);
                for (int j = j1; j <= j2; ++j)
                {
                    A(s2, size_t(j)) = (A(s1, size_t(j)) - A(s1, size_t(j - 1))) / NDU(size_t(pk + 1), size_t(rk + j));
                    d += A(s2, size_t(j)) * NDU(size_t(rk + j), size_t(pk));
                }
                if (int(r) <= pk)
                {
                    A(s2, kk) = -A(s1, kk - 1) / NDU(size_t(pk + 1), r);
                    d += A(s2, kk) * NDU(r, size_t(pk));
                }
                ders[kk * (p + 1) + r] = d;
                std::swap(s1, s2);
            }
        }

        // Multiply through by the falling-factorial factors p!/(p-k)!. Kept in T
        // (exact for these integer magnitudes) to avoid integer overflow for
        // high derivative orders at high degree.
        T fac = T(p);
        for (size_t kk = 1; kk <= n; ++kk)
        {
            for (size_t j = 0; j <= p; ++j)
                ders[kk * (p + 1) + j] *= fac;
            fac *= T(p - kk);
        }
    }

    /**
     * @brief Allocation-free d-th derivative of the non-zero basis functions.
     *
     * Fills Nd[0..p] with d^d N_{i-p+j,p}/du^d for j in [0, p]. Uses the cheap
     * A2.2 pass for d == 0 and the one-pass A2.3 otherwise. Requires
     * p <= bspline_stack_max_degree.
     */
    template <typename T>
    void basis_ders(size_t i, size_t p, size_t d, T u, const std::vector<T> &k, T *Nd)
    {
        if (d == 0)
        {
            basis_funcs(i, p, u, k, Nd);
            return;
        }
        T ders[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
        ders_basis_funs(i, p, d, u, k, ders);
        for (size_t j = 0; j <= p; ++j)
            Nd[j] = ders[d * (p + 1) + j];
    }

    template <typename T>
    auto find_span(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        // Piegl & Tiller, "The NURBS Book", Algorithm A2.1. Returns the knot
        // span s (as an iterator into k) such that k[s] <= u < k[s+1], with
        // u == k[n+1] mapped to the last non-empty span n. This matches the
        // half-open support convention (u < k[i+1]) used by basis_function,
        // which is what makes span reduction correct for DERIVATIVES at knots
        // of multiplicity > 1 (a plain lower_bound returns the left span there,
        // which is consistent for values but not for one-sided derivatives).
        if (u >= k[n + 1])
            return std::next(k.begin(), n);
        auto low = p, high = n + 1, mid = (low + high) / 2;
        while (u < k[mid] || u >= k[mid + 1])
        {
            if (u < k[mid]) high = mid;
            else            low  = mid;
            mid = (low + high) / 2;
        }
        return std::next(k.begin(), mid);
    }
    
    template <typename T>
    auto find_span2(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        if (u >= k[n + 1])
            return std::next(k.begin(), n);
        auto low = p, high = n + 1, mid = (low + high) / 2;
        while (u < k[mid] || u >= k[mid + 1])
        {
            if (u < k[mid])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }
        return std::next(k.begin(), mid);
    }

    template <typename T>
    auto find_span_range(size_t n, size_t p, T u, const std::vector<T> &k)
    {
        size_t i2 = find_span(n, p, u, k) - k.begin();
        const size_t i2_upper = (k.size() > p + 2) ? k.size() - p - 2 : 0;
        i2 = std::min(i2, i2_upper);
        size_t i1 = (i2 > p) ? i2 - p : 0;

        return std::make_pair(i1,i2);
    }

    // template <typename T>
    // auto eval_index_span(size_t n, size_t p, T u, const std::vector<T> &k, size_t d, bool use_span_reduction, size_t &i_min, size_t &i_max)
    // {
    //     if (use_span_reduction)//Reducing span for few pole makes things worst
    //     {
    //         i_max = find_span(n, p, u, k) - k.begin();
    //         i_max = std::min(i_max, k.size() - p - 2);
    //         i_min = std::max(int(0),int(i_max-p));
    //     }
    // }

    template <typename T>
    T basis_function(T u, size_t i, size_t p, size_t d, const std::vector<T> &k)
    {
        return basis_function(u, std::next(k.begin(), i), p, d, k.end());
    }

    /**
     * @brief BSpline curve evaluation using simple recursive basis functions
     * 
     * @tparam T 
     * @tparam dim 
     * @param u parameter on curve
     * @param k flat knots
     * @param poles poles
     * @param p degree
     * @param d derivative order
     * @param use_span_reduction 
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    auto eval_value_decasteljau(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> std::array<T, dim>
    {
        // One-pass allocation-free evaluation. Only the p+1 basis functions of
        // the knot span containing u are non-zero, so the sum is always span-
        // reduced regardless of `use_span_reduction` (the parameter is kept for
        // API compatibility): functions outside the span contribute exact zero.
        std::array<T, dim> pt;
        pt.fill(T(0));
        if (d > p)
            return pt; // derivative order above degree => identically zero

        const size_t n_poles = poles.size();
        size_t i = find_span(n_poles, p, u, k) - k.begin();
        const size_t i_upper = (k.size() > p + 2) ? k.size() - p - 2 : 0;
        i = std::min(i, i_upper);
        const size_t i_min = (i >= p) ? i - p : 0;

        if (p <= bspline_stack_max_degree)
        {
            T Nd[bspline_stack_max_degree + 1];
            basis_ders(i, p, d, u, k, Nd);
            for (size_t r = 0; r <= p; ++r)
            {
                const T N = Nd[r];
                const auto &pole = poles[i_min + r];
                for (size_t c = 0; c < dim; ++c)
                    pt[c] += N * pole[c];
            }
        }
        else
        {
            // Pathological degree: fall back to the recursive basis (correct,
            // exponential) so very high degrees still evaluate.
            for (size_t ip = i_min; ip <= i; ++ip)
            {
                const T N = basis_function(u, ip, p, d, k);
                const auto &pole = poles[ip];
                for (size_t c = 0; c < dim; ++c)
                    pt[c] += N * pole[c];
            }
        }

        return pt;
    }

    /**
     * @brief BSpline curve derivatives of all orders 0..n in a single pass.
     *
     * Computes C^(k)(u) for k in [0, n] and writes them into CK[0..n] (caller
     * provides storage for n+1 points). A single allocation-free A2.3 pass yields
     * every derivative order at once, so this replaces n+1 separate
     * eval_value_decasteljau calls. Orders above the degree are exactly zero.
     *
     * @param CK output buffer of at least n+1 points
     */
    template <typename T, size_t dim>
    void eval_ders_decasteljau(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p, size_t n, std::array<T, dim> *CK)
    {
        for (size_t kk = 0; kk <= n; ++kk)
            CK[kk].fill(T(0));

        const size_t du = std::min(n, p); // orders above the degree stay zero
        const size_t n_poles = poles.size();
        size_t i = find_span(n_poles, p, u, k) - k.begin();
        const size_t i_upper = (k.size() > p + 2) ? k.size() - p - 2 : 0;
        i = std::min(i, i_upper);
        const size_t i_min = (i >= p) ? i - p : 0;

        if (p <= bspline_stack_max_degree)
        {
            T ders[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
            ders_basis_funs(i, p, du, u, k, ders);
            for (size_t kk = 0; kk <= du; ++kk)
                for (size_t r = 0; r <= p; ++r)
                {
                    const T N = ders[kk * (p + 1) + r];
                    const auto &pole = poles[i_min + r];
                    for (size_t c = 0; c < dim; ++c)
                        CK[kk][c] += N * pole[c];
                }
        }
        else
        {
            // Pathological degree: recursive fallback (correct, exponential).
            for (size_t kk = 0; kk <= du; ++kk)
                for (size_t ip = i_min; ip <= i; ++ip)
                {
                    const T N = basis_function(u, ip, p, kk, k);
                    const auto &pole = poles[ip];
                    for (size_t c = 0; c < dim; ++c)
                        CK[kk][c] += N * pole[c];
                }
        }
    }

    // Forward declaration: defined later in this file, used by the n > max
    // fallback of eval_rational_ders below. Default arguments live here (on the
    // first declaration) so the later definition must not repeat them.
    template <typename T, size_t dim>
    auto eval_rational_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim + 1>> &poles, size_t p, size_t d = 0, bool use_span_reduction = true) -> std::array<T, dim>;

    /**
     * @brief Rational BSpline curve derivatives of all orders 0..n in one pass.
     *
     * Piegl & Tiller A4.2 ("The NURBS Book"): evaluate the homogeneous curve
     * derivatives Aders[0..n] (numerator) and the weight derivatives wders[0..n]
     * once via eval_ders_decasteljau on the weighted poles, then apply the
     * rational quotient rule. Writes the projected derivatives into CK[0..n].
     *
     * @param poles weighted poles (last coordinate is the weight)
     * @param CK    output buffer of at least n+1 points
     */
    template <typename T, size_t dim>
    void eval_rational_ders(T u, const std::vector<T> &k, const points_vector<T, dim + 1> &poles, size_t p, size_t n, std::array<T, dim> *CK)
    {
        if (n > bspline_stack_max_degree)
        {
            // Extremely high derivative order: fall back to the per-order path.
            for (size_t kk = 0; kk <= n; ++kk)
                CK[kk] = eval_rational_value_simple<T, dim>(u, k, poles, p, kk, false);
            return;
        }

        std::array<std::array<T, dim + 1>, bspline_stack_max_degree + 1> Aders;
        eval_ders_decasteljau<T, dim + 1>(u, k, poles, p, n, Aders.data());

        const T w0 = Aders[0].back();
        for (size_t kk = 0; kk <= n; ++kk)
        {
            std::array<T, dim> v;
            for (size_t c = 0; c < dim; ++c)
                v[c] = Aders[kk][c];
            for (size_t i = 1; i <= kk; ++i)
            {
                const T b = binomial_law<T>(kk, i) * Aders[i].back();
                for (size_t c = 0; c < dim; ++c)
                    v[c] -= b * CK[kk - i][c];
            }
            for (size_t c = 0; c < dim; ++c)
                CK[kk][c] = v[c] / w0;
        }
    }

    /**
     * @brief BSpline surface mixed derivatives, all orders up to (nu, nv), one pass.
     *
     * Computes S^(ku,kv)(u, v) for ku in [0, nu], kv in [0, nv] and writes them
     * row-major into SKL: SKL[ku*(nv+1)+kv]. A single A2.3 pass per parametric
     * direction yields every order, replacing (nu+1)*(nv+1) separate evaluations.
     * Mixed orders above the respective degree are exactly zero.
     *
     * @param SKL output buffer of at least (nu+1)*(nv+1) points
     */
    template <typename T, size_t dim>
    void eval_ders_decasteljau(T u, T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t nu, size_t nv, std::array<T, dim> *SKL)
    {
        for (size_t a = 0; a <= nu; ++a)
            for (size_t b = 0; b <= nv; ++b)
                SKL[a * (nv + 1) + b].fill(T(0));

        const size_t du = std::min(nu, p);
        const size_t dv = std::min(nv, q);
        const size_t n_polesU = ku.size() - p - 1;
        const size_t n_polesV = kv.size() - q - 1;

        size_t i = find_span(n_polesU, p, u, ku) - ku.begin();
        size_t j = find_span(n_polesV, q, v, kv) - kv.begin();
        const size_t i_upper = (ku.size() > p + 2) ? ku.size() - p - 2 : 0;
        const size_t j_upper = (kv.size() > q + 2) ? kv.size() - q - 2 : 0;
        i = std::min(i, i_upper);
        j = std::min(j, j_upper);
        const size_t i_min = (i >= p) ? i - p : 0;
        const size_t j_min = (j >= q) ? j - q : 0;

        if (p <= bspline_stack_max_degree && q <= bspline_stack_max_degree)
        {
            T dersU[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
            T dersV[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
            ders_basis_funs(i, p, du, u, ku, dersU);
            ders_basis_funs(j, q, dv, v, kv, dersV);
            for (size_t ku_ = 0; ku_ <= du; ++ku_)
                for (size_t kv_ = 0; kv_ <= dv; ++kv_)
                {
                    auto &dst = SKL[ku_ * (nv + 1) + kv_];
                    for (size_t rj = 0; rj <= q; ++rj)
                    {
                        const T nvw = dersV[kv_ * (q + 1) + rj];
                        const size_t row = n_polesU * (j_min + rj);
                        for (size_t ri = 0; ri <= p; ++ri)
                        {
                            const T w = dersU[ku_ * (p + 1) + ri] * nvw;
                            const auto &pole = poles[i_min + ri + row];
                            for (size_t c = 0; c < dim; ++c)
                                dst[c] += w * pole[c];
                        }
                    }
                }
        }
        else
        {
            // Pathological degree: recursive fallback (correct, exponential).
            for (size_t ku_ = 0; ku_ <= du; ++ku_)
                for (size_t kv_ = 0; kv_ <= dv; ++kv_)
                {
                    auto &dst = SKL[ku_ * (nv + 1) + kv_];
                    for (size_t jj = j_min; jj <= j; ++jj)
                    {
                        const T Nvv = basis_function(v, jj, q, kv_, kv);
                        for (size_t ii = i_min; ii <= i; ++ii)
                        {
                            const T Nuu = basis_function(u, ii, p, ku_, ku);
                            const auto &pole = poles[ii + n_polesU * jj];
                            for (size_t c = 0; c < dim; ++c)
                                dst[c] += Nuu * Nvv * pole[c];
                        }
                    }
                }
        }
    }


    // TODO Seems dead code
    /**
     * @brief BSpline curve evaluation using simple recursive basis functions
     *
     * @tparam T 
     * @tparam dim 
     * @param u parameter on curve
     * @param k flat knots
     * @param poles poles
     * @param p degree
     * @param d derivative order
     * @param use_span_reduction 
     * @return T
     */
    template <typename T, size_t dim>
    auto eval_value_decasteljau(T u, const std::vector<T> &k, const std::vector<T> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> T
    {
        T pt{};
        size_t n = poles.size();
        size_t i_max, i_min;
        eval_index_span(n, p, u, k, d, use_span_reduction, i_min, i_max);

        for (auto i = i_min; i <= i_max; i++)
        {
            auto N = basis_function(u, i, p, d, k);
            pt += N * poles[i];
        }
        
        return pt;
    }


    template <typename T, size_t dim>
    auto eval_value_par(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p, size_t d = 0,bool use_span_reduction =true) -> std::array<T, dim>
    {
        std::array<T, dim> pt;
        pt.fill(0);
        size_t n_poles = poles.size();
        size_t i_max{n_poles-1}, i_min{0};
        if (use_span_reduction )// correct for d != 0 too: find_span is half-open-consistent and basis derivatives share the basis support
        {
            i_max = find_span(n_poles, p, u, k) - k.begin();
            const size_t i_max_upper = (k.size() > p + 2) ? k.size() - p - 2 : 0;
            i_max = std::min(i_max, i_max_upper);
            i_min = (i_max > p) ? i_max - p : 0;
        }
        // The benifit of paralelization is not obvious
        std::vector<T> Ni(i_max+1-i_min);
        auto indexes = make_range(i_min,i_max);
        std::transform(
            GBS_PAR_EXEC
            indexes.begin(),
            indexes.end(),
            Ni.begin(),
            [u,p,d,&k](size_t i) {
                return basis_function(u, i, p, d, k);
            });
        
        points_vector<T, dim> poles_in_span(poles.begin()+i_min,poles.begin()+i_max+1);
        std::transform(
            GBS_PAR_EXEC
            Ni.begin(),
            Ni.end(),
            poles_in_span.begin(),
            poles_in_span.begin(),
            [](auto N_, const auto &pole_)
            {
                return pole_ * N_;
            }
        );

        // // return point<T,dim>{};        
        
        return std::reduce(
            GBS_PAR_EXEC
            poles_in_span.cbegin(),
            poles_in_span.cend()
        );

        // // auto indexes = make_range(i_min,i_max);

        // for (auto i = i_min; i <= i_max; i++)
        // {
        //     auto N = basis_function(u, i, p, d, k);
        //     // auto pole = poles[i];
        //     // for( size_t d {} ; d < dim ; d++ )
        //     // {
        //     //     pt[d] += N * pole[d];
        //     // }
        //     std::transform(
        //         // std::execution::par,
        //         poles[i].begin(),
        //         poles[i].end(),
        //         pt.begin(),
        //         pt.begin(),
        //         [&] (const auto val_, const auto tot_)
        //         {
        //             return tot_ + N * val_;
        //         }
        //     );
        // }
        
        // return pt;
    }

    template <typename T, size_t dim>
    auto eval_value_deboor_cox(T u, const std::vector<T> &k, const points_vector<T, dim> &poles, size_t p)
    {
        auto n_poles = poles.size();

        auto i =  std::min<size_t>( find_span2(n_poles, p, u, k) - k.begin(), k.size() - p - 2);

        std::vector<T> N(p+1);
        basis_funcs( i, p, u, k, N );

        point<T,dim> pt{0., 0., 0.};
        
        std::for_each(N.begin(), N.end(), [&] (const auto &N_)
        {
            std::transform(
                poles[i-p].begin(),
                poles[i-p].end(),
                pt.begin(),
                pt.begin(),
                [N_] (const auto pole_coordinate, const auto pt_coordinate)
                {
                    return pt_coordinate + N_ * pole_coordinate;
                }
            );
            i++;
        });
        return pt;
    }
    /**
     * @brief BSpline surface evaluation using simple recursive basis functions
     * 
     * @tparam T 
     * @tparam dim 
     * @param u first parameter on curve
     * @param v second parameter on curve
     * @param ku first flats knots
     * @param kv second flats knots
     * @param poles poles
     * @param p u degree
     * @param q v degree
     * @param du derivative order on u
     * @param dv derivative order on v
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    std::array<T, dim> eval_value_decasteljau(T u, T v, const std::vector<T> &ku, const std::vector<T> &kv, const std::vector<std::array<T, dim>> &poles, size_t p, size_t q, size_t du = 0, size_t dv = 0)
    {
        // One-pass allocation-free tensor-product evaluation: compute the
        // du-th derivative of the p+1 non-zero u-basis functions and the dv-th
        // derivative of the q+1 non-zero v-basis functions once each, then form
        // the outer product over the (p+1)x(q+1) non-zero poles.
        std::array<T, dim> pt;
        pt.fill(T(0));
        const size_t n_polesU = ku.size() - p - 1;
        const size_t n_polesV = kv.size() - q - 1;
        if (du > p || dv > q)
            return pt; // derivative order above degree => identically zero

        size_t i = find_span(n_polesU, p, u, ku) - ku.begin();
        size_t j = find_span(n_polesV, q, v, kv) - kv.begin();
        const size_t i_upper = (ku.size() > p + 2) ? ku.size() - p - 2 : 0;
        const size_t j_upper = (kv.size() > q + 2) ? kv.size() - q - 2 : 0;
        i = std::min(i, i_upper);
        j = std::min(j, j_upper);
        const size_t i_min = (i >= p) ? i - p : 0;
        const size_t j_min = (j >= q) ? j - q : 0;

        if (p <= bspline_stack_max_degree && q <= bspline_stack_max_degree)
        {
            T Nu[bspline_stack_max_degree + 1];
            T Nv[bspline_stack_max_degree + 1];
            basis_ders(i, p, du, u, ku, Nu);
            basis_ders(j, q, dv, v, kv, Nv);
            for (size_t rj = 0; rj <= q; ++rj)
            {
                const T nv = Nv[rj];
                const size_t row = n_polesU * (j_min + rj);
                for (size_t ri = 0; ri <= p; ++ri)
                {
                    const T w = Nu[ri] * nv;
                    const auto &pole = poles[i_min + ri + row];
                    for (size_t c = 0; c < dim; ++c)
                        pt[c] += w * pole[c];
                }
            }
        }
        else
        {
            // Pathological degree: recursive fallback (correct, exponential).
            for (size_t jj = j_min; jj <= j; ++jj)
            {
                const T Nvv = basis_function(v, jj, q, dv, kv);
                for (size_t ii = i_min; ii <= i; ++ii)
                {
                    const T Nuu = basis_function(u, ii, p, du, ku);
                    const auto &pole = poles[ii + n_polesU * jj];
                    for (size_t c = 0; c < dim; ++c)
                        pt[c] += Nuu * Nvv * pole[c];
                }
            }
        }
        return pt;
    }

    /**
     * @brief Utility function, add weights to poles array
     * 
     * @tparam T      : precision of curves
     * @tparam dim    : space dimension of curve (aka 1D, 2D, 3D,...)
     * @param poles   : non rational poles
     * @param weights : weights
     * @return std::vector<std::array<T, dim + 1>> : the unified poles
     */
    template <typename T, size_t dim>
    auto merge_weights(const std::vector<std::array<T, dim>> &poles, const std::vector<T> &weights) -> std::vector<std::array<T, dim + 1>>
    {
        std::vector<std::array<T, dim + 1>> p(poles.size());
        std::transform(
            GBS_PAR_EXEC
            poles.begin(), poles.end(), weights.begin(), p.begin(),
            [](const auto &v, const auto &c) {
                std::array<T, dim + 1> n;
                // std::copy(v.begin(), v.end(), n.begin());
                std::transform(v.begin(), v.end(), n.begin(),[&](const auto &v_){return v_*c;});
                n.back() = c;
                return n;
            });
        return p;
    }
    /**
     * @brief Project point with rational definition if last coord is null -> return origin point
     * 
     * @tparam T : precision of curves
     * @tparam dim : space dimension of curve (aka 1D, 2D, 3D,...)
     * @param pt : point coordinates to project
     * @return std::array<T, dim - 1> 
     */
    template <typename T, size_t dim>
    auto weight_projection(const std::array<T, dim> &pt) -> std::array<T, dim - 1>
    {
        std::array<T, dim - 1> r;
        if(pt.back()!=T{0}) 
        {
            std::transform(pt.begin(), std::next(pt.end(), -1), r.begin(), [&pt](const auto &pt_)
                           { return pt_ / pt.back(); });
        }
        else
        {
            r.fill(T{0.});
        }
        
        return r;
    }
    template <typename T, size_t dim, bool rational>
    inline auto weight_projected_pole(const point<T,dim+rational> &p) -> point<T,dim>
    {
        std::array<T, dim> r;
        if( rational )
        {
            if(p.back()!=T{0}) 
            {
                std::transform(p.begin(), std::next(p.end(), -1), r.begin(), [&p](const auto &p_)
                            { return p_ / p.back(); });
            }
            else
            {
                r.fill(T{0.});
            }
            
        }
        else
        {
            std::copy(p.begin(), p.end(), r.begin());
        }
        return r;
    }

    /**
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles_and_weights compact poles and weights storage
     * @param poles vector receiving projected poles 
     * @param weights vector receiving weights
     */
    template <typename T, size_t dim> 
    auto separate_weights(const std::vector<std::array<T, dim+1>> &poles_and_weights,std::vector<std::array<T, dim>> &poles, std::vector<T> &weights) ->void
    {
        poles.resize(poles_and_weights.size());
        weights.resize(poles_and_weights.size());

        std::transform(
            GBS_PAR_EXEC
            poles_and_weights.begin(),
            poles_and_weights.end(),
            poles.begin(),
            [](const auto &pw_) { return weight_projection<T, dim+1>(pw_); });

        std::transform(
            GBS_PAR_EXEC
            poles_and_weights.begin(),
            poles_and_weights.end(),
            weights.begin(),
            [](const auto &po_) {return po_.back(); } );
    }
    /**
     * @brief Add weight coordinates to a pole, the other coordinates are scaled
     * 
     * @tparam T 
     * @tparam dim 
     * @param pole 
     * @param w 
     * @return std::array<T, dim+1> 
     */
    template <typename T, size_t dim>
    auto add_weight(std::array<T, dim> pole, T w) -> std::array<T, dim+1>
    {
        std::array<T, dim+1> new_pole;
        std::transform(
            pole.begin(),
            pole.end(),
            new_pole.begin(),
            [&w](const auto &x_){return x_*w;}
        );
        new_pole.back()=w;
        return new_pole;
    }
    /**
     * @brief Add weight corrdinate colum, the other coordinates are scaled
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @return std::vector<std::array<T, dim+1>> 
     */
    template <typename T, size_t dim> 
    auto add_weights_coord(const std::vector<std::array<T, dim>> &poles) -> std::vector<std::array<T, dim+1>>
    {
        std::vector<std::array<T, dim+1>> poles_with_weights(poles.size());
        std::transform(
            GBS_PAR_EXEC
            poles.begin(),
            poles.end(),
            poles_with_weights.begin(),
            [](const auto p_)
            {

                return add_weight(p_,T(1.));
            }
        );
        return poles_with_weights;
    }

    /**
     * @brief Builds a pole vector containing weights
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param weights 
     * @return std::vector<std::array<T, dim+1>> 
     */
    template <typename T, size_t dim> 
    auto add_weights_coord(const std::vector<std::array<T, dim>> &poles, const std::vector<T> &weights) -> std::vector<std::array<T, dim+1>>
    {
        if (poles.size() != weights.size())
        {
            throw std::length_error("BSCurveGeneral: wrong pole vector length.");
        }
        std::vector<std::array<T, dim+1>> poles_with_weights(poles.size());
        std::transform(
            GBS_PAR_EXEC
            poles.begin(),
            poles.end(),
            weights.begin(),
            poles_with_weights.begin(),
            [](const auto p_, const auto w_)
            {

                return add_weight(p_,w_);
            }
        );
        return poles_with_weights;
    }

    /**
     * @brief Scale poles with origin {0. ... 0.}
     * 
     * @tparam T 
     * @tparam dim 
     * @param poles 
     * @param scale 
     */
    template <typename T, size_t dim>
    auto scale_poles(std::vector<std::array<T, dim>> &poles, T scale) -> void
    {
        std::transform(
            GBS_PAR_EXEC
            poles.begin(),
            poles.end(),
            poles.begin(),
            [&scale](const auto &p_) { return p_ * scale; });
    }

    template <typename T, size_t dim>
    auto scale_weights(std::vector<std::array<T, dim>> &poles, T mean_value = T(1.0)) -> void
    {
        // auto mean_value_ = std::reduce(
        //     std::execution::par,
        //     poles.begin(),
        //     poles.end(),
        //     T(0.),
        //     [](const auto &s_, const auto &v_) { s_.back() + v_.back(); });
        T mean_value_ = 0.;
        std::for_each(
            // std::execution::par,
            poles.begin(),
            poles.end(),
            [&mean_value_] (const auto &v_) mutable { mean_value_ += v_.back(); });
        mean_value_ /= poles.size();
        scale_poles(poles, mean_value / mean_value_);
    }

    template <typename T, size_t dim>
    auto scaled_poles(const std::vector<std::array<T, dim>> &poles, T scale) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_poles(new_poles,scale);
        return new_poles;
    }

    template <typename T, size_t dim>
    auto scaled_weights(const std::vector<std::array<T, dim>> &poles, T mean_value = T(1.0)) -> std::vector<std::array<T, dim>>
    {
        std::vector<std::array<T, dim>> new_poles(poles);
        scale_weights(new_poles);
        return new_poles;
    }
    /**
     * @brief Rational evaluation os a pol set and knots the last coordinate stands for weight and the others ares scales according to this weight
     * 
     * @tparam T 
     * @tparam dim 
     * @param u 
     * @param k 
     * @param poles 
     * @param p 
     * @param d 
     * @return std::array<T, dim> 
     */
    template <typename T, size_t dim>
    auto eval_rational_value_simple(T u, const std::vector<T> &k, const std::vector<std::array<T, dim+1>> &poles, size_t p, size_t d,bool use_span_reduction) -> std::array<T, dim>
    {
        if (d == 0)
        {
            return weight_projection(eval_value_decasteljau(u, k, poles , p, d, use_span_reduction));
        }
        else
        {
            auto wu = eval_value_decasteljau<T,dim+1>(u, k, poles, p, 0, false).back();
            auto Ckw = eval_value_decasteljau<T,dim+1>(u, k, poles, p, d, false);
            Ckw.back() = 1.;                  //
            auto Ak = weight_projection(Ckw); // not real projection just drop last coord
            std::array<T, dim> sum{Ak};
            for (size_t i = 1; i <= d; i++)
            {
                auto wi = eval_value_decasteljau<T,dim+1>(u, k, poles, p, i, false).back();
                auto C = eval_rational_value_simple<T,dim>(u, k, poles, p, d - i);
                sum = sum - binomial_law<T>(d, i) * wi * C;
            }
            sum = sum / wu;
            return sum;
        }
    }

} // namespace gbs