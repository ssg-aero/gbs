#pragma once

#include <gbs/execution.h>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <gbs/bscurve.h>
#include <limits>
#include <algorithm>

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

// interp_sparse_threshold is defined in gbs/gbsconstants.h.

/**
 * @brief Solve the collocation system N * X = B for all dim right-hand sides,
 *        choosing a dense or sparse LU by problem size.
 *
 * The collocation matrix is assembled banded (p+1 non-zeros per row). For large
 * systems an Eigen::SparseLU with a fill-reducing ordering exploits that sparsity
 * (~O(n*p^2)); for small systems a dense partialPivLu is faster and is kept to
 * avoid any regression on low point counts. The factorization is computed once
 * and reused across the dim columns of B.
 */
template <typename T>
auto solve_collocation(const MatrixX<T> &N, const MatrixX<T> &B) -> MatrixX<T>
{
    if (N.rows() < interp_sparse_threshold)
        return N.partialPivLu().solve(B);

    Eigen::SparseMatrix<T> Ns = N.sparseView();
    Ns.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
    solver.analyzePattern(Ns);
    solver.factorize(Ns);
    return solver.solve(B);
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

// The banded collocation matrix: each parameter contributes nc rows (derivative
// orders 0..nc-1), each with only deg+1 non-zeros, lying in a narrow band
// (bandwidth = deg for simple-multiplicity interpolation). The band math below
// mirrors fill_basis_band() in knotsfunctions.ixx — keep the two in sync. Pre-#96
// the matrix was built dense (n x n setZero) then rescanned via sparseView() — two
// O(n^2) passes that made the build scale super-linearly. We now assemble straight
// into the band (band_lu) and solve with a no-pivot band LU; the sparse path below
// is the fallback for the non-banded / small-pivot cases.

// Sparse assembly (fallback path). nc is a runtime argument so a single instance
// serves every constraint count.
template <typename T>
auto assemble_collocation_sparse(const std::vector<T> &k_flat, const std::vector<T> &u,
                                 size_t deg, size_t nc, Eigen::Index n_poles) -> Eigen::SparseMatrix<T>
{
    std::vector<Eigen::Triplet<T>> trip;
    trip.reserve(size_t(n_poles) * (deg + 1));
    const size_t n_basis = k_flat.size() - deg - 1;
    const size_t dd = std::min<size_t>(nc - 1, deg); // orders above degree are zero
    for (size_t i = 0; i < u.size(); ++i)
    {
        const T ui = u[i];
        const auto [span, i_min, r_max] = basis_band_extent(n_basis, deg, ui, k_flat);
        if (deg <= bspline_stack_max_degree)
        {
            T ders[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
            ders_basis_funs(span, deg, dd, ui, k_flat, ders);
            for (size_t d = 0; d <= dd; ++d)
                for (size_t r = 0; r <= r_max; ++r)
                    trip.emplace_back(Eigen::Index(nc * i + d), Eigen::Index(i_min + r), ders[d * (deg + 1) + r]);
        }
        else
        {
            for (size_t d = 0; d <= dd; ++d)
                for (size_t r = 0; r <= r_max; ++r)
                    trip.emplace_back(Eigen::Index(nc * i + d), Eigen::Index(i_min + r),
                                      basis_function(ui, i_min + r, deg, d, k_flat));
        }
    }
    Eigen::SparseMatrix<T> Ns(n_poles, n_poles);
    Ns.setFromTriplets(trip.begin(), trip.end());
    Ns.makeCompressed();
    return Ns;
}

// In-place banded LU, NO pivoting. The POINT-interpolation collocation matrix is
// totally positive (Schoenberg-Whitney), so Gaussian elimination without pivoting
// is stable (de Boor) and stays strictly inside the band — no fill, no symbolic
// analysis, no supernodes, unlike a general sparse LU (which spent ~60% of the
// interpolation build on its numeric factorization, issue #96). The matrix is
// assembled DIRECTLY into the band from the basis (no SparseMatrix at all).
// Derivative-constrained (Hermite, nc>=2) and pathological systems are NOT totally
// positive: no-pivot can hit a small/zero pivot or a large multiplier, so
// factorize_inplace returns false (small-pivot or growth guard) and the caller
// falls back to the PIVOTED sparse/dense solve — making the band path provably as
// accurate as partial pivoting, never silently wrong.
template <typename T>
struct band_lu
{
    Eigen::Index n = 0;
    int kl = 0, ku = 0, W = 0;
    std::vector<T> M; // band(i,j) stored at i*W + (j - i + kl)

    inline T &at(Eigen::Index i, Eigen::Index j) { return M[std::size_t(i) * std::size_t(W) + std::size_t(j - i + kl)]; }
    inline T at(Eigen::Index i, Eigen::Index j) const { return M[std::size_t(i) * std::size_t(W) + std::size_t(j - i + kl)]; }

    // Structured (point/Hermite) interpolation: row nc*i+d is derivative order d at
    // u[i]. Unlike assemble_and_factorize_constraints there is NO `W*2 > n` bail:
    // here the band width is structurally narrow (W = nc + 2*deg-ish, independent of
    // n) because the parameters are sorted and the knots match, so the band path is
    // always efficient. (The constraints overload takes arbitrary, possibly
    // non-banded, orderings, hence its guard.)
    auto assemble_and_factorize(const std::vector<T> &k_flat, const std::vector<T> &u,
                                size_t deg, size_t nc, Eigen::Index n_poles) -> bool
    {
        if (deg > bspline_stack_max_degree)
            return false; // high degree -> generic sparse path
        n = n_poles;
        const size_t n_basis = k_flat.size() - deg - 1;
        const size_t dd = std::min<size_t>(nc - 1, deg);
        const size_t np = u.size();

        // pass 1: spans + band extent
        std::vector<size_t> spans(np), imins(np), rmaxs(np);
        long kl_ = 0, ku_ = 0;
        for (size_t i = 0; i < np; ++i)
        {
            const auto [span, i_min, r_max] = basis_band_extent(n_basis, deg, u[i], k_flat);
            spans[i] = span; imins[i] = i_min; rmaxs[i] = r_max;
            kl_ = std::max(kl_, long(nc * i + dd) - long(i_min));
            ku_ = std::max(ku_, long(i_min + r_max) - long(nc * i));
        }
        kl = int(kl_ < 0 ? 0 : kl_);
        ku = int(ku_ < 0 ? 0 : ku_);
        W = kl + ku + 1;
        M.assign(std::size_t(n) * std::size_t(W), T(0));

        // pass 2: basis values straight into the band
        T scale = T(0);
        T ders[(bspline_stack_max_degree + 1) * (bspline_stack_max_degree + 1)];
        for (size_t i = 0; i < np; ++i)
        {
            ders_basis_funs(spans[i], deg, dd, u[i], k_flat, ders);
            const size_t i_min = imins[i], r_max = rmaxs[i];
            for (size_t d = 0; d <= dd; ++d)
                for (size_t r = 0; r <= r_max; ++r)
                {
                    const T v = ders[d * (deg + 1) + r];
                    at(Eigen::Index(nc * i + d), Eigen::Index(i_min + r)) = v;
                    scale = std::max(scale, std::abs(v));
                }
        }

        return factorize_inplace(scale);
    }

    // In-place no-pivot banded LU. Returns false (caller falls back to the pivoted
    // sparse/dense solve) on a near-zero pivot OR a large multiplier. The latter is
    // the growth guard for #1: no-pivot is only committed to when it is stable.
    // Point interpolation is totally positive so multipliers stay small (~1) and
    // the guard never trips; Hermite/ill-conditioned systems trip it (small pivot →
    // big multiplier, or a structural zero pivot) and defer to the pivoted solver,
    // making the band path provably as accurate as partial pivoting.
    auto factorize_inplace(T scale) -> bool
    {
        const T tol = scale * std::numeric_limits<T>::epsilon() * T(16);
        for (Eigen::Index k = 0; k < n; ++k)
        {
            const T akk = at(k, k);
            if (std::abs(akk) <= tol)
                return false;
            const Eigen::Index imax = std::min<Eigen::Index>(k + kl, n - 1);
            for (Eigen::Index i = k + 1; i <= imax; ++i)
            {
                const T lik = at(i, k) / akk;
                if (std::abs(lik) > interp_band_max_growth<T>)
                    return false; // no-pivot would be unstable here -> pivoted fallback
                at(i, k) = lik;
                const Eigen::Index jmax = std::min<Eigen::Index>(k + ku, n - 1);
                for (Eigen::Index j = k + 1; j <= jmax; ++j)
                    at(i, j) -= lik * at(k, j);
            }
        }
        return true;
    }

    // Assemble + factorize for GENERAL constraints (constrPoint): row i is the
    // ds[i]-th derivative of the basis at us[i]. The rows MUST be sorted by (u,d)
    // so the matrix is banded — sorting is a pure row permutation, so the poles
    // (columns) still come out in natural order. Returns false (caller falls back
    // to the dense/sparse solve) if a derivative order exceeds the degree, the
    // system is not actually banded (or is tiny), or a pivot is too small.
    auto assemble_and_factorize_constraints(const std::vector<T> &us, const std::vector<size_t> &ds,
                                            const std::vector<T> &k_flat, size_t deg) -> bool
    {
        if (deg > bspline_stack_max_degree)
            return false;
        n = Eigen::Index(us.size());
        const size_t n_basis = k_flat.size() - deg - 1;
        std::vector<size_t> spans(us.size()), imins(us.size()), rmaxs(us.size());
        long kl_ = 0, ku_ = 0;
        for (size_t i = 0; i < us.size(); ++i)
        {
            if (ds[i] > deg)
                return false; // derivative above degree -> degenerate (zero) row
            const auto [span, i_min, r_max] = basis_band_extent(n_basis, deg, us[i], k_flat);
            spans[i] = span; imins[i] = i_min; rmaxs[i] = r_max;
            kl_ = std::max(kl_, long(i) - long(i_min));
            ku_ = std::max(ku_, long(i_min + r_max) - long(i));
        }
        kl = int(kl_ < 0 ? 0 : kl_);
        ku = int(ku_ < 0 ? 0 : ku_);
        W = kl + ku + 1;
        if (std::size_t(W) * 2 > std::size_t(n)) // not really banded (or tiny) -> dense/sparse is better
            return false;
        M.assign(std::size_t(n) * std::size_t(W), T(0));
        T scale = T(0);
        T Nd[bspline_stack_max_degree + 1];
        for (size_t i = 0; i < us.size(); ++i)
        {
            basis_ders(spans[i], deg, ds[i], us[i], k_flat, Nd);
            for (size_t r = 0; r <= rmaxs[i]; ++r)
            {
                const T v = Nd[r];
                at(Eigen::Index(i), Eigen::Index(imins[i] + r)) = v;
                scale = std::max(scale, std::abs(v));
            }
        }
        return factorize_inplace(scale);
    }

    // Solve L U x = b (x and b are distinct length-n buffers).
    void solve(const T *b, T *x) const
    {
        for (Eigen::Index i = 0; i < n; ++i) // forward: unit-lower L
        {
            T s = b[i];
            const Eigen::Index kmin = std::max<Eigen::Index>(0, i - kl);
            for (Eigen::Index k = kmin; k < i; ++k)
                s -= at(i, k) * x[k];
            x[i] = s;
        }
        for (Eigen::Index i = n - 1; i >= 0; --i) // back: upper U
        {
            T s = x[i];
            const Eigen::Index jmax = std::min<Eigen::Index>(i + ku, n - 1);
            for (Eigen::Index j = i + 1; j <= jmax; ++j)
                s -= at(i, j) * x[j];
            x[i] = s / at(i, i);
        }
    }
};

// Collocation factorization with a band-LU fast path and a sparse-LU fallback,
// reused across all right-hand sides (issue #96).
template <typename T>
struct collocation_solver
{
    band_lu<T> band;
    bool banded = false;
    Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::NaturalOrdering<int>> sparse;

    void prepare(const std::vector<T> &k_flat, const std::vector<T> &u,
                 size_t deg, size_t nc, Eigen::Index n_poles)
    {
        banded = band.assemble_and_factorize(k_flat, u, deg, nc, n_poles);
        if (!banded)
        {
            Eigen::SparseMatrix<T> Ns = assemble_collocation_sparse<T>(k_flat, u, deg, nc, n_poles);
            sparse.analyzePattern(Ns);
            sparse.factorize(Ns);
        }
    }
    // Multi-RHS solve: every column of B is solved against the single
    // factorization. The band path writes each solution straight into the output
    // column (B is column-major, so .data() is contiguous; b and x are distinct
    // matrices) — no per-column temporary. The sparse path solves all RHS at once.
    auto solve(const MatrixX<T> &B) -> MatrixX<T>
    {
        if (!banded)
            return sparse.solve(B);
        MatrixX<T> X(B.rows(), B.cols());
        for (Eigen::Index c = 0; c < B.cols(); ++c)
            band.solve(B.col(c).data(), X.col(c).data());
        return X;
    }
};

template <typename T,size_t dim,size_t nc>
    auto build_poles(const std::vector<constrType<T,dim,nc> > &Q, const std::vector<T> &k_flat,const std::vector<T> &u, size_t deg) -> std::vector<std::array<T,dim> >
{
    static_assert(nc >= 1, "build_poles needs at least one constraint per point (nc >= 1)");
    const size_t n_pt = Q.size();
    const Eigen::Index n_poles = Eigen::Index(Q.size() * nc);

    // Banded collocation matrix assembled directly into band storage and solved by
    // a no-pivot band LU (issue #96), falling back to a sparse LU only on a small
    // pivot. Factorized once; the `dim` spatial components are solved as one
    // multi-RHS block.
    collocation_solver<T> solver;
    solver.prepare(k_flat, u, deg, nc, n_poles);

    MatrixX<T> B(n_poles, Eigen::Index(dim));
    for (size_t i = 0; i < n_pt; ++i)
        for (size_t deriv = 0; deriv < nc; ++deriv)
            for (size_t d = 0; d < dim; ++d)
                B(Eigen::Index(nc * i + deriv), Eigen::Index(d)) = Q[i][deriv][d];

    MatrixX<T> X = solver.solve(B);

    std::vector<std::array<T, dim>> poles(static_cast<std::size_t>(n_poles));
    for (Eigen::Index i = 0; i < n_poles; ++i)
        for (size_t d = 0; d < dim; ++d)
            poles[size_t(i)][d] = X(i, Eigen::Index(d));

    return poles;
}

// Batched, constrained pole solve. Every column in `Q_cols` is interpolated
// through the SAME (k_flat, u, deg) constraint system, so the banded collocation
// matrix is assembled and factorized ONCE and every column is solved together as
// extra right-hand sides. Equivalent to calling the single-column overload above
// once per column, but without re-building/re-factorizing the matrix each time
// (the lone real perf win of the spine loft, issue #40 PR2 / #44).
//
// The returned poles are the per-column results concatenated in order: column 0's
// `n_pt*nc` poles, then column 1's, ... — matching repeated single-column calls.
template <typename T, size_t dim, size_t nc>
auto build_poles(const std::vector<std::vector<constrType<T, dim, nc>>> &Q_cols,
                 const std::vector<T> &k_flat, const std::vector<T> &u, size_t deg)
    -> std::vector<std::array<T, dim>>
{
    static_assert(nc >= 1, "build_poles needs at least one constraint per point (nc >= 1)");
    if (Q_cols.empty())
        return {};

    const size_t n_cols = Q_cols.size();
    const size_t n_pt = Q_cols.front().size();
    const auto n_poles = Eigen::Index(n_pt * nc);

    // Banded collocation matrix assembled directly into band storage and solved by
    // a no-pivot band LU (issue #96), falling back to a sparse LU on a small pivot.
    // Factorized once; every (column, spatial component) is one RHS, solved as a
    // single multi-RHS block.
    collocation_solver<T> solver;
    solver.prepare(k_flat, u, deg, nc, n_poles);

    MatrixX<T> B(n_poles, Eigen::Index(n_cols * dim));
    for (size_t c = 0; c < n_cols; ++c)
        for (size_t i = 0; i < n_pt; ++i)
            for (size_t deriv = 0; deriv < nc; ++deriv)
                for (size_t d = 0; d < dim; ++d)
                    B(Eigen::Index(nc * i + deriv), Eigen::Index(c * dim + d)) = Q_cols[c][i][deriv][d];

    MatrixX<T> X = solver.solve(B);

    std::vector<std::array<T, dim>> poles(size_t(n_poles) * n_cols);
    for (size_t c = 0; c < n_cols; ++c)
        for (Eigen::Index i = 0; i < n_poles; ++i)
            for (size_t d = 0; d < dim; ++d)
                poles[c * size_t(n_poles) + size_t(i)][d] = X(i, Eigen::Index(c * dim + d));

    return poles;
}

template <typename T, size_t dim>
auto build_poles(const std::vector<constrPoint<T, dim>> &Q, const std::vector<T> &k_flat, size_t deg) -> std::vector<std::array<T, dim>>
{
    const size_t n_poles = Q.size();

    // General constraints can arrive in any order, so the collocation matrix is
    // not banded as given. Sorting the constraints by (u, d) is a pure ROW
    // permutation that makes it banded (the poles are the COLUMNS, so they stay in
    // natural order — no un-permutation needed). We then solve with the no-pivot
    // band LU (issue #96); the band assembler bails (-> dense/sparse fallback
    // below) for tiny, non-banded or ill-conditioned systems.
    // The common interpolation case already arrives sorted by (u, d), so only copy
    // + sort when it isn't (avoids an O(n) copy + O(n log n) sort that permutes
    // nothing).
    auto by_ud = [](const auto &a, const auto &b) { return a.u < b.u || (a.u == b.u && a.d < b.d); };
    const std::vector<constrPoint<T, dim>> *Qsrc = &Q;
    std::vector<constrPoint<T, dim>> Qsorted;
    if (!std::is_sorted(Q.begin(), Q.end(), by_ud))
    {
        Qsorted = Q;
        std::sort(Qsorted.begin(), Qsorted.end(), by_ud);
        Qsrc = &Qsorted;
    }
    const std::vector<constrPoint<T, dim>> &Qs = *Qsrc;
    std::vector<T> us(n_poles);
    std::vector<size_t> ds(n_poles);
    for (size_t i = 0; i < n_poles; ++i) { us[i] = Qs[i].u; ds[i] = Qs[i].d; }

    band_lu<T> band;
    if (band.assemble_and_factorize_constraints(us, ds, k_flat, deg))
    {
        std::vector<std::array<T, dim>> poles(n_poles);
        std::vector<T> b(n_poles), x(n_poles);
        for (size_t d = 0; d < dim; ++d)
        {
            for (size_t i = 0; i < n_poles; ++i) b[i] = Qs[i].v[d];
            band.solve(b.data(), x.data());
            for (size_t i = 0; i < n_poles; ++i) poles[i][d] = x[i];
        }
        return poles;
    }

    // Fallback: the original dense/sparse collocation solve (handles any order).
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_poles, n_poles);
    N.setZero();
    for (size_t i = 0; i < n_poles; ++i)
        fill_basis_row(N, Eigen::Index(i), Q[i].u, k_flat, deg, Q[i].d);

    MatrixX<T> B(n_poles, dim);
    for (size_t i = 0; i < n_poles; ++i)
        for (size_t d = 0; d < dim; ++d)
            B(Eigen::Index(i), Eigen::Index(d)) = Q[i].v[d];

    MatrixX<T> X = solve_collocation(N, B);

    std::vector<std::array<T, dim>> poles(n_poles);
    for (size_t i = 0; i < n_poles; ++i)
        for (size_t d = 0; d < dim; ++d)
            poles[i][d] = X(Eigen::Index(i), Eigen::Index(d));

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

// Batch interpolation (#91): build N INDEPENDENT curves from N point sets, each
// through the single-entity interpolate() above. The outer loop over the sets is
// threshold-guarded parallel via build_batch (gate gbs::parallel_batch_min_size):
// the builds are independent (each owns its solver/poles), so the result is
// bit-identical to a serial loop regardless of thread count. The single per-curve
// solve stays sequential — there is no useful intra-build parallel axis (#84).
template <typename T, size_t dim>
auto interpolate(const std::vector<points_vector<T, dim>> &Q_lst, size_t p, KnotsCalcMode mode)
    -> std::vector<BSCurve<T, dim>>
{
    return build_batch(Q_lst, [p, mode](const points_vector<T, dim> &Q)
                       { return interpolate<T, dim>(Q, p, mode); });
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