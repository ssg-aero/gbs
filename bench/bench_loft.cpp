// Performance audit harness for the spine-loft pole solve (issue #40 PR2 / #44).
//
// The spine loft interpolates, for every u-pole column, a v-direction constraint
// system (position + spine tangent). The v-collocation matrix is IDENTICAL for
// every column, yet the old code called the single-column build_poles() inside
// the loop over u-poles, so that matrix was assembled and factorized n_poles_u
// times. PR2 adds a batched build_poles() that assembles + factorizes ONCE and
// solves every column as multiple right-hand sides.
//
// This bench isolates exactly that kernel: it builds the same constraint columns
// the spine loft feeds to build_poles, then times
//   OLD: a loop of single-column build_poles()  (re-factorize per column)
//   NEW: one batched build_poles(Q_cols, ...)    (factorize once)
// and prints the speedup. Both overloads are exercised against the same data, so
// the comparison is self-contained (no baseline checkout needed) and the poles
// are verified identical to machine precision.
//
// We sweep:
//   - n_poles_u  (columns sharing one factorization) -> redundant factorizations
//   - n_poles_v  (sections) -> size of the v-collocation matrix

#include <gbs/bscinterp.h>
#include <gbs/bsctools.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using clk = std::chrono::steady_clock;
volatile double sink = 0.0;

template <typename F>
double best_ms(int reps, F &&f)
{
    double best = 1e300;
    for (int r = 0; r < reps; ++r)
    {
        auto t0 = clk::now();
        f();
        auto t1 = clk::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

// Build the (position + tangent) constraint columns and the v-system exactly as
// the spine loft does: n_poles_u columns, each with n_poles_v position+tangent
// constraints, sharing one degree-q v-collocation matrix.
struct loft_solve_input
{
    std::vector<std::vector<gbs::constrType<double, 3, 2>>> Q_cols;
    std::vector<double> kv_flat;
    std::vector<double> v;
    size_t q;
};

static loft_solve_input make_input(size_t n_poles_u, size_t n_poles_v)
{
    loft_solve_input in;
    in.q = std::max<size_t>(std::min<size_t>(3, n_poles_v - 1), 1);
    in.v.resize(n_poles_v);
    for (size_t j = 0; j < n_poles_v; ++j)
        in.v[j] = double(j) / double(n_poles_v - 1);
    in.kv_flat = gbs::build_simple_mult_flat_knots<double>(in.v.front(), in.v.back(), n_poles_v * 2, in.q);

    in.Q_cols.resize(n_poles_u, std::vector<gbs::constrType<double, 3, 2>>(n_poles_v));
    for (size_t i = 0; i < n_poles_u; ++i)
        for (size_t j = 0; j < n_poles_v; ++j)
        {
            double u = double(i) / double(n_poles_u - 1);
            double t = double(j) / double(n_poles_v - 1);
            std::array<double, 3> pos = {std::cos(3.0 * u) + t, std::sin(3.0 * u), 2.0 * t};
            std::array<double, 3> tg = {0.1, 0.2, 2.0};
            in.Q_cols[i][j] = gbs::constrType<double, 3, 2>{pos, tg};
        }
    return in;
}

// OLD path: re-factorize the v-system per column.
static std::vector<std::array<double, 3>> solve_per_column(const loft_solve_input &in)
{
    std::vector<std::array<double, 3>> poles;
    for (const auto &Q : in.Q_cols)
    {
        auto poles_v = gbs::build_poles(Q, in.kv_flat, in.v, in.q);
        poles.insert(poles.end(), poles_v.begin(), poles_v.end());
    }
    return poles;
}

// NEW path: factorize once, solve all columns as multiple RHS.
static std::vector<std::array<double, 3>> solve_batched(const loft_solve_input &in)
{
    return gbs::build_poles(in.Q_cols, in.kv_flat, in.v, in.q);
}

static double max_diff(const std::vector<std::array<double, 3>> &a,
                       const std::vector<std::array<double, 3>> &b)
{
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t d = 0; d < 3; ++d)
            m = std::max(m, std::abs(a[i][d] - b[i][d]));
    return m;
}

int main()
{
    std::printf("SPINE-loft pole solve: per-column re-factorize vs batched factorize-once\n\n");

    std::printf("3 sections (6x6 v-system), vs n_poles_u (columns)\n");
    std::printf("%-12s %12s %12s %10s %12s\n", "n_poles_u", "old ms", "new ms", "speedup", "max diff");
    std::printf("-----------------------------------------------------------------\n");
    for (size_t nu : {32u, 64u, 128u, 256u, 512u})
    {
        auto in = make_input(nu, 3);
        int reps = nu <= 128 ? 9 : 5;
        double old_ms = best_ms(reps, [&] { sink += solve_per_column(in)[0][0]; });
        double new_ms = best_ms(reps, [&] { sink += solve_batched(in)[0][0]; });
        double diff = max_diff(solve_per_column(in), solve_batched(in));
        std::printf("%-12zu %12.3f %12.3f %9.2fx %12.2e\n", nu, old_ms, new_ms, old_ms / new_ms, diff);
    }

    std::printf("\n256 columns, vs n_poles_v (v-system size = 2*n_poles_v)\n");
    std::printf("%-12s %12s %12s %10s %12s\n", "n_poles_v", "old ms", "new ms", "speedup", "max diff");
    std::printf("-----------------------------------------------------------------\n");
    for (size_t nv : {3u, 5u, 9u, 17u, 33u})
    {
        auto in = make_input(256, nv);
        double old_ms = best_ms(5, [&] { sink += solve_per_column(in)[0][0]; });
        double new_ms = best_ms(5, [&] { sink += solve_batched(in)[0][0]; });
        double diff = max_diff(solve_per_column(in), solve_batched(in));
        std::printf("%-12zu %12.3f %12.3f %9.2fx %12.2e\n", nv, old_ms, new_ms, old_ms / new_ms, diff);
    }

    std::printf("[sink=%g]\n", sink);
    return 0;
}
