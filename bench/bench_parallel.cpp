// Parallelization audit harness (#84).
//
// Measures serial vs parallel (std::execution::seq vs ::par) wall time for the
// embarrassingly-parallel multi-point evaluation kernels of the library, to:
//   1. test the standing claim at gbs/bscurve.h:67
//      ("Parallelization doesn't seems to be worth") empirically;
//   2. locate the crossover N where par overtakes seq;
//   3. report speedup at scale AND scaling vs core count;
//   4. confirm the parallel result is BIT-IDENTICAL to the serial one
//      (independent writes to distinct output slots -> deterministic).
//
// The kernels mirror the real ones (file:line in comments). This bench does NOT
// change any algorithm; it reproduces the kernel shape on a std::transform.
//
// Build/run:  scripts/build_bench.sh bench_parallel

#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>

#include <oneapi/tbb/global_control.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <execution>
#include <functional>
#include <vector>

using namespace gbs; // bring the gbs vecop operators (^, +, *, /) into scope
using clk = std::chrono::steady_clock;
volatile double sink = 0.0;

template <typename F>
static double best_ms(int reps, F &&f)
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

// A representative degree-3, 3D curve with a non-trivial pole count (so find_span
// + basis evaluation do realistic work).
static auto make_curve(size_t n_poles)
{
    gbs::points_vector<double, 3> Q(n_poles);
    for (size_t i = 0; i < n_poles; ++i)
    {
        double t = double(i) / double(n_poles - 1);
        Q[i] = {std::cos(6.0 * t), std::sin(6.0 * t), 2.0 * t};
    }
    return gbs::interpolate<double, 3>(Q, size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH);
}

static auto make_surface(size_t n_side)
{
    std::vector<std::array<double, 3>> Q(n_side * n_side);
    for (size_t j = 0; j < n_side; ++j)
        for (size_t i = 0; i < n_side; ++i)
        {
            double u = double(i) / double(n_side - 1);
            double v = double(j) / double(n_side - 1);
            Q[j * n_side + i] = {u, v, std::sin(3.0 * u) * std::cos(3.0 * v)};
        }
    return gbs::interpolate<double, 3>(Q, n_side, size_t{3}, size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH);
}

static std::vector<double> param_list(std::array<double, 2> b, size_t n)
{
    std::vector<double> u(n);
    for (size_t i = 0; i < n; ++i)
        u[i] = b[0] + (b[1] - b[0]) * (double(i) + 0.5) / double(n);
    return u;
}

template <typename Pts>
static bool bit_identical(const Pts &a, const Pts &b)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t d = 0; d < a[i].size(); ++d)
            if (a[i][d] != b[i][d])
                return false;
    return true;
}

// ---- kernels (shape-identical to the library functions cited) --------------

// gbs/bscurve.h:64 values() / bscanalysis.h:505 make_points(): value(u)
template <typename Pol, typename Crv>
static auto eval_curve(Pol pol, const Crv &crv, const std::vector<double> &u)
{
    gbs::points_vector<double, 3> out(u.size());
    std::transform(pol, u.begin(), u.end(), out.begin(),
                   [&crv](double u_) { return crv(u_); });
    return out;
}

// gbs/bscurve.h:159 d_dms(): d_dm(u) = C'(u)/|C'(u)|  (heavier per element)
template <typename Pol, typename Crv>
static auto eval_curve_ddm(Pol pol, const Crv &crv, const std::vector<double> &u)
{
    gbs::points_vector<double, 3> out(u.size());
    std::transform(pol, u.begin(), u.end(), out.begin(),
                   [&crv](double u_) { return crv.d_dm(u_); });
    return out;
}

// gbs/bssurf.h surface value(u,v)
template <typename Pol, typename Srf>
static auto eval_surface(Pol pol, const Srf &srf, const std::vector<std::array<double, 2>> &uv)
{
    gbs::points_vector<double, 3> out(uv.size());
    std::transform(pol, uv.begin(), uv.end(), out.begin(),
                   [&srf](const auto &p) { return srf(p[0], p[1]); });
    return out;
}

// gbs/offsets.h:7 offset_points kernel: value + du + dv + normal + offset
template <typename Pol, typename Srf>
static auto eval_offset(Pol pol, const Srf &srf, const std::vector<std::array<double, 2>> &uv)
{
    gbs::points_vector<double, 3> out(uv.size());
    std::transform(pol, uv.begin(), uv.end(), out.begin(),
                   [&srf](const auto &p) {
                       auto pt = srf(p[0], p[1]);
                       auto tu = srf(p[0], p[1], 1, 0);
                       auto tv = srf(p[0], p[1], 0, 1);
                       auto n = tu ^ tv;
                       n = n / gbs::norm(n);
                       return pt + n * 0.1;
                   });
    return out;
}

static std::vector<std::array<double, 2>> uv_list(std::array<double, 4> b, size_t n)
{
    std::vector<std::array<double, 2>> uv(n);
    // a pseudo-scattered but deterministic fill across the domain
    for (size_t i = 0; i < n; ++i)
    {
        double a = (double(i) + 0.5) / double(n);
        double c = std::fmod(a * 7.0, 1.0);
        uv[i] = {b[0] + (b[1] - b[0]) * a, b[2] + (b[3] - b[2]) * c};
    }
    return uv;
}

template <typename KernelSeq, typename KernelPar>
static void sweep(const char *title, KernelSeq kseq, KernelPar kpar,
                  const std::vector<size_t> &sizes,
                  std::function<size_t(size_t)> reps_of)
{
    std::printf("\n%s\n", title);
    std::printf("%-10s %12s %12s %10s %12s\n", "N", "seq ms", "par ms", "speedup", "determinism");
    std::printf("--------------------------------------------------------------\n");
    for (size_t N : sizes)
    {
        int reps = int(reps_of(N));
        auto rs = kseq(N);
        auto rp = kpar(N);
        bool same = bit_identical(rs.first, rp.first);
        sink += rs.first[0][0] + rp.first[0][0];
        std::printf("%-10zu %12.3f %12.3f %9.2fx %12s\n", N, rs.second, rp.second,
                    rs.second / rp.second, same ? "bit-identical" : "DIFFERS");
    }
}

int main()
{
    using namespace gbs;
    unsigned hw = tbb::this_task_arena::max_concurrency();
    std::printf("Parallel eval audit (#84)  hw concurrency = %u\n", hw);

    auto crv = make_curve(60);
    auto bc = crv.bounds();

    auto srf = make_surface(40);
    auto bs = srf.bounds();

    // ---- 1. curve value(u) -------------------------------------------------
    {
        std::vector<size_t> sizes{1000, 10000, 100000, 1000000, 4000000};
        auto kseq = [&](size_t N) {
            auto u = param_list(bc, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_curve(std::execution::seq, crv, u); });
            return std::make_pair(r, ms);
        };
        auto kpar = [&](size_t N) {
            auto u = param_list(bc, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_curve(std::execution::par, crv, u); });
            return std::make_pair(r, ms);
        };
        sweep("[1] CURVE value(u)  (make_points / values kernel, deg 3 / 60 poles)",
              kseq, kpar, sizes, [](size_t) { return 1; });
    }

    // ---- 2. curve d_dm(u) --------------------------------------------------
    {
        std::vector<size_t> sizes{1000, 10000, 100000, 1000000};
        auto kseq = [&](size_t N) {
            auto u = param_list(bc, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_curve_ddm(std::execution::seq, crv, u); });
            return std::make_pair(r, ms);
        };
        auto kpar = [&](size_t N) {
            auto u = param_list(bc, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_curve_ddm(std::execution::par, crv, u); });
            return std::make_pair(r, ms);
        };
        sweep("[2] CURVE d_dm(u)  (d_dms kernel: C'(u)/|C'(u)|)",
              kseq, kpar, sizes, [](size_t) { return 1; });
    }

    // ---- 3. surface value(u,v) --------------------------------------------
    {
        std::vector<size_t> sizes{1000, 10000, 100000, 1000000};
        auto kseq = [&](size_t N) {
            auto uv = uv_list(bs, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_surface(std::execution::seq, srf, uv); });
            return std::make_pair(r, ms);
        };
        auto kpar = [&](size_t N) {
            auto uv = uv_list(bs, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_surface(std::execution::par, srf, uv); });
            return std::make_pair(r, ms);
        };
        sweep("[3] SURFACE value(u,v)  (deg (3,3) / 40x40 poles)",
              kseq, kpar, sizes, [](size_t) { return 1; });
    }

    // ---- 4. surface offset point (heavy kernel) ---------------------------
    {
        std::vector<size_t> sizes{1000, 10000, 100000, 1000000};
        auto kseq = [&](size_t N) {
            auto uv = uv_list(bs, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_offset(std::execution::seq, srf, uv); });
            return std::make_pair(r, ms);
        };
        auto kpar = [&](size_t N) {
            auto uv = uv_list(bs, N);
            points_vector<double, 3> r;
            double ms = best_ms(N <= 100000 ? 7 : 3, [&] { r = eval_offset(std::execution::par, srf, uv); });
            return std::make_pair(r, ms);
        };
        sweep("[4] SURFACE offset point  (offsets.h kernel: value + du + dv + normal)",
              kseq, kpar, sizes, [](size_t) { return 1; });
    }

    // ---- 5. thread scaling on a fixed, large workload ----------------------
    {
        std::printf("\n[5] THREAD SCALING  (surface offset point, N = 1,000,000)\n");
        std::printf("%-10s %12s %10s %12s\n", "threads", "ms", "speedup", "efficiency");
        std::printf("--------------------------------------------------\n");
        const size_t N = 1000000;
        auto uv = uv_list(bs, N);
        double t1 = 0;
        for (unsigned nt : {1u, 2u, 4u, 8u, 16u, 32u, std::min(64u, hw)})
        {
            tbb::global_control gc(tbb::global_control::max_allowed_parallelism, nt);
            points_vector<double, 3> r;
            double ms = best_ms(3, [&] { r = eval_offset(std::execution::par, srf, uv); });
            sink += r[0][0];
            if (nt == 1u)
                t1 = ms;
            std::printf("%-10u %12.3f %9.2fx %11.0f%%\n", nt, ms, t1 / ms, 100.0 * (t1 / ms) / nt);
        }
    }

    std::printf("\n[sink=%g]\n", sink);
    return 0;
}
