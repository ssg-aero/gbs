// Performance study of the arc-length (curvilinear abscissa) derivative
// evaluators added on feat/vectorized-curvilinear-derivatives:
//
//   Curves   : d_dm, d_dm2          + vectorized d_dms,  d_dm2s
//   Surfaces : d_dmu, d_dmv,
//              d_dmu2, d_dmv2        + vectorized d_dmus, d_dmvs, d_dmu2s, d_dmv2s
//
// The study answers three questions:
//   1. Throughput (ns / evaluation, M evals / s) for each kernel.
//   2. Overhead of the vectorized wrapper vs a hand-written scalar loop.
//   3. Cost of the 2nd-order kernels relative to the 1st-order ones
//      (d_dm2 issues TWO De Boor passes: value(.,1) AND value(.,2)).
//
// Built without C++20 modules, so knot vectors / sample lists are built by
// hand (no make_range / build_simple_mult_flat_knots, which live in .ixx).

#include <gbs/bscurve.h>
#include <gbs/bssurf.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using clk = std::chrono::steady_clock;

namespace
{
    // Clamped, uniformly-spaced flat knot vector for `n_poles` control points
    // of degree `p` on [0, 1]. Size = n_poles + p + 1.
    std::vector<double> clamped_uniform_knots(std::size_t n_poles, std::size_t p)
    {
        const std::size_t n_int = n_poles - p - 1;     // interior knots
        std::vector<double> k;
        k.reserve(n_poles + p + 1);
        for (std::size_t i = 0; i <= p; ++i) k.push_back(0.0);
        for (std::size_t i = 1; i <= n_int; ++i)
            k.push_back(static_cast<double>(i) / static_cast<double>(n_int + 1));
        for (std::size_t i = 0; i <= p; ++i) k.push_back(1.0);
        return k;
    }

    // Deterministic, mildly wiggly pole coordinate so the curve/surface is
    // non-degenerate (non-zero, varying tangents) everywhere.
    double wiggle(std::size_t i, std::size_t axis)
    {
        const double t = 0.37 * static_cast<double>(i) + 1.3 * static_cast<double>(axis);
        // cheap pseudo-trig without <cmath> in the hot path of construction
        return std::sin(t) + 0.5 * std::cos(2.0 * t) + 0.1 * static_cast<double>(i);
    }

    gbs::BSCurve<double, 3> make_curve(std::size_t n_poles, std::size_t p)
    {
        gbs::points_vector<double, 3> poles(n_poles);
        for (std::size_t i = 0; i < n_poles; ++i)
            poles[i] = {wiggle(i, 0), wiggle(i, 1), wiggle(i, 2)};
        return gbs::BSCurve<double, 3>(poles, clamped_uniform_knots(n_poles, p), p);
    }

    gbs::BSSurface<double, 3> make_surface(std::size_t nu, std::size_t nv,
                                           std::size_t p, std::size_t q)
    {
        gbs::points_vector<double, 3> poles(nu * nv);
        // Library storage: U-first (inner), V-outer.
        for (std::size_t j = 0; j < nv; ++j)
            for (std::size_t i = 0; i < nu; ++i)
            {
                const std::size_t idx = j * nu + i;
                poles[idx] = {static_cast<double>(i) + 0.3 * wiggle(j, 0),
                              static_cast<double>(j) + 0.3 * wiggle(i, 1),
                              0.5 * wiggle(i + j, 2)};
            }
        return gbs::BSSurface<double, 3>(poles,
                                         clamped_uniform_knots(nu, p),
                                         clamped_uniform_knots(nv, q),
                                         p, q);
    }

    // Sample parameters strictly inside (0, 1) to stay clear of the clamped ends.
    std::vector<double> sample_u(std::size_t n)
    {
        std::vector<double> u(n);
        for (std::size_t i = 0; i < n; ++i)
            u[i] = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        return u;
    }

    std::vector<std::array<double, 2>> sample_uv(std::size_t n)
    {
        std::vector<std::array<double, 2>> uv(n);
        // Two coprime-ish strides so (u, v) sweeps the patch, not the diagonal.
        for (std::size_t i = 0; i < n; ++i)
        {
            const double u = (static_cast<double>((i * 7 + 3) % n) + 0.5) / static_cast<double>(n);
            const double v = (static_cast<double>((i * 13 + 5) % n) + 0.5) / static_cast<double>(n);
            uv[i] = {u, v};
        }
        return uv;
    }

    // sink to defeat dead-code elimination
    volatile double g_sink = 0.0;

    struct Result { double ns_per_eval; double meval_per_s; };

    template <typename F>
    Result bench(std::size_t n_evals, int repeats, F &&f)
    {
        double best_ns = 1e300;
        for (int r = 0; r < repeats; ++r)
        {
            const auto t0 = clk::now();
            const double acc = f();
            const auto t1 = clk::now();
            g_sink += acc;
            const double ns =
                std::chrono::duration<double, std::nano>(t1 - t0).count();
            best_ns = std::min(best_ns, ns);
        }
        const double per = best_ns / static_cast<double>(n_evals);
        return {per, 1000.0 / per};   // 1000/ns = evals per microsecond = M/s
    }

    void row(const std::string &name, std::size_t n, Result scalar, Result vec)
    {
        const double overhead = (vec.ns_per_eval / scalar.ns_per_eval - 1.0) * 100.0;
        std::cout << std::left << std::setw(10) << name
                  << std::right << std::setw(10) << n
                  << std::setw(12) << std::fixed << std::setprecision(2) << scalar.ns_per_eval
                  << std::setw(12) << scalar.meval_per_s
                  << std::setw(12) << vec.ns_per_eval
                  << std::setw(12) << vec.meval_per_s
                  << std::setw(11) << std::showpos << overhead << std::noshowpos << "%\n";
    }

    void header(const std::string &title)
    {
        std::cout << "\n" << title << "\n";
        std::cout << std::left << std::setw(10) << "kernel"
                  << std::right << std::setw(10) << "N"
                  << std::setw(12) << "scal ns/ev"
                  << std::setw(12) << "scal M/s"
                  << std::setw(12) << "vec ns/ev"
                  << std::setw(12) << "vec M/s"
                  << std::setw(12) << "vec ovhd" << "\n";
        std::cout << std::string(78, '-') << "\n";
    }
}

int main()
{
    using namespace gbs;
    const std::vector<std::size_t> sizes = {1000, 10000, 100000, 1000000};
    const int repeats = 7;

    std::cout << "Arc-length derivative performance study\n";
    std::cout << "clang-21 / -O3 / double / dim=3\n";

    // ---- Curves: degree-3, 60 poles --------------------------------------
    {
        const auto crv = make_curve(/*n_poles=*/60, /*p=*/3);
        header("CURVE  BSCurve<double,3>  degree=3, 60 poles");
        for (auto n : sizes)
        {
            const auto u = sample_u(n);

            auto s1 = bench(n, repeats, [&]{
                double a = 0; for (double x : u) a += crv.d_dm(x)[0]; return a; });
            auto v1 = bench(n, repeats, [&]{
                auto r = crv.d_dms(u); return r.empty() ? 0.0 : r.front()[0]; });
            row("d_dm", n, s1, v1);

            auto s2 = bench(n, repeats, [&]{
                double a = 0; for (double x : u) a += crv.d_dm2(x)[0]; return a; });
            auto v2 = bench(n, repeats, [&]{
                auto r = crv.d_dm2s(u); return r.empty() ? 0.0 : r.front()[0]; });
            row("d_dm2", n, s2, v2);
        }
    }

    // ---- Curve degree sweep (N fixed) ------------------------------------
    {
        const std::size_t n = 200000;
        const auto u = sample_u(n);
        std::cout << "\nCURVE degree sweep (N=" << n << ", 60 poles)\n";
        std::cout << std::left << std::setw(10) << "degree"
                  << std::right << std::setw(14) << "d_dm ns/ev"
                  << std::setw(14) << "d_dm2 ns/ev"
                  << std::setw(14) << "ratio 2/1" << "\n";
        std::cout << std::string(52, '-') << "\n";
        for (std::size_t p : {2u, 3u, 5u, 7u})
        {
            const auto crv = make_curve(60, p);
            auto a = bench(n, repeats, [&]{
                double s = 0; for (double x : u) s += crv.d_dm(x)[0]; return s; });
            auto b = bench(n, repeats, [&]{
                double s = 0; for (double x : u) s += crv.d_dm2(x)[0]; return s; });
            std::cout << std::left << std::setw(10) << p
                      << std::right << std::setw(14) << std::fixed << std::setprecision(2) << a.ns_per_eval
                      << std::setw(14) << b.ns_per_eval
                      << std::setw(14) << (b.ns_per_eval / a.ns_per_eval) << "\n";
        }
    }

    // ---- Surfaces: bidegree (3,3), 30x30 poles ---------------------------
    {
        const auto srf = make_surface(/*nu=*/30, /*nv=*/30, /*p=*/3, /*q=*/3);
        header("SURFACE  BSSurface<double,3>  bidegree=(3,3), 30x30 poles");
        for (auto n : sizes)
        {
            const auto uv = sample_uv(n);

            auto s1 = bench(n, repeats, [&]{
                double a = 0; for (auto &p : uv) a += srf.d_dmu(p[0], p[1])[0]; return a; });
            auto v1 = bench(n, repeats, [&]{
                auto r = srf.d_dmus(uv); return r.empty() ? 0.0 : r.front()[0]; });
            row("d_dmu", n, s1, v1);

            auto s2 = bench(n, repeats, [&]{
                double a = 0; for (auto &p : uv) a += srf.d_dmu2(p[0], p[1])[0]; return a; });
            auto v2 = bench(n, repeats, [&]{
                auto r = srf.d_dmu2s(uv); return r.empty() ? 0.0 : r.front()[0]; });
            row("d_dmu2", n, s2, v2);
        }
    }

    std::cout << "\n[sink=" << g_sink << "]\n";
    return 0;
}
