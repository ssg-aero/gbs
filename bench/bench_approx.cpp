// Performance audit harness for B-spline least-squares APPROXIMATION.
//
// Curve : gbs::approx(pts, p, n_poles, u, fix_bound)
//   fix_bound=false -> gbs::approx(...,k_flat)  : build_poles_matrix
//                      (BANDED one-pass assembly, ders) + dense colPivHouseholderQr
//   fix_bound=true  -> gbs::approx_bound_fixed  : also BANDED assembly via
//                      fill_basis_row since PR2 (#65) + dense colPivHouseholderQr
//   => post-PR2 both paths are banded: the A/B collapses to ~1.0 and is flat in
//      degree (it used to read off the recursive ~2^p blow-up, now removed).
//
// Surface: gbs::approx(grid, ku, kv, u, v, p, q)  (bssapprox.h)
//   fill_poles_matrix : BANDED tensor-product fill since PR2 (only the (p+1)(q+1)
//   non-zero columns per row) + dense colPivHouseholderQr on a
//   (nu*nv) x (np_u*np_v) matrix. The dense SOLVE is the remaining cost (PR3).
//
// We read off the scaling exponent (curve vs #points and vs degree; surface vs
// grid side N) and the recursive-vs-banded assembly gap (now closed for curves).

#include <gbs/bscapprox.h>
#include <gbs/bssapprox.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
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

static gbs::points_vector<double, 3> sample_curve(size_t n)
{
    gbs::points_vector<double, 3> Q(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = double(i) / double(n - 1);
        Q[i] = {std::cos(6.0 * t), std::sin(6.0 * t), 2.0 * t};
    }
    return Q;
}

static std::vector<double> uniform_params(size_t n)
{
    std::vector<double> u(n);
    for (size_t i = 0; i < n; ++i)
        u[i] = double(i) / double(n - 1);
    return u;
}

static std::vector<std::array<double, 3>> sample_grid(size_t nu, size_t nv)
{
    std::vector<std::array<double, 3>> Q(nu * nv); // U-first, V-outer
    for (size_t j = 0; j < nv; ++j)
        for (size_t i = 0; i < nu; ++i)
        {
            double u = double(i) / double(nu - 1);
            double v = double(j) / double(nv - 1);
            Q[j * nu + i] = {u, v, std::sin(3.0 * u) * std::cos(3.0 * v)};
        }
    return Q;
}

static std::string xprev(double ms, double prev)
{
    return prev > 0 ? (std::to_string(ms / prev) + "x") : std::string("-");
}

int main()
{
    using namespace gbs;

    // ---------------------------------------------------------------------
    // CURVE: banded one-pass (fix_bound=false) vs recursive (fix_bound=true),
    // degree 3, n_poles = n_points/4, vs #points.
    // ---------------------------------------------------------------------
    std::printf("CURVE approx, degree=3, n_poles=n/4, vs #points (total approx() ms)\n");
    std::printf("%-10s %10s %14s %14s %10s\n", "n_points", "n_poles", "banded ms", "recursive ms", "rec/band");
    std::printf("--------------------------------------------------------------------\n");
    for (size_t n : {50u, 100u, 200u, 400u, 800u})
    {
        auto Q = sample_curve(n);
        auto u = uniform_params(n);
        size_t p = 3;
        size_t n_poles = std::max<size_t>(n / 4, p + 1);
        int reps = n <= 200 ? 7 : 3;
        double band = best_ms(reps, [&] {
            auto c = approx<double, 3>(Q, p, n_poles, u, false);
            sink += c.poles()[0][0];
        });
        double rec = best_ms(reps, [&] {
            auto c = approx<double, 3>(Q, p, n_poles, u, true);
            sink += c.poles()[0][0];
        });
        std::printf("%-10zu %10zu %14.3f %14.3f %10s\n", n, n_poles, band, rec,
                    xprev(rec, band).c_str());
    }

    // ---------------------------------------------------------------------
    // CURVE: vs degree (assembly-dominated: n fixed, n_poles fixed).
    // Exposes the 2^p blow-up of the recursive path vs the flat banded path.
    // ---------------------------------------------------------------------
    std::printf("\nCURVE approx, n_points=400, n_poles=50, vs degree (total approx() ms)\n");
    std::printf("%-8s %14s %14s %10s\n", "degree", "banded ms", "recursive ms", "rec/band");
    std::printf("------------------------------------------------------\n");
    {
        auto Q = sample_curve(400);
        auto u = uniform_params(400);
        size_t n_poles = 50;
        for (size_t p : {1u, 2u, 3u, 5u, 7u, 9u})
        {
            double band = best_ms(5, [&] {
                auto c = approx<double, 3>(Q, p, n_poles, u, false);
                sink += c.poles()[0][0];
            });
            double rec = best_ms(5, [&] {
                auto c = approx<double, 3>(Q, p, n_poles, u, true);
                sink += c.poles()[0][0];
            });
            std::printf("%-8zu %14.3f %14.3f %10s\n", p, band, rec, xprev(rec, band).c_str());
        }
    }

    // ---------------------------------------------------------------------
    // SURFACE grid approx, bidegree (3,3), n_poles = (N/2)^2, vs grid side N.
    // Recursive dense Kronecker assembly + dense QR.
    // ---------------------------------------------------------------------
    std::printf("\nSURFACE approx, bidegree=(3,3), n_poles=(N/2)^2, vs grid side N (NxN)\n");
    std::printf("%-6s %10s %10s %12s %10s\n", "N", "pts", "poles", "time ms", "x prev");
    std::printf("------------------------------------------------------\n");
    double prev = 0;
    for (size_t N : {8u, 12u, 16u, 24u, 32u})
    {
        size_t np = std::max<size_t>(N / 2, 4u); // poles per direction
        auto Q = sample_grid(N, N);
        auto u = uniform_params(N);
        auto v = uniform_params(N);
        auto ku = build_simple_mult_flat_knots<double>(0.0, 1.0, np, 3);
        auto kv = build_simple_mult_flat_knots<double>(0.0, 1.0, np, 3);
        int reps = N <= 16 ? 4 : 2;
        double ms = best_ms(reps, [&] {
            auto poles = approx<double, 3>(Q, ku, kv, u, v, size_t{3}, size_t{3});
            sink += poles[0][0];
        });
        std::printf("%-6zu %10zu %10zu %12.3f %10s\n", N, N * N, np * np, ms, xprev(ms, prev).c_str());
        prev = ms;
    }

    std::printf("[sink=%g]\n", sink);
    return 0;
}
