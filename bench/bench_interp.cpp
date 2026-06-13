// Performance audit harness for B-spline interpolation.
//
// Curve : gbs::interpolate(points, p, CHORD_LENGTH)  -> build_poles
//           -> build_poles_matrix (dense, recursive basis_function, all n^2
//              entries) + dense partialPivLu (O(n^3)).
// Surface: gbs::interpolate(grid, n_v, p, q, CHORD_LENGTH) -> build_poles
//           -> fill_poles_matrix (dense Kronecker, (nu*nv)^2 entries) +
//              dense partialPivLu (O((nu*nv)^3)).
//
// We measure total interpolate() wall time and read off the scaling exponent
// (curve vs #points and vs degree; surface vs grid side N).

#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>

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

int main()
{
    using namespace gbs;
    const auto mode = KnotsCalcMode::CHORD_LENGTH;

    std::printf("CURVE interpolation, degree=3, vs #points (total interpolate() ms)\n");
    std::printf("%-10s %12s %12s\n", "n_points", "time ms", "x prev");
    std::printf("-----------------------------------------\n");
    double prev = 0;
    for (size_t n : {50u, 100u, 200u, 400u, 800u})
    {
        auto Q = sample_curve(n);
        int reps = n <= 200 ? 7 : 3;
        double ms = best_ms(reps, [&] {
            auto c = interpolate<double, 3>(Q, size_t{3}, mode);
            sink += c.poles()[0][0];
        });
        std::printf("%-10zu %12.3f %12s\n", n, ms,
                    prev > 0 ? (std::to_string(ms / prev) + "x").c_str() : "-");
        prev = ms;
    }

    std::printf("\nCURVE interpolation, n_points=200, vs degree (total ms)\n");
    std::printf("%-10s %12s\n", "degree", "time ms");
    std::printf("--------------------------\n");
    for (size_t p : {1u, 2u, 3u, 5u, 7u, 9u})
    {
        auto Q = sample_curve(200);
        double ms = best_ms(5, [&] {
            auto c = interpolate<double, 3>(Q, p, mode);
            sink += c.poles()[0][0];
        });
        std::printf("%-10zu %12.3f\n", p, ms);
    }

    std::printf("\nSURFACE interpolation, bidegree=(3,3), vs grid side N (NxN) (total ms)\n");
    std::printf("%-10s %12s %12s %12s\n", "N", "poles", "time ms", "x prev");
    std::printf("-------------------------------------------------------\n");
    prev = 0;
    for (size_t N : {8u, 12u, 16u, 24u, 32u})
    {
        auto Q = sample_grid(N, N);
        int reps = N <= 16 ? 5 : 2;
        double ms = best_ms(reps, [&] {
            auto s = interpolate<double, 3>(Q, size_t{N}, size_t{3}, size_t{3}, mode);
            sink += s.poles()[0][0];
        });
        std::printf("%-10zu %12zu %12.3f %12s\n", N, N * N, ms,
                    prev > 0 ? (std::to_string(ms / prev) + "x").c_str() : "-");
        prev = ms;
    }

    std::printf("[sink=%g]\n", sink);
    return 0;
}
