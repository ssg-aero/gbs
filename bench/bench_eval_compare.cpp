// Curve value evaluator: before/after for issue #30.
// Column "recursive (old)" is a local copy of the pre-#30 hot path that
// BSCurve::value used: span-reduced, but per-pole RECURSIVE basis_function
// (O(p . 2^p)). Column "A2.3 (new)" is the current eval_value_decasteljau
// (allocation-free Piegl & Tiller A2.3, O(p^2)). For reference we also show
// the library's eval_value_deboor_cox (O(p^2) but allocates a vector per call).

#include <gbs/bscurve.h>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

// Pre-#30 hot path, preserved here for an honest before/after measurement:
// span-reduced sum over the p+1 poles, each weight from the recursive Cox-de Boor
// basis_function (no memoization => ~2^p work per evaluation).
template <typename T, std::size_t dim>
static std::array<T, dim> eval_recursive_old(T u, const std::vector<T> &k,
                                             const gbs::points_vector<T, dim> &poles,
                                             std::size_t p, std::size_t d = 0)
{
    std::array<T, dim> pt; pt.fill(T(0));
    const std::size_t n_poles = poles.size();
    std::size_t i_max = gbs::find_span(n_poles, p, u, k) - k.begin();
    const std::size_t i_max_upper = (k.size() > p + 2) ? k.size() - p - 2 : 0;
    i_max = std::min(i_max, i_max_upper);
    const std::size_t i_min = (i_max > p) ? i_max - p : 0;
    for (std::size_t i = i_min; i <= i_max; ++i)
    {
        const T N = gbs::basis_function(u, i, p, d, k);
        for (std::size_t c = 0; c < dim; ++c) pt[c] += N * poles[i][c];
    }
    return pt;
}

using clk = std::chrono::steady_clock;
volatile double sink = 0.0;

static std::vector<double> clamped_uniform_knots(size_t n_poles, size_t p)
{
    size_t n_int = n_poles - p - 1;
    std::vector<double> k;
    for (size_t i = 0; i <= p; ++i) k.push_back(0.0);
    for (size_t i = 1; i <= n_int; ++i) k.push_back(double(i) / double(n_int + 1));
    for (size_t i = 0; i <= p; ++i) k.push_back(1.0);
    return k;
}

template <typename F>
double best_ns(size_t n, int reps, F &&f)
{
    double best = 1e300;
    for (int r = 0; r < reps; ++r)
    {
        auto t0 = clk::now();
        double a = f();
        auto t1 = clk::now();
        sink += a;
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count());
    }
    return best / double(n);
}

int main()
{
    using namespace gbs;
    const size_t N = 200000, n_poles = 60;
    const int reps = 7;
    std::vector<double> u(N);
    for (size_t i = 0; i < N; ++i) u[i] = (double(i) + 0.5) / double(N);

    std::printf("Curve VALUE (d=0), %zu poles, N=%zu, clang-21 -O3\n", n_poles, N);
    std::printf("%-8s %18s %16s %12s %16s\n", "degree", "recursive(old) ns", "A2.3(new) ns", "speedup", "deboor-cox ns");
    std::printf("----------------------------------------------------------------------------\n");
    for (size_t p : {2u, 3u, 5u, 7u})
    {
        auto k = clamped_uniform_knots(n_poles, p);
        points_vector<double, 3> poles(n_poles);
        for (size_t i = 0; i < n_poles; ++i)
            poles[i] = {std::sin(0.7 * i), std::cos(0.5 * i), 0.1 * i};

        double a = best_ns(N, reps, [&]{ double s = 0;
            for (double x : u) s += eval_recursive_old(x, k, poles, p, size_t{0})[0]; return s; });
        double b = best_ns(N, reps, [&]{ double s = 0;
            for (double x : u) s += eval_value_decasteljau(x, k, poles, p, 0)[0]; return s; });
        double c = best_ns(N, reps, [&]{ double s = 0;
            for (double x : u) s += eval_value_deboor_cox(x, k, poles, p)[0];     return s; });
        std::printf("%-8zu %18.2f %16.2f %11.1fx %16.2f\n", p, a, b, a / b, c);
    }
    std::printf("[sink=%g]\n", sink);
    return 0;
}
