// Performance and accuracy of the extrema (gbs/extrema.h), to compare NLopt (COBYLA, master
// before the removal of NLopt) and Newton's method (gbs/intersect/projection.h).
//
// The same source is compiled against both versions of the headers. For each case it prints:
//   - the time per call (best of a few repetitions),
//   - the evaluations of the geometry per call (value() calls, counted by wrappers) for the
//     local cases; the global cases run on the geometries themselves (a wrapper would hide a
//     B-spline from the subdivision, which uses its poles) and show "-",
//   - the true distance reached, recomputed from the returned parameters (the NLopt version
//     returned the square of the distance), against the exact distance when it is known.
// With an output file, it also writes the distance of every query, to compare the versions
// query by query (bench/python/compare_extrema.py).
//
// Usage: bench_extrema [out.csv]

#include <gbs/extrema.h>
#include <gbs/bscbuild.h>
#include <gbs/bscinterp.h>
#include <gbs/bselementary.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <vector>

using namespace gbs;
using T = double;
using P = point<T, 3>;
using clk = std::chrono::steady_clock;

static std::atomic<long> n_eval{0};

/// Counts the evaluations of a curve.
struct CountingCurve : Curve<T, 3>
{
    const Curve<T, 3> &c;
    explicit CountingCurve(const Curve<T, 3> &c_) : c{c_} {}
    auto value(T u, size_t d = 0) const -> P override
    {
        n_eval.fetch_add(1, std::memory_order_relaxed);
        return c.value(u, d);
    }
    auto bounds() const -> std::array<T, 2> override { return c.bounds(); }
};

/// Counts the evaluations of a surface.
struct CountingSurface : Surface<T, 3>
{
    const Surface<T, 3> &s;
    explicit CountingSurface(const Surface<T, 3> &s_) : s{s_} {}
    auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> P override
    {
        n_eval.fetch_add(1, std::memory_order_relaxed);
        return s.value(u, v, du, dv);
    }
    auto bounds() const -> std::array<T, 4> override { return s.bounds(); }
};

struct Case
{
    std::string name;
    int n;                                ///< queries
    std::function<T(int)> distance;       ///< runs query i, returns the true distance reached
    std::function<T(int)> exact{nullptr}; ///< exact distance of query i, if known
};

static FILE *out{};
static int n_throw{};

/// The true distance of query i, NaN when the extrema threw.
static T guarded(const std::function<T(int)> &f, int i)
{
    try
    {
        return f(i);
    }
    catch (const std::exception &)
    {
        ++n_throw;
        return std::numeric_limits<T>::quiet_NaN();
    }
}

static void run(const Case &c)
{
    // accuracy and evaluations (one pass), then the time (best of 3 passes)
    std::vector<T> d(c.n);
    n_eval = 0;
    n_throw = 0;
    for (int i = 0; i < c.n; ++i)
        d[i] = guarded(c.distance, i);
    const int thrown = n_throw;
    const double evals = double(n_eval.load()) / c.n;
    double best = 1e300;
    for (int r = 0; r < 3; ++r)
    {
        const auto t0 = clk::now();
        for (int i = 0; i < c.n; ++i)
            d[i] = guarded(c.distance, i);
        best = std::min(best, std::chrono::duration<double, std::micro>(clk::now() - t0).count() / c.n);
    }
    T err_max{}, err_sum{};
    int bad{};
    if (c.exact)
        for (int i = 0; i < c.n; ++i)
        {
            if (std::isnan(d[i]))
                continue;
            const T e = d[i] - c.exact(i);
            err_max = std::max(err_max, std::abs(e));
            err_sum += std::abs(e);
            bad += std::abs(e) > T(1e-6);
        }
    char ev[32] = "         -";
    if (evals > 0)
        std::snprintf(ev, sizeof ev, "%10.0f", evals);
    if (c.exact)
        std::printf("%-26s %10.1f us %s evals   err max %.2e  mean %.2e  > 1e-6: %d/%d\n", c.name.c_str(), best, ev, err_max, err_sum / c.n,
                    bad, c.n);
    else
        std::printf("%-26s %10.1f us %s evals   (no exact distance)\n", c.name.c_str(), best, ev);
    if (thrown)
        std::printf("%-26s %d/%d queries threw\n", "", thrown, c.n);
    if (out)
        for (int i = 0; i < c.n; ++i)
            std::fprintf(out, "%s,%d,%.17g,%.17g\n", c.name.c_str(), i, d[i], c.exact ? c.exact(i) : T(-1));
}

int main(int argc, char **argv)
{
    if (argc > 1)
        out = std::fopen(argv[1], "w");
    std::mt19937 rng{42};
    std::uniform_real_distribution<T> U{0, 1};
    const T pi = std::numbers::pi_v<T>;
    const T tol = T(1e-6);
    const ax2<T, 3> ax{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}};

    // geometries
    const auto circle_ = build_circle<T, 3>(T(1));
    const CountingCurve circle{circle_};
    const std::vector<constrType<T, 3, 1>> Q{{{0., 0., 0.}}, {{1., 0., 0.}}, {{1., 1., 0.}}, {{1., 1., 2.}}, {{0., 1., 1.}}, {{0., -1., 4.}}};
    const auto bspline_ = interpolate(Q, 3, KnotsCalcMode::CHORD_LENGTH);
    const CountingCurve bspline{bspline_};
    const std::vector<constrType<T, 3, 1>> Q2{{{-1., 0.5, 1.}}, {{0., 2., 0.5}}, {{1., 0.5, 2.}}, {{2., -1., 1.5}}, {{1.5, 0.5, 3.}}};
    const auto bspline2_ = interpolate(Q2, 3, KnotsCalcMode::CHORD_LENGTH);
    const CountingCurve bspline2{bspline2_};
    const auto sphere_ = build_sphere<T>(T(1), ax);
    const CountingSurface sphere{sphere_};
    const auto torus_ = build_torus<T>(T(3), T(1), ax);
    const CountingSurface torus{torus_};

    // queries
    const int N = 200;
    std::vector<P> pc(N), ps(N), pt(N);
    std::vector<T> pc_d(N), ps_d(N), pt_d(N);
    for (int i = 0; i < N; ++i)
    {
        const T a = 2 * pi * U(rng), r = T(0.5) + U(rng), z = T(0.5) * (U(rng) - T(0.5));
        pc[i] = {r * std::cos(a), r * std::sin(a), z};
        pc_d[i] = std::hypot(r - 1, z);
        const T b = pi * (U(rng) - T(0.5)) * T(0.98), rs = T(0.3) + T(1.4) * U(rng);
        ps[i] = {rs * std::cos(b) * std::cos(a), rs * std::cos(b) * std::sin(a), rs * std::sin(b)};
        ps_d[i] = std::abs(rs - 1);
        const T rt = T(0.2) + T(1.2) * U(rng); // distance to the core circle of the torus
        pt[i] = P{3 * std::cos(a), 3 * std::sin(a), 0} + rt * P{std::cos(b) * std::cos(a), std::cos(b) * std::sin(a), std::sin(b)};
        pt_d[i] = std::abs(rt - 1);
    }
    // nearest of 100 samples, moved by 0.02: start of the local overloads
    auto start_curve = [](const Curve<T, 3> &c, const P &p) {
        const auto [a, b] = c.bounds();
        T best = a, dmin = 1e300;
        for (int k = 0; k <= 100; ++k)
            if (const T u = a + (b - a) * k / 100, d = norm(c.value(u) - p); d < dmin)
                dmin = d, best = u;
        return std::clamp(best + T(0.02), a, b);
    };
    auto start_surface = [](const Surface<T, 3> &s, const P &p) {
        const auto [u1, u2, v1, v2] = s.bounds();
        std::array<T, 2> best{u1, v1};
        T dmin = 1e300;
        for (int i = 0; i <= 40; ++i)
            for (int j = 0; j <= 40; ++j)
                if (const T u = u1 + (u2 - u1) * i / 40, v = v1 + (v2 - v1) * j / 40, d = norm(s.value(u, v) - p); d < dmin)
                    dmin = d, best = {u, v};
        return std::array<T, 2>{std::clamp(best[0] + T(0.02), u1, u2), std::clamp(best[1] + T(0.02), v1, v2)};
    };
    std::vector<T> u0_circle(N);
    std::vector<std::array<T, 2>> uv0_sphere(N);
    for (int i = 0; i < N; ++i)
    {
        u0_circle[i] = start_curve(circle_, pc[i]);
        uv0_sphere[i] = start_surface(sphere_, ps[i]);
    }

    // segments outside the unit circle / sphere: exact distance to the sphere
    const int NS = 20;
    std::vector<std::array<P, 2>> seg(N);
    std::vector<T> seg_d(N);
    for (int i = 0; i < N; ++i)
    {
        const T a = 2 * pi * U(rng), b = pi * (U(rng) - T(0.5)) * T(0.8), h = T(1.1) + U(rng);
        const P n{std::cos(b) * std::cos(a), std::cos(b) * std::sin(a), std::sin(b)};
        P t = cross(n, P{0, 0, 1});
        t = t / norm(t);
        const T s0 = -T(1.5) * U(rng), s1 = T(1.5) * U(rng); // the foot of the centre lies inside
        seg[i] = {h * n + s0 * t, h * n + s1 * t};
        seg_d[i] = h - 1;
    }
    std::vector<BSCurve<T, 3>> seg_c;
    for (int i = 0; i < N; ++i)
        seg_c.push_back(build_segment<T, 3>(seg[i][0], seg[i][1]));
    std::vector<CountingCurve> segs;
    std::vector<std::array<T, 3>> tuv0_seg(N);
    for (int i = 0; i < N; ++i)
    {
        segs.emplace_back(seg_c[i]);
        const auto [a, b] = seg_c[i].bounds();
        const auto uv = start_surface(sphere_, seg_c[i]((a + b) / 2));
        tuv0_seg[i] = {(a + b) / 2 + T(0.02), uv[0], uv[1]};
    }

    std::printf("%-26s %13s %16s\n", "case", "time/call", "evals/call");
    run({"PC circle (global)", N, [&](int i) { const auto r = extrema_curve_point(circle_, pc[i], tol); return norm(circle_(r[0]) - pc[i]); },
         [&](int i) { return pc_d[i]; }});
    run({"PC circle (local)", N,
         [&](int i) { const auto r = extrema_curve_point(circle, pc[i], u0_circle[i], tol); return norm(circle_(r[0]) - pc[i]); },
         [&](int i) { return pc_d[i]; }});
    run({"PC B-spline (global)", N, [&](int i) {
             const auto p = P{pc[i][0] + T(0.5), pc[i][1] + T(0.5), T(4) * pc[i][2] + T(2)};
             const auto r = extrema_curve_point(bspline_, p, tol);
             return norm(bspline_(r[0]) - p);
         }});
    run({"PS sphere (global)", N, [&](int i) { const auto r = extrema_surf_pnt(sphere_, ps[i], tol); return norm(sphere_(r[0], r[1]) - ps[i]); },
         [&](int i) { return ps_d[i]; }});
    run({"PS sphere (local)", N,
         [&](int i) {
             const auto r = extrema_surf_pnt(sphere, ps[i], uv0_sphere[i][0], uv0_sphere[i][1], tol);
             return norm(sphere_(r[0], r[1]) - ps[i]);
         },
         [&](int i) { return ps_d[i]; }});
    run({"PS torus (global)", N, [&](int i) { const auto r = extrema_surf_pnt(torus_, pt[i], tol); return norm(torus_(r[0], r[1]) - pt[i]); },
         [&](int i) { return pt_d[i]; }});
    run({"CC circle-segment (global)", N, [&](int i) {
             const auto r = extrema_curve_curve(circle_, seg_c[i], tol);
             return norm(circle_(r[0]) - seg_c[i](r[1]));
         }});
    run({"CC B-splines (global)", 20, [&](int i) {
             const auto r = extrema_curve_curve(bspline_, bspline2_, tol, default_nlopt_algo, 30 + i, 30);
             return norm(bspline_(r[0]) - bspline2_(r[1]));
         }});
    run({"CS sphere-segment (global)", NS,
         [&](int i) { const auto r = extrema_surf_curve(sphere_, seg_c[i], tol); return norm(sphere_(r[0], r[1]) - seg_c[i](r[2])); },
         [&](int i) { return seg_d[i]; }});
    run({"CS sphere-segment (local)", N,
         [&](int i) {
             const auto &x0 = tuv0_seg[i];
             const auto r = extrema_surf_curve(sphere, segs[i], x0[0], x0[1], x0[2], tol);
             return norm(sphere_(r[0], r[1]) - seg_c[i](r[2]));
         },
         [&](int i) { return seg_d[i]; }});
    if (out)
        std::fclose(out);
}
