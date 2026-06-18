// Competitive benchmark — gbs vs OpenCASCADE (OCCT), C++ core front (#86).
//
// Measures gbs against OCCT (the primary C++ baseline, already bridged via
// `occt_utils`) on the same B-spline operations, on one machine, same flags,
// best-of-N wall time. The goal is to POSITION gbs honestly: every number is
// reported with its exact conditions, and the lags are in the same tables as
// the wins. No cherry-picking.
//
// Fairness rules applied here:
//   * Construction / gbs->OCCT conversion is done ONCE, OUTSIDE every timed
//     region (except where the operation IS a build: interpolation/approx).
//   * Same instances feed both libraries (same poles/knots/degree/param list).
//   * Evaluation is reported on TWO regimes, because the libraries expose
//     different surfaces:
//       - "per-call": a scalar C++ loop calling value(u) / Value(u). OCCT has
//         no batch evaluator, so this is the like-for-like evaluator quality
//         comparison (gbs compile-time double/dim vs OCCT runtime dispatch).
//       - "bulk": gbs's batched evaluator values(u_lst) — parallel since #88 —
//         vs OCCT's only option, the same scalar loop. This is the API-level
//         comparison ("how a user evaluates N points").
//   * In-place ops (knot insert, degree elevate) time a fresh COPY + the op on
//     BOTH sides (a user must copy to keep the original); copy is included on
//     both, stated.
//
// Structural caveat reported in the audit: gbs is templated on <double, dim=3>
// (compile-time), OCCT is runtime-typed double/3D. Both are -O3 Release.
//
// Build/run:  scripts/compare_bench.sh   (needs the gbs-occt8 micromamba env)

#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/bscapprox.h>

#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/surfacesbuild.h>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <Standard_Version.hxx>

#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/info.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

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

// ----------------------------------------------------------------------------
// sample geometry (shared by both libraries)
// ----------------------------------------------------------------------------
static gbs::points_vector<double, 3> sample_curve_pts(size_t n)
{
    gbs::points_vector<double, 3> Q(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = double(i) / double(n - 1);
        Q[i] = {std::cos(6.0 * t), std::sin(6.0 * t), 2.0 * t};
    }
    return Q;
}

static std::vector<std::array<double, 3>> sample_grid_pts(size_t nu, size_t nv)
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

static std::vector<double> param_list(std::array<double, 2> b, size_t n)
{
    std::vector<double> u(n);
    for (size_t i = 0; i < n; ++i)
        u[i] = b[0] + (b[1] - b[0]) * (double(i) + 0.5) / double(n);
    return u;
}

static std::vector<std::array<double, 2>> uv_list(std::array<double, 4> b, size_t n)
{
    std::vector<std::array<double, 2>> uv(n);
    for (size_t i = 0; i < n; ++i)
    {
        double a = (double(i) + 0.5) / double(n);
        double c = std::fmod(a * 7.0, 1.0);
        uv[i] = {b[0] + (b[1] - b[0]) * a, b[2] + (b[3] - b[2]) * c};
    }
    return uv;
}

// ----------------------------------------------------------------------------
static void hdr(const char *title)
{
    std::printf("\n%s\n", title);
    std::printf("%-12s %14s %14s %12s %10s\n", "N", "gbs ms", "occt ms", "ratio", "winner");
    std::printf("--------------------------------------------------------------------\n");
}
static void row(size_t N, double gbs_ms, double occ_ms)
{
    double ratio = occ_ms / gbs_ms;                  // >1 => gbs faster
    const char *who = ratio > 1.0 ? "gbs" : "occt";
    std::printf("%-12zu %14.4f %14.4f %11.2fx %10s\n", N, gbs_ms, occ_ms, ratio, who);
}

int main()
{
    using namespace gbs;
    const auto mode = KnotsCalcMode::CHORD_LENGTH;
    unsigned hw = oneapi::tbb::info::default_concurrency();
    std::printf("=== gbs vs OCCT %s — C++ core front (#86) ===\n", OCC_VERSION_COMPLETE);
    std::printf("hardware concurrency (TBB) = %u ; best-of-N wall time ; -O3 Release\n", hw);

    // ---- build the shared curve / surface ONCE (outside timed eval regions) -
    auto crv = interpolate<double, 3>(sample_curve_pts(60), size_t{3}, mode);
    auto bc = crv.bounds();
    auto occ_crv = occt_utils::BSplineCurve(crv); // Handle(Geom_BSplineCurve)

    auto srf = interpolate<double, 3>(sample_grid_pts(40, 40), size_t{40}, size_t{3}, size_t{3}, mode);
    auto bs = srf.bounds();
    auto occ_srf = occt_utils::BSplineSurface(srf); // Handle(Geom_BSplineSurface)

    // =====================================================================
    // 1. CURVE point evaluation — per-call (scalar loop, like-for-like)
    // =====================================================================
    hdr("[1] CURVE value(u) d=0  — PER-CALL scalar loop (gbs single value() vs OCCT Value())");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto u = param_list(bc, N);
        int reps = N <= 100000 ? 7 : 3;
        double g = best_ms(reps, [&] { double s = 0; for (double x : u) s += crv.value(x)[0]; sink += s; });
        double o = best_ms(reps, [&] { double s = 0; for (double x : u) s += occ_crv->Value(x).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 2. CURVE point evaluation — BULK (gbs parallel values() #88 vs OCCT loop)
    // =====================================================================
    hdr("[2] CURVE value(u) d=0  — BULK (gbs values(u_lst) parallel #88 vs OCCT scalar loop)");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto u = param_list(bc, N);
        int reps = N <= 100000 ? 7 : 3;
        points_vector<double, 3> r;
        double g = best_ms(reps, [&] { r = crv.values(u, 0); sink += r[0][0]; });
        double o = best_ms(reps, [&] { double s = 0; for (double x : u) s += occ_crv->Value(x).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 3. CURVE 1st derivative — BULK vs OCCT DN(u,1) loop
    // =====================================================================
    hdr("[3] CURVE 1st derivative  — BULK (gbs values(u_lst,1) vs OCCT DN(u,1) loop)");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto u = param_list(bc, N);
        int reps = N <= 100000 ? 7 : 3;
        points_vector<double, 3> r;
        double g = best_ms(reps, [&] { r = crv.values(u, 1); sink += r[0][0]; });
        double o = best_ms(reps, [&] { double s = 0; for (double x : u) s += occ_crv->DN(x, 1).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 4. CURVE 2nd derivative — BULK vs OCCT DN(u,2) loop
    // =====================================================================
    hdr("[4] CURVE 2nd derivative  — BULK (gbs values(u_lst,2) vs OCCT DN(u,2) loop)");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto u = param_list(bc, N);
        int reps = N <= 100000 ? 7 : 3;
        points_vector<double, 3> r;
        double g = best_ms(reps, [&] { r = crv.values(u, 2); sink += r[0][0]; });
        double o = best_ms(reps, [&] { double s = 0; for (double x : u) s += occ_crv->DN(x, 2).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 5. SURFACE point evaluation — BULK (gbs parallel vs OCCT loop)
    // =====================================================================
    hdr("[5] SURFACE value(u,v) d=0  — BULK (gbs values(uv_lst) parallel #88 vs OCCT Value loop)");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto uv = uv_list(bs, N);
        int reps = N <= 100000 ? 7 : 3;
        points_vector<double, 3> r;
        double g = best_ms(reps, [&] { r = srf.values(uv, 0, 0); sink += r[0][0]; });
        double o = best_ms(reps, [&] { double s = 0; for (auto &p : uv) s += occ_srf->Value(p[0], p[1]).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 6. SURFACE 1st derivative d/du — BULK vs OCCT DN(u,v,1,0) loop
    // =====================================================================
    hdr("[6] SURFACE d/du  — BULK (gbs values(uv_lst,1,0) vs OCCT DN(u,v,1,0) loop)");
    for (size_t N : {1000u, 10000u, 100000u, 1000000u})
    {
        auto uv = uv_list(bs, N);
        int reps = N <= 100000 ? 7 : 3;
        points_vector<double, 3> r;
        double g = best_ms(reps, [&] { r = srf.values(uv, 1, 0); sink += r[0][0]; });
        double o = best_ms(reps, [&] { double s = 0; for (auto &p : uv) s += occ_srf->DN(p[0], p[1], 1, 0).X(); sink += s; });
        row(N, g, o);
    }

    // =====================================================================
    // 7. INTERPOLATION — the operation IS a build; whole build timed.
    //    Includes the real 201-point case (#34). gbs CHORD_LENGTH vs OCCT
    //    GeomAPI_Interpolate (its own param scheme) — both interpolate, deg 3.
    // =====================================================================
    hdr("[7] CURVE interpolation deg 3  — full build (gbs interpolate vs OCCT GeomAPI_Interpolate)");
    for (size_t N : {50u, 100u, 201u, 400u, 800u})
    {
        auto Q = sample_curve_pts(N);
        int reps = N <= 201 ? 7 : 3;
        double g = best_ms(reps, [&] {
            auto c = interpolate<double, 3>(Q, size_t{3}, mode);
            sink += c.poles()[0][0];
        });
        // OCCT input array (built outside the timed region — matches gbs, whose
        // input vector is also pre-built; both time only the solve+assembly).
        Handle(TColgp_HArray1OfPnt) pts = new TColgp_HArray1OfPnt(1, int(N));
        for (size_t i = 0; i < N; ++i) pts->SetValue(int(i) + 1, gp_Pnt(Q[i][0], Q[i][1], Q[i][2]));
        double o = best_ms(reps, [&] {
            GeomAPI_Interpolate interp(pts, Standard_False, 1.0e-7);
            interp.Perform();
            auto c = interp.Curve();
            sink += c->Pole(1).X();
        });
        row(N, g, o);
    }

    // =====================================================================
    // 8. APPROXIMATION (least-squares fit) — STRUCTURAL caveat: gbs fits a
    //    FIXED pole count; OCCT GeomAPI_PointsToBSpline picks knots adaptively
    //    from a tolerance. Not bit-for-bit like-for-like; we fix deg=3 on both
    //    and give OCCT a loose tol so it stays near the same pole budget.
    //    Reported with the caveat (see audit).
    // =====================================================================
    hdr("[8] CURVE approximation deg 3 (~N/4 poles)  — gbs approx vs OCCT PointsToBSpline (tol=1e-3)");
    for (size_t N : {100u, 201u, 400u, 800u})
    {
        auto Q = sample_curve_pts(N);
        size_t n_poles = std::max<size_t>(8, N / 4);
        int reps = N <= 201 ? 7 : 3;
        double g = best_ms(reps, [&] {
            auto c = approx(Q, size_t{3}, n_poles, mode);
            sink += c.poles()[0][0];
        });
        TColgp_Array1OfPnt pts(1, int(N));
        for (size_t i = 0; i < N; ++i) pts.SetValue(int(i) + 1, gp_Pnt(Q[i][0], Q[i][1], Q[i][2]));
        double o = best_ms(reps, [&] {
            GeomAPI_PointsToBSpline fit(pts, 3, 3, GeomAbs_C2, 1.0e-3);
            auto c = fit.Curve();
            sink += c->Pole(1).X();
        });
        row(N, g, o);
    }

    // =====================================================================
    // 9. KNOT INSERTION — fresh copy + one insertion, both sides (copy incl.)
    // =====================================================================
    {
        std::printf("\n[9] KNOT insertion (fresh copy + insert 1 knot ; copy included both sides)\n");
        std::printf("%-12s %14s %14s %12s %10s\n", "reps", "gbs ms", "occt ms", "ratio", "winner");
        std::printf("--------------------------------------------------------------------\n");
        double um = 0.5 * (bc[0] + bc[1]);
        const int R = 20000;
        double g = best_ms(5, [&] {
            for (int i = 0; i < R; ++i) { auto c = crv; c.insertKnot(um, 1); sink += c.poles()[0][0]; }
        });
        double o = best_ms(5, [&] {
            for (int i = 0; i < R; ++i) {
                Handle(Geom_BSplineCurve) c = Handle(Geom_BSplineCurve)::DownCast(occ_crv->Copy());
                c->InsertKnot(um, 1, 1e-9); sink += c->Pole(1).X();
            }
        });
        std::printf("%-12d %14.4f %14.4f %11.2fx %10s\n", R, g, o, o / g, o / g > 1 ? "gbs" : "occt");
    }

    // =====================================================================
    // 10. DEGREE ELEVATION — fresh copy + elevate by 1, both sides
    // =====================================================================
    {
        std::printf("\n[10] DEGREE elevation (fresh copy + elevate +1 ; copy included both sides)\n");
        std::printf("%-12s %14s %14s %12s %10s\n", "reps", "gbs ms", "occt ms", "ratio", "winner");
        std::printf("--------------------------------------------------------------------\n");
        const int R = 5000;
        size_t deg = crv.degree();
        double g = best_ms(5, [&] {
            for (int i = 0; i < R; ++i) { auto c = crv; c.increaseDegree(1); sink += c.poles()[0][0]; }
        });
        double o = best_ms(5, [&] {
            for (int i = 0; i < R; ++i) {
                Handle(Geom_BSplineCurve) c = Handle(Geom_BSplineCurve)::DownCast(occ_crv->Copy());
                c->IncreaseDegree(int(deg) + 1); sink += c->Pole(1).X();
            }
        });
        std::printf("%-12d %14.4f %14.4f %11.2fx %10s\n", R, g, o, o / g, o / g > 1 ? "gbs" : "occt");
    }

    std::printf("\n[sink=%g]\n", sink);
    return 0;
}
