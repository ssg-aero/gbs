#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/execution.h>

#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Coverage for #88: bulk evaluators route std::transform through the size-gated
// gbs::transform_threshold helper. The gate only chooses the execution policy by
// input size; it must NEVER change the result. Every check below compares the
// bulk path to a per-element serial reference with EXACT equality (operator==,
// no tolerance): each output slot is written independently with no reduction, so
// the result is bit-identical regardless of policy or thread count.
//
// To exercise BOTH branches of the gate deterministically we drive the
// evaluators with input sizes below and above the threshold AND temporarily
// pin gbs::parallel_min_size to 0 (always parallel) / SIZE_MAX (always serial).

namespace
{
    // RAII override of the global gate so a failing assertion cannot leak a
    // mutated threshold into the next test case.
    struct scoped_min_size
    {
        std::size_t saved;
        explicit scoped_min_size(std::size_t v) : saved(gbs::parallel_min_size) { gbs::parallel_min_size = v; }
        ~scoped_min_size() { gbs::parallel_min_size = saved; }
    };

    // Degree-3 open curve with a non-trivial control net (the light "value"
    // kernel of the audit).
    gbs::BSCurve<double, 3> make_curve()
    {
        const std::vector<double> knots{0., 0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1., 1.};
        const gbs::points_vector<double, 3> poles{
            {0., 0., 0.}, {1., 2., 0.}, {2., -1., 1.}, {3., 1., 2.},
            {4., 0., 1.}, {5., 2., 0.}, {6., -1., 3.},
        };
        return gbs::BSCurve<double, 3>(poles, knots, 3);
    }

    // Bidegree-(3,3) open surface.
    gbs::BSSurface<double, 3> make_surface()
    {
        const std::vector<double> ku{0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
        const std::vector<double> kv{0., 0., 0., 0., 0.5, 1., 1., 1., 1.};
        gbs::points_vector<double, 3> poles;
        poles.reserve(5 * 5);
        for (int j = 0; j < 5; ++j)
            for (int i = 0; i < 5; ++i)
                poles.push_back({double(i), double(j), std::sin(0.7 * i) + std::cos(0.5 * j)});
        return gbs::BSSurface<double, 3>(poles, ku, kv, 3, 3);
    }

    std::vector<double> make_params(std::size_t n, std::array<double, 2> bnds)
    {
        std::vector<double> u(n);
        const double span = bnds[1] - bnds[0];
        for (std::size_t i = 0; i < n; ++i)
            u[i] = bnds[0] + span * (double(i) + 0.5) / double(n);
        return u;
    }

    std::vector<std::array<double, 2>> make_uv(std::size_t n, std::array<double, 4> bnds)
    {
        std::vector<std::array<double, 2>> uv(n);
        const double su = bnds[1] - bnds[0];
        const double sv = bnds[3] - bnds[2];
        for (std::size_t i = 0; i < n; ++i)
        {
            const double t = (double(i) + 0.5) / double(n);
            uv[i] = {bnds[0] + su * t, bnds[2] + sv * (1.0 - t)};
        }
        return uv;
    }

    // Sizes straddling the default gate (gbs::parallel_min_size == 1000).
    constexpr std::size_t N_small = 50;    // below the gate -> serial branch
    constexpr std::size_t N_large = 4096;  // above the gate -> parallel branch
}

// ---- the helper itself -----------------------------------------------------

TEST(tests_parallel_eval, transform_threshold_unary_matches_serial)
{
    std::vector<double> in(N_large);
    std::iota(in.begin(), in.end(), 0.0);
    const auto op = [](double x) { return std::sin(x) * 2.0 - x; };

    std::vector<double> ref(in.size());
    std::transform(in.begin(), in.end(), ref.begin(), op);

    for (std::size_t gate : {std::size_t(0), std::size_t(1000), std::numeric_limits<std::size_t>::max()})
    {
        scoped_min_size guard(gate);
        std::vector<double> out(in.size());
        gbs::transform_threshold(in.begin(), in.end(), out.begin(), op);
        EXPECT_EQ(out, ref) << "gate=" << gate; // exact, policy must not change result
    }
}

TEST(tests_parallel_eval, transform_threshold_binary_matches_serial)
{
    std::vector<double> a(N_large), b(N_large);
    std::iota(a.begin(), a.end(), 0.0);
    std::iota(b.begin(), b.end(), 1.0);
    const auto op = [](double x, double y) { return x * y - std::cos(x); };

    std::vector<double> ref(a.size());
    std::transform(a.begin(), a.end(), b.begin(), ref.begin(), op);

    for (std::size_t gate : {std::size_t(0), std::size_t(1000), std::numeric_limits<std::size_t>::max()})
    {
        scoped_min_size guard(gate);
        std::vector<double> out(a.size());
        gbs::transform_threshold(a.begin(), a.end(), b.begin(), out.begin(), op);
        EXPECT_EQ(out, ref) << "gate=" << gate;
    }
}

// ---- curve bulk evaluators -------------------------------------------------

namespace
{
    template <class BulkFn, class RefFn>
    void check_curve_eval(BulkFn bulk, RefFn ref_elem)
    {
        const auto crv = make_curve();
        for (std::size_t n : {N_small, N_large})
        {
            const auto u = make_params(n, crv.bounds());
            // Per-element serial reference.
            gbs::points_vector<double, 3> ref(n);
            std::transform(u.begin(), u.end(), ref.begin(),
                           [&](double uu) { return ref_elem(crv, uu); });
            // Bulk path with the gate forced to both branches.
            for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
            {
                scoped_min_size guard(gate);
                const auto got = bulk(crv, u);
                ASSERT_EQ(got.size(), ref.size());
                for (std::size_t i = 0; i < n; ++i)
                    EXPECT_EQ(got[i], ref[i]) << "n=" << n << " gate=" << gate << " i=" << i;
            }
        }
    }
}

TEST(tests_parallel_eval, curve_values_bit_identical)
{
    check_curve_eval(
        [](const auto &c, const auto &u) { return c.values(u, 0); },
        [](const auto &c, double u) { return c.value(u, 0); });
}

TEST(tests_parallel_eval, curve_values_first_derivative_bit_identical)
{
    check_curve_eval(
        [](const auto &c, const auto &u) { return c.values(u, 1); },
        [](const auto &c, double u) { return c.value(u, 1); });
}

TEST(tests_parallel_eval, curve_d_dms_bit_identical)
{
    check_curve_eval(
        [](const auto &c, const auto &u) { return c.d_dms(u); },
        [](const auto &c, double u) { return c.d_dm(u); });
}

TEST(tests_parallel_eval, curve_d_dm2s_bit_identical)
{
    check_curve_eval(
        [](const auto &c, const auto &u) { return c.d_dm2s(u); },
        [](const auto &c, double u) { return c.d_dm2(u); });
}

// ---- surface bulk evaluators -----------------------------------------------

namespace
{
    template <class BulkFn, class RefFn>
    void check_surface_eval(BulkFn bulk, RefFn ref_elem)
    {
        const auto srf = make_surface();
        for (std::size_t n : {N_small, N_large})
        {
            const auto uv = make_uv(n, srf.bounds());
            gbs::points_vector<double, 3> ref(n);
            std::transform(uv.begin(), uv.end(), ref.begin(),
                           [&](const auto &p) { return ref_elem(srf, p); });
            for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
            {
                scoped_min_size guard(gate);
                const auto got = bulk(srf, uv);
                ASSERT_EQ(got.size(), ref.size());
                for (std::size_t i = 0; i < n; ++i)
                    EXPECT_EQ(got[i], ref[i]) << "n=" << n << " gate=" << gate << " i=" << i;
            }
        }
    }
}

TEST(tests_parallel_eval, surface_values_bit_identical)
{
    check_surface_eval(
        [](const auto &s, const auto &uv) { return s.values(uv, 0, 0); },
        [](const auto &s, const auto &p) { return s.value(p[0], p[1], 0, 0); });
}

TEST(tests_parallel_eval, surface_d_dmus_bit_identical)
{
    check_surface_eval(
        [](const auto &s, const auto &uv) { return s.d_dmus(uv); },
        [](const auto &s, const auto &p) { return s.d_dmu(p[0], p[1]); });
}

TEST(tests_parallel_eval, surface_d_dmvs_bit_identical)
{
    check_surface_eval(
        [](const auto &s, const auto &uv) { return s.d_dmvs(uv); },
        [](const auto &s, const auto &p) { return s.d_dmv(p[0], p[1]); });
}

TEST(tests_parallel_eval, surface_d_dmu2s_bit_identical)
{
    check_surface_eval(
        [](const auto &s, const auto &uv) { return s.d_dmu2s(uv); },
        [](const auto &s, const auto &p) { return s.d_dmu2(p[0], p[1]); });
}

TEST(tests_parallel_eval, surface_d_dmv2s_bit_identical)
{
    check_surface_eval(
        [](const auto &s, const auto &uv) { return s.d_dmv2s(uv); },
        [](const auto &s, const auto &p) { return s.d_dmv2(p[0], p[1]); });
}

// The two-range surface values(u_lst, v_lst) overload.
TEST(tests_parallel_eval, surface_values_two_ranges_bit_identical)
{
    const auto srf = make_surface();
    const auto bnds = srf.bounds();
    for (std::size_t n : {N_small, N_large})
    {
        const auto u = make_params(n, {bnds[0], bnds[1]});
        const auto v = make_params(n, {bnds[2], bnds[3]});
        gbs::points_vector<double, 3> ref(n);
        for (std::size_t i = 0; i < n; ++i)
            ref[i] = srf.value(u[i], v[i], 0, 0);
        for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
        {
            scoped_min_size guard(gate);
            const auto got = srf.values(u, v, 0, 0);
            ASSERT_EQ(got.size(), ref.size());
            for (std::size_t i = 0; i < n; ++i)
                EXPECT_EQ(got[i], ref[i]) << "n=" << n << " gate=" << gate << " i=" << i;
        }
    }
}
