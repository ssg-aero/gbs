#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>
#include <gbs/bscinterp.h>
#include <gbs/bscapprox.h>
#include <gbs/execution.h>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

// Coverage for #91: the batch interp/approx builders run the outer loop over N
// INDEPENDENT entities through gbs::build_batch, which routes the loop to the
// parallel execution policy only above gbs::parallel_batch_min_size. The gate
// only chooses the policy by the build COUNT; it must NEVER change the result.
// Each build constructs its own object (own solver/poles, no shared state) and
// writes a distinct pre-indexed slot with no reduction, so every batch result is
// bit-identical to the per-entity serial loop — checked here with EXACT equality
// (operator==, no tolerance). Both gate branches are exercised by pinning
// gbs::parallel_batch_min_size to 0 (always parallel) / SIZE_MAX (always serial).

namespace
{
    // RAII override of the global batch gate so a failing assertion cannot leak a
    // mutated threshold into the next test case.
    struct scoped_batch_min_size
    {
        std::size_t saved;
        explicit scoped_batch_min_size(std::size_t v)
            : saved(gbs::parallel_batch_min_size) { gbs::parallel_batch_min_size = v; }
        ~scoped_batch_min_size() { gbs::parallel_batch_min_size = saved; }
    };

    // K distinct point sets: a helix whose phase shifts per entity, so every build
    // produces a genuinely different curve (not K copies of one solve).
    std::vector<gbs::points_vector<double, 3>> make_inputs(std::size_t K, std::size_t n_pts)
    {
        std::vector<gbs::points_vector<double, 3>> inputs(K);
        for (std::size_t k = 0; k < K; ++k)
        {
            auto &Q = inputs[k];
            Q.resize(n_pts);
            const double ph = double(k) / double(K == 0 ? 1 : K);
            for (std::size_t i = 0; i < n_pts; ++i)
            {
                const double t = double(i) / double(n_pts - 1);
                Q[i] = {std::cos(6.0 * t + ph), std::sin(6.0 * t + ph), 2.0 * t};
            }
        }
        return inputs;
    }

    // Exact, component-by-component equality of two curves (poles + knots + degree).
    void expect_same_curve(const gbs::BSCurve<double, 3> &got,
                           const gbs::BSCurve<double, 3> &ref, std::size_t k)
    {
        ASSERT_EQ(got.degree(), ref.degree());
        ASSERT_EQ(got.poles().size(), ref.poles().size());
        ASSERT_EQ(got.knotsFlats().size(), ref.knotsFlats().size());
        for (std::size_t i = 0; i < ref.poles().size(); ++i)
            EXPECT_EQ(got.poles()[i], ref.poles()[i]) << "k=" << k << " pole=" << i;
        for (std::size_t i = 0; i < ref.knotsFlats().size(); ++i)
            EXPECT_EQ(got.knotsFlats()[i], ref.knotsFlats()[i]) << "k=" << k << " knot=" << i;
    }

    // Sizes straddling any plausible gate: small (serial branch) and large enough
    // to exercise the real parallel branch when the gate is forced to 0.
    constexpr std::size_t K_small = 5;
    constexpr std::size_t K_large = 200;
}

// ---- the generic helper ----------------------------------------------------

TEST(tests_batch_build, build_batch_matches_serial)
{
    std::vector<int> in(K_large);
    for (std::size_t i = 0; i < in.size(); ++i) in[i] = int(i);
    const auto build = [](int x) { return x * x - 3 * x; };

    std::vector<int> ref(in.size());
    std::transform(in.begin(), in.end(), ref.begin(), build);

    for (std::size_t gate : {std::size_t(0), std::size_t(32), std::numeric_limits<std::size_t>::max()})
    {
        const auto got = gbs::build_batch(in, build, gate);
        EXPECT_EQ(got, ref) << "gate=" << gate;
    }
}

TEST(tests_batch_build, build_batch_empty_input)
{
    const std::vector<int> in;
    const auto got = gbs::build_batch(in, [](int x) { return x; });
    EXPECT_TRUE(got.empty());
}

// ---- batch interpolate -----------------------------------------------------

TEST(tests_batch_build, batch_interpolate_bit_identical)
{
    for (std::size_t K : {K_small, K_large})
    {
        const auto inputs = make_inputs(K, 80);

        // Per-entity serial reference (the loop the batch API replaces).
        std::vector<gbs::BSCurve<double, 3>> ref;
        ref.reserve(K);
        for (const auto &Q : inputs)
            ref.push_back(gbs::interpolate<double, 3>(Q, std::size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH));

        for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
        {
            scoped_batch_min_size guard(gate);
            const auto got = gbs::interpolate<double, 3>(inputs, std::size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH);
            ASSERT_EQ(got.size(), ref.size());
            for (std::size_t k = 0; k < K; ++k)
                expect_same_curve(got[k], ref[k], k);
        }
    }
}

// A malformed entity must surface as a CATCHABLE exception on either side of the
// gate — never std::terminate (which the parallel branch would call if a build's
// exception were allowed to escape the element function). With the gate forced to
// 0 (always parallel) the un-trapped version of this would crash the test process.
TEST(tests_batch_build, batch_interpolate_propagates_exception)
{
    auto inputs = make_inputs(K_large, 80);
    inputs[K_large / 2].resize(2); // too few points for a degree-3 interpolation

    // Wrap in a lambda so the comma in interpolate<double, 3> does not split the
    // EXPECT_ANY_THROW macro arguments.
    const auto run = [&] { return gbs::interpolate<double, 3>(inputs, std::size_t{3}, gbs::KnotsCalcMode::CHORD_LENGTH); };
    for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
    {
        scoped_batch_min_size guard(gate);
        EXPECT_ANY_THROW((void)run()) << "gate=" << gate;
    }
}

// ---- batch approx ----------------------------------------------------------

TEST(tests_batch_build, batch_approx_bit_identical)
{
    for (std::size_t K : {K_small, K_large})
    {
        const auto inputs = make_inputs(K, 400);

        std::vector<gbs::BSCurve<double, 3>> ref;
        ref.reserve(K);
        for (const auto &Q : inputs)
            ref.push_back(gbs::approx<double, 3>(Q, std::size_t{3}, std::size_t{50}, gbs::KnotsCalcMode::CHORD_LENGTH));

        for (std::size_t gate : {std::size_t(0), std::numeric_limits<std::size_t>::max()})
        {
            scoped_batch_min_size guard(gate);
            const auto got = gbs::approx<double, 3>(inputs, std::size_t{3}, std::size_t{50}, gbs::KnotsCalcMode::CHORD_LENGTH);
            ASSERT_EQ(got.size(), ref.size());
            for (std::size_t k = 0; k < K; ++k)
                expect_same_curve(got[k], ref[k], k);
        }
    }
}
