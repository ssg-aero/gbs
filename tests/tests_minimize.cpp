#include <doctest_gtest.hpp>
#include <gbs/minimize.h>
#include <gbs/curves>
#include <gbs/extrema.h>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <numbers>
#include <stdexcept>
#include <string>

namespace
{
    constexpr auto inf = std::numeric_limits<double>::infinity();
}

TEST(tests_minimize, rosenbrock)
{
    auto f = [](const std::array<double, 2> &x) { return 100. * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1. - x[0], 2); };
    auto r = gbs::nelder_mead(f, std::array{-1.2, 1.}, {-5., -5.}, {5., 5.}, 1e-10);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 1., 1e-7);
    ASSERT_NEAR(r.x[1], 1., 1e-7);
    ASSERT_LT(r.f, 1e-14);
    ASSERT_DOUBLE_EQ(r.f, f(r.x));
}

TEST(tests_minimize, quadratic_5d)
{
    // ill-scaled and coupled
    const std::array<double, 5> xs{0.1, -0.2, 0.3, 0.4, -0.5};
    auto f = [&](const std::array<double, 5> &x) {
        double s{};
        for (std::size_t i{}; i < 5; ++i)
            s += double(i + 1) * std::pow(x[i] - xs[i], 2);
        return s + (x[0] - xs[0]) * (x[1] - xs[1]);
    };
    auto r = gbs::nelder_mead(f, std::array<double, 5>{}, {-1., -1., -1., -1., -1.}, {1., 1., 1., 1., 1.}, 1e-10);
    ASSERT_TRUE(r.converged);
    for (std::size_t i{}; i < 5; ++i)
        ASSERT_NEAR(r.x[i], xs[i], 1e-7);
}

TEST(tests_minimize, active_bounds)
{
    // minimum outside the box: the corner (1, 0), reached although the simplex is clamped on the bounds
    auto f = [](const std::array<double, 2> &x) { return std::pow(x[0] - 2., 2) + std::pow(x[1] + 1., 2) + x[0] * x[1]; };
    auto r = gbs::nelder_mead(f, std::array{0.5, 0.5}, {0., 0.}, {1., 1.}, 1e-10);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 1., 1e-9);
    ASSERT_NEAR(r.x[1], 0., 1e-9);
}

TEST(tests_minimize, one_variable_not_smooth)
{
    // sums of absolute values (angle errors) have no derivative at their minimum
    auto f = [](const std::array<double, 1> &x) { return std::abs(x[0] - 0.3); };
    auto r = gbs::nelder_mead(f, std::array{0.9}, {0.}, {1.}, 1e-12);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 0.3, 1e-11);
}

TEST(tests_minimize, start_on_upper_bound)
{
    auto f = [](const std::array<double, 1> &x) { return std::pow(x[0] - 0.25, 2); };
    auto r = gbs::nelder_mead(f, std::array{3.}, {0.}, {1.}, 1e-10); // x0 clamped to 1
    ASSERT_NEAR(r.x[0], 0.25, 1e-8);
}

TEST(tests_minimize, unbounded_and_held_variables)
{
    auto f = [](const std::array<double, 3> &x) { return std::pow(x[0] - 40., 2) + std::pow(x[1] - 1., 2) + std::pow(x[2] - 1., 2); };
    // x[1] held by lo == hi
    auto r = gbs::nelder_mead(f, std::array{0., 2., 0.}, {-inf, 2., -inf}, {inf, 2., inf}, 1e-10);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 40., 1e-6);
    ASSERT_DOUBLE_EQ(r.x[1], 2.);
    ASSERT_NEAR(r.x[2], 1., 1e-8);
}

TEST(tests_minimize, nan_is_rejected)
{
    // f undefined below 0.2: the simplex steps there and is pushed back
    auto f = [](const std::array<double, 1> &x) { return x[0] < 0.2 ? std::nan("") : std::pow(x[0] - 0.5, 2); };
    auto r = gbs::nelder_mead(f, std::array{0.21}, {0.}, {1.}, 1e-10);
    ASSERT_NEAR(r.x[0], 0.5, 1e-8);
}

TEST(tests_minimize, nan_at_start)
{
    auto f = [](const std::array<double, 1> &x) { return x[0] < 0.2 ? std::nan("") : std::pow(x[0] - 0.5, 2); };
    auto r = gbs::nelder_mead(f, std::array{0.15}, {0.}, {1.}, 1e-10);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 0.5, 1e-8);
}

TEST(tests_minimize, nowhere_defined)
{
    auto f = [](const std::array<double, 2> &) { return std::nan(""); };
    auto r = gbs::nelder_mead(f, std::array{0.5, 0.5}, {0., 0.}, {1., 1.}, 1e-6);
    ASSERT_FALSE(r.converged);
    ASSERT_EQ(r.f, inf);
}

TEST(tests_minimize, corner_with_kink)
{
    // minimum (0.02, 1): x[1] on its bound, x[0] at a kink next to its bound; clamping the
    // trial points instead of changing variables lays the simplex on x[0] = 0 for good
    auto f = [](const std::array<double, 2> &x) {
        return std::pow(x[0] - 0.02, 2) + std::pow(x[1] - 1.1, 2) + std::abs(x[0] - 0.02) + 0.5 * (x[0] - 0.02) * (x[1] - 1.1);
    };
    auto r = gbs::nelder_mead(f, std::array{0.3, 0.5}, {0., 0.}, {1., 1.}, 1e-10);
    ASSERT_TRUE(r.converged);
    ASSERT_NEAR(r.x[0], 0.02, 1e-8);
    ASSERT_NEAR(r.x[1], 1., 1e-8);
}

TEST(tests_minimize, restarts)
{
    // separable sums of w (x - xs)² + |x - xs|, minimum at xs clamped into the box: without the
    // restarts, about one problem in thirteen stops short of it (up to 4e-4 away)
    std::mt19937 g(7);
    auto u = [&] { return double(g()) / 4294967296.; }; // the same on every standard library
    for (int k{}; k < 40; ++k)
    {
        std::array<double, 3> xs, x0, w;
        for (std::size_t i{}; i < 3; ++i)
        {
            const double a = u();
            xs[i] = a < 0.3 ? -0.2 * u() : (a < 0.6 ? 1. + 0.2 * u() : u());
            x0[i] = u();
            w[i] = std::pow(10., 2. * u() - 1.);
        }
        auto f = [&](const std::array<double, 3> &x) {
            double s{};
            for (std::size_t i{}; i < 3; ++i)
                s += w[i] * std::pow(x[i] - xs[i], 2) + std::abs(x[i] - xs[i]);
            return s;
        };
        auto r = gbs::nelder_mead(f, x0, {0., 0., 0.}, {1., 1., 1.}, 1e-10);
        ASSERT_TRUE(r.converged);
        for (std::size_t i{}; i < 3; ++i)
            ASSERT_NEAR(r.x[i], std::clamp(xs[i], 0., 1.), 1e-7);
    }
}

TEST(tests_minimize, initial_step)
{
    // a tilted double well: from x0 = 0.9, a tenth of the box stays in the right well, half of
    // it reaches the deeper left one
    auto f = [](const std::array<double, 1> &x) { return std::pow(x[0] * x[0] - 1., 2) + 0.3 * x[0]; };
    const auto near = gbs::nelder_mead(f, std::array{0.9}, {-2.}, {2.}, 1e-10);
    ASSERT_NEAR(near.x[0], 0.96015, 1e-5);
    const auto tenth = gbs::nelder_mead(f, std::array{0.9}, {-2.}, {2.}, 1e-10, 1000, 0.1);
    ASSERT_EQ(tenth.x[0], near.x[0]); // the default
    const auto far = gbs::nelder_mead(f, std::array{0.9}, {-2.}, {2.}, 1e-10, 1000, 0.5);
    ASSERT_NEAR(far.x[0], -1.03558, 1e-5);
    ASSERT_THROW(gbs::nelder_mead(f, std::array{0.9}, {-2.}, {2.}, 1e-10, 1000, 0.), std::invalid_argument);
    ASSERT_THROW(gbs::nelder_mead(f, std::array{0.9}, {-2.}, {2.}, 1e-10, 1000, std::nan("")), std::invalid_argument);
}

TEST(tests_minimize, max_evaluations)
{
    std::size_t calls{};
    auto f = [&](const std::array<double, 2> &x) {
        ++calls;
        return 100. * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1. - x[0], 2);
    };
    auto r = gbs::nelder_mead(f, std::array{-1.2, 1.}, {-5., -5.}, {5., 5.}, 1e-12, 30);
    ASSERT_FALSE(r.converged);
    ASSERT_EQ(r.evaluations, calls);
    ASSERT_LE(calls, 30u);
    ASSERT_LT(r.f, f(std::array{-1.2, 1.}));
}

TEST(tests_minimize, float)
{
    auto f = [](const std::array<float, 2> &x) { return (x[0] - 0.5f) * (x[0] - 0.5f) + (x[1] + 0.25f) * (x[1] + 0.25f); };
    auto r = gbs::nelder_mead(f, std::array{0.f, 0.f}, {-1.f, -1.f}, {1.f, 1.f}, 1e-5f);
    ASSERT_NEAR(r.x[0], 0.5f, 1e-4f);
    ASSERT_NEAR(r.x[1], -0.25f, 1e-4f);
}

TEST(tests_minimize, unwrap)
{
    // the review case: vertices at z = 0.1 and 0.1 + 2 pi give the same x, their centroid
    // 0.1 + pi is the mirror point of the box; the trial point is moved to the period of the
    // reference, x unchanged. pi - z (beyond the bound) is how the simplex straddles it: kept
    const gbs::minimize_detail::Box<double, 2> box{{0., -1.}, {1., std::numeric_limits<double>::infinity()}};
    constexpr double pi = std::numbers::pi;
    const auto z = box.unwrap({0.1 + 2. * pi, -0.3}, {0.1, 0.3});
    ASSERT_NEAR(z[0], 0.1, 1e-15);
    ASSERT_EQ(z[1], -0.3); // one bound: -z is left alone
    ASSERT_NEAR(box.x(z)[0], box.x({0.1 + 2. * pi, -0.3})[0], 1e-15);
    ASSERT_EQ(box.unwrap({pi - 1.4, 0.}, {1.4, 0.})[0], pi - 1.4);
}

TEST(tests_minimize, callables)
{
    // a tolerance of the type of a literal: T is deduced from the arrays only
    auto quad = [](const std::array<float, 1> &x) { return (x[0] - 0.25f) * (x[0] - 0.25f); };
    ASSERT_NEAR(gbs::nelder_mead(quad, std::array{0.f}, {-1.f}, {1.f}, 1e-5).x[0], 0.25f, 1e-4f);
    // a mutable lambda and a stateful functor, called non-const
    std::size_t calls{};
    auto counted = [&calls, last = 0.](const std::array<double, 1> &x) mutable {
        ++calls;
        last = x[0];
        return std::abs(last - 0.5);
    };
    const auto r = gbs::nelder_mead(counted, std::array{0.1}, {0.}, {1.}, 1e-10);
    ASSERT_EQ(r.evaluations, calls);
    ASSERT_NEAR(r.x[0], 0.5, 1e-9);
    struct Counter
    {
        std::size_t n{};
        auto operator()(const std::array<double, 1> &x) -> double { return ++n, x[0] * x[0]; }
    } counter;
    ASSERT_EQ(gbs::nelder_mead(counter, std::array{0.5}, {-1.}, {1.}, 1e-8).evaluations, counter.n);
}

// a callable whose value does not convert to T is refused at the call
template <typename F>
concept minimizable = requires(F f) { gbs::nelder_mead(f, std::array{0.5}, std::array{0.}, std::array{1.}, 1e-6); };
static_assert(minimizable<double (*)(const std::array<double, 1> &)>);
static_assert(!minimizable<std::string (*)(const std::array<double, 1> &)>);

TEST(tests_minimize, invalid_arguments)
{
    auto f = [](const std::array<double, 1> &x) { return x[0] * x[0]; };
    ASSERT_THROW(gbs::nelder_mead(f, std::array{0.5}, {0.}, {1.}, 0.), std::invalid_argument);
    ASSERT_THROW(gbs::nelder_mead(f, std::array{0.5}, {0.}, {1.}, 1e-6, 0), std::invalid_argument);
    ASSERT_THROW(gbs::nelder_mead(f, std::array{0.5}, {1.}, {0.}, 1e-6), std::invalid_argument);
    ASSERT_THROW(gbs::nelder_mead(f, std::array{std::nan("")}, {0.}, {1.}, 1e-6), std::invalid_argument);
}

TEST(tests_minimize, agrees_with_extrema)
{
    // the distance between two curves, without derivatives, against the Newton projection
    gbs::BSCurve<double, 2> c1{{{0., 0.}, {1., 1.}, {2., 0.}}, {0., 0., 0., 1., 1., 1.}, 2};
    gbs::BSCurve<double, 2> c2{{{0., 2.}, {1., 0.8}, {2., 2.}}, {0., 0., 0., 1., 1., 1.}, 2};
    auto d = [&](const std::array<double, 2> &u) { return gbs::distance(c1(u[0]), c2(u[1])); };
    auto r = gbs::nelder_mead(d, std::array{0.2, 0.8}, {0., 0.}, {1., 1.}, 1e-10);
    auto [u1, u2, dx] = gbs::extrema_curve_curve(c1, c2, 1e-10);
    ASSERT_NEAR(r.f, dx, 1e-12);
    ASSERT_NEAR(r.x[0], u1, 1e-6);
    ASSERT_NEAR(r.x[1], u2, 1e-6);
}
