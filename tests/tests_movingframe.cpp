#include <doctest_gtest.hpp>
#include <gbs/movingframe.h>
#include <gbs/bscinterp.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    /// x = (a cos t, a sin t, b t)
    class Helix : public Curve<T, 3>
    {
        T a_, b_, t_;

    public:
        Helix(T a, T b, T turns) : a_{a}, b_{b}, t_{2 * pi * turns} {}
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            const T c = a_ * std::cos(t + T(d) * pi / 2), s = a_ * std::sin(t + T(d) * pi / 2);
            return {c, s, d == 0 ? b_ * t : (d == 1 ? b_ : 0.)};
        }
        auto bounds() const -> std::array<T, 2> override { return {0., t_}; }
    };

    void expect_orthonormal(const Frame3<T> &f)
    {
        ASSERT_NEAR(norm(f.x), 1., 1e-12);
        ASSERT_NEAR(norm(f.y), 1., 1e-12);
        ASSERT_NEAR(norm(f.z), 1., 1e-12);
        ASSERT_NEAR(f.x * f.z, 0., 1e-12);
        ASSERT_NEAR(norm(cross(f.x, f.y) - f.z), 0., 1e-12);
    }

    /// A closed knot-like curve, not planar: (cos t + 2 cos 2t, sin t - 2 sin 2t, sin 3t) (a trefoil).
    class Trefoil : public Curve<T, 3>
    {
    public:
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            auto c = [d](T w, T x) { return std::pow(w, T(d)) * std::cos(w * x + T(d) * pi / 2); };
            auto s = [d](T w, T x) { return std::pow(w, T(d)) * std::sin(w * x + T(d) * pi / 2); };
            return {c(1, t) + 2 * c(2, t), s(1, t) - 2 * s(2, t), s(3, t)};
        }
        auto bounds() const -> std::array<T, 2> override { return {0., 2 * pi}; }
    };
}

TEST(tests_movingframe, rotation_minimizing_on_a_helix)
{
    // helix: the rotation-minimizing normal turns against the Frenet one at the torsion rate,
    // theta(t) = -tau |x'| t, tau = b / (a² + b²)
    const T a = 2., b = 0.5;
    auto h = std::make_shared<Helix>(a, b, 2.);
    const T tol = 1e-9;
    MovingFrame<T> f(h, 0., 4 * pi, FrameLaw::RotationMinimizing, tol);
    MovingFrame<T> fr(h, 0., 4 * pi, FrameLaw::Frenet);
    const T tau = b / (a * a + b * b), speed = std::sqrt(a * a + b * b);
    const auto f0 = f(0.), n0 = fr(0.);
    const T theta0 = std::atan2(f0.x * n0.y, f0.x * n0.x);
    for (int i = 0; i <= 37; ++i)
    {
        const T t = 4 * pi * i / 37.;
        const auto g = f(t), n = fr(t);
        expect_orthonormal(g);
        expect_orthonormal(n);
        ASSERT_NEAR(norm(g.o - h->value(t)), 0., 1e-14);
        ASSERT_NEAR(norm(g.z - n.z), 0., 1e-12);
        const T theta = std::atan2(g.x * n.y, g.x * n.x);
        ASSERT_NEAR(std::remainder(theta - theta0 + tau * speed * t, 2 * pi), 0., 1e-7);
    }
    ASSERT_EQ(f.closure_angle(), 0.);
}

TEST(tests_movingframe, rotation_minimizing_on_a_planar_curve)
{
    // on a planar curve the normal of the plane stays fixed in the frame (no torsion)
    std::vector<std::array<T, 3>> p{{0., 0., 0.}, {1., 2., 0.}, {3., -1., 0.}, {5., 1., 0.}, {6., 3., 0.}};
    std::vector<T> u{0., 0.25, 0.5, 0.75, 1.};
    std::vector<constrType<T, 3, 1>> q;
    for (const auto &x : p)
        q.push_back({x});
    auto c = std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(q, u, std::size_t{3}));
    MovingFrame<T> f(c, 0., 1., FrameLaw::RotationMinimizing, 1e-10);
    const point<T, 3> z{0., 0., 1.};
    const auto f0 = f(0.);
    for (int i = 0; i <= 50; ++i)
    {
        const auto g = f(i / 50.);
        expect_orthonormal(g);
        ASSERT_NEAR(g.x * z, f0.x * z, 1e-9);
        ASSERT_NEAR(g.y * z, f0.y * z, 1e-9);
    }
    // and so is the fixed binormal along z, the same up to a constant rotation
    MovingFrame<T> fb(c, 0., 1., FrameLaw::FixedBinormal, 1e-10, {}, false, {0., 0., 3.});
    for (int i = 0; i <= 10; ++i)
        ASSERT_NEAR(norm(fb(i / 10.).y - z), 0., 1e-12);
}

TEST(tests_movingframe, closed_path)
{
    // around a trefoil the frame comes back rotated by its holonomy: spread, it closes
    auto c = std::make_shared<Trefoil>();
    MovingFrame<T> f(c, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-9, {}, true);
    ASSERT_GT(std::abs(f.closure_angle()), 1e-3);
    const auto a = f(0.), b = f(2 * pi);
    ASSERT_NEAR(norm(a.x - b.x), 0., 1e-8);
    ASSERT_NEAR(norm(a.y - b.y), 0., 1e-8);
    // continuous: no jump between close parameters
    for (int i = 0; i < 400; ++i)
    {
        const T t = 2 * pi * i / 400.;
        expect_orthonormal(f(t));
        ASSERT_LT(norm(f(t + 1e-6).x - f(t).x), 1e-4);
    }
    // the same path open does not close
    MovingFrame<T> o(c, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-9);
    ASSERT_GT(norm(o(0.).x - o(2 * pi).x), 1e-3);
}

TEST(tests_movingframe, breaks_and_accuracy)
{
    // the angle tolerance is reached: a coarse frame against a fine one
    auto h = std::make_shared<Helix>(1., 1., 1.);
    MovingFrame<T> coarse(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-4);
    MovingFrame<T> fine(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-12);
    for (int i = 0; i <= 20; ++i)
        ASSERT_LT(norm(coarse(2 * pi * i / 20.).x - fine(2 * pi * i / 20.).x), 1e-4);
    // breaks inside the range change nothing on a smooth path
    std::vector<T> br{1., 2.5};
    MovingFrame<T> b(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-12, br);
    for (int i = 0; i <= 20; ++i)
        ASSERT_LT(norm(b(2 * pi * i / 20.).x - fine(2 * pi * i / 20.).x), 1e-10);
}

TEST(tests_movingframe, refusals)
{
    auto line = std::make_shared<BSCurve<T, 3>>(std::vector<std::array<T, 3>>{{0., 0., 0.}, {1., 1., 0.}}, std::vector<T>{0., 0., 1., 1.}, 1);
    ASSERT_THROW(MovingFrame<T>(line, 0., 1., FrameLaw::Frenet), std::domain_error);
    ASSERT_THROW(MovingFrame<T>(line, 0., 1., FrameLaw::FixedBinormal, 1e-8, {}, false, {1., 1., 0.}), std::domain_error);
    ASSERT_THROW(MovingFrame<T>(line, 0., 1., FrameLaw::FixedBinormal), std::invalid_argument);
    ASSERT_THROW(MovingFrame<T>(line, 1., 1., FrameLaw::RotationMinimizing), std::invalid_argument);
    // along a straight line the rotation-minimizing frame does not turn
    MovingFrame<T> f(line, 0., 1., FrameLaw::RotationMinimizing);
    ASSERT_NEAR(norm(f(0.).x - f(1.).x), 0., 1e-15);
}

namespace
{
    /// Two quarters of circle of radius 1 joined G1 at u = 1, turning opposite ways (an S in z = 0).
    class SBend : public Curve<T, 3>
    {
    public:
        auto value(T u, size_t d = 0) const -> std::array<T, 3> override
        {
            const T h = pi / 2;
            if (u <= 1.)
            {
                // center (1, 0): (1 - cos a, sin a), a = h u
                const T a = h * u, f = std::pow(h, T(d));
                if (d == 0)
                    return {1 - std::cos(a), std::sin(a), 0.};
                return {-f * std::cos(a + T(d) * pi / 2), f * std::sin(a + T(d) * pi / 2), 0.};
            }
            // center (1, 2): (1 + sin b, 2 - cos b), b = h (u - 1)
            const T b = h * (u - 1), f = std::pow(h, T(d));
            if (d == 0)
                return {1 + std::sin(b), 2 - std::cos(b), 0.};
            return {f * std::sin(b + T(d) * pi / 2), -f * std::cos(b + T(d) * pi / 2), 0.};
        }
        auto bounds() const -> std::array<T, 2> override { return {0., 2.}; }
    };

    /// A curve that gives no second derivatives (as an edge of an offset surface).
    class FirstOnly : public Curve<T, 3>
    {
        std::shared_ptr<Curve<T, 3>> c_;

    public:
        explicit FirstOnly(std::shared_ptr<Curve<T, 3>> c) : c_{std::move(c)} {}
        auto value(T u, size_t d = 0) const -> std::array<T, 3> override
        {
            if (d > 1)
                throw std::runtime_error("FirstOnly: only first derivatives");
            return c_->value(u, d);
        }
        auto bounds() const -> std::array<T, 2> override { return c_->bounds(); }
        auto has_second_derivatives() const -> bool override { return false; }
    };
}

TEST(tests_movingframe, reversals_between_samples)
{
    // a planar S whose inflection is between the 32 samples: Frenet refused, rotation-minimizing not
    std::vector<constrType<T, 3, 1>> q;
    std::vector<T> u;
    for (int i = 0; i <= 40; ++i)
    {
        u.push_back(i / 40.);
        q.push_back({point<T, 3>{4 * u.back(), std::sin(5.2 * u.back() - 0.37), 0.}});
    }
    auto s = std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(q, u, std::size_t{3}));
    ASSERT_THROW((void)MovingFrame<T>(s, 0., 1., FrameLaw::Frenet), std::domain_error);
    ASSERT_NO_THROW((void)MovingFrame<T>(s, 0., 1., FrameLaw::RotationMinimizing));
    // the S of two arcs: Frenet reverses at the break
    auto sb = std::make_shared<SBend>();
    ASSERT_NEAR(norm(sb->value(1. - 1e-12) - sb->value(1.)), 0., 1e-11);
    const std::vector<T> br{1.};
    ASSERT_THROW((void)MovingFrame<T>(sb, 0., 2., FrameLaw::Frenet, 1e-8, br), std::domain_error);
    MovingFrame<T> rm(sb, 0., 2., FrameLaw::RotationMinimizing, 1e-9, br);
    const point<T, 3> z{0., 0., 1.};
    ASSERT_NEAR(rm(2.).x * z, rm(0.).x * z, 1e-9);
    // an arc of the plane xz whose tangent passes along z between the samples: a fixed binormal z refused
    auto arc = std::make_shared<EllipseCurve<T>>(2., 2., ax2<T, 3>{{{0., 0., 0.}, {0., -1., 0.}, {1., 0., 0.}}}, -0.31, 0.5);
    ASSERT_THROW((void)MovingFrame<T>(arc, -0.31, 0.5, FrameLaw::FixedBinormal, 1e-8, {}, false, {0., 0., 1.}), std::domain_error);
    // a helix of pitch angle 45 degrees, ten turns in one piece: the Frenet binormal turns by 112
    // degrees between two samples, but continuously (the intervals halved)
    auto h = std::make_shared<Helix>(1., 1., 10.);
    MovingFrame<T> f(h, 0., 20 * pi, FrameLaw::Frenet);
    for (int i = 0; i <= 100; ++i)
        expect_orthonormal(f(20 * pi * i / 100.));
}

TEST(tests_movingframe, derivatives_and_accuracy)
{
    // without second derivatives: Frenet refused, rotation-minimizing from the first ones
    auto c = std::make_shared<FirstOnly>(std::make_shared<Helix>(2., 0.5, 1.));
    ASSERT_THROW((void)MovingFrame<T>(c, 0., 2 * pi, FrameLaw::Frenet), std::domain_error);
    MovingFrame<T> f(c, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-9);
    MovingFrame<T> g(std::make_shared<Helix>(2., 0.5, 1.), 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-9);
    for (int i = 0; i <= 10; ++i)
        ASSERT_LT(norm(f(2 * pi * i / 10.).x - g(2 * pi * i / 10.).x), 1e-12);
    // an angle_tol that is not positive; one that cannot be reached: a runtime error, not a refusal
    auto h = std::make_shared<Helix>(1., 1., 1.);
    ASSERT_THROW((void)MovingFrame<T>(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 0.), std::invalid_argument);
    ASSERT_THROW((void)MovingFrame<T>(h, 0., 2 * pi, FrameLaw::RotationMinimizing, std::nan("")), std::invalid_argument);
    ASSERT_THROW((void)MovingFrame<T>(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-300), std::runtime_error);
    // a break given twice counts once
    const std::vector<T> twice{1., 1., 2.5};
    MovingFrame<T> b(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-12, twice);
    MovingFrame<T> once(h, 0., 2 * pi, FrameLaw::RotationMinimizing, 1e-12, std::vector<T>{1., 2.5});
    for (int i = 0; i <= 20; ++i)
        ASSERT_LT(norm(b(2 * pi * i / 20.).x - once(2 * pi * i / 20.).x), 1e-14);
}
