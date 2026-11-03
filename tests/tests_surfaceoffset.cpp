#include <doctest_gtest.hpp>
#include <gbs/surfaceoffset.h>
#include <gbs/bssurf.h>
#include <gbs/primitives.h>

#include <numbers>

using namespace gbs;

TEST(tests_surf_offset, ctor)
{
    using T = double;
    const size_t dim = 3;

    // auto p_srf = std::make_share<BSSurface<T,dim>>(
    //     points_vector<T,dim>
    // )
}

namespace
{
    // a sphere of radius 2 about (1, -1, 0.5), its axis along (1, 1, 1)
    auto sphere_() -> SphericalSurface<double>
    {
        const Frame3<double> f{ax2<double, 3>{{{1., -1., 0.5}, {1., 1., 1.}, {1., -1., 0.}}}};
        return {2., f, {0., 2. * std::numbers::pi, -std::numbers::pi / 2., std::numbers::pi / 2.}};
    }
}

// The offset of a sphere is a sphere, its poles included.
TEST(tests_surf_offset, surface_offset_of_a_sphere)
{
    const auto s = std::make_shared<gbs::SphericalSurface<double>>(sphere_());
    auto f = [](double, double) { return 0.5; };
    const gbs::SurfaceOffset<double, 3, decltype(f)> offset{s, f};
    ASSERT_EQ(offset.bounds(), s->bounds());
    const gbs::point<double, 3> c{1., -1., 0.5};
    const double pi = std::numbers::pi;
    const double side = gbs::normal_direction(*s, 1., 0.) * (s->value(1., 0.) - c) > 0. ? 1. : -1.;
    for (double u : {0., 1., 2. * pi})
        for (double v : {-pi / 2., 0., 1., pi / 2.})
            ASSERT_NEAR(gbs::norm(offset.value(u, v) - c), 2. + side * 0.5, 1e-9);
    ASSERT_THROW(offset.value(0., 0., 1, 0), std::runtime_error);
}
