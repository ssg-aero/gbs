#include <doctest_gtest.hpp>
#include <gbs/curves>
#include <gbs-render/vtkGbsRender.h>
#include <gbs/bscbuild.h>
#include <gbs/bscshaping.h>
using namespace gbs;
using T = double;
const bool PLOT_ON = false;

TEST(curve_shaping, move_to_point)
{
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    T u = 2.3;

    auto crv = BSCurve<T,3>(poles,k,p);
    point<T,3> deltap{0.1,0.2,0.3}; 
    std::vector<std::array<T,3> > points{crv(u),crv(u)+deltap};

    auto new_crv = moved_to_point(crv, points[1], u);
    ASSERT_NEAR( distance(points[1], new_crv(u)), 0., 1e-9 );

    if(PLOT_ON)
    {
        auto seg = build_segment(points[0], points[1]);
        plot(
            gbs::crv_dsp<T, 3, false>{
                .c = &crv,
                .col_crv = {1., 0., 0.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },
            gbs::crv_dsp<T, 3, false>{
                .c = &new_crv,
                .col_crv = {0., 0., 1.},
                .poles_on = true,
                .col_poles = {0., 1., 0.},
                .col_ctrl = {0., 0., 0.},
                .show_curvature = false,
            },

            gbs::crv_dsp<T, 3, false>{
                .c = &seg,
                .col_crv = {0., 1., 0.},
            },
            points);
    }
}

TEST(curve_shaping, singular_constraints_throw)
{
    // two constraints the same make the constraint system singular: it throws instead of
    // moving the poles by meaningless values
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T, 3>> poles{{0., 0., 0.}, {0., 1., 0.}, {1., 1., 0.}, {1., 1., 1.}, {1., 1., 2.}, {3., 1., 1.}, {0., 4., 1.}};
    const auto crv = BSCurve<T, 3>(poles, k, 2);
    const point<T, 3> d{0.1, 0.2, 0.3};
    // fewer constraints than poles moved (minimal norm), as many (square)
    ASSERT_THROW(moved_to_constraints_delta(crv, std::vector<bsc_constraint<T, 3>>{{2.3, d, 0}, {2.3, d, 0}}, 1, 5), std::invalid_argument);
    ASSERT_THROW(moved_to_constraints_delta(crv, std::vector<bsc_constraint<T, 3>>{{2.3, d, 0}, {2.3, d, 0}, {3.1, d, 0}}, 2, 4), std::invalid_argument);
}

// Degree 30, above bspline_stack_max_degree: the constraint matrix is filled by
// basis_ders (A2.2/A2.3 on the heap), where the recursive basis_function took
// seconds per entry (exponential in the degree). The curve meets its
// constraints: a point and a first derivative at two parameters.
TEST(curve_shaping, move_to_constraints_high_degree)
{
    const size_t p = 30;
    std::vector<T> k(p + 1, 0.);
    for (T x : {0.5, 0.9, 1.4})
        k.push_back(x);
    k.insert(k.end(), p + 1, 2.);
    const size_t np = k.size() - p - 1;
    points_vector<T, 3> poles(np);
    for (size_t i = 0; i < np; ++i)
        poles[i] = {T(i) / T(np - 1), std::sin(0.4 * T(i)), 0.};
    const BSCurve<T, 3> crv(poles, k, p);

    const point<T, 3> D{0., 0., 0.3};
    const std::vector<bsc_constraint<T, 3>> constraints{
        {0.7, crv(0.7) + D, 0},
        {0.7, crv(0.7, 1), 1},
        {1.6, crv(1.6) - D, 0},
    };
    const auto shaped = moved_to_constraints(crv, constraints);
    for (const auto &[u, x, d] : constraints)
    {
        CAPTURE(u);
        CAPTURE(d);
        ASSERT_LT(norm(shaped(u, d) - x), 1e-8 * std::max(1., norm(x)));
    }
}
