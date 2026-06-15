#include <chrono>
#include <doctest_gtest.hpp>
#include <topology/baseIntersection.h>
#include <topology/robustPredicates.h>
#include <array>
#include <cmath>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include <gbs/vecop.ixx>
#endif
#include <random> 

using namespace gbs;

TEST(base_intersections, in_circle)
{
    {

        std::array<double,2> a{0.0, 0.0};
        std::array<double,2> b{0.1, 0.0};
        std::array<double,2> c{0.5, 0.5};


        ASSERT_TRUE( are_ccw(a, b, c) );
        ASSERT_FALSE(are_ccw(b, a, c) );

        ASSERT_TRUE( orient_2d(a, b, c) > 0. ); // 
        ASSERT_TRUE( orient_2d(b, a, c) < 0. );

        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, a ) , 0.);
        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, b ) , 0.);
        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, c ) , 0.);

        auto O = ( a + b + c ) / 3.;
        ASSERT_TRUE( in_circle<double>(a, b, c, O ) > 0.);

        {
            std::array<double,2> d{0.5,-0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0.);
        }

        {
            std::array<double,2> d{0.0,0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) > 0.);
        }

        {
            std::array<double,2> d{1.0,0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0.);
        }

        {
            std::array<double,2> d{0.5,0.3};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0.);
        }

    }

    {
        size_t n = 500;
        std::random_device rd;  
        std::mt19937 gen(rd()); 
        std::uniform_real_distribution<double> distrib(-2, 2);
        const double sqrt2 = std::sqrt(2);
        for (size_t i{}; i < n; i++)
        {
            std::array<double, 2> d{distrib(gen), distrib(gen)};

            auto test = in_circle<double>(
                {-1, -1},
                {1, -1},
                {0., sqrt2},
                d);
            auto sq_r = d[0] * d[0] + d[1] * d[1];
            ASSERT_TRUE(
                test >= 0. ? sq_r <= 2. : sq_r > 2.
            );
        }
    }
}

TEST(base_intersections, orient2d_robust_exact)
{
    using gbs::orient2d;
    using gbs::orient_2d;

    // Basic sign agreement on a well-separated, non-degenerate triangle.
    {
        std::array<double, 2> a{0.0, 0.0}, b{1.0, 0.0}, c{0.0, 1.0};
        ASSERT_GT(orient2d(a, b, c), 0.0); // CCW
        ASSERT_LT(orient2d(a, c, b), 0.0); // CW
        ASSERT_EQ(orient2d(a, b, a), 0.0); // degenerate -> exactly zero
    }

    // Hand-crafted near-degenerate configuration: the exact orientation is +1
    // (origin O w.r.t. the directed segment A->B) but the products A.x*B.y and
    // A.y*B.x both round to the SAME double, so the naive determinant cancels to
    // 0 (or, with FMA contraction, to some other value). The robust predicate
    // must still report the exact sign. This is exactly the failure mode that
    // made the cavity construction platform-dependent (cf. issue #48).
    {
        const double p = std::ldexp(1.0, 27); // 2^27, exactly representable
        std::array<double, 2> A{p + 1.0, p}, B{p + 2.0, p + 1.0}, O{0.0, 0.0};
        // exact: A.x*B.y - A.y*B.x = (2^27+1)(2^27+1) - (2^27)(2^27+2) = +1
        ASSERT_GT(orient2d(A, B, O), 0.0);
        ASSERT_LT(orient2d(B, A, O), 0.0);
    }

#ifdef __SIZEOF_INT128__
    // Torture test against an exact 128-bit integer reference, on near-parallel
    // configurations (heavy cancellation). The robust predicate's sign must
    // match the exact sign on every single case.
    {
        std::mt19937_64 gen(0xB0BCA7ULL);
        // k near 2^27 so the cross products k*(k+.) are near 2^54, where a double
        // can no longer represent consecutive integers: the exact determinant is
        // +/-1 but lives entirely in the rounding noise of the products.
        std::uniform_int_distribution<long long> kbig(1LL << 26, (1LL << 27) - 3);

        const size_t n = 100000;
        size_t naive_wrong = 0;
        for (size_t i = 0; i < n; ++i)
        {
            long long k = kbig(gen);
            long long ax = k + 1, ay = k, bx = k + 2, by = k + 1; // exact det = +1
            int sref = 1;
            if (gen() & 1ULL) // swap a<->b => exact det = -1
            {
                std::swap(ax, bx);
                std::swap(ay, by);
                sref = -1;
            }

            std::array<double, 2> A{(double)ax, (double)ay}, B{(double)bx, (double)by}, O{0.0, 0.0};
            double r = orient2d(A, B, O);
            int srob = (r > 0) - (r < 0);
            ASSERT_EQ(srob, sref);

            double nv = orient_2d(A, B, O);
            if (((nv > 0) - (nv < 0)) != sref)
                ++naive_wrong;
        }
        std::cout << "orient2d torture: naive determinant disagreed with the exact"
                  << " sign on " << naive_wrong << " / " << n
                  << " near-degenerate cases (robust predicate: 0 disagreements)\n";
    }
#endif
}

TEST(base_intersections, seg_seg_strict_intersection)
{
    std::array<double, 2> a{-1., 2}, b{1.5, -0.5}, c{-1, -1};
    {
        ASSERT_FALSE(seg_seg_strict_intersection(a, b, c, a));
        ASSERT_FALSE(seg_seg_strict_intersection(a, b, c, b));
        ASSERT_FALSE(seg_seg_strict_intersection(a, b, c, c));
    }
    {
        std::array<double, 2> d{1., 1.};
        ASSERT_TRUE(seg_seg_strict_intersection(a, b, c, d));
    }
    {
        std::array<double, 2> d{0., 1.};
        ASSERT_FALSE(seg_seg_strict_intersection(a, b, c, d));
    }
    {
        std::array<double, 2> d{2.,-1.};
        ASSERT_FALSE(seg_seg_strict_intersection(a, b, c, d));
    }
}

TEST(base_intersections, in_triangle)
{
    std::array<double, 2> a{1,1}, b{4,-2}, c{3,3};
    { // On corner
        ASSERT_TRUE(in_triangle(a, b, c, a));
        ASSERT_TRUE(in_triangle(a, b, c, b));
        ASSERT_TRUE(in_triangle(a, b, c, c));
    }
    {// inside
        std::array<double, 2> d{2., 1.};
        ASSERT_TRUE(in_triangle(a, b, c, d));
    }

    {// On side
        std::array<double, 2> d{3., -1.}; 
        ASSERT_TRUE(in_triangle(a, b, c, d)); 
    }

    {// Outside
        std::array<double, 2> d{1., -2.}; 
        ASSERT_FALSE(in_triangle(a, b, c, d));
        d = {5,-4};
        ASSERT_FALSE(in_triangle(a, b, c, d));
        d = {5,-1};
        ASSERT_FALSE(in_triangle(a, b, c, d));
        d = {4,2};
        ASSERT_FALSE(in_triangle(a, b, c, d));
        d = {1,3};
        ASSERT_FALSE(in_triangle(a, b, c, d));
    }
}

TEST(base_intersections, in_triangle_barycentric)
{
    std::array<double, 2> a{1,1}, b{4,-2}, c{3,3};
    { // On corner
        ASSERT_TRUE(in_triangle_barycentric(a, b, c, a));
        ASSERT_TRUE(in_triangle_barycentric(a, b, c, b));
        ASSERT_TRUE(in_triangle_barycentric(a, b, c, c));
    }
    {// inside
        std::array<double, 2> d{2., 1.};
        ASSERT_TRUE(in_triangle_barycentric(a, b, c, d));
    }

    {// On side
        std::array<double, 2> d{3., -1.}; 
        ASSERT_TRUE(in_triangle_barycentric(a, b, c, d)); 
    }

    {// Outside
        std::array<double, 2> d{1., -2.}; 
        ASSERT_FALSE(in_triangle_barycentric(a, b, c, d));
        d = {5,-4};
        ASSERT_FALSE(in_triangle_barycentric(a, b, c, d));
        d = {5,-1};
        ASSERT_FALSE(in_triangle_barycentric(a, b, c, d));
        d = {4,2};
        ASSERT_FALSE(in_triangle_barycentric(a, b, c, d));
        d = {1,3};
        ASSERT_FALSE(in_triangle_barycentric(a, b, c, d));
    }

    {
        auto generate_random_points = [](size_t num_points)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> dis(-10.0, 10.0);

            std::vector<std::array<double, 2>> points(num_points);
            for (auto &point : points)
            {
                point[0] = dis(gen);
                point[1] = dis(gen);
            }
            return points;
        };

        size_t num_points = 100000;
        size_t num_runs = 100;

        std::array<double, 2> a{1.0, 1.0};
        std::array<double, 2> b{3.0, 1.0};
        std::array<double, 2> c{2.0, 3.0};

        auto points = generate_random_points(num_points);

        auto start_time = std::chrono::high_resolution_clock::now();
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time;

        double in_triangle_time = 0;
        double in_triangle_barycentric_time = 0;

        for (size_t run = 0; run < num_runs; ++run)
        {
            start_time = std::chrono::high_resolution_clock::now();
            for (const auto &point : points)
            {
                auto test = in_triangle<double>(a, b, c, point);
            }
            end_time = std::chrono::high_resolution_clock::now();
            elapsed_time = end_time - start_time;
            in_triangle_time += elapsed_time.count();

            start_time = std::chrono::high_resolution_clock::now();
            for (const auto &point : points)
            {
                auto test = in_triangle_barycentric<double>(a, b, c, point);
            }
            end_time = std::chrono::high_resolution_clock::now();
            elapsed_time = end_time - start_time;
            in_triangle_barycentric_time += elapsed_time.count();
        }

        std::cout << "in_triangle average time: "
                  << in_triangle_time / num_runs
                  << " seconds" << std::endl;

        std::cout << "in_triangle_barycentric average time: "
                  << in_triangle_barycentric_time / num_runs
                  << " seconds" << std::endl;
    }
}