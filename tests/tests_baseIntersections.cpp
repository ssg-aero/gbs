#include <gtest/gtest.h>
#include <topology/baseIntersection.h>
#include <array>
#include <gbs/vecop.h>
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