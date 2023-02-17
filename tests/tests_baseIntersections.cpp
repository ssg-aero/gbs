#include <gtest/gtest.h>
#include <topology/baseIntersection.h>

using namespace gbs

TEST(tests_mesh_delauney2d, in_circle)
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