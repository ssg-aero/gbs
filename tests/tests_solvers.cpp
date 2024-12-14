#include <gtest/gtest.h>
#include <gbs/solvers.h>
#include <gbs/gbslib.h>
#ifdef GBS_USE_MODULES
    import vecop;
#endif
#include <numbers>
using std::cos;
using std::sin;
using std::exp;
using std::numbers::sqrt2;
using std::numbers::pi;
using gbs::operator-;

TEST(tests_solvers, solve_N_nlop)
{

        auto r = 0.2;
        gbs::point<double, 2> Oc{0., 1.};
        gbs::point<double, 2> Od{0., 0.};
        gbs::point<double, 2> D{1., 0.};

        auto f = [=, &Oc, &Od, &D](const std::vector<double> &x) {
            auto th = x[0];
            auto u = x[1];
            auto Pc = gbs::point<double, 2>{r * cos(th) + Oc[0], r * sin(th) + Oc[1]};
            auto Pd = gbs::point<double, 2>{Od[0] + u * D[0], Od[1] + u * D[1]};
            auto d = gbs::norm(Pc - Pd);
            return d;
        };

        std::vector<double> x{0., 0.};
        std::vector<double> lb{0., -1.};
        std::vector<double> hb{2 * pi, 1.};
        auto tol = 1.e-8;
        auto minf = gbs::solve_N_nlop(
            f,
            x, lb, hb,
            tol,
            // "LN_COBYLA"
            // "GN_ORIG_DIRECT"
            nlopt::GN_ORIG_DIRECT
            );

        ASSERT_NEAR(minf, 0.8, 2. * tol);
        ASSERT_NEAR(x[0], 3. * pi / 2., 2. * tol);
        ASSERT_NEAR(x[1], 0., 2. * tol);

    // solve with gradian
    

}

TEST(tests_solvers, solve_N_nlop_f)
{
    using T = float;
    T r = 0.2;
    gbs::point<T, 2> Oc{0., 1.};
    gbs::point<T, 2> Od{0., 0.};
    gbs::point<T, 2> D{1., 0.};

    auto f = [=, &Oc, &Od, &D](const std::vector<double> &x)
    {
        T th = x[0];
        T u = x[1];
        auto Pc = gbs::point<T, 2>{r * cos(th) + Oc[0], r * sin(th) + Oc[1]};
        auto Pd = gbs::point<T, 2>{Od[0] + u * D[0], Od[1] + u * D[1]};
        auto d = gbs::norm(Pc - Pd);
        return d;
    };

    std::vector<T> x{0., 0.};
    std::vector<T> lb{0., -1.};
    std::vector<T> hb{2 * pi, 1.};
    T tol = 1.e-5; // float requires la lower tol
    auto minf = gbs::solve_N_nlop(
        f,
        x, lb, hb,
        tol,
        // "LN_COBYLA"
        // "GN_ORIG_DIRECT"
        nlopt::GN_ORIG_DIRECT);

    ASSERT_NEAR(minf, T(0.8), 2. * tol);
    ASSERT_NEAR(x[0], 3. * pi / 2., 25. * tol);
    ASSERT_NEAR(x[1], 0., 2. * tol);

    // solve with gradian
    

}

TEST(tests_solvers, solve_D_nlop)
    {
        auto a = 2.;
        auto b = 1.;
        auto c = 3.;
        auto d = 1.;
        auto Y1= 2.;
        auto Y2= 6.;
        auto e = [](double x){return exp(x);};
        auto f = [a,b,c,d,Y1,Y2,&e](const std::vector<double> &X)
        {
            auto x = X[0];
            auto y = X[1];
            auto r1=a*x*x+b*y-Y1;
            auto r2=c*x+d*e(y)-Y2;

            return std::vector<double>{r1,r2};
        };
        auto g = [a,b,c,d,&e](const std::vector<double> &X,const std::vector<double> &r) {
            // std::vector<double> grad(X.size());
            auto x = X[0];
            auto y = X[1];
            auto df0dx = 2.*a*x ;
            auto df1dx = c ;
            auto df0dy = b ;
            auto df1dy = d * e(y) ;
            // f outo put is euclidian norm
            auto g1 = 2.* ( r[0] * df0dx + r[1] * df1dx );
            auto g2 = 2.* ( r[0] * df0dy + r[1] * df1dy );
            
            return std::vector<double>{g1,g2};
        };

        std::vector<double> x{0., 0.};
        std::vector<double> lb{-10., -10.};
        std::vector<double> hb{10., 10.};
        auto tol = 1.e-8;
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol,
            // nlopt::LN_COBYLA
            // nlopt::LD_TNEWTON_PRECOND
            nlopt::LD_CCSAQ
            );
        ASSERT_NEAR(minf,0.,tol);
    }

    TEST(tests_solvers, solve_D_nlop_f)
    {
        float a = 2.;
        float b = 1.;
        float c = 3.;
        float d = 1.;
        float Y1= 2.;
        float Y2= 6.;
        auto e = [](float x){return exp(x);};
        auto f = [a,b,c,d,Y1,Y2,&e](const std::vector<double> &X)
        {
            float x = X[0];
            float y = X[1];
            float r1=a*x*x+b*y-Y1;
            float r2=c*x+d*e(y)-Y2;

            return std::vector<double>{r1,r2};
        };
        auto g = [a,b,c,d,&e](const std::vector<double> &X,const std::vector<double> &r) {
            // std::vector<double> grad(X.size());
            float x = X[0];
            float y = X[1];
            float df0dx = 2.*a*x ;
            float df1dx = c ;
            float df0dy = b ;
            float df1dy = d * e(y) ;
            // f outo put is euclidian norm
            float g1 = 2.* ( r[0] * df0dx + r[1] * df1dx );
            float g2 = 2.* ( r[0] * df0dy + r[1] * df1dy );
            
            return std::vector<double>{g1,g2};
        };

        std::vector<float> x{0., 0.};
        std::vector<float> lb{-10., -10.};
        std::vector<float> hb{10., 10.};
        float tol = 1.e-5;
        float minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol,
            // nlopt::LN_COBYLA
            // nlopt::LD_TNEWTON_PRECOND
            nlopt::LD_CCSAQ
            );
        ASSERT_NEAR(minf,0.,tol);
    }