#include <gtest/gtest.h>
#include <gbs/solvers.h>
#include <gbs/gbslib.h>
#include <gbs/vecop.h>
#include <numbers>
TEST(tests_solvers, solve_N_nlop)
{
    using std::cos;
    using std::sin;
    using std::numbers::sqrt2;
    using std::numbers::pi;
    using gbs::operator-;

    auto r = 0.2;
    gbs::point<double,2> Oc{0.,1.};
    gbs::point<double,2> Od{0.,0.};
    gbs::point<double,2> D{1.,0.};

    auto f = [=,&Oc,&Od,&D](const std::vector<double> &x)
    {
        auto th = x[0];
        auto  u = x[1];
        auto Pc = gbs::point<double,2>{r*cos(th)+Oc[0],r*sin(th)+Oc[1]};
        auto Pd = gbs::point<double,2>{Od[0]+u*D[0],Od[1]+u*D[1]};
        auto d = gbs::norm(Pc-Pd);
        return d;
    };

    std::vector<double> x{0.,0.};
    std::vector<double> lb{0.,-1.};
    std::vector<double> hb{2*pi,1.};
    auto tol = 1.e-8;
    auto minf = gbs::solve_N_nlop(
        f,
        x,lb,hb,
        tol,
        // "LN_COBYLA"
        "GN_ORIG_DIRECT"
    );

    ASSERT_NEAR(minf,0.8,2.*tol);
    ASSERT_NEAR(x[0],3.*pi/2.,2.*tol);
    ASSERT_NEAR(x[1],0.,2.*tol);

}