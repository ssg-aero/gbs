
#include <iostream>
#include <vector>
#include <iomanip>
#include <nlopt.hpp>
#include <gtest/gtest.h>
#include <gbslib/bscurve.h>
#include <gbslib/vecop.h>
#include <gbslib/bssinterp.h>

#include <algorithm>
#include <functional>

using gbs::operator-;

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
//   if (!grad.empty()) {
//     grad[0] = 0.0;
//     grad[1] = 0.5 / sqrt(x[1]);
//   }
  return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
  double a = d->a, b = d->b;
//   if (!grad.empty()) {
//     grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
//     grad[1] = -1.0;
//   }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

TEST(tests_nlopt, example)
{

    //   nlopt::opt opt("LD_MMA", 2);
    nlopt::opt opt("LN_COBYLA", 2);
    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL;
    lb[1] = 0;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(myvfunc, NULL);
    my_constraint_data data[2] = {{2, 0}, {-1, 1}};
    opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
    opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234;
    x[1] = 5.678;
    double minf;

    try
    {
        opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
                  << std::setprecision(10) << minf << std::endl;
        // return EXIT_SUCCESS;
    }
    catch (std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
        // return EXIT_FAILURE;
    }
}

TEST(tests_nlopt, projC)
 {
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
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
    double u = 2.3;

    auto c = gbs::BSCurve(poles,k,p);


    std::array<double,3> pt = c.value(u);
    auto delta=0.0;
    pt[1]+=delta;

    class UserData
    {
        public:
        UserData(const gbs::BSCurve<double,3> &crv,const std::array<double,3> &p) : c{crv}, pt{p} {}
        gbs::BSCurve<double, 3> c;
        std::array<double, 3> pt;
    };
    UserData data(c,pt);

    auto f = [](const std::vector<double> &x, std::vector<double> &grad, void *user_data)
    {
        auto p_d = (UserData*)(user_data);
        auto c_u = p_d->c.value(x[0]);
        if (!grad.empty()) {
            auto dc_u = p_d->c.value(x[0],1);
               grad[0] = 0; 
            for(int i = 0 ; i < 3 ; i++)
            {
                grad[0] += 2*dc_u[i]*(c_u[i]-p_d->pt[i]);
            }
        }
        return gbs::sq_norm( c_u - p_d->pt );
    };

    nlopt::opt opt("LN_COBYLA", 1);
    // nlopt::opt opt("LD_MMA", 1);
    std::vector<double> lb(1),hb(1);
    lb[0] = k.front();
    hb[0] = k.back();

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(hb);
    opt.set_min_objective(f, &data);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(1);
    x[0] = 0.1234;
    double minf;

    try
    {
        opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] <<  ") = "
                  << std::setprecision(10) << minf << std::endl;
        ASSERT_NEAR(x[0],u,1e-4);
    }
    catch (std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
        ASSERT_TRUE(false);
    }
 }

 TEST(tests_nlopt, projS)
 {
    const std::vector<std::array<double,3> > points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };
    size_t p = 1;
    size_t q = 2;
    auto srf = gbs::interpolate(points,4,p,q,gbs::KnotsCalcMode::CHORD_LENGTH);
    auto u = 0.3;
    auto v = 0.7;


    std::array<double,3> pt = srf.value(u,v);
    auto delta=0.0;
    pt[1]+=delta;

    class ExtreamPS
    {
    public:
        ExtreamPS(const gbs::BSSurface<double, 3> &srf, const std::array<double, 3> &p) : s{srf}, pt{p} {}
        gbs::BSSurface<double, 3> s;
        std::array<double, 3> pt;
    };

    ExtreamPS data(srf, pt);


    auto f = [](const std::vector<double> &x, std::vector<double> &grad, void *user_data)
    {
        auto p_d = (ExtreamPS*)(user_data);
        // if (!grad.empty()) {
        //     auto dc_u = p_d->s.value(x[0],x[1],1);
        //        grad[0] = 0; 
        //     for(int i = 0 ; i < 3 ; i++)
        //     {
        //         grad[0] += 2*dc_u[i]*(c_u[i]-p_d->pt[i]);
        //     }
        // }
        return gbs::norm( p_d->s.value(x[0],x[1]) - p_d->pt );
    };

    nlopt::opt opt("LN_COBYLA", 2);
    std::vector<double> lb(2),hb(2);
    lb[0] = 0.;
    lb[1] = 0.;
    hb[0] = 1.;
    hb[1] = 1.;

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(hb);
    opt.set_min_objective(f, &data);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 0.1234;
    x[1] = 0.1234;

    double minf;

    try
    {
        opt.optimize(x, minf);
        std::cout << "found minimum at f(" << x[0] << " ; "<< x[1] <<  ") = "
                  << std::setprecision(10) << minf << std::endl;
        ASSERT_NEAR(x[0],u,1e-4);
        ASSERT_NEAR(x[1],v,1e-4);
    }
    catch (std::exception &e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
        ASSERT_TRUE(false);
    }    

 }