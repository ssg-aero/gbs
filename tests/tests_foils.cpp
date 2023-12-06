#include <gtest/gtest.h>
#include <gbs/bscbuild.h>

import vecop;
#include <gbs/bscinterp.h>
#include <gbs/bscapprox.h>
#include <gbs/bsctools.h>
#include <gbs/bssbuild.h>

#include <gbs-render/vtkGbsRender.h>

#include <Eigen/Dense>

const double PI = acos(-1.);

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

TEST(tests_foils, type1)
{
    auto r_ba = 5.;
    auto e_ba = 5.;

    auto r_bf = 0.25;
    auto e_bf = 1.;

    auto c = 25.;

    auto arc_ba = gbs::build_ellipse<double,3>(r_ba ,r_ba/ e_ba, {r_ba, 0., 0.});
    arc_ba.trim(2.5/8.,1./2.);
    arc_ba.reverse();

    auto arc_bf = gbs::build_ellipse<double,3>(r_bf ,r_bf/ e_bf, {c - r_bf, 0., 0.});
    arc_bf.trim(0.,0.5/8.);
    arc_bf.reverse();

    size_t p = 3;
    std::vector<double> u = {0.,1.};
    // std::vector<double> k = {0., 0., 0.,1./4.,1/2.,3./4, 1. ,1. ,1.};
    // std::vector<double> k = {0., 0., 0.,0.1,1/2.,0.9, 1. ,1. ,1.};
    std::vector<double> k = {0., 0., 0.,0., 0.33, 0.66 ,1., 1. ,1. ,1.};
    // std::vector<double> k = gbs::build_simple_mult_flat_knots<double>(u, 8, 5);
    auto pt0 = arc_ba.end();
    auto pt3 = arc_bf.begin();
    auto tg0 = arc_ba.end(1);
    auto tg3 = arc_bf.begin(1);
    auto cu0 = arc_ba.end(2);
    auto cu3 = arc_bf.begin(2);
    auto to0 = arc_ba.end(3);
    auto to3 = arc_bf.begin(3);

    std::vector<gbs::constrType<double, 3, 3>> Q =
        {
            {{pt0, tg0, cu0/*, to0*/}},
            {{pt3, tg3, cu3/*, to3*/}},
        };

    auto poles = gbs::build_poles(Q,k,u,p);


    auto arc = gbs::BSCurve<double,3>(poles,k,p);

}

auto thicken_foil = [](const auto &u, const auto camber_line, auto f_thickness, auto chord, bool side1, size_t n_poles , size_t p,auto x1 = 0.01,auto x2=0.99) {
    gbs::points_vector_2d_d points_side1(u.size());
    std::transform(u.begin(), u.end(), points_side1.begin(),
                   [&](const auto u_) {
                       auto t_ = camber_line.value(u_, 1);
                       std::swap(t_[0], t_[1]);
                       t_[0] = -t_[0];
                       t_ = t_ * (1 / gbs::norm(t_));
                       auto l_ = gbs::length(camber_line, 0., u_) / chord;
                       auto x_ = x1 * (1. - l_) + x2 * (l_);
                       auto ep = f_thickness(x_);
                       auto p_ = camber_line.value(u_);
                       auto side = side1 ? 1. : -1.;
                       return p_ + side * t_ * ep;
                   });
    return gbs::approx(points_side1, p, n_poles, u, true);
};

TEST(tests_foils, type2)
{
    using bsc2d = gbs::BSCurve<double, 2>;
    auto b1 = PI / 24;
    auto b2 = -PI / 12;
    std::vector<double> k_cl = {0., 0., 0., 0.5, 1., 1., 1.};
    auto t1 = 0.3;
    auto t2 = 0.3;
    gbs::points_vector_2d_d poles_cl =
        {
            {0., 0.},
            {t1 * cos(b1), t1 * sin(b1)},
            {1. - t2 * cos(b1), t2 * sin(b1)},
            {1., 0.},
        };
    bsc2d camber_line{poles_cl, k_cl, 2};

    auto chord = gbs::length(camber_line,0.,1.);

    auto u = gbs::make_range(0.,1.,20);

    auto t = 0.1;
    auto f_thickness_naca =[&t](const auto x_){return 5 * 0.1 *(0.2969*sqrt(x_)-0.1260*x_-0.3516*x_*x_+0.2843*x_*x_*x_-0.1015*x_*x_*x_*x_);};

    auto n_pole_side1 = 10;
    auto n_pole_side2 = 12;
    auto deg_side = 5;
    auto side1 = thicken_foil(u,camber_line,f_thickness_naca,chord,true,n_pole_side1,deg_side,0.01,0.95);
    auto side2 = thicken_foil(u,camber_line,f_thickness_naca,chord,false,n_pole_side2,deg_side,0.01,0.95);

    side1.reverse();
    auto le = gbs::c3_connect(side1,side2);
    auto te = gbs::c3_connect(side2,side1);

    auto foil_2d = gbs::join(side1,le);
    foil_2d = gbs::join(foil_2d,side2);
    foil_2d = gbs::join(foil_2d,te);
}

TEST(tests_foils, type2_blade)
{
    using bsc2d = gbs::BSCurve<double, 2>;
    // Camberline definition
    auto b1 = PI / 24;
    auto b2 = -PI / 12;
    std::vector<double> k_cl = {0., 0., 0., 0.5, 1., 1., 1.};
    auto t1 = 0.3;
    auto t2 = 0.3;
    gbs::points_vector_2d_d poles_cl =
        {
            {0., 0.},
            {t1 * cos(b1), t1 * sin(b1)},
            {1. - t2 * cos(b1), t2 * sin(b1)},
            {1., 0.},
        };
    bsc2d camber_line{poles_cl, k_cl, 2};
    // Thickness law definition
    auto chord = gbs::length(camber_line,0.,1.);
    auto u = gbs::make_range(0.,1.,100);
    auto t = 0.1;
    auto f_thickness_naca =[&t](const auto x_){return 5 * 0.1 *(0.2969*sqrt(x_)-0.1260*x_-0.3516*x_*x_+0.2843*x_*x_*x_-0.1015*x_*x_*x_*x_);};
    // Convertion from thickness law to NURBS
    auto n_pole_side1 = 10;
    auto n_pole_side2 = 12;
    auto deg_side = 5;
    auto side1 = thicken_foil(u,camber_line,f_thickness_naca,chord,true,n_pole_side1,deg_side,0.01,0.95);
    auto side2 = thicken_foil(u,camber_line,f_thickness_naca,chord,false,n_pole_side2,deg_side,0.01,0.95);
    // Definition of leading and trailling edges
    side1.reverse();
    auto le = gbs::c3_connect(side1,side2);
    auto te = gbs::c3_connect(side2,side1);
    // full foil definition
    auto foil_2d = gbs::join(side1,le);
    foil_2d = gbs::join(foil_2d,side2);
    foil_2d = gbs::join(foil_2d,te);
    // Spatial positioning
    auto foil1 = gbs::add_dimension(foil_2d,0.3);
    auto poles2 = foil_2d.poles();
    std::transform(
        poles2.begin(),
        poles2.end(),
        poles2.begin(),
        [](const auto &p_){
            auto p = gbs::translated(p_,std::array<double,2>{-0.5,0.});
            gbs::rotate(p,PI/12.);
            gbs::translate(p,std::array<double,2>{0.5,0.});
            return p;}
    );
    auto foil2 = gbs::add_dimension(bsc2d(poles2,foil_2d.knotsFlats(),foil_2d.degree()),0.6);
    auto poles3 = foil_2d.poles();
    std::transform(
        poles3.begin(),
        poles3.end(),
        poles3.begin(),
        [](const auto &p_){
            auto p = gbs::translated(p_,std::array<double,2>{0.0,0.2});
            return p;}
    );
    auto foil3 = gbs::add_dimension(bsc2d(poles3,foil_2d.knotsFlats(),foil_2d.degree()),0.9);


    std::list<gbs::BSCurve<double,3>> bs_lst {foil1,foil2,foil3};
    if(PLOT_ON)
        gbs::plot(
            gbs::loft( bs_lst ),
            gbs::crv_dsp<double,3,false>{
                .c =&foil1
                },
            
            gbs::crv_dsp<double,3,false>{
                .c =&foil2
                },
            gbs::crv_dsp<double,3,false>{
                .c =&foil3
                });


}