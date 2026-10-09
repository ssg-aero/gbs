#include <doctest_gtest.hpp>
#include <gbs-io/print.h>
#include <gbs-io/fromjson.h>
#include <gbs-io/fromjson2.h>
#include <gbs-io/iges.h>
#include <gbs-render/vtkGbsRender.h>
#include <gbs-render/vtkgridrender.h>
#include <gbs-mesh/mshedge.h>
#include <gbs-mesh/tfi.h>
#include <numbers>
#include <filesystem>
#include <fstream>
#include <sstream>
#ifdef GBS_USE_MODULES
   import knots_functions;
   import vecop;
#endif
#include <gbs/bssapprox.h>
#include "tests_helpers.h"


using namespace gbs;
using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

TEST(tests_io, json_bsc)
{

   rapidjson::Document document;
   std::string dir = get_directory(__FILE__);
   parse_file((dir+"/in/crv_bs.json").c_str(), document);

   auto curve1d = bscurve_direct<double, 1>(document["1d_entities"][0]);
   ASSERT_EQ(curve1d.degree(), 3);
   auto curve2d = bscurve_direct<double, 2>(document["2d_entities"][0]);
   ASSERT_EQ(curve2d.poles().size(), 4);
   ASSERT_EQ(curve2d.knotsFlats().size(), 8);
   auto curvei2d = bscurve_interp_cn<double, 2>(document["2d_entities"][1]);
   ASSERT_NEAR(curvei2d.value(0.5)[0], 0.0, 1e-6);
   ASSERT_NEAR(curvei2d.value(0.5)[1], 0.5, 1e-6);
   auto curve3d = bscurve_direct<double, 3>(document["3d_entities"][0]);
   ASSERT_NEAR(curve3d.poles()[2][2], 1.0, 1e-6);
}

TEST(tests_io, bscurve_direct_mults)
{
   rapidjson::Document document;
   std::string dir = get_directory(__FILE__);
   parse_file((dir+"/in/bscurve_direct_mults.json").c_str(), document);
   auto crv1 = bscurve_direct_mults<double,2>(document["curves"][0]);
   ASSERT_EQ(crv1.degree(), 3);
   ASSERT_EQ(crv1.mults()[0], 4);
   ASSERT_EQ(crv1.mults()[1], 4);
   ASSERT_DOUBLE_EQ(crv1.knots()[0], 0.);
   ASSERT_DOUBLE_EQ(crv1.knots()[1], 1.);
   ASSERT_DOUBLE_EQ(crv1.poles()[0][0], 0.);
   ASSERT_DOUBLE_EQ(crv1.poles()[0][1], 0.);
   ASSERT_DOUBLE_EQ(crv1.poles()[1][0], 0.3);
   ASSERT_DOUBLE_EQ(crv1.poles()[1][1], 0.1);
}

TEST(tests_io, meridian_channel)
{
   std::array<double, 3> col_crv = {1., 0., 0.};
   float line_width = 2.f;
   auto f_dspc = [line_width, &col_crv](const auto &crv_l) {
      std::vector<gbs::crv_dsp<double, 2, false>> crv_dsp;
      for (auto &crv : crv_l)
      {
         crv_dsp.push_back(
             gbs::crv_dsp<double, 2, false>{
                 .c = &(crv),
                 .col_crv = col_crv,
                 .poles_on = false,
                 .line_width = line_width,
             });
      }
      return crv_dsp;
   };

   rapidjson::Document document;
   std::string dir = get_directory(__FILE__);
   // parse_file("../tests/in/test_channel_solve_cax.json",document);
   parse_file((dir+"/in/test_channel_solve_roue_ct.json").c_str(), document);

   std::list<gbs::BSCurve2d_d> hub_curves;
   std::list<gbs::BSCurve2d_d> shr_curves;
   std::list<gbs::BSCurve2d_d> ml_curves;
   for (auto &val : document["hub_curves"].GetArray())
   {
      hub_curves.push_back(gbs::make_bscurve<double, 2>(val));
   }

   for (auto &val : document["shr_curves"].GetArray())
   {
      shr_curves.push_back(gbs::make_bscurve<double, 2>(val));
   }

   for (auto &val : document["mean_lines"].GetArray())
   {
      ml_curves.push_back(gbs::make_bscurve<double, 2>(val));
   }

   auto hub_pnts = gbs::make_point_vec<double, 2>(document["hub_corner_points"]);
   auto shr_pnts = gbs::make_point_vec<double, 2>(document["shr_corner_points"]);

   std::list<gbs::BSCurve2d_d> crv_l;
   crv_l.insert(crv_l.end(), hub_curves.begin(), hub_curves.end());
   crv_l.insert(crv_l.end(), shr_curves.begin(), shr_curves.end());

   std::list<gbs::BSCurve2d_d> crv_m;
   for (auto i = 0; i < hub_pnts.size(); i++)
   {
      gbs::points_vector_2d_d pts_(2);
      pts_[0] = hub_pnts[i];
      pts_[1] = shr_pnts[i];
      crv_m.push_back(gbs::interpolate(pts_, 1, gbs::KnotsCalcMode::CHORD_LENGTH));
   }

   if(PLOT_ON)
   {
      auto crv_dsp = f_dspc(crv_l);
      col_crv = {0., 1., 0.};
      auto crv_dsp_ml = f_dspc(ml_curves);
      col_crv = {0., 0., 0.};
      line_width = 4.f;
      gbs::plot(
          crv_dsp,
          crv_dsp_ml,
          f_dspc(crv_m));
   }
}

auto build_channel_curves(std::vector<gbs::BSCurve2d_d> &crv_m, std::vector<gbs::BSCurve2d_d> &crv_l, std::vector<double> &u_m, std::vector<double> &u_l)
{
   rapidjson::Document document;
   std::string dir = get_directory(__FILE__);
   // parse_file("../tests/in/test_channel_solve_cax.json",document);
   parse_file((dir+"/in/test_channel_solve_roue_ax.json").c_str(), document);

   auto ml_crv = gbs::make_bscurve<double, 2>(document["mean_lines"].GetArray()[0]);
   ml_crv.changeBounds(0., 1.);

   ASSERT_TRUE(document["hub_curves"].GetArray()[0].HasMember("constraints"));
   ASSERT_TRUE(document["shr_curves"].GetArray()[0].HasMember("constraints"));
   auto hub_cstr = make_constraints_vec<double, 2, 2>(document["hub_curves"].GetArray()[0]["constraints"]);
   auto shr_cstr = make_constraints_vec<double, 2, 2>(document["shr_curves"].GetArray()[0]["constraints"]);

   u_l = gbs::knots_and_mults(ml_crv.knotsFlats()).first;
   auto hub_crv = gbs::interpolate<double, 2, 2>(hub_cstr, u_l);
   auto shr_crv = gbs::interpolate<double, 2, 2>(shr_cstr, u_l);

   crv_m = {hub_crv, ml_crv, shr_crv};

   auto hub_pnts = gbs::make_point_vec<double, 2>(document["hub_corner_points"]);
   auto shr_pnts = gbs::make_point_vec<double, 2>(document["shr_corner_points"]);
   auto ml_pnts = gbs::make_point_vec<double, 2>(document["mean_lines"].GetArray()[0]["constraints"].GetArray()[0]);

   u_m = {0., 0.5, 1.};
   for (auto i = 0; i < hub_pnts.size(); i++)
   {
      gbs::points_vector_2d_d pts_(3);
      pts_[0] = hub_pnts[i];
      pts_[1] = ml_pnts[i];
      pts_[2] = shr_pnts[i];
      crv_l.push_back(gbs::interpolate(pts_, u_m, 1));
   }
}

auto f_dspc(const std::vector<gbs::BSCurve2d_d> &crv_l, std::array<double, 3> col_crv = {1., 0., 0.}, float line_width = 2.f)
{
   std::vector<gbs::crv_dsp<double, 2, false>> crv_dsp;
   for (auto &crv : crv_l)
   {
      crv_dsp.push_back(
          gbs::crv_dsp<double, 2, false>{
              .c = &(crv),
              .col_crv = col_crv,
              .poles_on = false,
              .line_width = line_width,
          });
   }
   return crv_dsp;
};

TEST(tests_io, meridian_channel_msh)
{
   std::vector<gbs::BSCurve2d_d> crv_m;
   std::vector<gbs::BSCurve2d_d> crv_l;
   std::vector<double> u_l;
   std::vector<double> u_m;
   build_channel_curves(crv_m, crv_l, u_m, u_l);
   //Buildind of blend functions
   auto alpha_i = build_tfi_blend_function(u_l, false);
   auto beta_j = build_tfi_blend_function(u_m, false);
   auto n_ksi = u_l.size();
   auto n_eth = u_m.size();

   auto X1 = [&](auto ksi, double eth) {
      auto X1_ = std::array<double, 2>{0., 0.};
      for (auto i = 0; i < n_ksi; i++)
      {
         X1_ += alpha_i[i](ksi) * crv_l[i].value(eth);
      }
      return X1_;
   };

   auto X2 = [&](auto ksi, double eth) {
      auto X2_ = X1(ksi, eth);
      for (auto j = 0; j < n_eth; j++)
      {
         X2_ += beta_j[j](eth) * (crv_m[j].value(ksi) - X1(ksi, u_m[j]));
      }
      return X2_;
   };

   gbs::points_vector<double, 2> pts;
   auto nj = 50;
   auto ni = 101;
   for (auto j = 0; j < nj; j++)
   {
      // auto j = 0;
      for (auto i = 0; i < ni; i++)
      {
         pts.push_back(
             X2(j / (nj - 1.), i / (ni - 1.)));
      }
   }

   if(PLOT_ON)
   {
      auto crv_l_dsp = f_dspc(crv_m);
      auto crv_m_dsp = f_dspc(crv_l, {0., 1., 0.}, 4.f);

      gbs::plot(
          crv_l_dsp,
          crv_m_dsp,
          pts);
   }
}

TEST(tests_io, meridian_channel_msh2)
{
   std::vector<gbs::BSCurve2d_d> crv_m;
   std::vector<gbs::BSCurve2d_d> crv_l;
   std::vector<double> u_l;
   std::vector<double> u_m;
   build_channel_curves(crv_m, crv_l, u_m, u_l);
   //Buildind of blend functions
   auto n_ksi = u_l.size();
   auto n_eth = u_m.size();

   const auto P = 2;
   const auto Q = 2;
   auto alpha_i = build_tfi_blend_function_with_derivatives<double, P>(u_l);
   auto beta_j = build_tfi_blend_function_with_derivatives<double, Q>(u_m);

   auto X1 = [&](auto ksi, double eth) {
      auto X1_ = std::array<double, 2>{0., 0.};
      for (auto i = 0; i < n_ksi; i++)
      {
         for (auto n = 0; n < P; n++)
         {
            X1_ += alpha_i[i][n](ksi) * crv_l[i].value(eth, n);
         }
      }
      return X1_;
   };

   auto X2 = [&](auto ksi, double eth) {
      auto X2_ = X1(ksi, eth);
      for (auto j = 0; j < n_eth; j++)
      {
         for (auto m = 0; m < Q; m++)
         {
            X2_ += beta_j[j][m](eth) * (crv_m[j].value(ksi, m) - X1(ksi, u_m[j])); // works only for the specific 2d planar case
         }
      }
      return X2_;
   };

   gbs::points_vector<double, 2> pts;
   auto nj = 50;
   auto ni = 101;
   for (auto j = 0; j < nj; j++)
   {
      // auto j = 0;
      for (auto i = 0; i < ni; i++)
      {
         pts.push_back(
             X2(j / (nj - 1.), i / (ni - 1.)));
      }
   }

   if(PLOT_ON)
   {
      auto crv_l_dsp = f_dspc(crv_m);
      auto crv_m_dsp = f_dspc(crv_l, {0., 1., 0.}, 4.f);

      gbs::plot(
          crv_l_dsp,
          crv_m_dsp,
          pts);
   }
}

TEST(tests_io, meridian_channel_msh3)
{
   std::vector<gbs::BSCurve2d_d> crv_m;
   std::vector<gbs::BSCurve2d_d> crv_l;
   std::vector<double> u_l;
   std::vector<double> u_m;
   build_channel_curves(crv_m, crv_l, u_m, u_l);
   //Buildind of blend functions
   auto alpha_i = build_tfi_blend_function(u_l, true);
   auto beta_j = build_tfi_blend_function(u_m, true);
   auto n_ksi = u_l.size();
   auto n_eth = u_m.size();

   auto X1 = [&](auto ksi, double eth) {
      auto X1_ = std::array<double, 2>{0., 0.};
      for (auto i = 0; i < n_ksi; i++)
      {
         X1_ += alpha_i[i](ksi) * crv_l[i].value(eth);
      }
      return X1_;
   };

   auto X2 = [&](auto ksi, double eth) {
      auto X2_ = X1(ksi, eth);
      for (auto j = 0; j < n_eth; j++)
      {
         X2_ += beta_j[j](eth) * (crv_m[j].value(ksi) - X1(ksi, u_m[j]));
      }
      return X2_;
   };

   gbs::points_vector<double, 2> pts;
   auto nj = 50;
   auto ni = 101;
   for (auto j = 0; j < nj; j++)
   {
      // auto j = 0;
      for (auto i = 0; i < ni; i++)
      {
         pts.push_back(
             X2(j / (nj - 1.), i / (ni - 1.)));
      }
   }
   if(PLOT_ON)
   {
      auto crv_l_dsp = f_dspc(crv_m);
      auto crv_m_dsp = f_dspc(crv_l, {0., 1., 0.}, 4.f);

      gbs::plot(
          crv_l_dsp,
          crv_m_dsp,
          pts);
   }
}

TEST(tests_io, meridian_channel_ed_msh)
{
   std::vector<gbs::BSCurve2d_d> crv_m;
   std::vector<gbs::BSCurve2d_d> crv_l;
   std::vector<double> u_l;
   std::vector<double> u_m;
   build_channel_curves(crv_m, crv_l, u_m, u_l);
   auto hub = std::make_shared<gbs::BSCurve2d_d>( crv_m.front() );
   auto inlet = std::make_shared<gbs::BSCurve2d_d>( crv_l.front() );
   auto shr = std::make_shared<gbs::BSCurve2d_d>( crv_m.back() );
   auto out = std::make_shared<gbs::BSCurve2d_d>( crv_l.back() );
   // auto p_hub = std::make_shared<gbs::BSCurve<double,2>>(&hub);

   auto ni = 11;
   auto nj = 80;
   auto ed_hub = msh_edge<double, 2>(hub);
   ed_hub.set_n_points(nj);
   ed_hub.compute_pnts();
   auto ed_inl = msh_edge<double, 2>(inlet);
   ed_inl.set_n_points(ni);
   ed_inl.update_law(0.01, 0.01, 1.2, 1.2);
   ed_inl.compute_pnts();
   auto ed_shr = msh_edge<double, 2>(shr);
   ed_shr.set_n_points(nj);
   ed_shr.compute_pnts();
   auto ed_out = msh_edge<double, 2>(out);
   ed_out.set_n_points(ni);
   ed_out.update_law(0.01, 0.01, 1.2, 1.2);
   ed_out.compute_pnts();
   // ed_msh.update_law(0.01,0.01);

   std::vector<gbs::points_vector<double, 2>> X_ksi{ed_hub.points(), ed_shr.points()}, X_eth{ed_inl.points(), ed_out.points()};
   auto n_ksi = X_ksi.size();
   auto n_eth = X_eth.size();

   std::vector<double> ksi_i{0., ni - 1.};
   std::vector<double> eth_j{0., nj - 1.};
   auto alpha_i = build_tfi_blend_function(ksi_i, true);
   auto beta_j = build_tfi_blend_function(eth_j, true);

   ASSERT_DOUBLE_EQ(alpha_i[0](0), 1.);
   ASSERT_DOUBLE_EQ(alpha_i[0](ni - 1.), 0.);
   ASSERT_DOUBLE_EQ(alpha_i[1](0), 0.);
   ASSERT_DOUBLE_EQ(alpha_i[1](ni - 1.), 1.);

   ASSERT_DOUBLE_EQ(beta_j[0](0), 1.);
   ASSERT_DOUBLE_EQ(beta_j[0](nj - 1.), 0.);
   ASSERT_DOUBLE_EQ(beta_j[1](0), 0.);
   ASSERT_DOUBLE_EQ(beta_j[1](nj - 1.), 1.);

   auto X1 = [&](auto ksi, double eth) {
      auto X1_ = std::array<double, 2>{0., 0.};
      for (auto i = 0; i < n_ksi; i++)
      {
         X1_ += alpha_i[i](ksi) * X_ksi[i][eth];
      }
      return X1_;
   };

   auto X2 = [&](auto ksi, double eth) {
      auto X2_ = X1(ksi, eth);
      for (auto j = 0; j < n_eth; j++)
      {
         X2_ += beta_j[j](eth) * (X_eth[j][ksi] - X1(ksi, eth_j[j]));
      }
      return X2_;
   };

   gbs::points_vector<double, 2> pts;
   for (auto j = 0; j < nj; j++)
   {
      // auto j = 0;
      for (auto i = 0; i < ni; i++)
      {
         pts.push_back(
             X2(i, j));
      }
   }
   if (PLOT_ON)
   {
      auto grid_actor = make_structuredgrid_actor(pts, ni, nj);

      std::array<double, 3> green_col{0., 1., 0.};
      gbs::plot(
          grid_actor,
          f_dspc({*hub}),
          make_actor(ed_hub.points(), 15., true, green_col),
          f_dspc({*inlet}),
          make_actor(ed_inl.points(), 15., true, green_col),
          f_dspc({*shr}),
          make_actor(ed_shr.points(), 15., true, green_col),
          f_dspc({*out}),
          make_actor(ed_out.points(), 15., true, green_col)
          // make_actor(pts,15.,true,std::array<double,3>{0.,1.,0.}.data())
      );
   }
}

TEST(tests_io, meridian_multi_channel)
{
   std::vector<gbs::BSCurve2d_d> crv_m;
   std::vector<gbs::BSCurve2d_d> crv_l;
   std::vector<double> u_l;
   std::vector<double> u_m;
   build_channel_curves(crv_m, crv_l, u_m, u_l);

   crv_m[0].trim(u_l[0], u_l[1]);
   crv_m[1].trim(u_l[0], u_l[1]);
   crv_m[2].trim(u_l[0], u_l[1]);

   gbs::print(crv_m[0]);

   if(PLOT_ON)
      gbs::plot(
         crv_m[0],
         crv_m[1],
         crv_m[2],
         crv_l[0],
         crv_l[1]);
}

TEST(tests_io, iges_curves)
{

   std::vector<double> k1 = {0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1.};
   std::vector<double> k2 = {0., 0., 0., 0.33, 0.66, 1., 1., 1.};
   std::vector<std::array<double, 4>> poles1 =
       {
           {0., 1., 0., 1},
           {1., 1., 0., 1.5},
           {1., 2., 0., 1.},
           {2., 3., 0., 0.5},
           {3., 3., 0., 1.2},
           {4., 2., 0., 1.},
       };

   std::vector<std::array<double, 3>> poles2 =
       {
           {4., 2., 0.},
           {5., 1., 1.},
           {5., 0., 1.},
           {4., -1., 1.},
           {3., -1., 0.},
       };
   size_t p1 = 2;
   size_t p2 = 2;

   gbs::BSCurveRational3d_d c1(poles1, k1, p1);
   gbs::BSCurve3d_d c2(poles2, k2, p2);

   // instantiate the IGES data object
   DLL_IGES model;
   gbs::add_geom(c1, model);
   gbs::add_geom(c2, model);
   // gbs::add_geom(srf, model);

   std::string dir = get_directory(__FILE__);
   model.Write((dir+"/out/tests_io_iges_curves.igs").c_str(), true);
}

TEST(tests_io, iges_surfaces)
{

   size_t nu = 30;
   size_t nv = 25;
   std::vector<std::array<double, 3>> points(nu * nv);
   std::vector<double> u(nu);
   std::vector<double> v(nv);
   auto x{0.};
   auto y{0.};
   auto z{0.};
   for (auto i = 0; i < nu; i++)
      u[i] = i / (nu - 1.);
   for (auto j = 0; j < nv; j++)
      v[j] = j / (nv - 1.);
   for (auto j = 0; j < nv; j++)
   {
      y = v[j] - 0.5;
      for (auto i = 0; i < nu; i++)
      {
         x = u[i] - 0.5;
         z = std::sin(std::numbers::pi * (std::sqrt(x * x + y * y)));
         points[i + j * nu] = {x, y, z};
      }
   }

   size_t n_poles_u = 15;
   size_t n_poles_v = 10;
   size_t deg_u = 5;
   size_t deg_v = 4;
   auto ku = gbs::build_simple_mult_flat_knots(0., 1., n_poles_u, deg_u);
   auto kv = gbs::build_simple_mult_flat_knots(0., 1., n_poles_v, deg_v);
   auto poles = gbs::approx(points, ku, kv, u, v, deg_u, deg_v);
   gbs::BSSurface<double,3> srf(poles, ku, kv, deg_u, deg_v);

    auto p_stream1 = std::make_shared<gbs::BSCurve<double,2>>(
        gbs::BSCurve<double,2>{
        {
            {0.,0.5},
            {0.3,0.5},
            {1.,0.5},
            {1.,1.0},
            {1.,1.5}
        }
        ,
        {0.,0.,0.,0.,0.3,1.,1.,1.,1.}
        ,
        3
    });

    gbs::ax2<double,3> ax {{
            {0.,0.,0.},
            {0.,0.,1.},
            {1.,0.,0.}
    }};

    auto p_stream_sheet1 = std::make_shared<gbs::SurfaceOfRevolution<double>>(
        gbs::SurfaceOfRevolution<double>{
            p_stream1,
            ax
        }
    );


   DLL_IGES model;
   gbs::add_geom(srf, model);
   gbs::add_geom(*p_stream_sheet1, model);
   std::string dir = get_directory(__FILE__);
   model.Write((dir+"/out/tests_io_iges_surfaces.igs").c_str(), true);
}
namespace
{
    // Parameter data of the single entity of an IGES file written by IgesWriter:
    // the P records (column 73) concatenated, numbers split on ',' up to ';'.
    std::vector<double> iges_parameters(const std::string &file)
    {
        std::ifstream in(file);
        std::string line, data;
        while (std::getline(in, line))
            if (line.size() >= 73 && line[72] == 'P')
                data += line.substr(0, 64);
        data = data.substr(0, data.find(';'));
        std::vector<double> values;
        std::stringstream ss(data);
        std::string tok;
        while (std::getline(ss, tok, ','))
            values.push_back(std::stod(tok));
        return values;
    }

    std::string iges_tmp(const std::string &name)
    {
        return (std::filesystem::temp_directory_path() / name).string();
    }
}

// IGES 126 / 128 store Cartesian control points and separate weights, always three
// coordinates; gbs stores rational poles in homogeneous form (w x, w y, [w z,] w).
TEST(tests_io, iges_rational_and_2d_poles)
{
    // rational 3D circle of radius 2: data = 126, K, M, 4 flags, K+M+2 knots, K+1 weights, 3(K+1) coordinates, V0, V1, normal
    {
        const auto c = build_circle<double, 3>(2., {1., 0., 0.});
        IgesWriter<double> w;
        w.add_geometry(c);
        const auto f = iges_tmp("gbs_writer_circle3d.igs");
        w.write(f);
        const auto d = iges_parameters(f);
        ASSERT_EQ(int(d[0]), 126);
        const int K = int(d[1]), M = int(d[2]);
        ASSERT_EQ(K + 1, int(c.poles().size()));
        ASSERT_EQ(d[5], 0.); // PROP3 = 0: rational
        const size_t iw = 7 + K + M + 2, ip = iw + K + 1;
        const auto cart = c.polesProjected();
        const auto wts = c.weights();
        for (int i = 0; i <= K; ++i)
        {
            ASSERT_NEAR(d[iw + i], wts[i], 1e-9);
            for (int k = 0; k < 3; ++k)
                ASSERT_NEAR(d[ip + 3 * i + k], cart[i][k], 1e-9);
        }
    }
    // rational 2D circle: z = 0 padding, Cartesian poles
    {
        const auto c = build_circle<double, 2>(1.5, {0.5, -1.});
        IgesWriter<double> w;
        w.add_geometry(c);
        const auto f = iges_tmp("gbs_writer_circle2d.igs");
        w.write(f);
        const auto d = iges_parameters(f);
        const int K = int(d[1]), M = int(d[2]);
        const size_t iw = 7 + K + M + 2, ip = iw + K + 1;
        const auto cart = c.polesProjected();
        for (int i = 0; i <= K; ++i)
        {
            ASSERT_NEAR(d[ip + 3 * i + 0], cart[i][0], 1e-9);
            ASSERT_NEAR(d[ip + 3 * i + 1], cart[i][1], 1e-9);
            ASSERT_NEAR(d[ip + 3 * i + 2], 0., 1e-15);
        }
    }
    // non-rational 2D segment: unit weights, z = 0
    {
        const auto s = build_segment<double, 2>({0., 1.}, {2., 3.});
        IgesWriter<double> w;
        w.add_geometry(s);
        const auto f = iges_tmp("gbs_writer_segment2d.igs");
        w.write(f);
        const auto d = iges_parameters(f);
        const int K = int(d[1]), M = int(d[2]);
        const size_t iw = 7 + K + M + 2, ip = iw + K + 1;
        ASSERT_EQ(K, 1);
        ASSERT_NEAR(d[iw], 1., 1e-15);
        ASSERT_NEAR(d[ip + 0], 0., 1e-15);
        ASSERT_NEAR(d[ip + 1], 1., 1e-15);
        ASSERT_NEAR(d[ip + 2], 0., 1e-15);
        ASSERT_NEAR(d[ip + 3], 2., 1e-15);
        ASSERT_NEAR(d[ip + 4], 3., 1e-15);
        ASSERT_NEAR(d[ip + 5], 0., 1e-15);
    }
    // rational surface: data = 128, K1, K2, M1, M2, 5 flags, knots u, knots v, weights, coordinates
    {
        const auto circle = build_circle<double, 3>(1.);
        points_vector<double, 4> poles;
        for (double z : {0., 2.})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z;
                poles.push_back(p);
            }
        const BSSurfaceRational<double, 3> srf{poles, circle.knotsFlats(), std::vector<double>{0., 0., 2., 2.}, circle.degree(), 1};
        IgesWriter<double> w;
        w.add_geometry(srf);
        const auto f = iges_tmp("gbs_writer_cylinder.igs");
        w.write(f);
        const auto d = iges_parameters(f);
        ASSERT_EQ(int(d[0]), 128);
        const int K1 = int(d[1]), K2 = int(d[2]), M1 = int(d[3]), M2 = int(d[4]);
        const size_t n = size_t(K1 + 1) * size_t(K2 + 1);
        ASSERT_EQ(n, srf.poles().size());
        const size_t iw = 10 + (K1 + M1 + 2) + (K2 + M2 + 2), ip = iw + n;
        const auto cart = srf.polesProjected();
        const auto wts = srf.weights();
        for (size_t i = 0; i < n; ++i)
        {
            ASSERT_NEAR(d[iw + i], wts[i], 1e-9);
            ASSERT_NEAR(std::hypot(d[ip + 3 * i], d[ip + 3 * i + 1]), std::hypot(cart[i][0], cart[i][1]), 1e-9);
            ASSERT_NEAR(d[ip + 3 * i + 2], cart[i][2], 1e-9);
        }
    }
}
