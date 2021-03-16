#include <vector>
#include <gbs/gbslib.h>
namespace gbs
{
template<typename T,size_t P>
auto build_tfi_blend_function_with_derivatives(const std::vector<T> &ksi_i) -> std::vector<std::array<gbs::BSCurve<T,1>,P>>
{
   auto n_ksi_i = ksi_i.size();
   std::vector<std::array<gbs::BSCurve<T, 1>, P>> alpha_i;

   for (int i = 0; i < n_ksi_i; i++)
   {
      alpha_i.push_back(std::array<gbs::BSCurve<double,1>,P>{});
      for (auto n = 0; n < P; n++)
      {
         std::vector<gbs::constrType<double, 1, P+1>> dji{n_ksi_i,{0.}};
         dji[i][n] = {1.};
         alpha_i.back()[n]=gbs::interpolate(dji, ksi_i);
      }
   }
   return alpha_i;
}

template<typename T>
auto build_tfi_blend_function(const std::vector<T> &ksi_i,bool slope_control) -> std::vector<gbs::BSCurve<T,1>>
{
   auto n_ksi_i = ksi_i.size();
   std::vector<gbs::BSCurve<double, 1>> alpha_i;

   for(int i = 0 ; i < n_ksi_i; i++)
   {
      if (slope_control)
      {
         std::vector<gbs::constrType<double, 1, 2>> dji{n_ksi_i,{0.}};
         dji[i][0] = {1.};
         alpha_i.push_back(gbs::interpolate(dji,ksi_i));
      }
      else
      {
         if ( n_ksi_i == 2)
         {
            std::vector<gbs::constrType<double, 1, 1>> dji{n_ksi_i,{0.}};
            dji[i][0] = {1.};
            alpha_i.push_back(gbs::interpolate(dji,ksi_i));
         }
         else
         {
            gbs::points_vector<double, 1> dji{n_ksi_i, {0.}};
            dji[i] = {1.};
            alpha_i.push_back(gbs::interpolate(dji, ksi_i, 2));
         }
      }
   }
   return alpha_i;
}
}