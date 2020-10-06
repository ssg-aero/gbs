#pragma once

#if defined(_WIN32)
#ifndef GBS_EXPORT
#define GBS_EXPORT __declspec( dllexport )
// For global variables :
#define GBS_EXPORTEXTERN __declspec( dllexport ) extern
#define GBS_EXPORTEXTERNC extern "C" __declspec( dllexport )
#endif  /* GBS_EXPORT */

#ifndef GBS_IMPORT
#define GBS_IMPORT __declspec( dllimport ) extern
#define GBS_IMPORTC extern "C" __declspec( dllimport )
#endif  /* GBS_IMPORT */
#endif  /*_WIN32 */

#include <vector>
#include <array>
namespace gbs
{
const double knot_eps = 1e-7;

template <typename T,size_t dim>
     using point = std::array<T,dim>;
template <typename T,size_t dim>
     using ax1 = std::array<point<T,dim>,2>;
template <typename T,size_t dim>
     using ax2 = std::array<point<T,dim>,3>;
template <typename T,size_t dim>
     using points_vector = std::vector< point<T,dim> >;

using points_vector_2d_f = points_vector<float, 2>;
using points_vector_2d_d = points_vector<double, 2>;
using points_vector_3d_f = points_vector<float, 3>;
using points_vector_3d_d = points_vector<double, 3>;

}