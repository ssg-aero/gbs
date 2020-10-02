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
     using PointArray = std::vector< std::array<T,dim> >;
}