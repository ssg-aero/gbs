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
#else
     #define GBS_EXPORT 
#endif  /*_WIN32 */

#include <vector>
#include <array>
#include <tuple>
namespace gbs
{
template< typename T>
constexpr double knot_eps = std::numeric_limits<T>::epsilon() * 100;

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
using points_vector_4d_f = points_vector<float, 4>;
using points_vector_4d_d = points_vector<double, 4>;

/**
 * @brief Applies a given functor to each element of a tuple at compile time.
 * 
 * This is a compile-time recursive function that iterates over each element of the tuple and applies
 * the given functor to it. The recursion ends when the index reaches the size of the tuple.
 *
 * @tparam Tuple The type of the tuple.
 * @tparam Functor The type of the functor to apply to each tuple element.
 * @tparam Index The current index in the tuple, defaulting to 0.
 * @param tpl The tuple to iterate over.
 * @param f The functor to apply to each element of the tuple.
 */
template <typename Tuple, typename Functor, size_t Index = 0>
auto tuple_for_each(const Tuple &tpl, const Functor &f) -> void
{
     constexpr auto tuple_size = std::tuple_size_v<Tuple>;
     if constexpr (Index < tuple_size)
     {
          f(std::get<Index>(tpl)); // Apply the functor to the current element.
          tuple_for_each<Tuple, Functor, Index + 1>(tpl, f); // Recursive call to the next element.
     }
}

/**
 * @brief Conditionally dereferences the given argument if it's a pointer; otherwise, returns it as is.
 * 
 * @tparam T The type of the argument, deduced automatically.
 * @param arg The argument to potentially dereference.
 * @return A reference to the dereferenced object if the argument is a pointer, or the argument itself otherwise.
 */
template<typename T>
decltype(auto) conditional_dereference(T&& arg) {
    if constexpr (std::is_pointer_v<std::remove_reference_t<T>>) {
        return *arg;  // Dereference if arg is a pointer
    } else {
        return arg;  // Return as-is if arg is not a pointer
    }
}

} // namespace gbs