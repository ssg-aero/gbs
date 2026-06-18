#pragma once
//
// Zero-copy-ish return for the bulk evaluators (#97).
//
// The bulk evaluators return a points_vector<T,dim> == std::vector<std::array<T,dim>>,
// whose storage is N*dim contiguous T (std::array<T,dim> is tightly packed).
// Going through py::cast(...) marshals that into a Python list of length-dim
// lists before numpy copies it back — O(N python objects), which dwarfs the
// (parallel, #88) evaluation itself (~500 ms vs ~11 ms for N=1e6). Instead we
// build the (N,dim) numpy array directly and memcpy the contiguous block once.
//
#include <pybind11/numpy.h>
#include <array>
#include <cstring>
#include <vector>

namespace py = pybind11;

// Build an (N, dim) numpy array from a contiguous vector<array<T,dim>> with a
// single memcpy. The returned array owns its buffer (a plain value return, no
// lifetime coupling to the temporary result).
template <typename T, std::size_t dim>
inline py::array_t<T> points_to_numpy(const std::vector<std::array<T, dim>> &pts)
{
    static_assert(sizeof(std::array<T, dim>) == dim * sizeof(T),
                  "array<T,dim> must be tightly packed for the bulk memcpy");
    const std::size_t n = pts.size();
    py::array_t<T> out({n, dim});
    if (n)
        std::memcpy(out.mutable_data(), pts.data(), n * dim * sizeof(T));
    return out;
}
