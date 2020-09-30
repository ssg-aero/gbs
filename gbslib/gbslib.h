#pragma once
namespace gbs
{
const double knot_eps = 1e-7;
template <typename T,size_t dim>
     using PointArray = std::vector< std::array<T,dim> >;
}