#pragma once
#include <vector>
#include <stdexcept>

namespace occt_utils
{
    // template <typename T>
    // std::vector<T> make_range(T v1, T v2, size_t n)
    // {
    //     if (n < 2)
    //     {
    //         throw std::length_error("2 points a required for a range");
    //     }
            
    //     std::vector<T> v(n);
    //     T step = ( v2 -v1 ) / T( n - 1);

    //     // T v2_check = v2-step;
    //     // std::generate(v.begin(),v.end(),[&,v_ = v1-step] () mutable 
    //     // { 
    //     //     return v_ >= v2_check ? v2 : v_+=step ; 
    //     // });
    //     // T v2_check = v2-step;
    //     //TODO check wich is faster
    //     v.front() = v1;
    //     std::generate(
    //         std::next(v.begin(),1),
    //         std::next(v.end()-1),
    //         [&,v_ = v1] () mutable 
    //             { 
    //                 return v_+=step ; 
    //             }
    //     );
    //     v.back() = v2;

    //     return v;
    // }

} // namespace occt_utils