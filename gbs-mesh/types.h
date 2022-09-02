#pragma once

    #if defined(__CUDACC__)
        #define _CPU_GPU_ __host__ __device__
        #define _GPU_ __device__
        #define _CPU_ __host__
        #include <thrust/device_vector.h>
        #include <thrust/tuple.h>
    #else
        #define _CPU_GPU_
        #define _GPU_
        #define _CPU_
        #include <vector>
        #include <tuple>
    #endif


namespace gbs
{
    using Real = float;
    #if defined(__CUDACC__)
        template<typename T>
        using vector = thrust::device_vector<T>;
        using Real3 = thrust::tuple<Real,Real,Real>;
        // template<size_t i>
        // using get = thrust::get<i>;
        using namespace thrust;
    #else
        template<typename T>
        using vector = std::vector<T>;
        using Real3 = std::tuple<Real,Real,Real>;
        // template<size_t i>
        // using get = std::get<i>;
        using namespace std;
    #endif
}
