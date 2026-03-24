#pragma once

// Apple libc++ does not ship parallel execution policies.
// Fall back to sequential algorithms on that platform.
#include <version>

#if defined(__cpp_lib_execution)
    #include <execution>
    #define GBS_HAS_EXECUTION_POLICIES 1
#else
    #define GBS_HAS_EXECUTION_POLICIES 0
#endif

// Use GBS_PAR_EXEC as the first argument to parallel STL algorithms.
// Expands to  std::execution::par,  when available, or nothing otherwise.
#if GBS_HAS_EXECUTION_POLICIES
    #define GBS_PAR_EXEC std::execution::par,
    #define GBS_SEQ_EXEC std::execution::seq,
#else
    #define GBS_PAR_EXEC
    #define GBS_SEQ_EXEC
#endif

namespace gbs
{
#if GBS_HAS_EXECUTION_POLICIES
    using parallel_execution_policy   = std::execution::parallel_policy;
    using sequenced_execution_policy  = std::execution::sequenced_policy;
    inline constexpr auto par_exec = std::execution::par;
    inline constexpr auto seq_exec = std::execution::seq;
#else
    struct parallel_execution_policy  {};
    struct sequenced_execution_policy {};
    inline constexpr parallel_execution_policy  par_exec{};
    inline constexpr sequenced_execution_policy seq_exec{};
#endif
}

