#pragma once

// Apple libc++ does not ship parallel execution policies.
// Fall back to sequential algorithms on that platform.
#include <version>
#include <numeric>  // std::reduce, std::transform_reduce (previously transitive via <execution>)
#include <algorithm>  // std::transform
#include <iterator>   // std::distance, std::begin/end/size
#include <cstddef>    // std::size_t
#include <vector>     // std::vector (build_batch output)
#include <type_traits> // std::decay_t (build_batch result type)

#include "gbsconstants.h"  // gbs::parallel_min_size

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

    // -------------------------------------------------------------------------
    // Size-gated parallel transform — the shared helper for bulk evaluators.
    // -------------------------------------------------------------------------
    // Dispatches std::transform to the parallel execution policy (GBS_PAR_EXEC)
    // only when the input range has at least `gbs::parallel_min_size` elements,
    // and runs serially below that. The parallel STL (libstdc++/MSVC) does not
    // cheaply fall back to serial for small ranges: it still spins up the TBB
    // task machinery (~9-15 us fixed cost), so an *unconditional* `par` regresses
    // small calls by up to 25x. The gate avoids that regression while keeping the
    // 2-21x win on large inputs. See bench/PARALLEL_AUDIT.md (#84/#92).
    //
    // Determinism: each iteration writes a distinct, pre-indexed output slot and
    // there is no reduction, so the result is BIT-IDENTICAL regardless of the
    // policy chosen or the thread count. The gate only selects the policy by
    // size; it never changes the result.
    //
    // Portability: on platforms without parallel policies (Apple libc++)
    // GBS_PAR_EXEC is empty, so both branches are serial and the gate is a
    // correctness-preserving no-op.

    /// Size-gated unary std::transform (single input range).
    template <class InputIt, class OutputIt, class UnaryOp>
    OutputIt transform_threshold(InputIt first, InputIt last, OutputIt out, UnaryOp op)
    {
        if (static_cast<std::size_t>(std::distance(first, last)) >= parallel_min_size)
            return std::transform(GBS_PAR_EXEC first, last, out, op);
        return std::transform(first, last, out, op);
    }

    /// Size-gated binary std::transform (two input ranges).
    template <class InputIt1, class InputIt2, class OutputIt, class BinaryOp>
    OutputIt transform_threshold(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt out, BinaryOp op)
    {
        if (static_cast<std::size_t>(std::distance(first1, last1)) >= parallel_min_size)
            return std::transform(GBS_PAR_EXEC first1, last1, first2, out, op);
        return std::transform(first1, last1, first2, out, op);
    }

    // -------------------------------------------------------------------------
    // Size-gated parallel BATCH build — the shared helper for batch builders.
    // -------------------------------------------------------------------------
    // Applies `build` to each input entity and returns the results in input
    // order. The sibling of transform_threshold for the interp/approx BATCH axis
    // (#91): building N INDEPENDENT curves/surfaces (loft-section sweeps,
    // multi-patch models, fitting batches) is embarrassingly parallel — each
    // build owns its solver/matrix/poles. Unlike transform_threshold, the gate is
    // the COUNT of builds (gbs::parallel_batch_min_size), not a per-element count:
    // each "element" here is a whole banded solve, so the crossover sits far below
    // parallel_min_size (see gbsconstants.h / PARALLEL_AUDIT.md [7]/[9]).
    //
    // Reentrancy / determinism: every build constructs its own object with NO
    // shared mutable state (the band-LU/sparse solver of #96 is a per-call stack
    // local) and writes a distinct, pre-indexed output slot, with no reduction.
    // The result is therefore BIT-IDENTICAL to the serial loop regardless of the
    // policy chosen or the thread count (#91, cf. #48). Below the gate, and on
    // platforms without parallel policies (Apple libc++ -> GBS_PAR_EXEC empty),
    // both branches are serial and the gate is a correctness-preserving no-op.
    //
    // `Result` (the builder's return type) must be default-constructible, since
    // the output is pre-sized for the parallel slot-by-slot write; BSCurve /
    // BSSurface are.
    template <class InputRange, class Builder>
    auto build_batch(const InputRange &inputs, Builder build,
                     std::size_t min_size = parallel_batch_min_size)
        -> std::vector<std::decay_t<decltype(build(*std::begin(inputs)))>>
    {
        using Result = std::decay_t<decltype(build(*std::begin(inputs)))>;
        std::vector<Result> out(static_cast<std::size_t>(std::size(inputs)));
        if (out.size() >= min_size)
            std::transform(GBS_PAR_EXEC std::begin(inputs), std::end(inputs), out.begin(), build);
        else
            std::transform(std::begin(inputs), std::end(inputs), out.begin(), build);
        return out;
    }
}

