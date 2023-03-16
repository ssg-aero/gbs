#pragma once
#include <chrono>

template<typename T, typename F, typename... Args>
auto measure_execution_time(F func, Args&&... args) {
    auto start_time = std::chrono::steady_clock::now();
    func(std::forward<Args>(args)...);
    auto end_time = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<T>(end_time - start_time).count();
}