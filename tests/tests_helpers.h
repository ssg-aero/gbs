#pragma once
#include <chrono>

template<typename T, typename F, typename... Args>
auto measure_execution_time(F func, Args&&... args) {
    auto start_time = std::chrono::steady_clock::now();
    func(std::forward<Args>(args)...);
    auto end_time = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<T>(end_time - start_time).count();
}

// Function to extract the directory from a file path
std::string get_directory(const std::string& file_path) {
    size_t found = file_path.find_last_of("/\\");
    return found != std::string::npos ? file_path.substr(0, found) : "";
}