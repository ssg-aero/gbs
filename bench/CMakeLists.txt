cmake_minimum_required(VERSION 3.28)
project(gbs_bench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED)

# gbs core is header-only: just point at the repo root + inc/.
get_filename_component(GBS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)

add_executable(bench_arc_length_derivatives bench_arc_length_derivatives.cpp)
target_include_directories(bench_arc_length_derivatives PRIVATE ${GBS_ROOT} ${GBS_ROOT}/inc)
target_link_libraries(bench_arc_length_derivatives PRIVATE Eigen3::Eigen)



add_executable(bench_eval_compare bench_eval_compare.cpp)
target_include_directories(bench_eval_compare PRIVATE ${GBS_ROOT} ${GBS_ROOT}/inc)
target_link_libraries(bench_eval_compare PRIVATE Eigen3::Eigen)

add_executable(bench_interp bench_interp.cpp)
target_include_directories(bench_interp PRIVATE ${GBS_ROOT} ${GBS_ROOT}/inc)
target_link_libraries(bench_interp PRIVATE Eigen3::Eigen)

add_executable(bench_loft bench_loft.cpp)
target_include_directories(bench_loft PRIVATE ${GBS_ROOT} ${GBS_ROOT}/inc)
target_link_libraries(bench_loft PRIVATE Eigen3::Eigen)
