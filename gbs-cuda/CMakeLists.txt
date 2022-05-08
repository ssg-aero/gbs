# cmake_minimum_required(VERSION 3.18)

# project(gbs-cuda LANGUAGES CXX CUDA)

# set( CMAKE_INCLUDE_CURRENT_DIR ON )

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 61) # Max 61 for SB2
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# add_executable(trinter src/trinter.cu)

add_subdirectory(tests)

