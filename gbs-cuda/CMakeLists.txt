# cmake_minimum_required(VERSION 3.18)

project(gbs-cuda LANGUAGES CXX)

#set(CMAKE_CUDA_STANDARD 17)
#set(CMAKE_CUDA_ARCHITECTURES 61) # Max 61 for SB2
#set(CMAKE_CUDA_STANDARD_REQUIRED ON)

option(BUILD_GBS_CUDA_TESTS "Enable the build of the GBS CUDA extension tests" OFF)

# add_executable(trinter src/trinter.cu)

set(GBS_CUDA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_CUDA_HEADERS
    ${GBS_CUDA_INCLUDE_DIR}/gbs/cuda/basisfunctions.cuh
)
add_library(gbs-cuda INTERFACE ${GBS_CUDA_HEADERS})

target_include_directories(
    gbs-cuda
    INTERFACE
        $<BUILD_INTERFACE:${GBS_CUDA_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)

install(FILES ${GBS_CUDA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs/cuda)

if(${BUILD_GBS_CUDA_TESTS})
    add_subdirectory(tests)
endif()
