# cmake_minimum_required(VERSION 3.17.0)

project(gbs-mesh LANGUAGES CXX)

INSTALL_HEADERS(${PROJECT_NAME})

if(${USE_CUDA})
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    if(%{CUDA_ARCHITECTURES})
        set(CMAKE_CUDA_ARCHITECTURES %{CUDA_ARCHITECTURES})
    else()
        set(CMAKE_CUDA_ARCHITECTURES 61) # Max 61 for SB2
    endif()
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif()