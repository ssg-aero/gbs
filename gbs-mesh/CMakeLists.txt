project(gbs-mesh LANGUAGES CXX)

INSTALL_HEADERS(${PROJECT_NAME})

if(${GBS_BUILD_TESTS} AND ${USE_CUDA})
    add_subdirectory(tests)
endif()