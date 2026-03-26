project(gbs-mesh LANGUAGES CXX)

INSTALL_HEADERS(${PROJECT_NAME})

if(GBS_BUILD_TESTS AND GBS_USE_CUDA)
    add_subdirectory(tests)
endif()