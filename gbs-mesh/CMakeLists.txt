project(gbs-mesh LANGUAGES CXX)

INSTALL_HEADERS(${PROJECT_NAME})

if(GBS_BUILD_TESTS)
    add_subdirectory(tests)
endif()