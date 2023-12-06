
project(gbs-render LANGUAGES CXX)

# gbslib
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB GBS_RENDER_SRC_LIST "*.cpp")
add_library(${PROJECT_NAME} SHARED ${GBS_RENDER_SRC_LIST})
target_link_libraries(${PROJECT_NAME} 
    ${VTK_LIBRARIES}
    TBB::tbb
    vecop
)

# tests
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(${GBS_BUILD_TESTS})
    add_subdirectory(tests)
endif(${GBS_BUILD_TESTS})

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)

install(TARGETS ${PROJECT_NAME})
INSTALL_HEADERS(${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )