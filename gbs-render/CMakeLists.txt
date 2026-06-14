
project(gbs-render LANGUAGES CXX)

# gbslib
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB GBS_RENDER_SRC_LIST "*.cpp")
add_library(${PROJECT_NAME} SHARED ${GBS_RENDER_SRC_LIST})
target_link_libraries(${PROJECT_NAME}
    gbs
    ${VTK_LIBRARIES}
    Eigen3::Eigen
    TBB::tbb
    ${GBS_MODULES}
)

# examples (scenes VTK interactives, ex-tests sans assertion).
# Opt-in : code VTK interactif, non executable en CI et sensible a la version de
# VTK (certains headers varient). On ne le construit que sur demande explicite.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(GBS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

vtk_module_autoinit(
        TARGETS ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_NAME render)
install(TARGETS ${PROJECT_NAME} EXPORT GbsTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
INSTALL_HEADERS(${PROJECT_NAME})
if ( MSVC )
        SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )