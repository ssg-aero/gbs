cmake_minimum_required(VERSION 3.18)

project(gbs-core LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(gbs_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GBS_CORE_HEADERS
    ${gbs_INCLUDE_DIR}/gbs/basisfunctions.h
    ${gbs_INCLUDE_DIR}/gbs/bscanalysis.h
    ${gbs_INCLUDE_DIR}/gbs/bscapprox.h
    ${gbs_INCLUDE_DIR}/gbs/bscbuild.h
    ${gbs_INCLUDE_DIR}/gbs/bscinterp.h
    ${gbs_INCLUDE_DIR}/gbs/bsctools.h
    ${gbs_INCLUDE_DIR}/gbs/bscurve.h
    ${gbs_INCLUDE_DIR}/gbs/bssanalysis.h
    ${gbs_INCLUDE_DIR}/gbs/bssapprox.h
    ${gbs_INCLUDE_DIR}/gbs/bssbuild.h
    ${gbs_INCLUDE_DIR}/gbs/bssinterp.h
    ${gbs_INCLUDE_DIR}/gbs/bsstools.h
    ${gbs_INCLUDE_DIR}/gbs/bssurf.h
    ${gbs_INCLUDE_DIR}/gbs/curvecomposite.h
    ${gbs_INCLUDE_DIR}/gbs/curveline.h
    ${gbs_INCLUDE_DIR}/gbs/curveoffset.h
    ${gbs_INCLUDE_DIR}/gbs/curveonsurface.h
    ${gbs_INCLUDE_DIR}/gbs/curvereparametrized.h
    ${gbs_INCLUDE_DIR}/gbs/curves
    ${gbs_INCLUDE_DIR}/gbs/curvescheck.h
    ${gbs_INCLUDE_DIR}/gbs/curvetrimmed.h
    ${gbs_INCLUDE_DIR}/gbs/exceptions.h
    ${gbs_INCLUDE_DIR}/gbs/extrema.h
    ${gbs_INCLUDE_DIR}/gbs/gbslib.h
    ${gbs_INCLUDE_DIR}/gbs/knotsfunctions.h
    ${gbs_INCLUDE_DIR}/gbs/maths.h
    ${gbs_INCLUDE_DIR}/gbs/offsets.h
    ${gbs_INCLUDE_DIR}/gbs/solvers.h
    ${gbs_INCLUDE_DIR}/gbs/surfaceoffset.h
    ${gbs_INCLUDE_DIR}/gbs/surfaceofrevolution.h
    ${gbs_INCLUDE_DIR}/gbs/surfaces
    ${gbs_INCLUDE_DIR}/gbs/transform.h
    ${gbs_INCLUDE_DIR}/gbs/transformpoints.h
    ${gbs_INCLUDE_DIR}/gbs/vecop.h
)
add_library(gbs-core INTERFACE ${gbs_INCLUDE_DIR})
add_library(gbs::core ALIAS gbs-core)
add_library(gbs ALIAS gbs-core)

# Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# TBB
find_package(TBB REQUIRED)
# Boost
find_package(Boost REQUIRED)
# NLOPT
find_package(NLopt REQUIRED)

target_include_directories(
    gbs-core
    INTERFACE
        $<BUILD_INTERFACE:${gbs_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
        ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(
    gbs-core
    INTERFACE
        TBB::tbb
        NLopt::nlopt
        ${Boost_LIBRARIES}
)
target_compile_features(gbs-core INTERFACE cxx_std_20)
set_property(TARGET gbs-core PROPERTY EXPORT_NAME gbs::core)

install(TARGETS gbs-core EXPORT gbs-targets)
install(FILES ${GBS_CORE_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gbs)

# CMake config
# Add helpers for cmake configuration
include(CMakePackageConfigHelpers)

# Definition
set(ConfigFileInstallDir ${CMAKECONFIG_INSTALL_DIR})
configure_package_config_file(gbsConfig.cmake.in
                            "${CMAKE_CURRENT_BINARY_DIR}/gbsConfig.cmake"
                            INSTALL_DESTINATION ${ConfigFileInstallDir})

# Makes the project importable from the build directory
export(EXPORT gbs-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/gbsTargets.cmake")

# Installation
install(EXPORT gbs-targets 
        FILE gbsTargets.cmake
        DESTINATION ${ConfigFileInstallDir})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gbsConfig.cmake
        DESTINATION ${ConfigFileInstallDir})
