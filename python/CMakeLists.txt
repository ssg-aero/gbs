project(pygbs LANGUAGES CXX)

# gslib
include_directories(${CMAKE_SOURCE_DIR})

find_package(libIGES REQUIRED)

find_package (Python3 COMPONENTS Interpreter Development)

find_package(pybind11 REQUIRED)
pybind11_add_module(
        ${PROJECT_NAME}
        gbsBindCurves.cpp
        gbsBindSurfaces.cpp
        gbsbindrender.cpp 
        gbsbindmesh.cpp
        gbsbindinterp.cpp
        gbsbindapprox.cpp
        gbsbindbuildcurve.cpp
        gbsBindBuildSurfaces.cpp
        gbsBindCurveTools.cpp
        gbsBindSurfaceTools.cpp
        gbsBindShaping.cpp
        gbsBindDiscretize.cpp
        gbsBindBrep.cpp
        gbsBindBrepOps.cpp
        gbsbind.cpp 
        # ${SRC_LIST}
)
# the bindings of the BREP kernel: up to 4 GB to compile (GBS_HEAVY_COMPILE_JOBS)
if(GBS_HEAVY_COMPILE_JOBS GREATER 0)
        set_property(TARGET ${PROJECT_NAME} PROPERTY JOB_POOL_LINK gbs_heavy)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.2)
                # only the sources that include the BREP kernel: the others (up to
                # 1.8 GB) need not wait for them (a pool for a source file: CMake 4.2)
                get_target_property(pygbs_sources ${PROJECT_NAME} SOURCES)
                foreach(source ${pygbs_sources})
                        file(STRINGS ${source} brep_includes REGEX "#include <gbs-brep/")
                        if(brep_includes)
                                set_property(SOURCE ${source} PROPERTY JOB_POOL_COMPILE gbs_heavy)
                        endif()
                endforeach()
        else()
                set_property(TARGET ${PROJECT_NAME} PROPERTY JOB_POOL_COMPILE gbs_heavy)
        endif()
endif()
target_link_libraries(${PROJECT_NAME}
        PRIVATE
                gbs
                Eigen3::Eigen
                ${VTK_LIBRARIES}
                gbs-render
                TBB::tbb
                libIGES::iges
                ${GBS_MODULES}
)

# Exposed as pygbs.gbs.__version__ (single source: the top-level project() VERSION)
target_compile_definitions(${PROJECT_NAME} PRIVATE GBS_VERSION_STR="${CMAKE_PROJECT_VERSION}")

vtk_module_autoinit(
        TARGETS  ${PROJECT_NAME}
        MODULES ${VTK_LIBRARIES}
        )

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "gbs")

file(GLOB PY_SRC_LIST "*.py")

if($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
        DESTINATION $ENV{SP_DIR}/${PROJECT_NAME}
        )
        install(FILES 
                        ${PY_SRC_LIST}
                DESTINATION 
                        $ENV{SP_DIR}/${PROJECT_NAME}
                )
else($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
                DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
        )
        install(FILES 
                        ${PY_SRC_LIST}
                DESTINATION 
                        ${Python3_SITELIB}/${PROJECT_NAME}
        )
endif()
