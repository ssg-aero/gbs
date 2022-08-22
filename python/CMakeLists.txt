cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)

# gslib
include_directories(${CMAKE_SOURCE_DIR})

find_package (Python3 COMPONENTS Interpreter Development)

find_package(pybind11 REQUIRED)
pybind11_add_module(
        ${PROJECT_NAME}
        gbsbindrender.cpp 
        gbsbindmesh.cpp
        gbsbindinterp.cpp
        gbsbind.cpp 
        # ${SRC_LIST}
)
target_link_libraries(
        ${PROJECT_NAME} PRIVATE 
        NLopt::nlopt 
        ${VTK_LIBRARIES} 
        gbs-render 
        TBB::tbb 
)

vtk_module_autoinit(
        TARGETS  ${PROJECT_NAME} 
        MODULES ${VTK_LIBRARIES}
        )

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "gbs")

file(GLOB SRC_LIST "*.py")

if($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
        DESTINATION $ENV{SP_DIR}/${PROJECT_NAME}
        )
        install(FILES 
                        # vistaplot.py
                        ${SRC_LIST}
                DESTINATION 
                        $ENV{SP_DIR}/${PROJECT_NAME}
                )
else($ENV{CONDA_BUILD})
        install(TARGETS ${PROJECT_NAME} 
                DESTINATION ${Python3_SITELIB}/${PROJECT_NAME}
        )
        install(FILES 
                        # vistaplot.py
                        ${SRC_LIST}
                DESTINATION 
                        ${Python3_SITELIB}/${PROJECT_NAME}
        )
endif()

