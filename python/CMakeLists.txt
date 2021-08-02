cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)

# gslib
include_directories(${CMAKE_SOURCE_DIR})

# Set installation directories
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; import os;print(get_python_lib())"
                OUTPUT_VARIABLE PY_SP_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(SP_DIR "${PY_SP_DIR}" CACHE PATH "Site-package directory")
message("  -> site-packages: " ${SP_DIR})

# file(GLOB SRC_LIST "pybind_*")

find_package(pybind11 REQUIRED)
pybind11_add_module(
        ${PROJECT_NAME} 
        gbsbind.cpp 
        # ${SRC_LIST}
)
target_link_libraries(${PROJECT_NAME} PRIVATE nlopt ${VTK_LIBRARIES} gbs-render)

vtk_module_autoinit(
        TARGETS  ${PROJECT_NAME} 
        MODULES ${VTK_LIBRARIES}
        )


install(TARGETS ${PROJECT_NAME} 
        DESTINATION ${SP_DIR})
