cmake_minimum_required(VERSION 3.17.0)

project(pygbs LANGUAGES CXX)
#
# gslib
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
# site-package dir
if(EXISTS $ENV{SP_DIR})
    file(TO_CMAKE_PATH $ENV{SP_DIR} SP_DIR)
    set(${PROJECT_NAME}_INSTALL_DIRECTORY ${SP_DIR}/${PROJECT_NAME})
else()
    if(UNIX)
        set(${PROJECT_NAME}_INSTALL_DIRECTORY ${CONDA_PREFIX}/lib/python${python_version_major}.${python_version_minor}/site-packages/${PROJECT_NAME})
    else()
        set(${PROJECT_NAME}_INSTALL_DIRECTORY ${CONDA_PREFIX}/Lib/site-packages/${PROJECT_NAME})
    endif()
endif()
find_package(pybind11 REQUIRED)
pybind11_add_module(${PROJECT_NAME} gbsbind.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE nlopt)

if ( MSVC )
    SET_OUT_DIR(${PROJECT_NAME} ${CMAKE_BINARY_DIR})
endif ( MSVC )
install(TARGETS ${PROJECT_NAME} 
    COMPONENT python
    LIBRARY  DESTINATION ${LIBRARY_LIB}
    RUNTIME  DESTINATION ${LIBRARY_BIN}
    INCLUDES DESTINATION ${LIBRARY_INC}
    ARCHIVE  DESTINATION ${LIBRARY_LIB}
)