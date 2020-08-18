cmake_minimum_required(VERSION 3.17.0)

project(gbslib_tests LANGUAGES CXX)
# c++ conf
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# gslib
message(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)

# Locate GTest
enable_testing()
set(GTEST_ROOT "C:/Program Files/googletest-distribution")
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
#
# occt
# 
set(HOME_DIR $ENV{USERPROFILE})
file(TO_CMAKE_PATH ${HOME_DIR} HOME_DIR)
if(EXISTS ${HOME_DIR}/Anaconda3/envs/dev-occt-7.4.0/Library/lib/cmake/opencascade)
    set(OpenCASCADE_DIR ${HOME_DIR}/Anaconda3/envs/dev-occt-7.4.0/Library/lib/cmake/opencascade)
else()
    set(OpenCASCADE_DIR ${HOME_DIR}/Miniconda3/envs/dev-occt-7.4.0/Library/lib/cmake/opencascade)
endif()
# message( ${HOME_DIR})
# message(${OpenCASCADE_DIR})
# set(OpenCASCADE_DIR "C:/Users/sebastien/workspace/3rdParty/opencascade-7.4.0/cmake")
find_package(OpenCASCADE REQUIRED)
link_directories(${OpenCASCADE_LIBRARY_DIR})
include_directories(${OpenCASCADE_INCLUDE_DIR})
add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING) # occt interators seems to be old

# set(CMAKE_PREFIX_PATH "C:/Users/sebastien/workspace/3rdParty/eigen/cmake")
# find_package (Eigen3 REQUIRED NO_MODULE)
include_directories("C:/Users/sebastien/workspace/3rdParty/eigen")
# include_directories("C:/Users/sebastien/Anaconda3/envs/dev-occt-7.4.0/Library/include/eigen3")

link_directories("C:/Users/sebastien/workspace2/occt-utils/build")
include_directories("C:/Users/sebastien/workspace2/occt-utils")

file(GLOB SRC_LIST "tests_*")
message(${SRC_LIST})
add_executable(tests ${SRC_LIST})

target_link_libraries(tests 
                            GTest::GTest 
                            GTest::Main 
                            ${OpenCASCADE_LIBRARIES} 
                            occt-utils 
                    )
