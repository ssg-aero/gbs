# Values of reference of the Catmull-Clark subdivision computed by OpenSubdiv (see reference.cpp).
# Not part of the build of gbs.
cmake_minimum_required(VERSION 3.20)
project(subd_reference LANGUAGES CXX)

find_package(TBB REQUIRED) # OpenSubdivConfig.cmake (conda-forge) links TBB::tbb without looking for it
find_package(OpenSubdiv REQUIRED)
find_package(RapidJSON REQUIRED)

add_executable(subd_reference reference.cpp)
target_compile_features(subd_reference PRIVATE cxx_std_20)
target_link_libraries(subd_reference PRIVATE OpenSubdiv::osdCPU)
target_include_directories(subd_reference PRIVATE ${RapidJSON_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})
set_target_properties(subd_reference PROPERTIES BUILD_RPATH "${CMAKE_PREFIX_PATH}/lib")
