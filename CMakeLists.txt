cmake_minimum_required (VERSION 2.8)
project (Tutorial)

SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(pimd pimd.cpp)
target_link_libraries(pimd tools fftw3 ${ARMADILLO_LIBRARIES})

include_directories("./tools" "./potentials" ${ARMADILLO_INCLUDE_DIRS})
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
