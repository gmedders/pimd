cmake_minimum_required (VERSION 2.8)
project (Tutorial)

SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

SET(GCC_COVERAGE_COMPILE_FLAGS "-g")
SET(GCC_COVERAGE_LINK_FLAGS    "-g")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(pimd pimd.cpp)
add_executable(rpmd rpmd.cpp)
target_link_libraries(pimd tools fftw3 ${ARMADILLO_LIBRARIES})
target_link_libraries(rpmd tools ${ARMADILLO_LIBRARIES})

include_directories("./tools" "./potentials" ${ARMADILLO_INCLUDE_DIRS})
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
