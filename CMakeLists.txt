cmake_minimum_required (VERSION 2.8)
project (Tutorial)

#SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-7)

SET(GCC_COVERAGE_COMPILE_FLAGS "-g -O3 --std=c++11")
SET(GCC_COVERAGE_LINK_FLAGS    "-g -O3 --std=c++11")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

file(GLOB tools_sources ./*/*.cpp)
#file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(pimd pimd.cpp)
target_link_libraries(pimd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rpmd rpmd.cpp)
target_link_libraries(rpmd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(ensemble_tcf_rpmd ensemble_tcf_rpmd.cpp)
target_link_libraries(ensemble_tcf_rpmd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(scan_potential scan_potential.cpp)
target_link_libraries(scan_potential tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rdf rdf.cpp)
target_link_libraries(rdf tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rdf-partial-traj rdf-partial-traj.cpp)
target_link_libraries(rdf-partial-traj tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(test-initial-conditions test-initial-conditions.cpp)
target_link_libraries(test-initial-conditions tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(test-randn test-randn.cpp)
target_link_libraries(test-randn tools)

include_directories("./tools" "./potentials" ${ARMADILLO_INCLUDE_DIRS})
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
