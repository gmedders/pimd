cmake_minimum_required (VERSION 2.8)
set(PROJECT_NAME_STR pimd)
project(${PROJECT_NAME_STR})

SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 --std=c++11")
SET(GCC_COVERAGE_LINK_FLAGS    "-O3 --std=c++11")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES})

file(GLOB tools_sources ./lib/*.cpp)
#file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(pimd ./src/pimd.cpp)
target_link_libraries(pimd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rpmd ./src/rpmd.cpp)
target_link_libraries(rpmd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(classical-md ./src/classical-md.cpp)
target_link_libraries(classical-md tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(ensemble_tcf_rpmd ./src/ensemble_tcf_rpmd.cpp)
target_link_libraries(ensemble_tcf_rpmd tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(scan_potential ./src/scan_potential.cpp)
target_link_libraries(scan_potential tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rdf ./src/rdf.cpp)
target_link_libraries(rdf tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(rdf-partial-traj ./src/rdf-partial-traj.cpp)
target_link_libraries(rdf-partial-traj tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(test-initial-conditions ./src/test-initial-conditions.cpp)
target_link_libraries(test-initial-conditions tools fftw3 ${ARMADILLO_LIBRARIES})

add_executable(test-randn ./src/test-randn.cpp)
target_link_libraries(test-randn tools)

# include_directories("./tools" "./potentials" ${ARMADILLO_INCLUDE_DIRS})
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
