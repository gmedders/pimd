cmake_minimum_required (VERSION 2.8)
project (Tutorial)

SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-6)

file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

add_executable(pimd pimd.cpp)
target_link_libraries(pimd tools fftw3)

include_directories("./tools" "./potentials")
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
