cmake_minimum_required (VERSION 2.8)
project (Tutorial)

file(GLOB tools_sources ./*.cpp ./*/*.cpp)
set( tools ${tools_sources})

add_library( tools ${tools})

add_executable(pimd pimd.cpp)
target_link_libraries(pimd tools)

include_directories("./tools")
#add_subdirectory(tools)
#
#target_link_libraries(pimd tools)
#
