cmake_minimum_required(VERSION 2.8.8)
set(PROJECT_NAME_STR pimd)
project(${PROJECT_NAME_STR})

#set( CMAKE_VERBOSE_MAKEFILE on )

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 --std=c++11")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3 --std=c++11")
endif()

#-------------------
# set common include folder for module
#-------------------
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES})
set(TCLAP_INCLUDES ${PROJECT_SOURCE_DIR}/include/tclap)
include_directories(${TCLAP_INCLUDES})

file(GLOB tools_sources ./lib/*.cpp)
set( tools ${tools_sources})
add_library( tools ${tools})

find_package(Armadillo REQUIRED)

#-------------------
# Main source files
#-------------------

file(GLOB MAIN_SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Set up each executable target in the src directory
foreach(_main_file ${MAIN_SRC_FILES})
    get_filename_component(_main_name ${_main_file} NAME_WE)
    add_executable(${_main_name} ${_main_file})
    target_link_libraries(${_main_name} tools fftw3 ${ARMADILLO_LIBRARIES})
endforeach()

#-------------------
# Download and build the googletest external dependency
#-------------------

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
get_filename_component(DEPS_ROOT "${PROJECT_BINARY_DIR}/deps" ABSOLUTE)

include(ExtProjectUtils)

find_package(Threads REQUIRED)

# Will download external CMakeable project from git repo release-1.8.0 and install it in $DEPS_ROOT
# This also will create "googletest.git" target, which we'll use as dependency for our test project
ExtProjectGit("https://github.com/google/googletest.git" "release-1.8.0" ${DEPS_ROOT} CMAKE_ARGS "${_OPT_CMAKE_ARGS}")

include_directories("${DEPS_ROOT}/include")
link_directories("${DEPS_ROOT}/lib")

#-------------------
# Test
#-------------------

enable_testing()

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)

# from list of files we'll create tests test_name.cpp -> test_name
foreach(_test_file ${TEST_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    add_dependencies(${_test_name} "googletest.git")
    target_link_libraries(${_test_name} tools gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} fftw3 ${ARMADILLO_LIBRARIES})
    add_test(${_test_name} ${_test_name})
    set_tests_properties(${_test_name} PROPERTIES TIMEOUT 5)
endforeach()
