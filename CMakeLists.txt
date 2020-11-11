project(bilevel)
cmake_minimum_required(VERSION 3.7)

set(CMAKE_CXX_STANDARD 11)

include(UserConfig.cmake)
set(DEPEND_DIR $ENV{PWD})

message(STATUS "The system information: ${CMAKE_SYSTEM}.${CMAKE_SYSTEM_PROCESSOR}")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

if (NOT ${CPLEX_LIB_DIR} STREQUAL "")
	find_library(CPLEXLIB cplex PATHS ${CPLEX_LIB_DIR})
	if(CPLEXLIB)
		message(STATUS "Found ${CPLEXLIB}")
	else()
		message(FATAL_ERROR "Could NOT find cplex library.")
	endif()
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)
option(UNIT_TESTING "Enable unit testing" OFF)

add_subdirectory(extra)
add_subdirectory(src)

enable_testing()
add_subdirectory(test)