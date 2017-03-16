cmake_minimum_required(VERSION 3.4)
project(CarbonDiffusion)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.cpp)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
add_library(InfoMatrix src/InfoMatrix.cpp)
add_library(AssembleMatrix src/AssembleMatrix.cpp)
add_library(solver src/CarbonDiffusion.cpp)
target_link_libraries(${PROJECT_NAME} InfoMatrix AssembleMatrix solver)

# Gtest
add_subdirectory(lib/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_executable( runUnitTests test/test_main.cpp)
target_link_libraries(runUnitTests gtest gtest_main)
target_link_libraries(runUnitTests InfoMatrix AssembleMatrix solver)
add_test( runUnitTests runUnitTests )

