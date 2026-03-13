cmake_minimum_required(VERSION 3.15)
project(QuantLibResultsAnalysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Locate the QuantLib build via its exported CMake config.
# The parent repo's build/ directory contains cmake/QuantLibConfig.cmake.
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../build/cmake")

find_package(QuantLib CONFIG REQUIRED)

message(STATUS "Found QuantLib: ${QuantLib_DIR}")

add_executable(generate_data generate_data.cpp)
target_link_libraries(generate_data PRIVATE QuantLib::QuantLib)

# Ensure the data and figures output directories exist at build time.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/figures")
