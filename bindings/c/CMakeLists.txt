cmake_minimum_required(VERSION 3.10)
project(gelfgren_c_examples C)

set(CMAKE_C_STANDARD 11)

# Find the Gelfgren library
set(GELFGREN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
set(GELFGREN_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../target/release")

# Add include directory
include_directories(${GELFGREN_INCLUDE_DIR})

# Link directories
link_directories(${GELFGREN_LIB_DIR})

# Basic usage example
add_executable(basic_usage basic_usage.c)
target_link_libraries(basic_usage gelfgren m pthread dl)

# Instructions
message(STATUS "To build the C examples:")
message(STATUS "  1. First build the Rust library: cargo build --release")
message(STATUS "  2. Then build this project: mkdir build && cd build && cmake .. && make")
message(STATUS "  3. Run: ./basic_usage")
