cmake_minimum_required(VERSION 3.10)
project(gelfgren_cpp_bindings CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the Gelfgren library
set(GELFGREN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../include")
set(GELFGREN_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../target/release")

# Add include directories
include_directories(${GELFGREN_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Link directories
link_directories(${GELFGREN_LIB_DIR})

# C++ example
add_executable(cpp_example example.cpp)
target_link_libraries(cpp_example gelfgren m pthread dl)

# Enable warnings
target_compile_options(cpp_example PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3>
)

# Instructions
message(STATUS "")
message(STATUS "Gelfgren C++ Bindings")
message(STATUS "=====================")
message(STATUS "")
message(STATUS "To build the C++ example:")
message(STATUS "  1. First build the Rust library: cargo build --release")
message(STATUS "  2. Then build this project:")
message(STATUS "     mkdir build && cd build")
message(STATUS "     cmake .. && make")
message(STATUS "  3. Run: ./cpp_example")
message(STATUS "")
