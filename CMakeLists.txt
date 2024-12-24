cmake_minimum_required(VERSION 3.10)
project(fluid_new)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)

add_executable(fluid_new ${CMAKE_CURRENT_SOURCE_DIR}/fluid_new.cpp)
target_link_libraries(fluid_new PRIVATE Types Simulation)