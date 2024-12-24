cmake_minimum_required(VERSION 3.26)

add_library(Simulation INTERFACE)
target_include_directories(Simulation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})