cmake_minimum_required(VERSION 3.26)

add_library(ThreadPool INTERFACE)
target_include_directories(ThreadPool INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})