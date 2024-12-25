add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE Types Simulation Selector ThreadPool)
