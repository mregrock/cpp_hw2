set(HEADERS
    fixed.hpp
    fast_fixed.hpp
    fixed_operators.hpp
)

set(SOURCES
    ../fluid_new.cpp
)

add_executable(fluid_new ${SOURCES} ${HEADERS})

target_include_directories(fluid_new PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})