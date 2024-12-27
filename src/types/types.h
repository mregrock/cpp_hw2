#pragma once

#include <vector>
#include <array>
#include <string>
#include <string_view>
#include <optional>
#include <variant>
#include <type_traits>
#include <limits>

constexpr size_t DEFAULT_N = 36, DEFAULT_M = 84;
constexpr size_t dynamic_size = std::numeric_limits<size_t>::max();

template<typename T, size_t N = dynamic_size, size_t M = dynamic_size>
using SimulationMatrix = std::conditional_t<N == dynamic_size || M == dynamic_size,
    std::vector<std::vector<T>>,
    std::array<std::array<T, M + (std::is_same_v<T, char> ? 1 : 0)>, N>
>;

std::vector<std::string> initial_field = {
    "####################################################################################",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                       .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #           .........                                  #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#            #                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............#                                                      #",
    "#..............#............################                     #                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "#...........................#....................................#                 #",
    "##################################################################                 #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "#                                                                                  #",
    "####################################################################################",
};

#ifndef TYPES
#define TYPES FLOAT,FIXED(31,17),FAST_FIXED(25, 11),FIXED(32, 16),DOUBLE,FAST_FIXED(32, 16)
#endif

#ifndef SIZES
#define SIZES /*S(36,84), S(18,42)*/
#endif

struct SizePair {
    size_t n, m;
    constexpr SizePair(size_t n_, size_t m_) : n(n_), m(m_) {}
};

template<size_t N, size_t M>
struct SizeType {
    static constexpr size_t n = N;
    static constexpr size_t m = M;
};
