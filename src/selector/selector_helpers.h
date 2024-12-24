#pragma once

#include <string>
#include <type_traits>

#include "../types/defs.h"
#include "../selector/selector_lists.hpp"

template<typename T>
struct is_fixed : std::false_type {};

template<size_t N, size_t K>
struct is_fixed<Fixed<N,K>> : std::true_type {};

template<typename T>
inline constexpr bool is_fixed_v = is_fixed<T>::value;

template<typename T>
struct is_fast_fixed : std::false_type {};

template<size_t N, size_t K>
struct is_fast_fixed<FastFixed<N,K>> : std::true_type {};

template<typename T>
inline constexpr bool is_fast_fixed_v = is_fast_fixed<T>::value;

std::pair<size_t, size_t> parse_fixed_params(const std::string& type) {
    size_t start = type.find('(') + 1;
    size_t comma = type.find(',', start);
    size_t end = type.find(')', comma);
    
    size_t N = std::stoul(type.substr(start, comma - start));
    size_t K = std::stoul(type.substr(comma + 1, end - comma - 1));
    return {N, K};
}

template<typename T>
static bool matches_type(const std::string& type) {
    if constexpr (std::is_same_v<T, float>) {
        return type == "FLOAT";
    } else if constexpr (std::is_same_v<T, double>) {
        return type == "DOUBLE";
    } else if constexpr (is_fixed_v<T>) {
        if (!type.starts_with("FIXED(")) return false;
        auto [bits, frac] = parse_fixed_params(type);
        return bits == T::bits && frac == T::frac;
    } else if constexpr (is_fast_fixed_v<T>) {
        if (!type.starts_with("FAST_FIXED(")) return false;
        auto [bits, frac] = parse_fixed_params(type);
        return bits == T::bits && frac == T::frac;
    }
    return false;
}