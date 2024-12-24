#pragma once

#include "defs.h"

constexpr SizePair parse_size(const char* s) {
    s += 2;
    size_t n = 0;
    while (*s >= '0' && *s <= '9') {
        n = n * 10 + (*s - '0');
        s++;
    }
    s++;
    size_t m = 0;
    while (*s >= '0' && *s <= '9') {
        m = m * 10 + (*s - '0');
        s++;
    }
    return SizePair(n, m);
}

template<typename T>
std::string get_pretty_type_name() {
    if constexpr (std::is_same_v<T, float>) {
        return "float";
    } else if constexpr (std::is_same_v<T, double>) {
        return "double";
    } else if constexpr (is_fixed_v<T>) {
        return "Fixed<" + std::to_string(T::bits) + "," + std::to_string(T::frac) + ">";
    } else if constexpr (is_fast_fixed_v<T>) {
        return "FastFixed<" + std::to_string(T::bits) + "," + std::to_string(T::frac) + ">";
    } else {
        return "unknown";
    }
}

std::string get_arg(std::string_view arg_name, int argc, char** argv, std::string_view default_value) {
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == arg_name) {
            return argv[i + 1];
        }
    }
    return std::string(default_value);
}

bool is_valid_type(const std::string& type) {
    if (type == "FLOAT" || type == "DOUBLE") return true;
    
    if (type.starts_with("FIXED(") || type.starts_with("FAST_FIXED(")) {
        size_t pos = type.find(',');
        if (pos == std::string::npos) return false;
        
        size_t end_pos = type.find(')', pos);
        if (end_pos == std::string::npos) return false;
        
        return true;
    }
    
    return false;
}