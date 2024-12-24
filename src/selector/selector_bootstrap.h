#pragma once

#include "selector_lists.hpp"
#include "selector_helpers.h"
#include "selector.hpp"

template<typename... Types>
bool try_all_type_combinations(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                             size_t n, size_t m, const std::string& filename) {
    #define S(N, M) SizeType<N, M>
    return TypeSelector<TypesList<Types...>, SizesList<SIZES>, TypesList<>>::try_combinations(p_type, v_type, v_flow_type, n, m, filename);
}

bool create_and_run_simulation(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type, 
                             size_t n, size_t m, const std::string& filename) {
    try {
        std::cerr << "\nTrying to create simulation with types:" << std::endl;
        std::cerr << "p_type: " << p_type << std::endl;
        std::cerr << "v_type: " << v_type << std::endl;
        std::cerr << "v_flow_type: " << v_flow_type << std::endl;
        #define FLOAT float
        #define DOUBLE double
        #define FIXED(N, K) Fixed<N, K>
        #define FAST_FIXED(N, K) FastFixed<N, K>
        if (!try_all_type_combinations<TYPES>(p_type, v_type, v_flow_type, n, m, filename)) {
            std::cerr << "Error: No matching type combination found" << std::endl;
            return false;
        }
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }
}