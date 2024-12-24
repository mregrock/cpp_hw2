#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <assert.h>
#include <string_view>
#include <optional>
#include <variant>
#include "src/types/defs.h"
#include "src/simulation/fluid_simulator.hpp"
#include "src/simulation/vector_field.hpp"
#include "src/selector/selector_bootstrap.h"
#include "src/simulation/simulation_helpers.h"
#include "src/types/parse_helpers.h"

int main(int argc, char** argv) {
    std::string p_type = get_arg("--p-type", argc, argv, "FAST_FIXED(32,16)");
    std::string v_type = get_arg("--v-type", argc, argv, "FIXED(31,17)");
    std::string v_flow_type = get_arg("--v-flow-type", argc, argv, "DOUBLE");
    std::string size_str = get_arg("--size", argc, argv, "S(36,84)");
    std::string filename = get_arg("--file", argc, argv, "field.txt");
    
    if (!is_valid_type(p_type)) {
        std::cerr << "Invalid p_type: " << p_type << std::endl;
        return 1;
    }
    if (!is_valid_type(v_type)) {
        std::cerr << "Invalid v_type: " << v_type << std::endl;
        return 1;
    }
    if (!is_valid_type(v_flow_type)) {
        std::cerr << "Invalid v_flow_type: " << v_flow_type << std::endl;
        return 1;
    }
    auto size = parse_size(size_str.c_str());
    
    if (!create_and_run_simulation(p_type, v_type, v_flow_type, size.n, size.m, filename)) {
        std::cerr << "Failed to create simulation with types: " << p_type << ", " << v_type << ", " << v_flow_type << std::endl;
        return 1;
    }
    
    return 0;
}
