#pragma once

#include "simulation_state.hpp"
#include "fluid_simulator.hpp"
#include "../selector/selector_helpers.h"

bool is_state_file(const std::string& filename) {
    return filename.ends_with(".bin");
}

std::vector<std::string> read_field(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Can't open file: " + filename);
    }
    std::vector<std::string> field;
    std::string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (!field.empty() && line.length() != field[0].length()) {
            throw std::runtime_error("Inconsistent line length in field");
        }
        field.push_back(line);
    }
    return field;
}

template<typename P, typename V, typename VF, typename Size=SizeType<dynamic_size, dynamic_size>>
void run_simulation(size_t n = DEFAULT_N, size_t m = DEFAULT_M, const std::string& filename="field.txt") {  
    std::cerr << "Running simulation with n=" << n << ", m=" << m << ", filename=" << filename << std::endl;
    try {
        if (is_state_file(filename)) {
            std::cerr << "Loading full state " << filename << std::endl;
            auto state = SimulationState<P, V, VF, Size::n, Size::m>::load(filename);
            FluidSimulator<P, V, VF, Size::n, Size::m> simulator(state);
            simulator.run();
        }
        else {
            std::cerr << "Loading initial field " << filename << std::endl;
            auto initial_field = read_field(filename);
            SimulationState<P, V, VF, Size::n, Size::m> state(initial_field);
            FluidSimulator<P, V, VF, Size::n, Size::m> simulator(state);
            simulator.run();
        }
        std::cerr << "Simulation completed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception in run_simulation: " << e.what() << std::endl;
        throw;
    }
}