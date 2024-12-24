#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "../types/types.h"
#include "vector_field.hpp"

template<
    typename PressureType,
    typename VelocityType,
    typename VFlowType,
    size_t N = dynamic_size,
    size_t M = dynamic_size
>
struct SimulationState {
    PressureType rho[256];
    SimulationMatrix<PressureType, N, M> p{};
    SimulationMatrix<PressureType, N, M> old_p{};

    SimulationMatrix<char, N, M> field{};
    SimulationMatrix<int, N, M> last_use{};
    int UT{0};
    std::mt19937 rnd{1337};
    SimulationMatrix<int, N, M> dirs{};

    VectorField<VelocityType, N, M> velocity{};
    VectorField<VFlowType, N, M> velocity_flow{};

    void resize_all(size_t n, size_t m) {
        field.resize(n);
        for (auto& row : field) {
            row.resize(m + 1);
        }
        p.resize(n);
        for (auto& row : p) {
            row.resize(m);
        }
        old_p.resize(n);
        for (auto& row : old_p) {
            row.resize(m);
        }
        last_use.resize(n);
        for (auto& row : last_use) {
            row.resize(m);
        }
        dirs.resize(n);
        for (auto& row : dirs) {
            row.resize(m);
        }
        velocity.resize(n, m);
        velocity_flow.resize(n, m);
    }

    SimulationState(const std::vector<std::string>& initial_field = {}) {
        if constexpr (N == dynamic_size || M == dynamic_size) {
            if (!initial_field.empty()) {
                size_t n = initial_field.size();
                size_t m = initial_field[0].size();
                resize_all(n, m);
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = 0; j < m; ++j) {
                        field[i][j] = initial_field[i][j];
                    }
                }
            }
        } else {
            if (!initial_field.empty()) {
                if (initial_field.size() != N || initial_field[0].size() != M) {
                    throw std::runtime_error("Field size mismatch");
                }
                for (size_t i = 0; i < N; ++i) {
                    for (size_t j = 0; j < M; ++j) {
                        field[i][j] = initial_field[i][j];
                    }
                }
            }
        }
    }

    size_t get_n() const {
        if constexpr (N == dynamic_size) {
            return field.size();
        } else {
            return N;
        }
    }

    size_t get_m() const {
        if constexpr (M == dynamic_size) {
            return field[0].size() - 1;
        } else {
            return M;
        }
    }

    void save(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open file for writing");
        }
        size_t n, m;
        if constexpr (N == dynamic_size || M == dynamic_size) {
            n = field.size();
            m = field[0].size() - 1;
        } else {
            n = N;
            m = M;
        }
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        out.write(reinterpret_cast<const char*>(&m), sizeof(m));
        
        for (size_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(field[i].data()), m);
        }

        for (size_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(p[i].data()), m);
        }
        for (size_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(old_p[i].data()), m);
        }
        for (size_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(last_use[i].data()), m);
        }
        out.write(reinterpret_cast<const char*>(&UT), sizeof(UT));
        out.write(reinterpret_cast<const char*>(rho), sizeof(rho));
        for (size_t i = 0; i < n; ++i) {
            out.write(reinterpret_cast<const char*>(dirs[i].data()), m);
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                out.write(reinterpret_cast<const char*>(velocity.v[i][j].data()), 4 * sizeof(VelocityType));
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                out.write(reinterpret_cast<const char*>(velocity_flow.v[i][j].data()), 4 * sizeof(VFlowType));
            }
        }
        out.close();
    }

    static SimulationState load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Cannot open file for reading: " + filename);
        }
        size_t n, m;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        in.read(reinterpret_cast<char*>(&m), sizeof(m));

        SimulationState<PressureType, VelocityType, VFlowType, N, M> state;
        state.resize_all(n, m);

        for (size_t i = 0; i < n; ++i) {
            in.read(reinterpret_cast<char*>(state.field[i].data()), m);
        }

        for (size_t i = 0; i < n; ++i) {
            in.read(reinterpret_cast<char*>(state.p[i].data()), m * sizeof(PressureType));
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                in.read(reinterpret_cast<char*>(state.velocity.v[i][j].data()), 
                        4 * sizeof(VelocityType));
            }
        }

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                in.read(reinterpret_cast<char*>(state.velocity_flow.v[i][j].data()), 
                        4 * sizeof(VFlowType));
            }
        }

        if (!in) {
            throw std::runtime_error("Error reading from file: " + filename);
        }
        return state;
    }
};
