#pragma once

#include "../types/types.h"
#include "vector_field.hpp"
#include "simulation_state.hpp"
#include "fixed_operators.hpp"
#include <utility>
#include <array>
#include <ranges>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <cassert>
#include <string_view>
#include <optional>
#include <variant>
#include <iostream>
#include "../thread_pool/thread_pool.hpp"

template<typename PressureType, typename VelocityType, typename VFlowType, size_t N = DEFAULT_N, size_t M = DEFAULT_M>
class FluidSimulator {
private:
    size_t runtime_n, runtime_m;

    static constexpr int autosave_interval = 10;

    static constexpr size_t T = 1'000'000;
    static constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    
    SimulationMatrix<char, N, M> field{};
    SimulationMatrix<PressureType, N, M> p{}, old_p{};
    SimulationMatrix<int, N, M> last_use{};
    int UT = 0;
    std::mt19937 rnd{1337};
    PressureType rho[256];

    VectorField<VelocityType, N, M> velocity{};
    VectorField<VFlowType, N, M> velocity_flow{};
    SimulationMatrix<int, N, M> dirs{};

    ThreadPool thread_pool;

    struct ParticleParams {
        char type;
        PressureType cur_p;
        std::array<VelocityType, 4> v;

        void swap_with(FluidSimulator* sim, int x, int y) {
            std::swap(sim->field[x][y], type);
            std::swap(sim->p[x][y], cur_p);
            std::swap(sim->velocity.v[x][y], v);
        }
    };

    template<typename T>
    PressureType to_pressure(T value) {
        if constexpr (std::is_same_v<T, PressureType>) {
            return value;
        } else {
            return PressureType(value);
        }
    }

    template<typename T>
    VelocityType to_velocity(T value) {
        if constexpr (std::is_same_v<T, VelocityType>) {
            return value;
        } else {
            return VelocityType(value);
        }
    }

    template<typename T>
    VFlowType to_flow(T value) {
        if constexpr (std::is_same_v<T, VFlowType>) {
            return value;
        } else {
            return VFlowType(value);
        }
    }

    size_t get_n() const {
        if constexpr (N == dynamic_size) {
            return runtime_n;
        } else {
            return N;
        }
    }

    size_t get_m() const {
        if constexpr (M == dynamic_size) {
            return runtime_m;
        } else {
            return M;
        }
    }

    template<typename Func>
    void parallel_for(size_t start, size_t end, size_t chunk_size, Func f) {
        std::vector<std::future<void>> futures;
        for (size_t i = start; i < end; i += chunk_size) {
            size_t local_end = std::min(i + chunk_size, end);
            futures.push_back(
                thread_pool.enqueue([=]() {
                    for (size_t j = i; j < local_end; ++j) {
                        f(j);
                    }
                })
            );
        }
        for (auto& future : futures) {
            future.get();
        }
    }

    void init_state(const SimulationState<PressureType, VelocityType, VFlowType, N, M>& state) {
        field = state.field;
        p = state.p;
        old_p = state.old_p;
        last_use = state.last_use;
        UT = state.UT;
        dirs = state.dirs;
        velocity = state.velocity;
        velocity_flow = state.velocity_flow;
    }

    std::tuple<VelocityType, bool, std::pair<int, int>> propagate_flow(int x, int y, VelocityType lim) {
        last_use[x][y] = UT - 1;
        VelocityType ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) continue;
                
                auto vp = std::min(lim, to_velocity(cap - flow));
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, to_flow(vp));
                    last_use[x][y] = UT;
                    return {vp, true, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, to_flow(t));
                    last_use[x][y] = UT;
                    return {t, prop && end != std::pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, false, {0, 0}};
    }

    VelocityType random01() {
        return VelocityType(static_cast<double>(rnd() & ((1 << 16) - 1)) / (1 << 16));
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > VelocityType(0)) {
                    stop = false;
                    break;
                }
            }
            if (!stop) return;
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > VelocityType(0)) continue;
            propagate_stop(nx, ny);
        }
    }

    VelocityType move_prob(int x, int y) {
        VelocityType sum = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) continue;
            auto v = velocity.get(x, y, dx, dy);
            if (v < VelocityType(0)) continue;
            sum += v;
        }
        return sum;
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<VelocityType, 4> tres;
            VelocityType sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < VelocityType(0)) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == VelocityType(0)) break;

            VelocityType p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > VelocityType(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);

        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < VelocityType(0)) {
                propagate_stop(nx, ny);
            }
        }
        if (ret && !is_first) {
            ParticleParams pp{};
            pp.swap_with(this, x, y);
            pp.swap_with(this, nx, ny);
            pp.swap_with(this, x, y);
        }
        return ret;
    }

    SimulationState<PressureType, VelocityType, VFlowType, N, M> create_state() const {
        SimulationState<PressureType, VelocityType, VFlowType, N, M> state;
        state.field = field;
        state.p = p;
        state.old_p = old_p;
        state.last_use = last_use;
        state.UT = UT;
        state.dirs = dirs;
        state.velocity = velocity;
        state.velocity_flow = velocity_flow;
        return state;
    }

public:
    FluidSimulator(const SimulationState<PressureType, VelocityType, VFlowType, N, M>& state) {
        field = state.field;
        rho[' '] = PressureType(0.01);
        rho['.'] = PressureType(1000);
        if constexpr (N == dynamic_size || M == dynamic_size) {
            runtime_n = field.size();
            runtime_m = field[0].size() - 1;
            init_state(state);
        }
        else {
            runtime_n = N;
            runtime_m = M;
        }
        for (size_t x = 0; x < get_n(); ++x) {
            for (size_t y = 0; y < get_m(); ++y) {
                if (field[x][y] == '#') continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
    }

    void run(bool benchmark = false, size_t ticks = 1000) {
        PressureType g = PressureType(0.1);
        for (size_t x = 0; x < get_n(); ++x) {
            field[x][get_m()] = '\0';
            if (!benchmark) {
                std::cout << field[x].data() << "\n";
            }
        }
        for (size_t i = 0; i < T; ++i) {
            if (!benchmark && i % autosave_interval == 0) {
                auto state = create_state();
                state.save("autosave.bin");
            }
            if (i == ticks) break;
            PressureType total_delta_p = 0;
            // Apply external forces
            size_t chunk_size = get_n() / thread_pool.get_thread_count();
            parallel_for(0, get_n(), chunk_size, [&](size_t x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] == '#') continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, VelocityType(g));
                }
            });

            // Apply forces from p
            old_p = p;
            parallel_for(0, get_n(), chunk_size, [&](size_t x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = to_pressure(delta_p);
                            auto& contr = velocity.get(nx, ny, -dx, -dy);
                            if (to_pressure(contr) * rho[(int)field[nx][ny]] >= force) {
                                contr -= to_velocity(force / rho[(int)field[nx][ny]]);
                                continue;
                            }
                            force -= to_pressure(contr) * rho[(int)field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, to_velocity(force / rho[(int)field[x][y]]));
                            std::atomic<PressureType>& atomic_p = reinterpret_cast<std::atomic<PressureType>&>(p[x][y]);
                            atomic_p.store(atomic_p.load() - force / dirs[x][y]);
                            
                            std::atomic<PressureType>& atomic_total = reinterpret_cast<std::atomic<PressureType>&>(total_delta_p);
                            atomic_total.store(atomic_total.load() - force / dirs[x][y]);
                        }
                    }
                }
            });
            // Make flow from velocities
            if constexpr (N == dynamic_size || M == dynamic_size) {
                velocity_flow.resize(runtime_n, runtime_m);
            }
            else {
                velocity_flow = {};
            }
            bool prop = false;
            do {
                UT += 2;
                prop = false;
                for (size_t x = 0; x < get_n(); ++x) {
                    for (size_t y = 0; y < get_m(); ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, VelocityType(1));
                            if (t > VelocityType(0)) prop = true;
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            parallel_for(0, get_n(), chunk_size, [&](size_t x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > VelocityType(0)) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = to_velocity(new_v);
                            auto force = to_pressure(old_v - new_v) * rho[(int)field[x][y]];
                            if (field[x][y] == '.') force *= PressureType(0.8);
                            if (field[x + dx][y + dy] == '#') {
                                std::atomic<PressureType>& atomic_p = reinterpret_cast<std::atomic<PressureType>&>(p[x][y]);
                                atomic_p.store(atomic_p.load() + force / dirs[x][y]);
                                std::atomic<PressureType>& atomic_total = reinterpret_cast<std::atomic<PressureType>&>(total_delta_p);
                                atomic_total.store(atomic_total.load() + force / dirs[x][y]);
                            } else {
                                std::atomic<PressureType>& atomic_p = reinterpret_cast<std::atomic<PressureType>&>(p[x + dx][y + dy]);
                                atomic_p.store(atomic_p.load() + force / dirs[x + dx][y + dy]);
                                std::atomic<PressureType>& atomic_total = reinterpret_cast<std::atomic<PressureType>&>(total_delta_p);
                                atomic_total.store(atomic_total.load() + force / dirs[x + dx][y + dy]);
                            }
                        }
                    }
                }
            });

            UT += 2;
            prop = false;
            for (size_t x = 0; x < get_n(); ++x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop && !benchmark) {
                std::cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < get_n(); ++x) {
                    field[x][get_m()] = '\0';
                    std::cout << field[x].data() << "\n";
                }
            }
        }
    }
};
