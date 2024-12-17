#include <iostream>
#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <assert.h>

using namespace std;

constexpr size_t DEFAULT_N = 36, DEFAULT_M = 84;
constexpr char initial_field[DEFAULT_N][DEFAULT_M + 1] = {
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

template<size_t N, size_t K>
struct Fixed {
    static_assert(N > K, "N must be greater than K");
    static_assert(N <= 64, "N must be less than or equal to 64");
    
    using StorageType = typename std::conditional<N <= 32, int32_t, int64_t>::type;
    
    constexpr Fixed(int v = 0): v(static_cast<StorageType>(v) << K) {}
    constexpr Fixed(float f): v(f * (StorageType(1) << K)) {}
    constexpr Fixed(double f): v(f * (StorageType(1) << K)) {}

    static constexpr Fixed from_raw(StorageType x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    StorageType v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;

    explicit operator double() const {
        return v / static_cast<double>(StorageType(1) << K);
    }

    friend Fixed operator/(Fixed a, int b) {
        return Fixed::from_raw(a.v / b);
    }

    friend Fixed operator*(Fixed a, int b) {
        return Fixed::from_raw(a.v * b);
    }

    friend Fixed operator*(int a, Fixed b) {
        return b * a;
    }
};

template<size_t N, size_t K>
Fixed<N,K> operator+(Fixed<N,K> a, Fixed<N,K> b) {
    return Fixed<N,K>::from_raw(a.v + b.v);
}

template<size_t N, size_t K>
Fixed<N,K> operator-(Fixed<N,K> a, Fixed<N,K> b) {
    return Fixed<N,K>::from_raw(a.v - b.v);
}

template<size_t N, size_t K>
Fixed<N,K> operator*(Fixed<N,K> a, Fixed<N,K> b) {
    using ST = typename Fixed<N,K>::StorageType;
    if constexpr (N <= 32) {
        return Fixed<N,K>::from_raw((static_cast<int64_t>(a.v) * b.v) >> K);
    } else {
        ST high = (a.v >> K) * b.v;
        ST low = (a.v & ((ST(1) << K) - 1)) * b.v >> K;
        return Fixed<N,K>::from_raw(high + low);
    }
}

template<size_t N, size_t K>
Fixed<N,K> operator/(Fixed<N,K> a, Fixed<N,K> b) {
    using ST = typename Fixed<N,K>::StorageType;
    if constexpr (N <= 32) {
        return Fixed<N,K>::from_raw((static_cast<int64_t>(a.v) << K) / b.v);
    } else {
        return Fixed<N,K>::from_raw((a.v << K) / b.v);
    }
}

template<size_t N, size_t K>
Fixed<N,K> &operator+=(Fixed<N,K> &a, Fixed<N,K> b) { return a = a + b; }

template<size_t N, size_t K>
Fixed<N,K> &operator-=(Fixed<N,K> &a, Fixed<N,K> b) { return a = a - b; }

template<size_t N, size_t K>
Fixed<N,K> &operator*=(Fixed<N,K> &a, Fixed<N,K> b) { return a = a * b; }

template<size_t N, size_t K>
Fixed<N,K> &operator/=(Fixed<N,K> &a, Fixed<N,K> b) { return a = a / b; }

template<size_t N, size_t K>
Fixed<N,K> operator-(Fixed<N,K> x) { return Fixed<N,K>::from_raw(-x.v); }

template<size_t N, size_t K>
Fixed<N,K> abs(Fixed<N,K> x) {
    Fixed<N,K> ret = x;
    if (ret.v < 0) ret.v = -ret.v;
    return ret;
}

template<size_t N, size_t K>
ostream &operator<<(ostream &out, Fixed<N,K> x) {
    return out << static_cast<double>(x);
}

template<typename NumericType, size_t N, size_t M>
struct VectorField {
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    array<NumericType, 4> v[N][M];

    NumericType &add(int x, int y, int dx, int dy, NumericType dv) {
        return get(x, y, dx, dy) += dv;
    }

    NumericType &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

template<typename NumericType, size_t N, size_t M>
struct SimulationState {
    NumericType rho[256];
    NumericType p[N][M]{};
    NumericType old_p[N][M]{};
    char field[N][M + 1];
    int last_use[N][M]{};
    int UT{0};
    mt19937 rnd{1337};
    int dirs[N][M]{};
    VectorField<NumericType, N, M> velocity{};
    VectorField<NumericType, N, M> velocity_flow{};

    SimulationState() {
        memcpy(field, initial_field, sizeof(field));
    }
};

template<
    typename NumericType = Fixed<32, 16>,
    size_t SizeN = DEFAULT_N,
    size_t SizeM = DEFAULT_M
>
class FluidSimulator {
private:
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

    SimulationState<NumericType, SizeN, SizeM>& state;

    struct ParticleParams {
        char type;
        NumericType cur_p;
        array<NumericType, 4> v;
        char (&field)[SizeN][SizeM + 1];
        NumericType (&p)[SizeN][SizeM];
        VectorField<NumericType, SizeN, SizeM>& velocity;

        ParticleParams(
            char (&f)[SizeN][SizeM + 1],
            NumericType (&p_arr)[SizeN][SizeM],
            VectorField<NumericType, SizeN, SizeM>& vel
        ) : field(f), p(p_arr), velocity(vel) {}

        void swap_with(int x, int y) {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };

    using FlowResult = tuple<NumericType, bool, pair<int, int>>;
    
    FlowResult do_propagate_flow(int x, int y, NumericType lim) {
        state.last_use[x][y] = state.UT - 1;
        NumericType ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (state.field[nx][ny] != '#' && state.last_use[nx][ny] < state.UT) {
                auto cap = state.velocity.get(x, y, dx, dy);
                auto flow = state.velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                auto vp = min(lim, cap - flow);
                if (state.last_use[nx][ny] == state.UT - 1) {
                    state.velocity_flow.add(x, y, dx, dy, vp);
                    state.last_use[x][y] = state.UT;
                    return FlowResult{vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = do_propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    state.velocity_flow.add(x, y, dx, dy, t);
                    state.last_use[x][y] = state.UT;
                    return FlowResult{t, prop && end != pair(x, y), end};
                }
            }
        }
        state.last_use[x][y] = state.UT;
        return FlowResult{ret, 0, {0, 0}};
    }

public:
    explicit FluidSimulator(SimulationState<NumericType, SizeN, SizeM>& simulation_state)
        : state(simulation_state)
    {
        initialize();
    }

    void initialize() {
        state.rho[' '] = NumericType(0.01);
        state.rho['.'] = NumericType(1000);

        for (size_t x = 0; x < SizeN; ++x) {
            for (size_t y = 0; y < SizeM; ++y) {
                if (state.field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    state.dirs[x][y] += (state.field[x + dx][y + dy] != '#');
                }
            }
        }
    }

    NumericType get_random01() {
        return NumericType::from_raw((state.rnd() & ((1 << 16) - 1)));
    }

    NumericType get_move_prob(int x, int y) {
        NumericType sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (state.field[nx][ny] == '#' || state.last_use[nx][ny] == state.UT) {
                continue;
            }
            auto v = state.velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    bool do_propagate_move(int x, int y, bool is_first) {
        state.last_use[x][y] = state.UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<NumericType, deltas.size()> tres;
            NumericType sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (state.field[nx][ny] == '#' || state.last_use[nx][ny] == state.UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = state.velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            NumericType p = get_random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(state.velocity.get(x, y, dx, dy) > 0 && state.field[nx][ny] != '#' && state.last_use[nx][ny] < state.UT);

            ret = (state.last_use[nx][ny] == state.UT - 1 || do_propagate_move(nx, ny, false));
        } while (!ret);
        
        state.last_use[x][y] = state.UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (state.field[nx][ny] != '#' && state.last_use[nx][ny] < state.UT - 1 && state.velocity.get(x, y, dx, dy) < 0) {
                do_propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp(state.field, state.p, state.velocity);
                pp.swap_with(x, y);
                pp.swap_with(nx, ny);
                pp.swap_with(x, y);
            }
        }
        return ret;
    }

    void do_propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (state.field[nx][ny] != '#' && state.last_use[nx][ny] < state.UT - 1 && state.velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        state.last_use[x][y] = state.UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (state.field[nx][ny] == '#' || state.last_use[nx][ny] == state.UT || state.velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            do_propagate_stop(nx, ny);
        }
    }

    void simulate_step(NumericType g, size_t iteration) {
        NumericType total_delta_p = 0;
        
        // Apply external forces
        for (size_t x = 0; x < SizeN; ++x) {
            for (size_t y = 0; y < SizeM; ++y) {
                if (state.field[x][y] == '#')
                    continue;
                if (state.field[x + 1][y] != '#')
                    state.velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(state.old_p, state.p, sizeof(state.p));
        for (size_t x = 0; x < SizeN; ++x) {
            for (size_t y = 0; y < SizeM; ++y) {
                if (state.field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (state.field[nx][ny] != '#' && state.old_p[nx][ny] < state.old_p[x][y]) {
                        auto delta_p = state.old_p[x][y] - state.old_p[nx][ny];
                        auto force = delta_p;
                        auto &contr = state.velocity.get(nx, ny, -dx, -dy);
                        if (contr * state.rho[(int) state.field[nx][ny]] >= force) {
                            contr -= force / state.rho[(int) state.field[nx][ny]];
                            continue;
                        }
                        force -= contr * state.rho[(int) state.field[nx][ny]];
                        contr = 0;
                        state.velocity.add(x, y, dx, dy, force / state.rho[(int) state.field[x][y]]);
                        state.p[x][y] -= force / NumericType(state.dirs[x][y]);
                        total_delta_p -= force / NumericType(state.dirs[x][y]);
                    }
                }
            }
        }

        // Make flow from velocities
        state.velocity_flow = {};
        bool prop = false;
        do {
            state.UT += 2;
            prop = 0;
            for (size_t x = 0; x < SizeN; ++x) {
                for (size_t y = 0; y < SizeM; ++y) {
                    if (state.field[x][y] != '#' && state.last_use[x][y] != state.UT) {
                        auto [t, local_prop, _] = do_propagate_flow(x, y, 1);
                        if (t > 0) {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < SizeN; ++x) {
            for (size_t y = 0; y < SizeM; ++y) {
                if (state.field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    auto old_v = state.velocity.get(x, y, dx, dy);
                    auto new_v = state.velocity_flow.get(x, y, dx, dy);
                    if (old_v > 0) {
                        assert(new_v <= old_v);
                        state.velocity.get(x, y, dx, dy) = new_v;
                        auto force = (old_v - new_v) * state.rho[(int) state.field[x][y]];
                        if (state.field[x][y] == '.')
                            force *= NumericType(0.8);
                        if (state.field[x + dx][y + dy] == '#') {
                            state.p[x][y] += force / NumericType(state.dirs[x][y]);
                            total_delta_p += force / NumericType(state.dirs[x][y]);
                        } else {
                            state.p[x + dx][y + dy] += force / NumericType(state.dirs[x + dx][y + dy]);
                            total_delta_p += force / NumericType(state.dirs[x + dx][y + dy]);
                        }
                    }
                }
            }
        }

        state.UT += 2;
        prop = false;
        for (size_t x = 0; x < SizeN; ++x) {
            for (size_t y = 0; y < SizeM; ++y) {
                if (state.field[x][y] != '#' && state.last_use[x][y] != state.UT) {
                    if (get_random01() < get_move_prob(x, y)) {
                        prop = true;
                        do_propagate_move(x, y, true);
                    } else {
                        do_propagate_stop(x, y, true);
                    }
                }
            }
        }

        if (prop) {
            cout << "Tick " << iteration << ":\n";
            for (size_t x = 0; x < SizeN; ++x) {
                cout << state.field[x] << "\n";
            }
        }
    }

    void run(NumericType g = NumericType(0.1)) {
        for (size_t i = 0; i < 1'000'000; ++i) {
            simulate_step(g, i);
        }
    }

    friend int main();
};

int main() {
    SimulationState<Fixed<32, 16>, DEFAULT_N, DEFAULT_M> state;
    FluidSimulator<Fixed<32, 16>> simulator(state);
    simulator.run();
    return 0;
}
