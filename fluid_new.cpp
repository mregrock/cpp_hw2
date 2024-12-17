#include <iostream>
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

#ifndef DTYPES
#define DTYPES FLOAT,FIXED(32,16),DOUBLE
#endif

#ifndef DSIZES
#define DSIZES S(36,84)
#endif

struct Size {
    size_t n, m;
    constexpr Size(size_t n_, size_t m_) : n(n_), m(m_) {}
};

constexpr Size parse_size(const char* s) {
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
    return Size(n, m);
}

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

    explicit operator float() const { return v / float(StorageType(1) << K); }
    explicit operator double() const { return v / double(StorageType(1) << K); }

    friend Fixed operator/(Fixed a, int b) {
        return Fixed::from_raw(a.v / b);
    }

    friend Fixed operator*(Fixed a, int b) {
        return Fixed::from_raw(a.v * b);
    }

    friend Fixed operator*(int a, Fixed b) {
        return b * a;
    }

    template<size_t N2, size_t K2>
    explicit operator Fixed<N2,K2>() const {
        if constexpr (K2 >= K)
            return Fixed<N2,K2>::from_raw(static_cast<typename Fixed<N2,K2>::StorageType>(v) << (K2 - K));
        else
            return Fixed<N2,K2>::from_raw(static_cast<typename Fixed<N2,K2>::StorageType>(v) >> (K - K2));
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

template<
    typename PressureType,
    typename VelocityType,
    typename VFlowType,
    size_t N = DEFAULT_N,
    size_t M = DEFAULT_M
>
struct SimulationState {
    PressureType rho[256];
    PressureType p[N][M]{};
    PressureType old_p[N][M]{};

    char field[N][M + 1];
    int last_use[N][M]{};
    int UT{0};
    mt19937 rnd{1337};
    int dirs[N][M]{};

    VectorField<VelocityType, N, M> velocity{};
    VectorField<VFlowType, N, M> velocity_flow{};

    SimulationState() {
        memcpy(field, initial_field, sizeof(field));
    }
};

template<typename T>
struct NumericTraits {
    static T from_raw(int32_t x) { return T(x) / T(1 << 16); }
};

template<size_t N, size_t K>
struct NumericTraits<Fixed<N,K>> {
    static Fixed<N,K> from_raw(typename Fixed<N,K>::StorageType x) {
        return Fixed<N,K>::from_raw(x);
    }
};

template<typename T>
struct is_fixed : std::false_type {};

template<size_t N, size_t K>
struct is_fixed<Fixed<N,K>> : std::true_type {};

template<typename T>
inline constexpr bool is_fixed_v = is_fixed<T>::value;

template<typename To, typename From>
To convert_to(const From& value) {
    if constexpr (std::is_same_v<To, From>) {
        return value;
    } else if constexpr (std::is_floating_point_v<To>) {
        if constexpr (std::is_floating_point_v<From>) {
            return static_cast<To>(value);
        } else {
            return static_cast<To>(static_cast<double>(value));
        }
    } else if constexpr (is_fixed_v<To>) {
        if constexpr (std::is_floating_point_v<From>) {
            return To(value);
        } else if constexpr (is_fixed_v<From>) {
            return static_cast<To>(value);
        } else {
            return To(value);
        }
    }
    return To(value);
}

template<typename T1, typename T2>
auto add(const T1& a, const T2& b) {
    using CommonType = std::conditional_t<
        std::is_floating_point_v<T1> || std::is_floating_point_v<T2>,
        double,
        std::conditional_t<is_fixed_v<T1>, T1,
        std::conditional_t<is_fixed_v<T2>, T2, T1>>>;
    return convert_to<CommonType>(a) + convert_to<CommonType>(b);
}

template<typename T1, typename T2>
auto subtract(const T1& a, const T2& b) {
    using CommonType = std::conditional_t<
        std::is_floating_point_v<T1> || std::is_floating_point_v<T2>,
        double,
        std::conditional_t<is_fixed_v<T1>, T1,
        std::conditional_t<is_fixed_v<T2>, T2, T1>>>;
    return convert_to<CommonType>(a) - convert_to<CommonType>(b);
}

template<typename T1, typename T2>
auto multiply(const T1& a, const T2& b) {
    using CommonType = std::conditional_t<
        std::is_floating_point_v<T1> || std::is_floating_point_v<T2>,
        double,
        std::conditional_t<is_fixed_v<T1>, T1,
        std::conditional_t<is_fixed_v<T2>, T2, T1>>>;
    return convert_to<CommonType>(a) * convert_to<CommonType>(b);
}

template<typename T1, typename T2>
auto divide(const T1& a, const T2& b) {
    using CommonType = std::conditional_t<
        std::is_floating_point_v<T1> || std::is_floating_point_v<T2>,
        double,
        std::conditional_t<is_fixed_v<T1>, T1,
        std::conditional_t<is_fixed_v<T2>, T2, T1>>>;
    return convert_to<CommonType>(a) / convert_to<CommonType>(b);
}

template<typename T1, typename T2>
auto min_value(const T1& a, const T2& b) {
    using CommonType = std::conditional_t<
        std::is_floating_point_v<T1> || std::is_floating_point_v<T2>,
        double,
        std::conditional_t<is_fixed_v<T1>, T1,
        std::conditional_t<is_fixed_v<T2>, T2, T1>>>;
    auto converted_a = convert_to<CommonType>(a);
    auto converted_b = convert_to<CommonType>(b);
    return converted_a < converted_b ? converted_a : converted_b;
}

template<
    typename PressureType = Fixed<32, 16>,
    typename VelocityType = Fixed<32, 16>,
    typename VFlowType = Fixed<32, 16>,
    size_t SizeN = DEFAULT_N,
    size_t SizeM = DEFAULT_M
>
class FluidSimulator {
private:
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

    SimulationState<PressureType, VelocityType, VFlowType, SizeN, SizeM>& state;

    struct ParticleParams {
        char type;
        PressureType cur_p;
        array<VelocityType, 4> v;
        char (&field)[SizeN][SizeM + 1];
        PressureType (&p)[SizeN][SizeM];
        VectorField<VelocityType, SizeN, SizeM>& velocity;

        ParticleParams(
            char (&f)[SizeN][SizeM + 1],
            PressureType (&p_arr)[SizeN][SizeM],
            VectorField<VelocityType, SizeN, SizeM>& vel
        ) : field(f), p(p_arr), velocity(vel) {}

        void swap_with(int x, int y) {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };

    using FlowResult = tuple<PressureType, bool, pair<int, int>>;
    
    FlowResult do_propagate_flow(int x, int y, PressureType lim) {
        state.last_use[x][y] = state.UT - 1;
        PressureType ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (state.field[nx][ny] != '#' && state.last_use[nx][ny] < state.UT) {
                auto cap = state.velocity.get(x, y, dx, dy);
                auto flow = state.velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                auto vp = min_value(lim, subtract(cap, flow));
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
    explicit FluidSimulator(SimulationState<PressureType, VelocityType, VFlowType, SizeN, SizeM>& simulation_state)
        : state(simulation_state)
    {
        initialize();
    }

    void initialize() {
        state.rho[' '] = PressureType(0.01);
        state.rho['.'] = PressureType(1000);

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

    PressureType get_random01() {
        return NumericTraits<PressureType>::from_raw(state.rnd() & ((1 << 16) - 1));
    }

    PressureType get_move_prob(int x, int y) {
        PressureType sum = 0;
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
            sum = add(sum, v);
        }
        return sum;
    }

    bool do_propagate_move(int x, int y, bool is_first) {
        state.last_use[x][y] = state.UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<PressureType, deltas.size()> tres;
            PressureType sum = 0;
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

            PressureType p = get_random01() * sum;
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

    void simulate_step(PressureType g, size_t iteration) {
        PressureType total_delta_p = 0;
        
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
                        state.p[x][y] -= force / PressureType(state.dirs[x][y]);
                        total_delta_p -= force / PressureType(state.dirs[x][y]);
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
                            force *= PressureType(0.8);
                        if (state.field[x + dx][y + dy] == '#') {
                            state.p[x][y] += force / PressureType(state.dirs[x][y]);
                            total_delta_p += force / PressureType(state.dirs[x][y]);
                        } else {
                            state.p[x + dx][y + dy] += force / PressureType(state.dirs[x + dx][y + dy]);
                            total_delta_p += force / PressureType(state.dirs[x + dx][y + dy]);
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

    void run(PressureType g = PressureType(0.1)) {
        for (size_t i = 0; i < 1'000'000; ++i) {
            simulate_step(g, i);
        }
    }
};

template<typename P, typename V, typename VF>
void run_simulation(size_t n, size_t m) {
    if (n != DEFAULT_N || m != DEFAULT_M) {
        cerr << "Only default size " << DEFAULT_N << "x" << DEFAULT_M << " is supported\n";
        return;
    }
    SimulationState<P, V, VF, DEFAULT_N, DEFAULT_M> state;
    FluidSimulator<P, V, VF, DEFAULT_N, DEFAULT_M> simulator(state);
    simulator.run();
}

#define FLOAT_TYPE if (type == "FLOAT") { run_simulation<float>(size.n, size.m); found = true; }
#define DOUBLE_TYPE if (type == "DOUBLE") { run_simulation<double>(size.n, size.m); found = true; }
#define FIXED_TYPE(N, K) if (type == "FIXED(" #N "," #K ")") { run_simulation<Fixed<N,K>>(size.n, size.m); found = true; }

string get_arg(string_view arg_name, int argc, char** argv, string_view default_value) {
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == arg_name) {
            return argv[i + 1];
        }
    }
    return string(default_value);
}

bool is_valid_type(const string& type) {
    if (type == "FLOAT" || type == "DOUBLE") return true;
    
    if (type.starts_with("FIXED(") || type.starts_with("FAST_FIXED(")) {
        size_t pos = type.find(',');
        if (pos == string::npos) return false;
        
        size_t end_pos = type.find(')', pos);
        if (end_pos == string::npos) return false;
        
        return true;
    }
    
    return false;
}

int main(int argc, char** argv) {
    string p_type = get_arg("--p-type", argc, argv, "FIXED(32,16)");
    string v_type = get_arg("--v-type", argc, argv, "FIXED(32,16)");
    string v_flow_type = get_arg("--v-flow-type", argc, argv, "FIXED(32,16)");
    string size_str = get_arg("--size", argc, argv, "S(36,84)");
    
    if (!is_valid_type(p_type)) {
        cerr << "Invalid p_type: " << p_type << endl;
        return 1;
    }
    if (!is_valid_type(v_type)) {
        cerr << "Invalid v_type: " << v_type << endl;
        return 1;
    }
    if (!is_valid_type(v_flow_type)) {
        cerr << "Invalid v_flow_type: " << v_flow_type << endl;
        return 1;
    }

    bool found = false;
    auto size = parse_size(size_str.c_str());

    if (p_type == "FLOAT" && v_type == "FIXED(32,5)" && v_flow_type == "DOUBLE") {
        run_simulation<float, Fixed<32,5>, double>(size.n, size.m);
        found = true;
    }

    if (!found) {
        cerr << "Unsupported type combination" << endl;
        return 1;
    }
    
    return 0;
}

template<typename T, size_t N, size_t K>
auto operator+(T a, Fixed<N,K> b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(a) + T(static_cast<T>(b));
    else
        return Fixed<N,K>(a) + b;
}

template<typename T, size_t N, size_t K>
auto operator-(T a, Fixed<N,K> b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(a) - T(static_cast<T>(b));
    else
        return Fixed<N,K>(a) - b;
}

template<typename T, size_t N, size_t K>
auto operator*(T a, Fixed<N,K> b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(a) * T(static_cast<T>(b));
    else
        return Fixed<N,K>(a) * b;
}

template<typename T, size_t N, size_t K>
auto operator/(T a, Fixed<N,K> b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(a) / T(static_cast<T>(b));
    else
        return Fixed<N,K>(a) / b;
}

template<typename T, size_t N, size_t K>
auto operator+(Fixed<N,K> a, T b) { return b + a; }

template<typename T, size_t N, size_t K>
auto operator-(Fixed<N,K> a, T b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(static_cast<T>(a)) - T(b);
    else
        return a - Fixed<N,K>(b);
}

template<typename T, size_t N, size_t K>
auto operator*(Fixed<N,K> a, T b) { return b * a; }

template<typename T, size_t N, size_t K>
auto operator/(Fixed<N,K> a, T b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(static_cast<T>(a)) / T(b);
    else
        return a / Fixed<N,K>(b);
}

template<typename T, size_t N, size_t K>
bool operator<(Fixed<N,K> a, T b) {
    if constexpr (std::is_floating_point_v<T>)
        return T(static_cast<T>(a)) < b;
    else
        return a < Fixed<N,K>(b);
}

template<typename T, size_t N, size_t K>
bool operator<(T a, Fixed<N,K> b) { return !(b < a); }

template<typename T1, typename T2>
auto& operator+=(T1& a, const T2& b) {
    a = add(a, b);
    return a;
}

template<typename T1, typename T2>
auto& operator-=(T1& a, const T2& b) {
    a = subtract(a, b);
    return a;
}
