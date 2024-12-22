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
#include <variant>
#include "src/fixed.hpp"
#include "src/fast_fixed.hpp"
#include "src/fixed_operators.hpp"

using namespace std;

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

template<typename NumericType, size_t N, size_t M>
struct VectorField {
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    SimulationMatrix<std::array<NumericType, 4>, N, M> v{};

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
    mt19937 rnd{1337};
    SimulationMatrix<int, N, M> dirs{};

    VectorField<VelocityType, N, M> velocity{};
    VectorField<VFlowType, N, M> velocity_flow{};

    SimulationState(const std::vector<std::string>& initial_field = {}) {
        std::cerr << "Creating SimulationState with N=" << N << ", M=" << M << std::endl;
        std::cerr << "Initial field empty: " << initial_field.empty() << std::endl;
        
        if constexpr (N == dynamic_size || M == dynamic_size) {
            std::cerr << "Dynamic size branch" << std::endl;
            if (!initial_field.empty()) {
                size_t n = initial_field.size();
                size_t m = initial_field[0].size();
                std::cerr << "Resizing to " << n << "x" << m << std::endl;
                
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
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = 0; j < m; ++j) {
                        field[i][j] = initial_field[i][j];
                    }
                }
            }
        } else {
            std::cerr << "Static size branch" << std::endl;
            if (!initial_field.empty()) {
                if (initial_field.size() != N || initial_field[0].size() != M) {
                    std::cerr << "Size mismatch: got " << initial_field.size() << "x" 
                                << initial_field[0].size() << " instead of " << N << "x" << M << std::endl;
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

template<typename T>
struct is_fast_fixed : std::false_type {};

template<size_t N, size_t K>
struct is_fast_fixed<FastFixed<N,K>> : std::true_type {};

template<typename T>
inline constexpr bool is_fast_fixed_v = is_fast_fixed<T>::value;

template<typename PressureType, typename VelocityType, typename VFlowType, size_t N = DEFAULT_N, size_t M = DEFAULT_M>
class FluidSimulator {
private:
    size_t runtime_n, runtime_m;

    static constexpr size_t T = 1'000'000;
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    
    SimulationMatrix<char, N, M> field{};
    SimulationMatrix<PressureType, N, M> p{}, old_p{};
    SimulationMatrix<int, N, M> last_use{};
    int UT = 0;
    mt19937 rnd{1337};
    PressureType rho[256];

    struct VectorField {
        SimulationMatrix<array<VelocityType, 4>, N, M> v{};
        
        VelocityType& get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }

        VelocityType& add(int x, int y, int dx, int dy, VelocityType dv) {
            return get(x, y, dx, dy) += dv;
        }

        void resize(size_t n, size_t m) {
            v.resize(n);
            for (auto& row : v) {
                row.assign(m, std::array<VelocityType, 4>{});
            }
        }
    };

    VectorField velocity{}, velocity_flow{};
    SimulationMatrix<int, N, M> dirs{};

    struct ParticleParams {
        char type;
        PressureType cur_p;
        array<VelocityType, 4> v;

        void swap_with(FluidSimulator* sim, int x, int y) {
            swap(sim->field[x][y], type);
            swap(sim->p[x][y], cur_p);
            swap(sim->velocity.v[x][y], v);
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

    void init_state(const SimulationState<PressureType, VelocityType, VFlowType, N, M>& state) {
        field = state.field;
        p = state.p;
        old_p = state.old_p;
        last_use = state.last_use;
        UT = state.UT;
        dirs = state.dirs;
        velocity.resize(runtime_n, runtime_m);
        velocity_flow.resize(runtime_n, runtime_m);
    }

    tuple<VelocityType, bool, pair<int, int>> propagate_flow(int x, int y, VelocityType lim) {
        last_use[x][y] = UT - 1;
        VelocityType ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) continue;
                
                auto vp = min(lim, cap - flow);
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    return {vp, true, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    return {t, prop && end != pair(x, y), end};
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
            array<VelocityType, 4> tres;
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
            size_t d = ranges::upper_bound(tres, p) - tres.begin();

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

    void run() {
        PressureType g = PressureType(0.1);
        for (size_t x = 0; x < get_n(); ++x) {
            field[x][get_m()] = '\0';
            cout << field[x].data() << "\n";
        }
        for (size_t i = 0; i < T; ++i) {
            PressureType total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < get_n(); ++x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] == '#') continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, VelocityType(g));
                }
            }

            // Apply forces from p
            old_p = p;
            for (size_t x = 0; x < get_n(); ++x) {
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
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }
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
            for (size_t x = 0; x < get_n(); ++x) {
                for (size_t y = 0; y < get_m(); ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > VelocityType(0)) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = to_pressure(old_v - new_v) * rho[(int)field[x][y]];
                            if (field[x][y] == '.') force *= PressureType(0.8);
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

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

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < get_n(); ++x) {
                    field[x][get_m()] = '\0';
                    cout << field[x].data() << "\n";
                }
            }
        }
    }
};

template<typename P, typename V, typename VF, typename Size=SizeType<dynamic_size, dynamic_size>>
void run_simulation(size_t n = DEFAULT_N, size_t m = DEFAULT_M) {    
    try {
        SimulationState<P, V, VF, Size::n, Size::m> state(initial_field);
        FluidSimulator<P, V, VF, Size::n, Size::m> simulator(state);
        
        simulator.run();
        std::cerr << "Simulation completed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception in run_simulation: " << e.what() << std::endl;
        throw;
    }
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

pair<size_t, size_t> parse_fixed_params(const string& type) {
    size_t start = type.find('(') + 1;
    size_t comma = type.find(',', start);
    size_t end = type.find(')', comma);
    
    size_t N = stoul(type.substr(start, comma - start));
    size_t K = stoul(type.substr(comma + 1, end - comma - 1));
    return {N, K};
}

template<typename T>
static bool matches_type(const string& type) {
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

template<typename... Types>
struct TypesList {
    static constexpr size_t size = sizeof...(Types);
    template<size_t I>
    using type_at = typename std::tuple_element<I, std::tuple<Types...>>::type;
};

template<typename... Sizes>
struct SizesList {
    static constexpr size_t size = sizeof...(Sizes);
    template<size_t I>
    using size_at = typename std::tuple_element<I, std::tuple<Sizes...>>::type;
};

template<typename AllowedTypes, typename AllowedSizes, typename SelectedTypes>
struct TypeSelector {
    template<typename... Selected>
    static bool try_combinations(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
        std::cerr << "Entering try_combinations with sizes n=" << n << ", m=" << m << std::endl;
        return try_all_p_types<0>(p_type, v_type, v_flow_type, n, m);
    }

private:
    template<size_t I>
    static bool try_all_p_types(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using P = typename AllowedTypes::template type_at<I>;
            return try_with_p_type<P>(p_type, v_type, v_flow_type, n, m) ||
                   try_all_p_types<I + 1>(p_type, v_type, v_flow_type, n, m);
        }
    }

    template<typename P>
    static bool try_with_p_type(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
        if (!matches_type<P>(p_type)) return false;
        return try_all_v_types<P, 0>(p_type, v_type, v_flow_type, n, m);
    }

    template<typename P, size_t I>
    static bool try_all_v_types(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using V = typename AllowedTypes::template type_at<I>;
            return try_with_v_type<P, V>(p_type, v_type, v_flow_type, n, m) ||
                   try_all_v_types<P, I + 1>(p_type, v_type, v_flow_type, n, m);
        }
    }

    template<typename P, typename V>
    static bool try_with_v_type(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
        if (!matches_type<V>(v_type)) return false;
        return try_all_vf_types<P, V, 0>(p_type, v_type, v_flow_type, n, m);
    }

    template<typename P, typename V, size_t I>
    static bool try_all_vf_types(const string& p_type, const string& v_type, const string& v_flow_type,
                                size_t n, size_t m) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using VF = typename AllowedTypes::template type_at<I>;
            return try_with_vf_type<P, V, VF>(p_type, v_type, v_flow_type, n, m) ||
                   try_all_vf_types<P, V, I + 1>(p_type, v_type, v_flow_type, n, m);
        }
    }

    template<typename P, typename V, typename VF>
    static bool try_with_vf_type(const string& p_type, const string& v_type, const string& v_flow_type,
                                size_t n, size_t m) {
        if (!matches_type<VF>(v_flow_type)) return false;

        std::cerr << "Found matching types and about to create simulation" << std::endl;
        
        if constexpr (AllowedSizes::size == 0) {
            run_simulation<P, V, VF>(n, m);
            return true;
        } else {
            return try_all_sizes<P, V, VF, 0>(p_type, v_type, v_flow_type, n, m);
        }
    }
    template<typename P, typename V, typename VF, size_t I>
    static bool try_all_sizes(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m) {
        if constexpr (I >= AllowedSizes::size) {
            return false;
        } else {
            using CurrentSize = typename AllowedSizes::template size_at<I>;
            return try_with_size<P, V, VF, CurrentSize>(p_type, v_type, v_flow_type, n, m) ||
                   try_all_sizes<P, V, VF, I + 1>(p_type, v_type, v_flow_type, n, m);
        }
    }

    template<typename P, typename V, typename VF, typename CurrentSize>
    static bool try_with_size(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m) {
        if (CurrentSize::n != n || CurrentSize::m != m) return false;
        
        std::cerr << "Found matching size: " << CurrentSize::n << "x" << CurrentSize::m << std::endl;
        run_simulation<P, V, VF, CurrentSize>(CurrentSize::n, CurrentSize::m);
        return true;
    }
};

template<typename... Types>
bool try_all_type_combinations(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m) {
    #define S(N, M) SizeType<N, M>
    return TypeSelector<TypesList<Types...>, SizesList<SIZES>, TypesList<>>::try_combinations(p_type, v_type, v_flow_type, n, m);
}

bool create_and_run_simulation(const string& p_type, const string& v_type, const string& v_flow_type, 
                             size_t n, size_t m) {
    try {
        cerr << "\nTrying to create simulation with types:" << endl;
        cerr << "p_type: " << p_type << endl;
        cerr << "v_type: " << v_type << endl;
        cerr << "v_flow_type: " << v_flow_type << endl;
        #define FLOAT float
        #define DOUBLE double
        #define FIXED(N, K) Fixed<N, K>
        #define FAST_FIXED(N, K) FastFixed<N, K>
        if (!try_all_type_combinations<TYPES>(p_type, v_type, v_flow_type, n, m)) {
            cerr << "Error: No matching type combination found" << endl;
            return false;
        }
        return true;
    }
    catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return false;
    }
}

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
    string p_type = get_arg("--p-type", argc, argv, "FAST_FIXED(32,16)");
    string v_type = get_arg("--v-type", argc, argv, "FIXED(31,17)");
    string v_flow_type = get_arg("--v-flow-type", argc, argv, "DOUBLE");
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

    auto size = parse_size(size_str.c_str());
    
    if (!create_and_run_simulation(p_type, v_type, v_flow_type, size.n, size.m)) {
        cerr << "Failed to create simulation with types: " << p_type << ", " << v_type << ", " << v_flow_type << endl;
        return 1;
    }
    
    return 0;
}
