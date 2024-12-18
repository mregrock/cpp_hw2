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

#ifndef TYPES
#define TYPES FLOAT,FIXED(31,17),FIXED(25, 11),FIXED(32, 16),DOUBLE
#endif

#ifndef SIZES
#define SIZES S(36,84)
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
    
    static constexpr size_t bits = N;
    static constexpr size_t frac = K;
    
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

template<typename PressureType, typename VelocityType, typename VFlowType, size_t N = DEFAULT_N, size_t M = DEFAULT_M>
class FluidSimulator {
private:
    static constexpr size_t T = 1'000'000;
    static constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    
    char field[N][M + 1];
    PressureType p[N][M]{}, old_p[N][M];
    int last_use[N][M]{};
    int UT = 0;
    mt19937 rnd{1337};
    PressureType rho[256];

    struct VectorField {
        array<VelocityType, 4> v[N][M]{};
        
        VelocityType& get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }

        VelocityType& add(int x, int y, int dx, int dy, VelocityType dv) {
            return get(x, y, dx, dy) += dv;
        }
    };

    VectorField velocity{}, velocity_flow{};
    int dirs[N][M]{};

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
        memcpy(field, state.field, sizeof(field));
        rho[' '] = PressureType(0.01);
        rho['.'] = PressureType(1000);
        
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#') continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
    }

    void run() {
        PressureType g = PressureType(0.1);
        
        for (size_t i = 0; i < T; ++i) {
            PressureType total_delta_p = 0;
            
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#') continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, VelocityType(g));
                }
            }

            // Apply forces from p
            memcpy(old_p, p, sizeof(p));
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = delta_p;
                            auto& contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int)field[nx][ny]] >= force) {
                                contr -= force / rho[(int)field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int)field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int)field[x][y]]);
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = false;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, VelocityType(1));
                            if (t > VelocityType(0)) prop = true;
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#') continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > VelocityType(0)) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int)field[x][y]];
                            if (field[x][y] == '.') force *= PressureType(0.8);
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / to_pressure(dirs[x][y]);
                                total_delta_p += force / to_pressure(dirs[x][y]);
                            } else {
                                p[x + dx][y + dy] += force / to_pressure(dirs[x + dx][y + dy]);
                                total_delta_p += force / to_pressure(dirs[x + dx][y + dy]);
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
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
                for (size_t x = 0; x < N; ++x) {
                    cout << field[x] << "\n";
                }
            }
        }
    }
};

template<typename P, typename V, typename VF>
void run_simulation(size_t n = DEFAULT_N, size_t m = DEFAULT_M) {
    static_assert(!std::is_same_v<P, void>, "PressureType cannot be void");
    static_assert(!std::is_same_v<V, void>, "VelocityType cannot be void");
    static_assert(!std::is_same_v<VF, void>, "VFlowType cannot be void");
    
    SimulationState<P, V, VF, DEFAULT_N, DEFAULT_M> state;
    FluidSimulator<P, V, VF, DEFAULT_N, DEFAULT_M> simulator(state);
    simulator.run();
}

template<typename T>
std::string get_pretty_type_name() {
    if constexpr (std::is_same_v<T, float>) {
        return "float";
    } else if constexpr (std::is_same_v<T, double>) {
        return "double";
    } else if constexpr (is_fixed_v<T>) {
        return "Fixed<" + std::to_string(T::bits) + "," + std::to_string(T::frac) + ">";
    } else {
        return "unknown";
    }
}

template<typename T = void>
struct TypeWrapper {
    using type = T;
    T value;
    
    TypeWrapper() = default;
    TypeWrapper(const T& v) : value(v) {}
    
    TypeWrapper& operator+=(const TypeWrapper& other) { value += other.value; return *this; }
    TypeWrapper& operator-=(const TypeWrapper& other) { value -= other.value; return *this; }
    TypeWrapper& operator*=(const TypeWrapper& other) { value *= other.value; return *this; }
    TypeWrapper& operator/=(const TypeWrapper& other) { value /= other.value; return *this; }
    
    template<typename U>
    TypeWrapper& operator+=(const U& other) { value += other; return *this; }
    template<typename U>
    TypeWrapper& operator-=(const U& other) { value -= other; return *this; }
    template<typename U>
    TypeWrapper& operator*=(const U& other) { value *= other; return *this; }
    template<typename U>
    TypeWrapper& operator/=(const U& other) { value /= other; return *this; }
    
    bool operator<(const TypeWrapper& other) const { return value < other.value; }
    bool operator>(const TypeWrapper& other) const { return value > other.value; }
    bool operator==(const TypeWrapper& other) const { return value == other.value; }
    bool operator<=(const TypeWrapper& other) const { return value <= other.value; }
    bool operator>=(const TypeWrapper& other) const { return value >= other.value; }
    
    template<typename U>
    bool operator<(const U& other) const { return value < other; }
    template<typename U>
    bool operator>(const U& other) const { return value > other; }
    template<typename U>
    bool operator<=(const U& other) const { return value <= other; }
    template<typename U>
    bool operator>=(const U& other) const { return value >= other; }
    template<typename U>
    bool operator==(const U& other) const { return value == other; }
    
    TypeWrapper operator+(const TypeWrapper& other) const { return TypeWrapper(value + other.value); }
    TypeWrapper operator-(const TypeWrapper& other) const { return TypeWrapper(value - other.value); }
    TypeWrapper operator*(const TypeWrapper& other) const { return TypeWrapper(value * other.value); }
    TypeWrapper operator/(const TypeWrapper& other) const { return TypeWrapper(value / other.value); }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic_v<U>>>
    TypeWrapper(U v) : value(T(v)) {}

    std::string get_type_string() const {
        return get_pretty_type_name<T>();
    }
};

template<>
struct TypeWrapper<void> {
    using type = void;
    
    TypeWrapper() = default;
    
    template<typename U>
    TypeWrapper(const U&) {}
    
    template<typename U>
    TypeWrapper& operator=(const U&) { return *this; }
    
    TypeWrapper& operator+=(const TypeWrapper&) { return *this; }
    TypeWrapper& operator-=(const TypeWrapper&) { return *this; }
    TypeWrapper& operator*=(const TypeWrapper&) { return *this; }
    TypeWrapper& operator/=(const TypeWrapper&) { return *this; }
    
    template<typename U>
    TypeWrapper& operator+=(const U&) { return *this; }
    template<typename U>
    TypeWrapper& operator-=(const U&) { return *this; }
    template<typename U>
    TypeWrapper& operator*=(const U&) { return *this; }
    template<typename U>
    TypeWrapper& operator/=(const U&) { return *this; }
    
    bool operator<(const TypeWrapper&) const { return false; }
    bool operator>(const TypeWrapper&) const { return false; }
    bool operator==(const TypeWrapper&) const { return true; }
    bool operator<=(const TypeWrapper&) const { return true; }
    bool operator>=(const TypeWrapper&) const { return true; }
    
    template<typename U>
    bool operator<(const U&) const { return false; }
    template<typename U>
    bool operator>(const U&) const { return false; }
    template<typename U>
    bool operator<=(const U&) const { return true; }
    template<typename U>
    bool operator>=(const U&) const { return true; }
    template<typename U>
    bool operator==(const U&) const { return true; }
    
    TypeWrapper operator+(const TypeWrapper&) const { return TypeWrapper{}; }
    TypeWrapper operator-(const TypeWrapper&) const { return TypeWrapper{}; }
    TypeWrapper operator*(const TypeWrapper&) const { return TypeWrapper{}; }
    TypeWrapper operator/(const TypeWrapper&) const { return TypeWrapper{}; }

    std::string get_type_string() const {
        return "void";
    }
};

template<typename T, typename U>
bool operator<(const U& a, const TypeWrapper<T>& b) { return b > a; }
template<typename T, typename U>
bool operator>(const U& a, const TypeWrapper<T>& b) { return b < a; }
template<typename T, typename U>
bool operator<=(const U& a, const TypeWrapper<T>& b) { return b >= a; }
template<typename T, typename U>
bool operator>=(const U& a, const TypeWrapper<T>& b) { return b <= a; }
template<typename T, typename U>
bool operator==(const U& a, const TypeWrapper<T>& b) { return b == a; }

pair<size_t, size_t> parse_fixed_params(const string& type) {
    size_t start = type.find('(') + 1;
    size_t comma = type.find(',', start);
    size_t end = type.find(')', comma);
    
    size_t N = stoul(type.substr(start, comma - start));
    size_t K = stoul(type.substr(comma + 1, end - comma - 1));
    return {N, K};
}

template<typename T>
struct TypeHolder {
    using type = T;
    
    TypeHolder() {
        cerr << "Creating TypeHolder<" << get_pretty_type_name<T>() << ">" << endl;
    }
};

template<std::size_t I, typename... Types>
struct TypesPackGetAt;

template<typename First, typename... Rest>
struct TypesPackGetAt<0, First, Rest...> {
    using type = First;
};

template<std::size_t I, typename First, typename... Rest>
struct TypesPackGetAt<I, First, Rest...> {
    using type = typename TypesPackGetAt<I-1, Rest...>::type;
};

template<typename... Types>
struct TypeParser {
    using ResultType = std::variant<TypeHolder<Types>...>;
    
    static ResultType parse_type(const string& type) {
        cerr << "\nTypeParser::parse_type called with type: " << type << endl;
        
        ResultType result{TypeHolder<typename TypesPackGetAt<0, Types...>::type>{}};
        bool found = false;
        
        auto try_match = [&](auto dummy) -> bool {
            using T = decltype(dummy);
            if (type_matches<T>(type)) {
                std::cout << "Match found! Type: " << get_pretty_type_name<T>() << std::endl;
                result = TypeHolder<T>{};
                found = true;
                return true;
            }
            std::cout << "Trying type: " << get_pretty_type_name<T>() << std::endl;
            return false;
        };
        
        (try_match(Types{}) || ...);
        
        if (!found) {
            cerr << "No matching type found for: " << type << endl;
            throw std::runtime_error("Type not found");
        }
        
        return result;
    }

private:
    template<typename T>
    static bool type_matches(const string& type) {
        if constexpr (std::is_same_v<T, float>) {
            return type == "FLOAT";
        } else if constexpr (std::is_same_v<T, double>) {
            return type == "DOUBLE";
        } else if constexpr (is_fixed_v<T>) {
            if (!type.starts_with("FIXED(")) {
                return false;
            }
            
            auto params = parse_fixed_params(type);
            return (params.first == T::bits && params.second == T::frac);
        }
        return false;
    }
};

template <typename... Ts>
struct WrapTypes {
    using type = std::variant<TypeWrapper<Ts>...>;
};

auto get_type_wrapper(const string& type) {
    #define FIXED(N,K) Fixed<N,K>
    #define FLOAT float
    #define DOUBLE double
    
    using Parser = TypeParser<TYPES>;
    auto result = Parser::parse_type(type);
    using TypeWrapperVariant = typename WrapTypes<TYPES>::type;
    #undef FIXED
    #undef FLOAT
    #undef DOUBLE
    return std::visit([](auto&& x) -> TypeWrapperVariant {
        using ActualType = typename std::decay_t<decltype(x)>::type;
        cerr << "Creating TypeWrapper<" << get_pretty_type_name<ActualType>() << ">" << endl;
        return TypeWrapper<ActualType>{};
    }, result);
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
    }
    return false;
}

template<typename... Types>
struct TypesList {
    static constexpr size_t size = sizeof...(Types);
    template<size_t I>
    using type_at = typename std::tuple_element<I, std::tuple<Types...>>::type;
};

template<typename AllowedTypes, typename SelectedTypes>
struct TypeSelector {
    template<typename... Selected>
    static bool try_combinations(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m) {
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

        cerr << "Found matching types:\n"
             << "P: " << get_pretty_type_name<P>() << "\n"
             << "V: " << get_pretty_type_name<V>() << "\n"
             << "VF: " << get_pretty_type_name<VF>() << "\n";

        run_simulation<P, V, VF>(n, m);
        return true;
    }
};

template<typename... Types>
bool try_all_type_combinations(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m) {
    return TypeSelector<TypesList<Types...>, TypesList<>>::try_combinations(p_type, v_type, v_flow_type, n, m);
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

template<size_t N, size_t K>
Fixed<N,K> operator+(float a, Fixed<N,K> b) {
    return Fixed<N,K>(a) + b;
}

template<size_t N, size_t K>
Fixed<N,K> operator+(Fixed<N,K> a, float b) {
    return a + Fixed<N,K>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator-(float a, Fixed<N,K> b) {
    return Fixed<N,K>(a) - b;
}

template<size_t N, size_t K>
Fixed<N,K> operator-(Fixed<N,K> a, float b) {
    return a - Fixed<N,K>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator*(float a, Fixed<N,K> b) {
    return Fixed<N,K>(a) * b;
}

template<size_t N, size_t K>
Fixed<N,K> operator*(Fixed<N,K> a, float b) {
    return a * Fixed<N,K>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator/(float a, Fixed<N,K> b) {
    return Fixed<N,K>(a) / b;
}

template<size_t N, size_t K>
Fixed<N,K> operator/(Fixed<N,K> a, float b) {
    return a / Fixed<N,K>(b);
}

template<size_t N, size_t K>
Fixed<N,K>& operator+=(float& a, Fixed<N,K> b) {
    a = static_cast<float>(Fixed<N,K>(a) + b);
    return a;
}

template<size_t N, size_t K>
float& operator-=(float& a, Fixed<N,K> b) {
    a = static_cast<float>(Fixed<N,K>(a) - b);
    return a;
}

template<size_t N, size_t K>
Fixed<N,K> operator+(double a, Fixed<N,K> b) {
    return Fixed<N,K>(a) + b;
}

template<size_t N, size_t K>
Fixed<N,K>& operator-=(Fixed<N,K>& a, float b) { 
    return a = a - Fixed<N,K>(b); 
}

template<size_t N, size_t K>
Fixed<N,K>& operator+=(Fixed<N,K>& a, float b) { 
    return a = a + Fixed<N,K>(b); 
}

template<size_t N, size_t K>
Fixed<N,K>& operator*=(Fixed<N,K>& a, float b) { 
    return a = a * Fixed<N,K>(b); 
}

int main(int argc, char** argv) {
    string p_type = get_arg("--p-type", argc, argv, "FIXED(32,16)");
    string v_type = get_arg("--v-type", argc, argv, "DOUBLE");
    string v_flow_type = get_arg("--v-flow-type", argc, argv, "FLOAT");
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
