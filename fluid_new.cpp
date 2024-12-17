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

constexpr size_t N = 36, M = 84;
constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

char field[N][M + 1] = {
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

struct Fixed {
    constexpr Fixed(int v): v(v << 16) {}
    constexpr Fixed(float f): v(f * (1 << 16)) {}
    constexpr Fixed(double f): v(f * (1 << 16)) {}
    constexpr Fixed(): v(0) {}

    static constexpr Fixed from_raw(int32_t x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    int32_t v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
};

Fixed operator+(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v + b.v);
}

Fixed operator-(Fixed a, Fixed b) {
    return Fixed::from_raw(a.v - b.v);
}

Fixed operator*(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v * b.v) >> 16);
}

Fixed operator/(Fixed a, Fixed b) {
    return Fixed::from_raw(((int64_t) a.v << 16) / b.v);
}

Fixed &operator+=(Fixed &a, Fixed b) {
    return a = a + b;
}

Fixed &operator-=(Fixed &a, Fixed b) {
    return a = a - b;
}

Fixed &operator*=(Fixed &a, Fixed b) {
    return a = a * b;
}

Fixed &operator/=(Fixed &a, Fixed b) {
    return a = a / b;
}

Fixed operator-(Fixed x) {
    return Fixed::from_raw(-x.v);
}

Fixed abs(Fixed x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}

ostream &operator<<(ostream &out, Fixed x) {
    return out << x.v / (double) (1 << 16);
}

struct VectorField {
    array<Fixed, deltas.size()> v[N][M];
    Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
        return get(x, y, dx, dy) += dv;
    }

    Fixed &get(int x, int y, int dx, int dy) {
        size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }
};

Fixed rho[256];
Fixed p[N][M]{}, old_p[N][M];
VectorField velocity{}, velocity_flow{};
int last_use[N][M]{};
int UT = 0;
mt19937 rnd(1337);
int dirs[N][M]{};

class FluidSimulator {
private:
    Fixed (&rho_ref)[256];
    Fixed (&p_ref)[N][M];
    Fixed (&old_p_ref)[N][M];
    VectorField& velocity_ref;
    VectorField& velocity_flow_ref;
    int (&last_use_ref)[N][M];
    int& UT_ref;
    mt19937& rnd_ref;
    int (&dirs_ref)[N][M];
    char (&field_ref)[N][M + 1];

    Fixed random01() {
        return Fixed::from_raw((rnd_ref() & ((1 << 16) - 1)));
    }

    Fixed move_prob(int x, int y) {
        Fixed sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field_ref[nx][ny] == '#' || last_use_ref[nx][ny] == UT_ref) {
                continue;
            }
            auto v = velocity_ref.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use_ref[x][y] = UT_ref - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<Fixed, deltas.size()> tres;
            Fixed sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field_ref[nx][ny] == '#' || last_use_ref[nx][ny] == UT_ref) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity_ref.get(x, y, dx, dy);
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

            Fixed p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity_ref.get(x, y, dx, dy) > 0 && field_ref[nx][ny] != '#' && last_use_ref[nx][ny] < UT_ref);

            ret = (last_use_ref[nx][ny] == UT_ref - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        
        last_use_ref[x][y] = UT_ref;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field_ref[nx][ny] != '#' && last_use_ref[nx][ny] < UT_ref - 1 && velocity_ref.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp(field_ref, p_ref, velocity_ref);
                pp.swap_with(x, y);
                pp.swap_with(nx, ny);
                pp.swap_with(x, y);
            }
        }
        return ret;
    }

    struct ParticleParams {
        char type;
        Fixed cur_p;
        array<Fixed, deltas.size()> v;
        char (&field)[N][M + 1];
        Fixed (&p)[N][M];
        VectorField& velocity;

        ParticleParams(char (&f)[N][M + 1], Fixed (&p_arr)[N][M], VectorField& vel)
            : field(f), p(p_arr), velocity(vel) {}

        void swap_with(int x, int y) {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field_ref[nx][ny] != '#' && last_use_ref[nx][ny] < UT_ref - 1 && velocity_ref.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use_ref[x][y] = UT_ref;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field_ref[nx][ny] == '#' || last_use_ref[nx][ny] == UT_ref || velocity_ref.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    tuple<Fixed, bool, pair<int, int>> propagate_flow(int x, int y, Fixed lim) {
        last_use_ref[x][y] = UT_ref - 1;
        Fixed ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field_ref[nx][ny] != '#' && last_use_ref[nx][ny] < UT_ref) {
                auto cap = velocity_ref.get(x, y, dx, dy);
                auto flow = velocity_flow_ref.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                auto vp = min(lim, cap - flow);
                if (last_use_ref[nx][ny] == UT_ref - 1) {
                    velocity_flow_ref.add(x, y, dx, dy, vp);
                    last_use_ref[x][y] = UT_ref;
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow_ref.add(x, y, dx, dy, t);
                    last_use_ref[x][y] = UT_ref;
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use_ref[x][y] = UT_ref;
        return {ret, 0, {0, 0}};
    }

public:
    FluidSimulator(
        Fixed (&rho_arr)[256],
        Fixed (&p_arr)[N][M],
        Fixed (&old_p_arr)[N][M],
        VectorField& velocity,
        VectorField& velocity_flow,
        int (&last_use_arr)[N][M],
        int& UT,
        mt19937& rnd,
        int (&dirs_arr)[N][M],
        char (&field_arr)[N][M + 1]
    ) : rho_ref(rho_arr),
        p_ref(p_arr),
        old_p_ref(old_p_arr),
        velocity_ref(velocity),
        velocity_flow_ref(velocity_flow),
        last_use_ref(last_use_arr),
        UT_ref(UT),
        rnd_ref(rnd),
        dirs_ref(dirs_arr),
        field_ref(field_arr)
    {
        initialize();
    }

    void initialize() {
        rho_ref[' '] = Fixed(0.01);
        rho_ref['.'] = Fixed(1000);

        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field_ref[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs_ref[x][y] += (field_ref[x + dx][y + dy] != '#');
                }
            }
        }
    }

    Fixed get_random01() {
        return random01();
    }

    Fixed get_move_prob(int x, int y) {
        return move_prob(x, y);
    }

    void do_propagate_move(int x, int y, bool is_first) {
        propagate_move(x, y, is_first);
    }

    void do_propagate_stop(int x, int y, bool force = false) {
        propagate_stop(x, y, force);
    }

    auto do_propagate_flow(int x, int y, Fixed lim) {
        return propagate_flow(x, y, lim);
    }

    friend int main();
};

int main() {
    FluidSimulator simulator(rho, p, old_p, velocity, velocity_flow, last_use, UT, rnd, dirs, field);
    
    Fixed g = 0.1;

    for (size_t i = 0; i < T; ++i) {
        
        Fixed total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                if (field[x + 1][y] != '#')
                    velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(old_p, p, sizeof(p));
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    int nx = x + dx, ny = y + dy;
                    if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                        auto delta_p = old_p[x][y] - old_p[nx][ny];
                        auto force = delta_p;
                        auto &contr = velocity.get(nx, ny, -dx, -dy);
                        if (contr * rho[(int) field[nx][ny]] >= force) {
                            contr -= force / rho[(int) field[nx][ny]];
                            continue;
                        }
                        force -= contr * rho[(int) field[nx][ny]];
                        contr = 0;
                        velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
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
            prop = 0;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        auto [t, local_prop, _] = simulator.do_propagate_flow(x, y, 1);
                        if (t > 0) {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    auto old_v = velocity.get(x, y, dx, dy);
                    auto new_v = velocity_flow.get(x, y, dx, dy);
                    if (old_v > 0) {
                        assert(new_v <= old_v);
                        velocity.get(x, y, dx, dy) = new_v;
                        auto force = (old_v - new_v) * rho[(int) field[x][y]];
                        if (field[x][y] == '.')
                            force *= 0.8;
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
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] != '#' && last_use[x][y] != UT) {
                    if (simulator.get_random01() < simulator.get_move_prob(x, y)) {
                        prop = true;
                        simulator.do_propagate_move(x, y, true);
                    } else {
                        simulator.do_propagate_stop(x, y, true);
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
