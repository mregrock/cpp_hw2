#pragma once

#include "../types/types.h"
#include <ranges>
#include <cassert>
#include <array>
#include <utility>

template<typename NumericType, size_t N, size_t M>
struct VectorField {
    static constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};
    SimulationMatrix<std::array<NumericType, 4>, N, M> v{};

    NumericType &add(int x, int y, int dx, int dy, NumericType dv) {
        return get(x, y, dx, dy) += dv;
    }

    NumericType &get(int x, int y, int dx, int dy) {
        size_t i = std::ranges::find(deltas, std::pair(dx, dy)) - deltas.begin();
        assert(i < deltas.size());
        return v[x][y][i];
    }

    void resize(size_t n, size_t m) {
        v.resize(n);
        for (auto& row : v) {
            row.assign(m, std::array<NumericType, 4>{});
        }
    }
};