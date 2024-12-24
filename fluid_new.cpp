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
#include "src/types/fixed.hpp"
#include "src/types/fast_fixed.hpp"
#include "src/types/fixed_operators.hpp"
#include "src/types/types.h"
#include "src/simulation/fluid_simulator.hpp"
#include "src/simulation/vector_field.hpp"

using namespace std;

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

vector<string> read_field(const string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Can't open file: " + filename);
    }
    vector<string> field;
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (!field.empty() && line.length() != field[0].length()) {
            throw std::runtime_error("Inconsistent line length in field");
        }
        field.push_back(line);
    }
    return field;
}

bool is_state_file(const string& filename) {
    return filename.ends_with(".bin");
}

template<typename P, typename V, typename VF, typename Size=SizeType<dynamic_size, dynamic_size>>
void run_simulation(size_t n = DEFAULT_N, size_t m = DEFAULT_M, const string& filename="field.txt") {  
    cerr << "Running simulation with n=" << n << ", m=" << m << ", filename=" << filename << endl;
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
                               size_t n, size_t m, const string& filename) {
        std::cerr << "Entering try_combinations with sizes n=" << n << ", m=" << m << ", filename=" << filename << std::endl;
        return try_all_p_types<0>(p_type, v_type, v_flow_type, n, m, filename);
    }

private:
    template<size_t I>
    static bool try_all_p_types(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m, const string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using P = typename AllowedTypes::template type_at<I>;
            return try_with_p_type<P>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_p_types<I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P>
    static bool try_with_p_type(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m, const string& filename) {
        if (!matches_type<P>(p_type)) return false;
        return try_all_v_types<P, 0>(p_type, v_type, v_flow_type, n, m, filename);
    }

    template<typename P, size_t I>
    static bool try_all_v_types(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m, const string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using V = typename AllowedTypes::template type_at<I>;
            return try_with_v_type<P, V>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_v_types<P, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V>
    static bool try_with_v_type(const string& p_type, const string& v_type, const string& v_flow_type,
                               size_t n, size_t m, const string& filename) {
        if (!matches_type<V>(v_type)) return false;
        return try_all_vf_types<P, V, 0>(p_type, v_type, v_flow_type, n, m, filename);
    }

    template<typename P, typename V, size_t I>
    static bool try_all_vf_types(const string& p_type, const string& v_type, const string& v_flow_type,
                                size_t n, size_t m, const string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using VF = typename AllowedTypes::template type_at<I>;
            return try_with_vf_type<P, V, VF>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_vf_types<P, V, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V, typename VF>
    static bool try_with_vf_type(const string& p_type, const string& v_type, const string& v_flow_type,
                                size_t n, size_t m, const string& filename) {
        if (!matches_type<VF>(v_flow_type)) return false;

        std::cerr << "Found matching types and about to create simulation" << std::endl;
        
        if constexpr (AllowedSizes::size == 0) {
            run_simulation<P, V, VF>(n, m, filename);
            return true;
        } else {
            return try_all_sizes<P, V, VF, 0>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }
    template<typename P, typename V, typename VF, size_t I>
    static bool try_all_sizes(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m, const string& filename) {
        if constexpr (I >= AllowedSizes::size) {
            return false;
        } else {
            using CurrentSize = typename AllowedSizes::template size_at<I>;
            return try_with_size<P, V, VF, CurrentSize>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_sizes<P, V, VF, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V, typename VF, typename CurrentSize>
    static bool try_with_size(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m, const string& filename) {
        if (CurrentSize::n != n || CurrentSize::m != m) return false;
        
        std::cerr << "Found matching size: " << CurrentSize::n << "x" << CurrentSize::m << std::endl;
        run_simulation<P, V, VF, CurrentSize>(CurrentSize::n, CurrentSize::m, filename);
        return true;
    }
};

template<typename... Types>
bool try_all_type_combinations(const string& p_type, const string& v_type, const string& v_flow_type,
                             size_t n, size_t m, const string& filename) {
    #define S(N, M) SizeType<N, M>
    return TypeSelector<TypesList<Types...>, SizesList<SIZES>, TypesList<>>::try_combinations(p_type, v_type, v_flow_type, n, m, filename);
}

bool create_and_run_simulation(const string& p_type, const string& v_type, const string& v_flow_type, 
                             size_t n, size_t m, const string& filename) {
    try {
        cerr << "\nTrying to create simulation with types:" << endl;
        cerr << "p_type: " << p_type << endl;
        cerr << "v_type: " << v_type << endl;
        cerr << "v_flow_type: " << v_flow_type << endl;
        #define FLOAT float
        #define DOUBLE double
        #define FIXED(N, K) Fixed<N, K>
        #define FAST_FIXED(N, K) FastFixed<N, K>
        if (!try_all_type_combinations<TYPES>(p_type, v_type, v_flow_type, n, m, filename)) {
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
    string filename = get_arg("--file", argc, argv, "field.txt");
    
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
    
    if (!create_and_run_simulation(p_type, v_type, v_flow_type, size.n, size.m, filename)) {
        cerr << "Failed to create simulation with types: " << p_type << ", " << v_type << ", " << v_flow_type << endl;
        return 1;
    }
    
    return 0;
}
