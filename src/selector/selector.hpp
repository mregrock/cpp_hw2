#pragma once

#include <string>
#include <type_traits>

#include "../types/types.h"
#include "selector_helpers.h"
#include "simulation_helpers.h"

template<typename AllowedTypes, typename AllowedSizes, typename SelectedTypes>
struct TypeSelector {
    template<typename... Selected>
    static bool try_combinations(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                               size_t n, size_t m, const std::string& filename) {
        std::cerr << "Entering try_combinations with sizes n=" << n << ", m=" << m << ", filename=" << filename << std::endl;
        return try_all_p_types<0>(p_type, v_type, v_flow_type, n, m, filename);
    }

private:
    template<size_t I>
    static bool try_all_p_types(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                               size_t n, size_t m, const std::string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using P = typename AllowedTypes::template type_at<I>;
            return try_with_p_type<P>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_p_types<I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P>
    static bool try_with_p_type(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                               size_t n, size_t m, const std::string& filename) {
        if (!matches_type<P>(p_type)) return false;
        return try_all_v_types<P, 0>(p_type, v_type, v_flow_type, n, m, filename);
    }

    template<typename P, size_t I>
    static bool try_all_v_types(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                               size_t n, size_t m, const std::string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using V = typename AllowedTypes::template type_at<I>;
            return try_with_v_type<P, V>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_v_types<P, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V>
    static bool try_with_v_type(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                               size_t n, size_t m, const std::string& filename) {
        if (!matches_type<V>(v_type)) return false;
        return try_all_vf_types<P, V, 0>(p_type, v_type, v_flow_type, n, m, filename);
    }

    template<typename P, typename V, size_t I>
    static bool try_all_vf_types(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                                size_t n, size_t m, const std::string& filename) {
        if constexpr (I >= AllowedTypes::size) {
            return false;
        } else {
            using VF = typename AllowedTypes::template type_at<I>;
            return try_with_vf_type<P, V, VF>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_vf_types<P, V, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V, typename VF>
    static bool try_with_vf_type(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                                size_t n, size_t m, const std::string& filename) {
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
    static bool try_all_sizes(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                             size_t n, size_t m, const std::string& filename) {
        if constexpr (I >= AllowedSizes::size) {
            return false;
        } else {
            using CurrentSize = typename AllowedSizes::template size_at<I>;
            return try_with_size<P, V, VF, CurrentSize>(p_type, v_type, v_flow_type, n, m, filename) ||
                   try_all_sizes<P, V, VF, I + 1>(p_type, v_type, v_flow_type, n, m, filename);
        }
    }

    template<typename P, typename V, typename VF, typename CurrentSize>
    static bool try_with_size(const std::string& p_type, const std::string& v_type, const std::string& v_flow_type,
                             size_t n, size_t m, const std::string& filename) {
        if (CurrentSize::n != n || CurrentSize::m != m) return false;
        
        std::cerr << "Found matching size: " << CurrentSize::n << "x" << CurrentSize::m << std::endl;
        run_simulation<P, V, VF, CurrentSize>(CurrentSize::n, CurrentSize::m, filename);
        return true;
    }
};