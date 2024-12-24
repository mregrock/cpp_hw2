#pragma once

#include <tuple>

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