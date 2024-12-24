#pragma once

#include "fixed.hpp"
#include "fast_fixed.hpp"

template<typename T1, typename T2>
struct MorePrecise;

template<size_t N1, size_t K1, size_t N2, size_t K2>
struct MorePrecise<Fixed<N1,K1>, Fixed<N2,K2>> {
    using type = std::conditional_t<
        (N1 > N2) || (N1 == N2 && K1 >= K2),
        Fixed<N1,K1>,
        Fixed<N2,K2>
    >;
};

template<size_t N1, size_t K1, size_t N2, size_t K2>
struct MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>> {
    using type = std::conditional_t<
        (N1 > N2) || (N1 == N2 && K1 >= K2),
        Fixed<N1,K1>,
        FastFixed<N2,K2>
    >;
};

template<size_t N1, size_t K1, size_t N2, size_t K2>
struct MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>> {
    using type = std::conditional_t<
        (N1 > N2) || (N1 == N2 && K1 >= K2),
        FastFixed<N1,K1>,
        Fixed<N2,K2>
    >;
};

template<size_t N1, size_t K1, size_t N2, size_t K2>
struct MorePrecise<FastFixed<N1,K1>, FastFixed<N2,K2>> {
    using type = std::conditional_t<
        (N1 > N2) || (N1 == N2 && K1 >= K2),
        FastFixed<N1,K1>,
        FastFixed<N2,K2>
    >;
};

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator-(const Fixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) - ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator-(const FastFixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) - ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator+(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) + ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator+(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) + ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator-(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) - ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator-(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) - ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator*(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) * ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator*(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) * ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator/(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) / ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator/(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) / ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) < ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) < ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<=(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) <= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<=(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) <= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<=(const FastFixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) <= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator<=(const Fixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) <= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator>(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) > ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator>(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) > ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator>=(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) >= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator>=(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) >= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator==(const Fixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator==(const FastFixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator==(const Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
bool operator==(const FastFixed<N1,K1>& a, const Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
Fixed<N1,K1>& operator+=(Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return a += ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
Fixed<N1,K1>& operator-=(Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return a -= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
Fixed<N1,K1>& operator*=(Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return a *= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
Fixed<N1,K1>& operator/=(Fixed<N1,K1>& a, const FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return a /= ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator-(Fixed<N1,K1>& a, FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) - ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator==(Fixed<N1,K1>& a, FastFixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<Fixed<N1,K1>, FastFixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}

template<size_t N1, size_t K1, size_t N2, size_t K2>
auto operator==(FastFixed<N1,K1>& a, Fixed<N2,K2>& b) {
    using ResultType = typename MorePrecise<FastFixed<N1,K1>, Fixed<N2,K2>>::type;
    return ResultType(a) == ResultType(b);
}