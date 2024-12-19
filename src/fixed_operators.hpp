#pragma once

#include "fixed.hpp"
#include "fast_fixed.hpp"

template<size_t N, size_t K>
Fixed<N,K> operator+(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a + static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator+(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) + b;
}

template<size_t N, size_t K>
Fixed<N,K> operator-(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a - static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator-(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) - b;
}

template<size_t N, size_t K>
Fixed<N,K> operator*(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a * static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator*(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) * b;
}

template<size_t N, size_t K>
Fixed<N,K> operator/(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a / static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K> operator/(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) / b;
}

template<size_t N, size_t K>
bool operator<(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a < static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
bool operator<(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) < b;
}

template<size_t N, size_t K>
bool operator<=(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a <= static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
bool operator<=(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) <= b;
}

template<size_t N, size_t K>
bool operator>(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a > static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
bool operator>(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) > b;
}

template<size_t N, size_t K>
bool operator>=(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a >= static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
bool operator>=(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) >= b;
}

template<size_t N, size_t K>
bool operator==(const Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a == static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
bool operator==(const FastFixed<N,K>& a, const Fixed<N,K>& b) {
    return static_cast<Fixed<N,K>>(a) == b;
}

template<size_t N, size_t K>
Fixed<N,K>& operator+=(Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a += static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K>& operator-=(Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a -= static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K>& operator*=(Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a *= static_cast<Fixed<N,K>>(b);
}

template<size_t N, size_t K>
Fixed<N,K>& operator/=(Fixed<N,K>& a, const FastFixed<N,K>& b) {
    return a /= static_cast<Fixed<N,K>>(b);
}