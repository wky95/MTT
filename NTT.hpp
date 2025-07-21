#pragma once

#include "ModInt.hpp"

template<uint32_t mod>
struct NTT {
    using i32 = int32_t;
    using i64 = int64_t;
    using u32 = uint32_t;
    using u64 = uint64_t;
    using mint = MontgomeryModInt<mod>;
    const mint g = 3;
    constexpr void ntt(vector<mint> &a, int n, int op) {
        vector<int> R(n);
        for (int i = 0; i < n; i++) R[i] = R[i/2]/2 + (i&1)*(n/2);
        for (int i = 0; i < n; i++)
            if (i < R[i]) swap(a[i], a[R[i]]);
        for (int m = 2; m <= n; m <<= 1) {
            mint w1 = g.pow((mod-1)/m * op + mod-1);
            for (int i = 0; i < n; i += m) {
                mint wk = 1;
                for (int k = 0; k < m / 2; k++) {
                    mint x = a[i+k], y = a[i+k+m/2] * wk;
                    a[i+k] = x + y;
                    a[i+k+m/2] = x - y;
                    wk *= w1;
                }
            }
        }
        if (op == -1) {
            for (int i = 0; i < n; i++) a[i] /= n;
        }
    }
    template<class RetType>
    constexpr vector<RetType> mult(const auto &a0, const auto &b0) {
        vector<mint> a(begin(a0), end(a0)), b(begin(b0), end(b0));
        int n = ssize(a), m = ssize(b), len = n+m==2 ? 1 : 2<<__lg(n+m-2);
        a.resize(len), b.resize(len);
        vector<mint> c(len);
        ntt(a, len, 1); ntt(b, len, 1);
        for (int i = 0; i < len; i++) c[i] = a[i] * b[i];
        ntt(c, len, -1);
        c.resize(n + m - 1);
        return vector<RetType>(begin(c), end(c));
    }
};