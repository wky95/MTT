#pragma once

#include "ModInt.hpp"
#include "NTT.hpp"

template<const uint32_t mod>
struct MTT {
    using i32 = int32_t;
    using i64 = int64_t;
    using u32 = uint32_t;
    using u64 = uint64_t;
    inline static constexpr u64 M0 = 1004535809;
    inline static constexpr u64 M1 = 104857601;
    inline static constexpr u64 M2 = 998244353;
    using mint0 = MontgomeryModInt<M0>;
    using mint1 = MontgomeryModInt<M1>;
    using mint2 = MontgomeryModInt<M2>;
    using mint = MontgomeryModInt<mod>;
    inline static NTT<M0> ntt0;
    inline static NTT<M1> ntt1;
    inline static NTT<M2> ntt2;
    static constexpr pair<i64, i64> exgcd(i64 a, i64 b) {
        if (b == 0) return {1, 0};
        auto [x1, y1] = exgcd(b, a % b);
        return {y1, x1 - a / b * y1};
    }
    static constexpr auto mult(const auto &a, const auto &b) {
        if (mod == M0) return ntt0.mult<mint>(a, b);
        if (mod == M1) return ntt1.mult<mint>(a, b);
        if (mod == M2) return ntt2.mult<mint>(a, b);
        int n = ssize(a), m = ssize(b);
        vector<mint0> c0 = ntt0.mult<mint0>(a, b);
        vector<mint1> c1 = ntt1.mult<mint1>(a, b);
        vector<mint2> c2 = ntt2.mult<mint2>(a, b);
        vector<mint> res(n + m - 1);
        i64 P1 = M0 * M1;
        for (int i = 0; i < n + m - 1; i++) {
            auto [k0, k1] = exgcd(M0, M1);
            k0 = ((c1[i] - i64(c0[i].get())) * k0).get();
            i64 x1 = (i64(c0[i].get()) + k0 * M0 % P1) % P1;
            tie(k0, k1) = exgcd(P1, M2);
            k0 = ((c2[i] - x1) * k0).get();
            res[i] = (x1 + k0 % mod * (P1 % mod)) % mod;
        }
        return res;
    }
};
