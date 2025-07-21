#pragma once

#include "ModInt.hpp"

template<const uint32_t mod>
struct Binomial {
    using mint = MontgomeryModInt<mod>;
    vector<mint> _fac, _facInv;
    Binomial(int size) : _fac(size), _facInv(size) {
        _fac[0] = 1;
        for(int i = 1; i < size; i++)
            _fac[i] = _fac[i - 1] * i;
        if (size > 0)
            _facInv.back() = mint(1) / _fac.back();
        for(int i = size - 2; i >= 0; i--)
            _facInv[i] = _facInv[i + 1] * (i + 1);
    }
    mint fac(int i) { return i < 0 ? 0 : _fac[i]; }
    mint faci(int i) { return i < 0 ? 0 : _facInv[i]; }
    mint inv(int i) { return _facInv[i] * _fac[i - 1]; }
    mint binom(int n, int r) { return r < 0 or n < r ? 0 : fac(n) * faci(r) * faci(n - r); }
    mint catalan(int i) { return binom(2 * i, i) - binom(2 * i, i + 1); }
    mint excatalan(int n, int m, int k) { //(+1) * n, (-1) * m, prefix sum > -k
        if (k > m) return binom(n + m, m);
        else if (k > m - n) return binom(n + m, m) - binom(n + m, m - k);
        else return mint(0);
    }
};
