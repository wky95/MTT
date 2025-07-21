template<const uint32_t mod>
struct MontgomeryModInt {
    using mint = MontgomeryModInt;
    using i32 = int32_t;
    using u32 = uint32_t;
    using u64 = uint64_t;
    static constexpr u32 get_r() {
        u32 ret = mod;
        for (i32 i = 0; i < 4; ++i) ret *= 2 - mod * ret;
        return ret;
    }
    static constexpr u32 r = get_r();
    static constexpr u32 n2 = -u64(mod) % mod;
    static_assert(mod < (1 << 30), "invalid, mod >= 2 ^ 30");
    static_assert((mod & 1) == 1, "invalid, mod % 2 == 0");
    static_assert(r * mod == 1, "this code has bugs.");
    u32 a;
    constexpr MontgomeryModInt() : a(0) {}
    constexpr MontgomeryModInt(const int64_t &b)
        : a(reduce(u64(b % mod + mod) * n2)){};
    static constexpr u32 reduce(const u64 &b) {
        return (b + u64(u32(b) * u32(-r)) * mod) >> 32;
    }
    constexpr mint &operator+=(const mint &b) {
        if (i32(a += b.a - 2 * mod) < 0) a += 2 * mod;
        return *this;
    }
    constexpr mint &operator-=(const mint &b) {
        if (i32(a -= b.a) < 0) a += 2 * mod;
        return *this;
    }
    constexpr mint &operator*=(const mint &b) {
        a = reduce(u64(a) * b.a);
        return *this;
    }
    constexpr mint &operator/=(const mint &b) {
        *this *= b.inverse();
        return *this;
    }
    constexpr mint operator+(const mint &b) const { return mint(*this) += b; }
    constexpr mint operator-(const mint &b) const { return mint(*this) -= b; }
    constexpr mint operator*(const mint &b) const { return mint(*this) *= b; }
    constexpr mint operator/(const mint &b) const { return mint(*this) /= b; }
    constexpr bool operator==(const mint &b) const {
        return (a >= mod ? a - mod : a) == (b.a >= mod ? b.a - mod : b.a);
    }
    constexpr bool operator!=(const mint &b) const {
        return (a >= mod ? a - mod : a) != (b.a >= mod ? b.a - mod : b.a);
    }
    constexpr mint operator-() const { return mint() - mint(*this); }
    constexpr mint operator+() const { return mint(*this); }
    constexpr mint pow(u64 n) const {
        mint ret(1), mul(*this);
        for (; n > 0; mul = mul * mul, n >>= 1) 
            if (n & 1) ret *= mul;
        return ret;
    }
    constexpr mint inverse() const {
        int x = get(), y = mod, u = 1, v = 0, t = 0, tmp = 0;
        while (y > 0) {
            t = x / y;
            x -= t * y, u -= t * v;
            tmp = x, x = y, y = tmp;
            tmp = u, u = v, v = tmp;
        }
        return mint{u};
    }
    friend ostream &operator<<(ostream &os, const mint &b) {
        return os << b.get();
    }
    friend istream &operator>>(istream &is, mint &b) {
        int64_t t;
        is >> t;
        b = MontgomeryModInt<mod>(t);
        return (is);
    }
    constexpr u32 get() const {
        u32 ret = reduce(a);
        return ret >= mod ? ret - mod : ret;
    }
    static constexpr u32 get_mod() { return mod; }

    template<const uint32_t NewMod>
    constexpr MontgomeryModInt(const MontgomeryModInt<NewMod> &x)
      : a(reduce(u64(x.get()) * n2))
    {}
};

template<const uint32_t mod>
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

template<const uint32_t mod>
struct MPoly {
    using mint = MontgomeryModInt<mod>;
    MTT<mod> mtt;
    vector<mint> a;
    inline mint& operator [](const int x) { return a[x]; }
    void resize(const int x) { a.resize(x); }
    MPoly() {}
    MPoly(const int x) { resize(x); }
    MPoly(const MPoly &b) { a = b.a; }
    MPoly(const MPoly &b, const int x) { a = b.a; resize(x); }
    MPoly(const vector<mint> b) { a = b; }
    MPoly(const vector<mint> &b, const int x) { a = b; resize(x); }
    inline void extend(const MPoly &b){ a.insert(a.end(), b.a.begin(), b.a.end()); }
    MPoly& operator+=(const int x) { a[0] += x; return *this; }
    MPoly& operator-=(const int x) { a[0] -= x; return *this; }
    MPoly& operator+=(const MPoly &b){
        if (ssize(a) < ssize(b.a)) resize(ssize(b.a));
        for (int i = 0; i < ssize(b.a); i++) a[i] += b.a[i];
        return *this;
    }
    MPoly& operator-=(const MPoly &b){
		if (ssize(a) < ssize(b.a)) resize(ssize(b.a));
		for (int i = 0; i < ssize(b.a); i++) a[i] -= b.a[i];
		return *this;
	}
    MPoly operator - () const {
		MPoly res(*this);
		for (auto &i : res.a) i = -i;
		return res;
	}
    MPoly& operator*=(const int x){
		for (auto &i : a) i *= x;
		return *this;
	}
    MPoly& operator*=(const MPoly &b) {
		return *this = mtt.mult(a, b.a);
	}
    MPoly& operator/=(const int x){
		return *this *= mint(x).inverse();
	}
    MPoly& mulFix(const MPoly &b) {
        int sz = ssize(a);
		this->a = mtt.mult(a, b.a); resize(sz);
		return *this;
	}
    MPoly mulFixconst(const MPoly &b) const {
        MPoly res;
        res.a = mtt.mult(a, b.a);
        res.resize(ssize(a));
		return res;
	}
#define make_mpoly_operator(op) inline MPoly operator op (const MPoly &b) const { MPoly res(*this); return (res op##=b); }
#define make_int_operator(op)   inline MPoly operator op (const int x)    const { MPoly res(*this); return (res op##=x); }
	make_mpoly_operator(+) make_mpoly_operator(-) make_mpoly_operator(*)
	make_int_operator(+)   make_int_operator(-)   make_int_operator(*)   make_int_operator(/)
#undef make_mpoly_operator
#undef make_int_operator
    MPoly diff() {
		MPoly res(ssize(a) - 1);
        for (int i = 1; i < ssize(a); i++) res.a[i-1] = a[i] * i;
		return res;
	}
    MPoly integ() {
		MPoly res(ssize(a) + 1);
		for (int i = 1; i <= ssize(a); i++) res.a[i] = a[i-1] / i;
		return res;
	}
    MPoly inv() {
        MPoly g(1); g[0] = a[0].inverse();
        if (ssize(a) > 1) {
            for (int l = 1, len = 2 << __lg(ssize(a) - 1); l < len; l <<= 1) 
                (g *= 2) -= MPoly(a, l << 1).mulFix(g * g);
            g.resize(ssize(a));
        }
        return g;
    }
    MPoly exp() {
        MPoly g(1); g[0] = 1;
        if (ssize(a) > 1) {
            for (int l = 1, len = 2 << __lg(ssize(a) - 1); l < len; l <<= 1) {
                g.resize(l << 1);
                g = ((MPoly(a, l << 1) + 1) - g.log()).mulFixconst(g);
            }
            g.resize(ssize(a));
        }
        return g;
    }
    MPoly log() { return (this->inv().mulFix(this->diff())).integ(); }
    MPoly pow(int64_t n) {
        int i = find_if_not(begin(a), end(a), [](auto x){ return x == mint(0); }) - begin(a);
        if (i and (n >= ssize(a) or n * i >= ssize(a))) return MPoly<mod>(ssize(a));
        if (i == ssize(a)) return MPoly<mod>(ssize(a)) + 1;
        MPoly b(ssize(a) - i);
        for (int j = 0; j < ssize(b.a); j++) b.a[j] = a[j + i] / a[i];
        MPoly res1 = (b.log() * (int)(mint(n % mod).get())).exp() * (int)(a[i].pow(n).get());
        MPoly res2(ssize(a));
        for (int j = 0; j < min<int32_t>(ssize(res1.a), (ssize(a) - n * i)); j++) 
            res2.a[j + n * i] = res1.a[j];
        return res2;
    }
    pair<MPoly, MPoly> divide(MPoly b) {
        // a = bQ + R, O(NlogN), b.back() != 0
        int n = ssize(a), m = ssize(b.a), k = n - m + 1;
        MPoly Q, R(*this);
        if (n >= m) {
            MPoly ra(vector<mint>(rbegin(a), rend(a)), k); 
            MPoly rb(vector<mint>(rbegin(b.a), rend(b.a)), k);
            Q = ra.mulFixconst(rb.inv()); reverse(begin(Q.a), end(Q.a));
            R -= b * Q;
        }
        while (ssize(Q.a) and Q.a.back() == 0) Q.a.pop_back();
        while (ssize(R.a) and R.a.back() == 0) R.a.pop_back();
        return {Q, R};
    }
    MPoly evaluate(MPoly x) { // calculate for each i y_i = a(x_i)
        if (!ssize(x.a)) return {};
        int n = ssize(x.a);
        vector<MPoly> up(n * 2);
        for (int i = 0; i < n; i++) up[i + n] = MPoly(vector<mint>{-x[i], 1});
        for (int i = n - 1; i > 0; i--) up[i] = up[i*2] * up[i*2+1];
        vector<MPoly> down(n * 2);
        down[1] = this->divide(up[1]).second;
        for (int i = 2; i < n * 2; i++) down[i] = down[i >> 1].divide(up[i]).second;
        MPoly y(n);
        for (int i = 0; i < n; i++) y[i] = down[i + n].a.empty() ? 0 : down[i + n][0];
        return y;
    }
    MPoly interpolate(MPoly y) { // calculate f -> f(*this) = y
        int n = a.size();
        vector<MPoly> up(n * 2);
        for (int i = 0; i < n; i++) up[i + n] = MPoly(vector<mint>{-a[i], 1});
        for (int i = n - 1; i > 0; i--) up[i] = up[i*2] * up[i*2+1];
        auto b = up[1].diff().evaluate(a);
        for (int i = 0; i < n; i++) b[i] = y[i] / b[i];
        vector<MPoly> down(n * 2);
        for (int i = 0; i < n; i++) down[i + n] = MPoly(vector<mint>{b[i]});
        for (int i = n - 1; i > 0; i--) {
            down[i] = down[i*2]*up[i*2+1] + down[i*2+1]*up[i*2];
        }
        return down[1];
    }
};

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
