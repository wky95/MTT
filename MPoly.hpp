#pragma once

#include "MTT.hpp"

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
    inline void extend(const MPoly &b){
		a.insert(a.end(), b.a.begin(), b.a.end());
	}
    MPoly& operator+=(const MPoly &b){
		if (ssize(a) < ssize(b.a)) resize(ssize(b.a));
		for (int i = 0; i < ssize(b.a); i++) a[i] += b.a[i];
		return *this;
	}
    MPoly& operator+=(const int x) {
        a[0] += x;
		return *this;
	}
    MPoly operator - () const {
		MPoly res(*this);
		for (auto &i : res.a) i = -i;
		return res;
	}
    MPoly& operator-=(const MPoly &b){
		if (ssize(a) < ssize(b.a)) resize(ssize(b.a));
		for (int i = 0; i < ssize(b.a); i++) a[i] -= b.a[i];
		return *this;
	}
    MPoly& operator-=(const int x) {
        a[0] -= x;
		return *this;
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
	make_int_operator(+) make_int_operator(-) make_int_operator(*) make_int_operator(/)
#undef make_mpoly_operator
#undef make_int_operator
    MPoly inv() {
        int sz = ssize(a);
        if (sz == 1) return MPoly({a[0].inverse()});
        int len = 2 << __lg(sz - 1);
        MPoly g(1); g.a[0] = a[0].inverse();
        for (int l = 1; l < len; l <<= 1) {
            (g *= 2) -= MPoly(a, l << 1).mulFix(g * g);
        }
        g.resize(sz);
        return g;
    }
    pair<MPoly, MPoly> divide(MPoly b) {
        // a = bQ + R, O(NlogN), b.back() != 0
        int n = ssize(a), m = ssize(b.a), k = n - m + 1;
        MPoly Q, R(*this);
        if (n >= m) {
            MPoly ra(vector<mint>(rbegin(a), rend(a)), k); 
            MPoly rb(vector<mint>(rbegin(b.a), rend(b.a)), k);
            Q = ra.mulFixconst(rb.inv()); reverse(all(Q.a));
            R -= b * Q;
        }
        while (ssize(Q.a) and Q.a.back() == 0) Q.a.pop_back();
        while (ssize(R.a) and R.a.back() == 0) R.a.pop_back();
        return {Q, R};
    }
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
    MPoly log() { return (this->inv().mulFix(this->diff())).integ(); }
    MPoly exp() {
        int len = (ssize(a) == 1 ? 1 : 2 << __lg(ssize(a) - 1));
        MPoly g(1); g[0] = 1;
        for (int l = 1; l < len; l <<= 1) {
            g.resize(l << 1);
            g = ((MPoly(a, l << 1) + 1) - g.log()).mulFixconst(g);
        }
        g.resize(ssize(a));
        return g;
    }
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