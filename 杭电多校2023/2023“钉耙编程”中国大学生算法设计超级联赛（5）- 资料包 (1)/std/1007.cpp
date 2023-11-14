#include <bits/stdc++.h>
using i64 = int64_t;
using u32 = uint32_t;
using u64 = uint64_t;
const int N = 1 << 20;
const int P = 998244353;
struct mint {
    int x;
    constexpr mint(int x = 0) : x(x) {}
    mint operator-() const { return x > 0 ? P - x : 0; }
    mint operator+(mint o) const { return x + o.x < P ? x + o.x : x + o.x - P; }
    mint operator-(mint o) const { return x - o.x < 0 ? x - o.x + P : x - o.x; }
    mint operator*(mint o) const { return int(u64(x) * o.x % P); }
    mint &operator+=(mint o) { return *this = *this + o; }
    mint &operator-=(mint o) { return *this = *this - o; }
    mint &operator*=(mint o) { return *this = *this * o; }
    mint inv() const { return pow(P - 2); }
    mint cbrt() const { return pow((P + P - 1) / 3); }
    mint pow(auto k) const {
        mint a = x;
        mint b = 1;
        for (; k; k >>= 1) {
            if (k & 1)
                b *= a;
            a *= a;
        }
        return b;
    }
    mint sqrt() const {
        if (pow(P >> 1).x != 1) return 0;
        mint a = pow(60);
        mint b = pow(119);
        for (int k = 21; k >= 0; --k)
            if (b.pow(1 << k).x != 1) {
                a *= mint(3).pow(P >> (k + 2));
                b *= mint(3).pow(P >> (k + 1));
            }
        return std::min(a.x, P - a.x);
    }
};
mint w[N];
mint fac[N];
mint inv[N];
mint ivv[N];

void __attribute__((constructor)) init() {
    w[N / 2] = 1;
    mint g = mint(3).pow(P / N);
    for (int i = N / 2 + 1; i < N; ++i) w[i] = w[i - 1] * g;
    for (int i = N / 2 - 1; i > 0; --i) w[i] = w[i << 1];
    fac[0] = fac[1] = 1;
    inv[0] = inv[1] = 1;
    ivv[0] = ivv[1] = 1;
    for (int i = 2; i < N; ++i) fac[i] = fac[i - 1] * i;
    for (int i = 2; i < N; ++i) inv[i] = inv[P % i] * (P - P / i);
    for (int i = 2; i < N; ++i) ivv[i] = ivv[i - 1] * inv[i];
}
int main() {
#ifdef LOCAL
    auto flush = [&]() {};
    auto ii = [&]() {
        int x;
        std::cin >> x;
        return x;
    };
    auto oo = [&](auto x, char c = 10) {
        std::cout << x << c << std::flush;
    };
#else
    char bufI[1 << 19], *ptrI = bufI, *endI = bufI + sizeof(bufI);
    char bufO[1 << 19], *ptrO = bufO, *endO = bufO + sizeof(bufO);
    fread(bufI, 1, sizeof(bufI), stdin);
    auto load = [&]() {
        memcpy(bufI, ptrI, endI - ptrI);
        fread(endI - ptrI + bufI, 1, ptrI - bufI, stdin);
        ptrI = bufI;
    };
    auto flush = [&]() {
        fwrite(bufO, 1, ptrO - bufO, stdout);
        ptrO = bufO;
    };
    auto ii = [&]() {
        if (endI - ptrI < 32) load();
        int x{};
        int n{};
        for (; *ptrI < 48; ++ptrI) n = *ptrI == 45;
        for (; *ptrI > 47; ++ptrI) x = x * 10 + *ptrI - 48;
        return n ? -x : +x;
    };
    auto oo = [&](auto x, char c = 10) {
        if (endO - ptrO < 32) flush();
        if (x < 0) x = -x, *ptrO++ = '-';
        char buf[20];
        char *end = buf + 20;
        char *ptr = buf + 20;
        *--ptr = c;
        for (; x >= 10; x /= 10)
            *--ptr = char(48 + x % 10);
        *--ptr = char(48 + x);
        memcpy(ptrO, ptr, end - ptr);
        ptrO += end - ptr;
    };
#endif
    int T = ii();
    while(T--) {
        int n = ii();
        int m = ii();
        mint a = ii();
        mint b = ii();
        mint p = a * b.inv();
        mint q = -p + 1;
        mint r = p * q.inv();
        mint ans;
        mint sum;
        mint tmp = q.pow(n);
        for (int i = 1; i <= n; ++i) {
            sum += mint(i).pow(m);
            tmp *= r * inv[i] * (n - i + 1);
            ans += tmp * sum;
        }
        oo(ans.x);
    }
    flush();
}
