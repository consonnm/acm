#include <bits/stdc++.h>
using i64 = int64_t;
using u32 = uint32_t;
using u64 = uint64_t;
const int N = 1 << 20;

int ch[N];
int dp[N];
int node[N];
int head[N];
int next[N];

void solve(int u, int w) {
    node[u] = 1;
    for (int v = head[u]; v; v = next[v]) {
        if (v == ch[u]) {
            solve(v, w);
            node[u] = node[v] + 2;
        } else {
            solve(v, v);
            dp[w] = std::max(dp[w], dp[v] + std::__lg(node[v]) + 1);
        }
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    auto ii = [&]() {
        int x;
        std::cin >> x;
        return x;
    };
    auto oo = [&](auto x, char c = 10) {
        std::cout << x << c << std::flush;
    };
    int T = ii();
    while(T--) {
        int n = ii();
        for (int i = 1; i <= n; ++i) {
            int f = ii();
            next[i] = head[f];
            head[f] = i;
        }
        for (int i = 1; i <= n; ++i) ch[i] = ii();
        solve(1, 1);
        oo(dp[1] + std::__lg(node[1]) + 1);
        for (int i = 0; i <= n; ++i) ch[i] = dp[i] = node[i] = head[i] = next[i] = 0;
    }
}
