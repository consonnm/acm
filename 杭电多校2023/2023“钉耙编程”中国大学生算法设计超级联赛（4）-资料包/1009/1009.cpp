#include <bits/stdc++.h>
using namespace std;
 
typedef long long i64;
 
constexpr int mod = 998244353;
 
i64 inv[200005] = {0, 1};
 
void solve() {
   int n;
   cin >> n;
   vector<vector<pair<int, int>>> g(n);
   for (int i = 0, u, v, w; i < n - 1; i++) {
      cin >> u >> v >> w;
      u--;
      v--;
      w--;
      g[u].emplace_back(v, w);
      g[v].emplace_back(u, w);
   }
   vector<int> sz(n, 1), dp(2*n), up(n);
   vector<vector<int>> stk(n);
   for (int i = n; i < 2*n; i++) {
      dp[i] = n;
   }
   function<void(int, int)> DFS = [&](int u, int f) {
      for (auto &[v, w] : g[u]) {
         if (v == f) continue;
         stk[w].push_back(v);
         DFS(v, u);
         stk[w].pop_back();
         sz[u] += sz[v];
         
         up[v] = stk[w].size() ? stk[w].back() : n + w;
         
         dp[up[v]] -= sz[v];
      }
      dp[u] += sz[u];
   };
   DFS(0, -1);
   
   i64 ans = 0;
   for (int i = 1; i < n; i++) {
      ans += 1ll * dp[i] * dp[up[i]] % mod;
      ans %= mod;
   }
   i64 res = 0;
   for (int k = 1; k <= n; k++) {
      res ^= (ans * k % mod * (k - 1) % mod * inv[n] % mod * inv[n - 1] % mod);
   }
   cout << res << '\n';
}
 
int main() {
   ios::sync_with_stdio(false);
   cin.tie(0);
   for (int i = 2; i <= 200000; i++) {
      inv[i] = (mod - mod / i) * inv[mod % i] % mod;
   }
   int t;
   cin >> t;
   while (t--) {
      solve();
   }
   return 0;
}
