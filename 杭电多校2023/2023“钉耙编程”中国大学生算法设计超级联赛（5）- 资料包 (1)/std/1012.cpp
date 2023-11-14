#include <bits/stdc++.h>
using namespace std;

inline int rd() {
  int x = 0;
  bool f = 0;
  char c = getchar();
  for (; !isdigit(c); c = getchar()) f |= (c == '-');
  for (; isdigit(c); c = getchar()) x = x * 10 + (c ^ 48);
  return f ? -x : x;
}

#define N 1000007
#define mod 1000000007

int fpow(int x, int t = mod - 2) {
	int res = 1;
	for (; t; t >>= 1, x = 1ll * x * x % mod)
		if (t & 1) res = 1ll * res * x % mod;
	return res;
}

int fac[N], ifac[N], deg[N], ans[N];

int C(int n, int m) {
	if (n < m) return 0;
	return 1ll * fac[n] * ifac[m] % mod * ifac[n - m] % mod;
}

void work() {
	int n, m; cin >> n >> m;
	for (int i = 1; i <= n; ++i) deg[i] = ans[i] = 0;
	for (int i = 1; i <= m; ++i) {
      	int x; cin >> x; ++deg[x];
      cin >> x; ++deg[x];
	}
	for (int i = 1; i <= n; ++i) 
		for (int j = 2; j <= deg[i]; ++j) 
			ans[j] = (ans[j] + C(deg[i], j)) % mod;
	int res = 0;
	for (int i = 2; i <= n; ++i) res ^= ans[i];
	printf("%d\n", res);
}

int main() {
  ios::sync_with_stdio(false);
  	cin.tie(0);
  	
	fac[0] = ifac[0] = 1;
	for (int i = 1; i < N; ++i) fac[i] = 1ll * fac[i - 1] * i % mod;
	ifac[N - 1] = fpow(fac[N - 1]);
	for (int i = N - 2; i; --i) ifac[i] = 1ll * ifac[i + 1] * (i + 1) % mod;

	int t; cin >> t;
	while (t--) work();
	return 0;
}