#include <bits/stdc++.h>
typedef long long int64;
typedef __int128_t int128;
std::mt19937_64 rnd(time(0));
constexpr int Pr[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
int64 power(int64 x, int64 y, int64 m)
{
	int64 z = 1;
	for (; y; y >>= 1, x = static_cast<int128>(x) * x % m)
		if (y & 1)
			z = static_cast<int128>(z) * x % m;
	return z;
}
bool MillerRabin(int64 n)
{
	for (int p : Pr)
	{
		if (n == p)
			return true;
		if (n % p == 0 or power(p % n, n - 1, n) != 1)
			return false;
		for (int64 x = n - 1; !(x & 1);)
		{
			int64 r = power(p % n, x >>= 1, n);
			if (r != 1 and r != n - 1)
				return false;
			if (r == n - 1)
				break;
		}
	}
	return true;
}
int64 find(int64 n)
{
	int64 s = 0, t = 0, c = rnd() % (n - 1) + 1;
	for (int step = 1;; step <<= 1, s = t)
	{
		int64 r = 1;
		for (int i = 1; i <= step; i++)
		{
			t = (static_cast<int128>(t) * t + c) % n;
			r = static_cast<int128>(r) * std::abs(s - t) % n;
			if (i % 127 == 0)
				if (int64 d = std::gcd(r, n); d > 1)
					return d;
		}
		if (int64 d = std::gcd(r, n); d > 1)
			return d;
	}
	return -1;
}
void PollardRho(int64 n, std::vector<int64> &pr)
{
	if (n == 1)
		return;
	if (MillerRabin(n))
	{
		pr.emplace_back(n);
		return;
	}
	int64 m = n;
	while (m == n)
		m = find(n);
	while (n % m == 0)
		n /= m;
	PollardRho(m, pr);
	PollardRho(n, pr);
}
constexpr int mod = 998244353, inv2 = (mod + 1) / 2;
struct Power
{
	int BS[(1 << 15) + 1], GS[(1 << 15) + 1];
	Power(int c)
	{
		BS[0] = 1;
		for (int i = 1; i <= (1 << 15); i++)
			BS[i] = static_cast<int64>(BS[i - 1]) * c % mod;
		GS[0] = 1;
		for (int i = 1; i <= (1 << 15); i++)
			GS[i] = static_cast<int64>(GS[i - 1]) * BS[1 << 15] % mod;
	}
	int operator[](long long n)
	{
		n %= mod - 1;
		return static_cast<int64>(GS[n >> 15]) * BS[n & ((1 << 15) - 1)] % mod;
	}
};
void solve(void)
{
	int64 n, m, ans = 0;
	scanf("%lld%lld", &n, &m);
	int c = (m % mod - 1 + mod) % mod;
	Power P(c);
	std::vector<int64> pr;
	PollardRho(n, pr);
	std::sort(pr.begin(), pr.end());
	pr.resize(std::unique(pr.begin(), pr.end()) - pr.begin());
	std::vector<int> vpn(pr.size()), vpm(pr.size());
	int64 _n = n, _m = m;
	for (int i = 0; i < (int)pr.size(); i++)
		while (_n % pr[i] == 0)
			_n /= pr[i], vpn[i]++;
	for (int i = 0; i < (int)pr.size(); i++)
		while (_m % pr[i] == 0)
			_m /= pr[i], vpm[i]++;
	auto F = [&](int64 n)
	{ return (P[n] + (n & 1 ? 1 : mod - 1)) % mod; };
	std::function<void(int, int64, int64, int64)> dfs = [&](int u, int64 g, int64 phi, int64 gm)
	{
		if (u == (int)pr.size())
		{
			phi %= mod;
			gm = (gm % mod - 1 + mod) % mod;
			ans = (ans + static_cast<int64>(F(n / g)) * phi % mod * gm + static_cast<int64>(F(n / g - 1)) * c % mod * phi) % mod;
			return;
		}
		for (int i = 0; i <= vpn[u]; i++)
			dfs(u + 1, g, phi, gm), g *= pr[u], phi *= i ? pr[u] : (pr[u] - 1), gm *= i < vpm[u] ? pr[u] : 1;
	};
	if (!(n & 1))
	{
		int n_div_2 = (n >> 1) % mod;
		ans = (ans + m % mod * P[n >> 1] % mod * n_div_2) % mod;
		if (!(m & 1))
			ans = (ans + m % mod * P[(n >> 1) - 1] % mod * n_div_2) % mod;
	}
	dfs(0, 1, 1, 1);
	printf("%lld\n", static_cast<int64>(ans) * power(m % mod, mod - 2, mod) % mod * power(n % mod, mod - 2, mod) % mod * inv2 % mod);
}
signed main()
{
	int T;
	for (scanf("%d", &T); T--;)
		solve();
	return 0;
}