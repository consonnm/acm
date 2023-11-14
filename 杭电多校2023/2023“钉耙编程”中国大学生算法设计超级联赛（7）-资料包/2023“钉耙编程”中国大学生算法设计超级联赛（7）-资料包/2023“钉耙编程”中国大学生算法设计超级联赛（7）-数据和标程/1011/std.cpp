#include <bits/stdc++.h>
#define int long long
using namespace std;
inline int read()
{
	int t = 0;
	char v = getchar();
	while (v < '0')
		v = getchar();
	while (v >= '0')
		t = (t << 3) + (t << 1) + v - 48, v = getchar();
	return t;
}
int t, X, A, B, ans;
inline int Sqrt(int x)
{
	int o = sqrt(x);
	while ((o + 1) * (o + 1) <= x)
		++o;
	while (o * o > x)
		--o;
	return o;
}
signed main()
{
	t = read();
	while (t--)
	{
		X = read(), A = read(), B = read(), ans = 0;
		while (X && ans <= 100)
			X = min(X - 1, min(X + A >> 1, Sqrt(X + B))), ++ans;
		printf("%lld\n", ans + X);
	}
}
