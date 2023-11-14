#include <bits/stdc++.h>
int read(void)
{
	int x = 0;
	char c = getchar();
	while (c < '0' or c > '9')
		c = getchar();
	while (c >= '0' and c <= '9')
		x = x * 10 + (c ^ 48), c = getchar();
	return x;
}
const int maxn = 1 << 20;
int n, a[maxn];
void solve(void)
{
	int f = 1, c = 0;
	for (int n = read(); n--;)
		read() > 1 ? (f = 0) : c++;
	puts(!f ? "499122177" : (c & 1 ? "1" : "0"));
}
signed main()
{
	for (int T = read(); T--;)
		solve();
	return 0;
}