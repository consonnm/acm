#include <bits/stdc++.h>
using namespace std;
const int N = 100010;
int t, n, a[N], b[N];
void modify(int x, int v)
{
	for (int i = x; i <= n; i += (i & (-i)))
	{
		b[i] += v;
	}
}
int query(int x)
{
	int res = 0;
	for (int i = x; i; i -= (i & (-i)))
	{
		res += b[i];
	}
	return res;
}
int main()
{
	scanf("%d", &t);
	for (int ii = 1; ii <= t; ii++)
	{
		scanf("%d", &n);
		if (n == 1)
		{
			scanf("%d", &a[1]);
			printf("0 1\n");
			continue;
		}
		int ans1 = 0, ans2 = 1, inv = 0;
		while (ans2 <= n)
		{
			ans2 <<= 1;
		}
		ans2--;
		ans2 ^= 2;
		for (int i = 1; i <= n; i++)
		{
			b[i] = 0;
		}
		for (int i = 1; i <= n; i++)
		{
			scanf("%d", &a[i]);
			int tmp = query(n) - query(a[i]);
			if (tmp & 1)
			{
				inv ^= 1;
			}
			modify(a[i], 1);
		}
		if (inv)
		{
			ans1 ^= 2, ans2 ^= 2;
		}
		printf("%d %d\n", ans1, ans2);
	}
	return 0;
}
