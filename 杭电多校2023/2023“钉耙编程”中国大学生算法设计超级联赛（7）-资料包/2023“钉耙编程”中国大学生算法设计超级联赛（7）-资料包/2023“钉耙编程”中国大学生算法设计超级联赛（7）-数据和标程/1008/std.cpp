#include <bits/stdc++.h>
using namespace std;
const int N = 100010;
int t, n, u, z, x, y, cnt, f[N], cha[N], vis[N], len[N], lg[N], st1[N][20], st2[N][20];
vector<int> v[N];
void dfs(int x, int fa)
{
	f[x] = fa;
	for (auto y : v[x])
	{
		if (y != fa)
		{
			dfs(y, x);
		}
	}
}
int dfs2(int x, int fa)
{
	int res = 0;
	for (auto y : v[x])
	{
		if (!vis[y] && y != fa)
		{
			res = max(res, dfs2(y, x) + 1);
		}
	}
	return res;
}
int query1(int l, int r)
{
	int k = lg[r - l + 1];
	return max(st1[l][k], st1[r - (1 << k) + 1][k]);
}
int query2(int l, int r)
{
	int k = lg[r - l + 1];
	return max(st2[l][k], st2[r - (1 << k) + 1][k]);
}
int main()
{
	scanf("%d", &t);
	for (int ii = 1; ii <= t; ii++)
	{
		scanf("%d%d%d", &n, &x, &y);
		lg[0] = -1, cnt = 0;
		for (int i = 1; i <= n; i++)
		{
			v[i].clear();
			vis[i] = 0;
			lg[i] = lg[i >> 1] + 1;
		}
		for (int i = 1; i <= n - 1; i++)
		{
			scanf("%d%d", &u, &z);
			v[u].push_back(z), v[z].push_back(u);
		}
		dfs(x, 0);
		while (y)
		{
			cha[++cnt] = y;
			vis[y] = 1;
			y = f[y];
		}
		for (int i = 1; i <= cnt; i++)
		{
			len[i] = dfs2(cha[i], 0);
			st1[i][0] = len[i] + i, st2[i][0] = len[i] - i;
		}
		for (int i = 1; i <= 19; i++)
		{
			for (int j = 1; j + (1 << i) - 1 <= n; j++)
			{
				st1[j][i] = max(st1[j][i - 1], st1[j + (1 << (i - 1))][i - 1]);
				st2[j][i] = max(st2[j][i - 1], st2[j + (1 << (i - 1))][i - 1]);
			}
		}
		int ca = cnt, cb = 1, cur = 1;
		while (1)
		{
			if (cur == 1)
			{
				if (ca - 1 == cb)
				{
					if (len[ca] > len[cb])
					{
						printf("1\n");
						break;
					}
					else
					{
						printf("0\n");
						break;
					}
				}
				if (len[ca] > query1(cb, ca - 1) - cb)
				{
					printf("1\n");
					break;
				}
				ca--;
			}
			else
			{
				if (cb + 1 == ca)
				{
					if (len[cb] > len[ca])
					{
						printf("0\n");
						break;
					}
					else
					{
						printf("1\n");
						break;
					}
				}
				if (len[cb] > query2(cb + 1, ca) + ca)
				{
					printf("0\n");
					break;
				}
				cb++;
			}
			cur ^= 1;
		}
	}
	return 0;
}
