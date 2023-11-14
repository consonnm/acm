#include <bits/stdc++.h>
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
int n, m, fa[700002], Siz[700002], bl[100002], N[500002], c[100002], V[500002], col[700002], rt[100002], ls[15000002], rs[15000002], siz[15000002], pos[15000002], tot, num, A[100002], B[100002], C[100002], O[500002], X[500002], Y[500002], Z[500002];
vector<int> G[100002], W[100002], S[100002], Q[500002];
inline void add(int x, int y, int z) { G[x].push_back(y), W[x].push_back(z); }
inline int root(int x) { return x == fa[x] ? x : fa[x] = root(fa[x]); }
inline void ins(int &p, int l, int r, int x, int y)
{
	if (!p)
		p = ++tot;
	siz[p] += Siz[y];
	if (l == r)
	{
		if (!pos[p])
			pos[p] = ++num, fa[num] = num, col[num] = l;
		fa[y] = pos[p], Siz[pos[p]] += Siz[y];
		return;
	}
	int mid = l + r >> 1;
	if (x <= mid)
		ins(ls[p], l, mid, x, y);
	else
		ins(rs[p], mid + 1, r, x, y);
}
vector<int> AA;
inline void Get(int p, int l, int r, int x, int y)
{
	if (!siz[p])
		return;
	if (l == r)
	{
		AA.push_back(pos[p]);
		pos[p] = siz[p] = 0;
		return;
	}
	int mid = l + r >> 1;
	if (x <= mid)
		Get(ls[p], l, mid, x, y);
	if (y > mid)
		Get(rs[p], mid + 1, r, x, y);
	siz[p] = siz[ls[p]] + siz[rs[p]];
}
inline int ask(int p, int l, int r, int x, int y)
{
	if (!siz[p] || (l >= x && r <= y))
		return siz[p];
	int mid = l + r >> 1, s = 0;
	if (x <= mid)
		s += ask(ls[p], l, mid, x, y);
	if (y > mid)
		s += ask(rs[p], mid + 1, r, x, y);
	return s;
}
inline void chk(int p, int l, int r, int x)
{
	if (l == r)
	{
		siz[p] = Siz[pos[p]];
		return;
	}
	int mid = l + r >> 1;
	if (x <= mid)
		chk(ls[p], l, mid, x);
	else
		chk(rs[p], mid + 1, r, x);
	siz[p] = siz[ls[p]] + siz[rs[p]];
}
inline void dfs(int x, int y)
{
	for (int i = 0; i < G[x].size(); ++i)
		if (G[x][i] ^ y)
		{
			if (C[W[x][i]] == m + 1)
				bl[G[x][i]] = bl[x];
			else
				bl[G[x][i]] = G[x][i];
			dfs(G[x][i], x);
		}
}
inline int Ask(int x) { return col[root(x)]; }
int main()
{
	n = read(), m = read();
	for (int i = 1; i <= n; ++i)
		c[i] = read();
	for (int i = 1; i < n; ++i)
		A[i] = read(), B[i] = read(), C[i] = m + 1, add(A[i], B[i], i), add(B[i], A[i], i);
	for (int i = 1; i <= m; ++i)
	{
		O[i] = read();
		if (O[i] == 1)
			X[i] = read(), C[X[i]] = i;
		else if (O[i] == 2)
			V[i] = read(), X[i] = read(), Y[i] = read(), Z[i] = read();
		else
			V[i] = read(), X[i] = read(), Y[i] = read();
	}
	bl[1] = 1, dfs(1, 1), num = n;
	for (int i = 1; i <= n; ++i)
		S[bl[i]].push_back(i), Siz[i] = 1, fa[i] = i, ins(rt[1], 1, n, c[i], i);
	for (int i = m; i; --i)
		if (O[i] == 1)
		{
			int x = bl[A[X[i]]], y = bl[B[X[i]]];
			if (S[x].size() > S[y].size())
				swap(x, y);
			N[i] = x;
			for (auto z : S[x])
				S[y].push_back(z), Q[i].push_back(z), bl[z] = y;
		}
	for (int i = 1; i <= n; ++i)
		bl[i] = 1;
	for (int i = 1; i <= m; ++i)
	{
		if (O[i] == 1)
		{
			for (auto z : Q[i])
			{
				--Siz[root(z)], chk(rt[bl[z]], 1, n, c[z] = Ask(z));
				bl[z] = N[i], ins(rt[bl[z]], 1, n, c[z], z);
			}
		}
		else if (O[i] == 2)
		{
			AA.clear(), Get(rt[bl[V[i]]], 1, n, X[i], Y[i]);
			++num, fa[num] = num;
			for (auto z : AA)
				fa[z] = num, Siz[num] += Siz[z];
			ins(rt[bl[V[i]]], 1, n, Z[i], num);
		}
		else
			printf("%d\n", ask(rt[bl[V[i]]], 1, n, X[i], Y[i]));
	}
}
