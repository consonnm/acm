// xtqqwq
#include <bits/stdc++.h>

#define fi first
#define se second
#define mp make_pair
#define pb push_back

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

template <typename T>
bool chkmax(T &x, T y) { return x < y ? x = y, 1 : 0; }
template <typename T>
bool chkmin(T &x, T y) { return x > y ? x = y, 1 : 0; }

int readint()
{
	int x = 0, f = 1;
	char ch = getchar();
	while (ch < '0' || ch > '9')
	{
		if (ch == '-')
			f = -1;
		ch = getchar();
	}
	while (ch >= '0' && ch <= '9')
	{
		x = x * 10 + ch - '0';
		ch = getchar();
	}
	return x * f;
}

int n, m, q, bcnt;
int occ[800005], rt[2][800005];
ll ans[1000005];
char s[400005];
vector<int> vec[2][800005];

namespace bit
{
	int timer;
	int mark[400005], val[400005];
	void add(int x, int c)
	{
		for (; x <= n; x += (x & (-x)))
		{
			if (mark[x] != timer)
				mark[x] = timer, val[x] = 0;
			val[x] += c;
		}
	}
	int ask(int x)
	{
		int ret = 0;
		for (; x; x -= (x & (-x)))
			ret += mark[x] == timer ? val[x] : 0;
		return ret;
	}
}

struct sam
{
	int lst, cnt, ncnt;
	int ch[800005][26], fa[800005], len[800005], siz[800005], tmp[800005], A[800005], rg[800005], tag[800005], ed[800005], sz[800005], son[800005], dfn[800005], rnk[800005], top[800005];
	vector<int> adj[800005];
	void clear()
	{
		for (int i = 1; i <= cnt; i++)
		{
			memset(ch[i], 0, sizeof(ch[i]));
			fa[i] = len[i] = siz[i] = tmp[i] = A[i] = rg[i] = tag[i] = sz[i] = son[i] = 0;
			adj[i].clear();
		}
		cnt = lst = ncnt = 0;
	}
	void ins(int c)
	{
		int p = lst, np = ++cnt;
		lst = np, len[np] = len[p] + 1, siz[np] = 1;
		for (; p && !ch[p][c]; p = fa[p])
			ch[p][c] = np;
		if (!p)
			return (void)(fa[np] = 1);
		int q = ch[p][c];
		if (len[q] == len[p] + 1)
			return (void)(fa[np] = q);
		int nq = ++cnt;
		len[nq] = len[p] + 1, rg[nq] = rg[q];
		memcpy(ch[nq], ch[q], sizeof(ch[nq]));
		fa[nq] = fa[q], fa[q] = fa[np] = nq;
		for (; ch[p][c] == q; p = fa[p])
			ch[p][c] = nq;
	}
	void dfs1(int u)
	{
		sz[u] = 1;
		for (auto v : adj[u])
		{
			dfs1(v);
			sz[u] += sz[v];
			if (sz[v] > sz[son[u]])
				son[u] = v;
		}
	}
	void dfs2(int u)
	{
		dfn[u] = ++ncnt, rnk[ncnt] = u;
		if (son[u])
			top[son[u]] = top[u], dfs2(son[u]);
		for (auto v : adj[u])
		{
			if (v == son[u])
				continue;
			top[v] = v;
			dfs2(v);
		}
	}
	void init()
	{
		for (int i = 1; i <= cnt; i++)
			tmp[len[i]]++;
		for (int i = 1; i <= n; i++)
			tmp[i] += tmp[i - 1];
		for (int i = cnt; i >= 1; i--)
			A[tmp[len[i]]--] = i;
		for (int i = cnt; i >= 2; i--)
		{
			int u = A[i];
			adj[fa[u]].pb(u);
			siz[fa[u]] += siz[u];
		}
		dfs1(1);
		top[1] = 1;
		dfs2(1);
	}
	int find(int x, int d)
	{
		while (len[fa[top[x]]] >= d)
			x = fa[top[x]];
		int L = dfn[top[x]], R = dfn[x], res = 0;
		while (L <= R)
		{
			int mid = (L + R) / 2;
			if (len[rnk[mid]] >= d)
				res = mid, R = mid - 1;
			else
				L = mid + 1;
		}
		return rnk[res];
	}
} T[2];

void init()
{
	T[0].cnt = T[0].lst = T[1].cnt = T[1].lst = 1;
	for (int i = 1; i <= n; i++)
		T[0].ins(s[i] - 'a'), T[0].rg[T[0].lst] = i, T[0].ed[i] = T[0].lst;
	for (int i = n; i >= 1; i--)
		T[1].ins(s[i] - 'a'), T[1].rg[T[1].lst] = i, T[1].ed[i] = T[1].lst;
	T[0].init();
	T[1].init();
	for (int i = 2; i <= T[0].cnt; i++)
	{
		int u = T[0].A[i];
		int L = T[0].rg[u] - T[0].len[u] + 1, R = T[0].rg[u];
		int id = T[1].find(T[1].ed[L], R - L + 1);
		if (T[1].len[id] == R - L + 1)
		{
			T[0].tag[u] = T[1].tag[id] = ++bcnt;
			rt[0][bcnt] = u, rt[1][bcnt] = id;
			occ[bcnt] = T[0].siz[u];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = T[i].cnt; j >= 2; j--)
		{
			int u = T[i].A[j];
			if (T[i].tag[u])
				continue;
			for (int k = 0; k < 26; k++)
				if (T[i].ch[u][k])
					T[i].tag[u] = T[i].tag[T[i].ch[u][k]];
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 2; j <= T[i].cnt; j++)
		{
			int u = T[i].A[j];
			vec[i][T[i].tag[u]].pb(u);
		}
	}
}

namespace part2
{
	ll val[800005], res[800005];
	vector<ll> pre[800005];
	vector<pii> chn[800005];
	vector<pair<pii, int>> qry[800005];
	void clear()
	{
		for (int i = 1; i <= T[0].cnt; i++)
			val[i] = res[i] = 0;
		for (int i = 1; i <= bcnt; i++)
			pre[i].clear(), chn[i].clear(), qry[i].clear();
	}
	void solve1()
	{
		for (int i = 1; i <= bcnt; i++)
			for (auto r : vec[0][i])
				val[r] += val[T[0].fa[r]];
		for (int i = 1; i <= bcnt; i++)
		{
			int pl = vec[0][i].size();
			ll now = 0;
			sort(chn[i].begin(), chn[i].end());
			for (auto r : chn[i])
				val[vec[0][i][vec[0][i].size() - r.se - 1]]--;
			for (auto r : vec[1][i])
			{
				int d = T[1].len[r] - T[1].len[T[1].fa[r]];
				while (pl > vec[0][i].size() - d)
					now += val[vec[0][i][--pl]];
				res[r] = now + res[T[1].fa[r]];
			}
			pl = chn[i].size();
			for (int j = 0; j < vec[1][i].size(); j++)
			{
				int r = vec[1][i][j];
				while (pl && chn[i][pl - 1].fi == vec[1][i].size() - j - 1)
					pl--;
				res[r] += chn[i].size() - pl;
			}
			pre[i].resize(vec[0][i].size());
			pre[i][0] = val[vec[0][i][0]];
			for (int j = 1; j < vec[0][i].size(); j++)
				pre[i][j] = pre[i][j - 1] + val[vec[0][i][j]];
		}
	}
	void addq(int x, int y, int i, int c)
	{
		int u = T[1].find(T[1].ed[x], y - x + 1);
		int bel = T[1].tag[u];
		ans[i] += c * res[T[1].fa[u]];
		int px = T[1].rg[u] - T[1].rg[rt[1][bel]];
		int py = T[1].rg[rt[1][bel]] + T[1].len[rt[1][bel]] - 1 - (T[1].rg[u] + y - x);
		int L = vec[0][bel].size() - T[1].len[u] + T[1].len[T[1].fa[u]], R = vec[0][bel].size() - py - 1;
		ans[i] += c * (pre[bel][R] - (L ? pre[bel][L - 1] : 0));
		qry[bel].pb(mp(mp(px, py), c > 0 ? i : -i));
	}
	void solve2()
	{
		for (int i = 1; i <= bcnt; i++)
		{
			sort(qry[i].begin(), qry[i].end());
			reverse(qry[i].begin(), qry[i].end());
			bit::timer++;
			int pl = chn[i].size();
			for (auto r : qry[i])
			{
				while (pl && chn[i][pl - 1].fi >= r.fi.fi)
					bit::add(chn[i][--pl].se + 1, 1);
				ll tmp = chn[i].size() - pl - bit::ask(r.fi.se);
				if (r.se > 0)
					ans[r.se] += tmp;
				else
					ans[-r.se] -= tmp;
			}
		}
	}
}

namespace part1
{
	ll val[800005];
	vector<pii> qry1[400005], qry2[400005];
	namespace bit
	{
		ll val[400005];
		void add(int x, ll c)
		{
			for (; x <= n; x += (x & (-x)))
				val[x] += c;
		}
		ll query(int x)
		{
			ll ret = 0;
			for (; x; x -= (x & (-x)))
				ret += val[x];
			return ret;
		}
		void change(int l, int r, ll x)
		{
			add(l, x);
			add(r + 1, -x);
		}
	}
	namespace lct
	{
		int ch[800005][2], fa[800005], tag[800005];
		bool rev[800005];
		bool nroot(int x) { return ch[fa[x]][0] == x || ch[fa[x]][1] == x; }
		bool son(int x) { return ch[fa[x]][1] == x; }
		void reverse(int x) { swap(ch[x][0], ch[x][1]), rev[x] ^= 1; }
		void pushdown(int x)
		{
			if (rev[x])
				reverse(ch[x][0]), reverse(ch[x][1]), rev[x] = 0;
			if (tag[x])
			{
				if (ch[x][0])
					tag[ch[x][0]] = tag[x];
				if (ch[x][1])
					tag[ch[x][1]] = tag[x];
			}
		}
		void pushall(int x)
		{
			if (nroot(x))
				pushall(fa[x]);
			pushdown(x);
		}
		void rotate(int x)
		{
			int y = fa[x], z = fa[y], k = son(x), w = ch[x][!k];
			if (nroot(y))
				ch[z][son(y)] = x;
			ch[x][!k] = y, ch[y][k] = w;
			if (w)
				fa[w] = y;
			fa[y] = x, fa[x] = z;
		}
		void splay(int x)
		{
			pushall(x);
			while (nroot(x))
			{
				int y = fa[x];
				if (nroot(y))
					rotate(son(x) == son(y) ? y : x);
				rotate(x);
			}
		}
		void access(int x, int r)
		{
			int y = 0, lst = 0;
			for (; x; x = fa[y = x])
			{
				splay(x);
				if (tag[x])
				{
					bit::change(max(1, lst - T[0].len[x] + 1), tag[x] - T[0].len[x], val[x]);
					lst = tag[x];
				}
				ch[x][1] = y;
			}
			tag[y] = r;
		}
		vector<pii> access2(int x, int r)
		{
			int y = 0;
			vector<pii> vec;
			for (; x; x = fa[y = x])
			{
				splay(x);
				if (tag[x])
					vec.pb(mp(tag[x], T[1].len[x]));
				ch[x][1] = y;
			}
			tag[y] = r;
			reverse(vec.begin(), vec.end());
			return vec;
		}
	}
	void clear()
	{
		for (int i = 1; i <= T[0].cnt; i++)
			val[i] = 0;
		for (int i = 1; i <= n; i++)
			bit::val[i] = 0;
		for (int i = 1; i <= n; i++)
			qry1[i].clear(), qry2[i].clear();
	}
	void solve()
	{
		for (int i = 1; i <= bcnt; i++)
			for (auto r : vec[0][i])
				val[r] += val[T[0].fa[r]];
		for (int i = 1; i <= T[1].cnt; i++)
			lct::fa[i] = T[1].fa[i], lct::ch[i][0] = lct::ch[i][1] = 0, lct::tag[i] = lct::rev[i] = 0;
		for (int i = n; i >= 1; i--)
		{
			vector<pii> vec = lct::access2(T[1].ed[i], i);
			sort(qry2[i].begin(), qry2[i].end());
			int pl = 0;
			for (auto r : qry2[i])
			{
				int res = 0;
				while (pl < vec.size() && vec[pl].fi + vec[pl].se - 1 <= r.fi)
					pl++;
				if (pl == vec.size())
					res = i + (vec.size() ? vec.back().se : 0);
				else
				{
					int tmp = r.fi - vec[pl].fi + 1;
					if (pl && tmp <= vec[pl - 1].se)
						res = i + vec[pl - 1].se;
					else
						res = i + tmp;
				}
				part2::addq(i, r.fi, r.se, 1);
				if (res > i)
					part2::addq(i, res - 1, r.se, -1);
			}
		}
		for (int i = 1; i <= T[0].cnt; i++)
			lct::fa[i] = T[0].fa[i], lct::ch[i][0] = lct::ch[i][1] = lct::tag[i] = 0;
		for (int i = 1; i <= n; i++)
		{
			lct::access(T[0].ed[i], i);
			for (auto r : qry1[i])
				ans[r.se] -= bit::query(r.fi);
		}
	}
}

void clear()
{
	part1::clear();
	part2::clear();
	T[0].clear();
	T[1].clear();
	for (int i = 1; i <= bcnt; i++)
	{
		vec[0][i].clear(), vec[1][i].clear();
		occ[i] = rt[0][i] = rt[1][i] = 0;
	}
	bcnt = 0;
	for (int i = 1; i <= q; i++)
		ans[i] = 0;
}

void solve()
{
	clear();
	n = readint();
	m = readint();
	q = readint();
	scanf("%s", s + 1);
	init();
	int x, y;
	for (int i = 1; i <= m; i++)
	{
		x = readint();
		y = readint();
		int u = T[0].find(T[0].ed[y], y - x + 1);
		int bel = T[0].tag[u];
		part2::chn[bel].pb(mp(T[0].rg[u] - y + x - (T[0].rg[rt[0][bel]] - T[0].len[rt[0][bel]] + 1), T[0].rg[rt[0][bel]] - T[0].rg[u]));
		part2::val[u]++;
		part1::val[u]++;
	}
	for (int i = 1; i <= q; i++)
	{
		x = readint();
		y = readint();
		part1::qry1[y].pb(mp(x, i));
		part1::qry2[x].pb(mp(y, i));
	}
	part2::solve1();
	part1::solve();
	part2::solve2();
	for (int i = 1; i <= q; i++)
		printf("%lld\n", ans[i]);
}

int main()
{
	int T = readint();
	while (T--)
		solve();
}