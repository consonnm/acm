#include <bits/stdc++.h>

#define pb push_back
#define fi first
#define se second
#define mp make_pair

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef long double ld;

template <typename T>
bool chkmin(T &x, T y) { return x > y ? x = y, 1 : 0; }
template <typename T>
bool chkmax(T &x, T y) { return x < y ? x = y, 1 : 0; }

ll readint()
{
	ll x = 0, f = 1;
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

int n, tot, fs, mn, rt, cur;
int v[200005], nxt[200005], h[100005], siz[100005], dep[100005], pl[100005], f[100005], dis[100005];
ll c[200005], val[20][100005], r[100005];
bool hv[100005], vis[100005];
vector<pll> vec[100005];
queue<int> q;

void addedge(int x, int y, ll z)
{
	v[++tot] = y;
	c[tot] = z;
	nxt[tot] = h[x];
	h[x] = tot;
	v[++tot] = x;
	c[tot] = z;
	nxt[tot] = h[y];
	h[y] = tot;
}

void getsize(int u, int fa)
{
	siz[u] = 1;
	for (int p = h[u]; p; p = nxt[p])
	{
		if (v[p] == fa || vis[v[p]])
			continue;
		getsize(v[p], u);
		siz[u] += siz[v[p]];
	}
}

void getroot(int u, int fa)
{
	int mx = 0;
	for (int p = h[u]; p; p = nxt[p])
	{
		if (v[p] == fa || vis[v[p]])
			continue;
		getroot(v[p], u);
		chkmax(mx, siz[v[p]]);
	}
	if (chkmin(mn, max(mx, fs - siz[u])))
		rt = u;
}

void dfs(int u, int fa)
{
	vec[rt].pb(mp(r[u] - val[cur][u], u));
	for (int p = h[u]; p; p = nxt[p])
	{
		if (v[p] == fa || vis[v[p]])
			continue;
		val[cur][v[p]] = val[cur][u] + c[p];
		dfs(v[p], u);
	}
}

void solve(int u)
{
	vis[u] = 1;
	getsize(u, -1);
	cur = dep[u];
	val[cur][u] = 0;
	vec[u].clear();
	dfs(u, -1);
	sort(vec[u].begin(), vec[u].end());
	for (int p = h[u]; p; p = nxt[p])
	{
		if (vis[v[p]])
			continue;
		fs = siz[v[p]], mn = n, rt = 0;
		getroot(v[p], -1);
		dep[rt] = dep[u] + 1;
		// adj[u].pb(rt);
		f[rt] = u;
		solve(rt);
	}
}

void push(int x)
{
	int bef = x;
	// cout<<"push "<<x<<' '<<dis[x]<<endl;
	while (x)
	{
		while (pl[x] >= 0 && vec[x][pl[x]].fi >= val[dep[x]][bef])
		{
			if (!hv[vec[x][pl[x]].se])
			{
				hv[vec[x][pl[x]].se] = 1;
				dis[vec[x][pl[x]].se] = dis[bef] + 1;
				q.push(vec[x][pl[x]].se);
			}
			pl[x]--;
		}
		x = f[x];
	}
}

int sum = 0;
void solve()
{
	n = readint();
	sum += n;
	// cerr<<"test "<<n<<' '<<sum<<endl;
	for (int i = 2; i <= n; i++)
		r[i] = readint();
	for (int i = 1; i <= n; i++)
		h[i] = f[i] = vis[i] = hv[i] = 0;
	tot = 0;
	int x, y;
	ll z;
	for (int i = 1; i < n; i++)
	{
		x = readint();
		y = readint();
		z = readint();
		addedge(x, y, z);
	}
	getsize(1, -1);
	fs = n, mn = n, rt = 0;
	getroot(1, -1);
	dep[rt] = 0;
	solve(rt);
	for (int i = 1; i <= n; i++)
		dis[i] = 1 << 30;
	for (int i = 1; i <= n; i++)
		hv[i] = 0;
	dis[1] = 0, hv[1] = 1;
	q.push(1);
	for (int i = 1; i <= n; i++)
		pl[i] = vec[i].size() - 1;
	while (!q.empty())
	{
		int t = q.front();
		q.pop();
		push(t);
	}
	for (int i = 2; i <= n; i++)
		printf("%d ", dis[i] == (1 << 30) ? -1 : dis[i]);
	printf("\n");
}

int main()
{
	int T = readint();
	while (T--)
		solve();
}