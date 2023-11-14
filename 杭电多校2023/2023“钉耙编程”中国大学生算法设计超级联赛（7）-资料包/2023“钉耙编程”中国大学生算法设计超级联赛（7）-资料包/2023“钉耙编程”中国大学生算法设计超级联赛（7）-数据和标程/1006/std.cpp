#include <bits/stdc++.h>
using namespace std;

int read();
typedef long long ll;
#define fr(i, l, r) for (int i = (l); i <= (r); ++i)
#define rf(i, l, r) for (int i = (l); i >= (r); --i)
#define fo(i, l, r) for (int i = (l); i < (r); ++i)
#define foredge(i, u, v) for (int i = fir[u], v; v = to[i], i; i = nxt[i])
#define filein(file) freopen(file ".in", "r", stdin)
#define fileout(file) freopen(file ".out", "w", stdout)

const int N = 1e6 + 5, M = 5005, MOD = 1e9 + 7, PHI = MOD - 1, B = 31624, dn1 = MOD >> 1, g = 5;
int n, m, a[N], b[N];

ll qpow(ll a, int x)
{
	ll z = 1;
	for (; x; x >>= 1, a = a * a % MOD)
		if (x & 1)
			z = z * a % MOD;
	return z;
}
int f[B + 5], p[4000];
bool np[B + 5];
namespace BSGS
{
	const int P = 2e6 + 3;
	int hd[P], nxt[N], ky[N], val[N], nds;
	void insert(int x, int v)
	{
		int hsh = x % P;
		ky[++nds] = x;
		val[nds] = v;
		nxt[nds] = hd[hsh];
		hd[hsh] = nds;
	}
	int query(int x)
	{
		for (int u = hd[x % P]; u; u = nxt[u])
		{
			if (ky[u] == x)
				return val[u];
		}
		return -1;
	}
	void init()
	{
		ll pw = 1, tmp = 1;
		int lim = int(1e3);
		fr(i, 1, lim) pw = pw * 5 % MOD;
		fr(i, 1, int(1e6) + 1)
			tmp = tmp * pw % MOD,
			insert(tmp, i * lim);
	}
	int dlog(ll x)
	{
		for (int i = 0;; x = x * g % MOD, ++i)
		{
			int res = query(x);
			if (~res)
				return res - i;
		}
	}
}
void prepare()
{
	BSGS::init();
	f[1] = 0;
	fr(i, 2, B)
	{
		if (!np[i])
			p[++p[0]] = i, f[i] = BSGS::dlog(i);
		for (int j = 1, ni; j <= p[0] && (ni = i * p[j]) <= B; ++j)
		{
			np[ni] = 1;
			f[ni] = (f[i] + f[p[j]]) % PHI;
			if (i % p[j] == 0)
				break;
		}
	}
}
ll dlog(int x)
{
	if (x <= B)
		return f[x];
	int q = MOD / x, r = MOD % x;
	return r * 2 <= x ? dlog(r) + dn1 - f[q] : dlog(x - r) - f[q + 1];
	// p=qx+r,x=(-r)/q
	// p=(q+1)x+(r-x),x=(x-r)/(q+1)
}

struct FPOW
{
	int bs[1 << 15 | 2], gs[1 << 15 | 2];
	void init()
	{
		*bs = 1;
		fr(i, 1, 1 << 15) bs[i] = 1ll * bs[i - 1] * g % MOD;
		*gs = 1;
		fr(i, 1, 1 << 15) gs[i] = 1ll * gs[i - 1] * bs[1 << 15] % MOD;
	}
	ll operator()(const int &i)
	{
		// assert(i>=0);
		return 1ll * bs[i & (1 << 15) - 1] * gs[i >> 15] % MOD;
	}
} G;
struct Query
{
	int l, r, id;
} q[M];
ll ans[M];

typedef pair<int, int> pii;
#define fi first
#define se second
ll now;
vector<pii> all;
int pre[N], nxt[N];
void rebuild(int L, int R)
{
	fr(i, L, R) all.push_back(pii(a[i], i));
	auto ptr = all.end() - (R - L + 1);
	sort(ptr, all.end());
	inplace_merge(all.begin(), ptr, all.end());
	int lst = 0;
	fo(i, 0, all.size())
	{
		pre[all[i].se] = lst;
		nxt[lst] = all[i].se;
		lst = all[i].se;
	}
	pre[0] = all.back().se;
	nxt[all.back().se] = 0;

	now = 1;
	fo(i, 1, all.size())
	{
		now = now * G(1ll * a[all[i].se] * b[all[i - 1].se] % PHI) % MOD;
		// assert(G(1ll*a[all[i].se]*b[all[i-1].se]%PHI) == qpow(a[all[i-1].se],a[all[i].se]));
	}
}
void del(int i)
{
	if (pre[i])
	{
		if (nxt[i])
			now = now * G(1ll * b[pre[i]] * (a[nxt[i]] - a[i]) % PHI) % MOD;
		else
			now = now * G(PHI - 1ll * b[pre[i]] * a[i] % PHI) % MOD;
	}
	if (nxt[i])
		now = now * G(PHI - 1ll * b[i] * a[nxt[i]] % PHI) % MOD;
	pre[nxt[i]] = pre[i];
	nxt[pre[i]] = nxt[i];
}
void add(int i) { pre[nxt[i]] = nxt[pre[i]] = i; }
int main()
{
	prepare();
	scanf("%d%d", &n, &m);
	fr(i, 1, n)
	{
		scanf("%d", a + i);
		b[i] = (dlog(a[i]) % PHI + PHI) % PHI;
	}
	fr(i, 1, m) scanf("%d%d", &q[i].l, &q[i].r), q[i].id = i;
	int bls = max(1, int(n / sqrt(m) * 1.2));
	sort(q + 1, q + m + 1, [&](const Query &a, const Query &b)
		 {
		if(a.l/bls!=b.l/bls) return a.l/bls>b.l/bls;
		return a.r>b.r; });
	G.init();
	for (int L = n / bls * bls, R = n, pt = 1; R >= 1; R = L - 1, L -= bls)
	{
		L = max(L, 1);
		rebuild(L, R);
		int l = L, r = n;
		for (; pt <= m && q[pt].l >= L; ++pt)
		{
			while (r > q[pt].r)
				del(r--);
			ll tmp = now;
			while (l < q[pt].l)
				del(l++);
			ans[q[pt].id] = now;
			while (l > L)
				add(--l);
			now = tmp;
		}
	}
	fr(i, 1, m) printf("%lld\n", ans[i]);
	return 0;
}
int read()
{
	static int x, c, f;
	x = 0, f = 1;
	do
		c = getchar(), (c == '-' && (f = -1));
	while (!isdigit(c));
	do
		x = x * 10 + (c & 15), c = getchar();
	while (isdigit(c));
	return x * f;
}
