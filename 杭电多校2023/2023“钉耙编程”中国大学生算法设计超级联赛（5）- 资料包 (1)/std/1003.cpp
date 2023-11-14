#include <bits/stdc++.h>
using namespace std ;
#define int long long
#define rep(i, a, b) for (ll i = (a); i <= (b); ++i)
#define per(i, a, b) for (ll i = (a); i >= (b); --i)
#define loop(it, v) for (auto it = v.begin(); it != v.end(); it++)
#define cont(i, x) for (ll i = head[x]; i; i = edge[i].nex)
#define clr(a) memset(a, 0, sizeof(a))
#define ass(a, cnt) memset(a, cnt, sizeof(a))
#define cop(a, b) memcpy(a, b, sizeof(a))
#define lowbit(x) (x & -x)
#define all(x) x.begin(), x.end()
#define SC(t, x) static_cast <t> (x)
#define ub upper_bound
#define lb lower_bound
#define pqueue priority_queue
#define mp make_pair
#define pb push_back
#define pof pop_front
#define pob pop_back
#define fi first
#define se second
#define y1 y1_
#define Pi acos(-1.0)
#define iv inline void
#define enter putchar('\n')
#define siz(x) ((ll)x.size())
#define file(x) freopen(x".in", "r", stdin),freopen(x".out", "w", stdout)
typedef double db ;
typedef long long ll ;
typedef unsigned long long ull ;
typedef pair <ll, ll> pii ;
typedef vector <ll> vi ;
typedef vector <pii> vii ;
typedef queue <ll> qi ;
typedef queue <pii> qii ;
typedef set <ll> si ;
typedef map <ll, ll> mii ;
typedef map <string, ll> msi ;
const ll maxn = 5e5 + 100 ;
const ll inf = 0x3f3f3f3f ;
const ll iinf = 1 << 30 ;
const ll linf = 2e18 ;
const ll mod = 1e9 + 7 ;
const double eps = 1e-7 ;
template <class T = ll> T chmin(T &a, T b) { return a = min(a, b);}
template <class T = ll> T chmax(T &a, T b) { return a = max(a, b);}
template <class T = ll> iv red(T &x) { x -= mod, x += x >> 31 & mod;}
template <class T = ll> T read()
{
	T f = 1, a = 0;
	char ch = getchar() ;
	while (!isdigit(ch)) { if (ch == '-') f = -1 ; ch = getchar() ; }
	while (isdigit(ch)) { a =  (a << 3) + (a << 1) + ch - '0' ; ch = getchar() ; }
	return a * f ;
}

ll t;

ll n, ans;

char s[maxn];

ll d[2][maxn];

void init()
{
	int l = 0, r = 0;
	rep(i, 1, n)
	{
		int k = 1;
		if(i <= r) k = min(d[1][l + r - i], r - i + 1);
		while(s[i + k] == s[i - k]) ++ k;
		d[1][i] = k;
		if(i + k - 1 > r) r = i + k - 1, l = i - k + 1;
	}
	l = r = 0;
	rep(i, 1, n)
	{
		int k = 0;
		if(i < r) k = min(r - i, d[0][l + r - i - 1]);
		while(s[i + k + 1] == s[i - k]) ++ k;
		d[0][i] = k;
		if(i + k > r) r = i + k, l = i - k + 1;
	}
	// rep(i, 1, n) printf("%lld %lld\n", d[0][i], d[1][i]);
}

ll tot, rt[2][maxn];

struct node
{
	ll ls, rs, sum;
}
tr[20000005];

ll newnode(int now)
{
	++ tot;
	tr[tot] = tr[now];
	return tot;
}

ll modify(int now, int l, int r, int x)
{
	ll ret = newnode(now);
	tr[ret].sum ++;
	if(l == r) return ret;
	int mid = (l + r) >> 1;
	if(x <= mid) tr[ret].ls = modify(tr[ret].ls, l, mid, x);
	else tr[ret].rs = modify(tr[ret].rs, mid + 1, r, x);
	return ret;
}

ll query(ll rt1, ll rt2, ll l, ll r, ll ql)
{
	if(ql <= l) return tr[rt2].sum - tr[rt1].sum;
	ll ret = 0, mid = (l + r) >> 1;
	if(ql <= mid) ret += query(tr[rt1].ls, tr[rt2].ls, l, mid, ql);
	ret += query(tr[rt1].rs, tr[rt2].rs, mid + 1, r, ql);
	return ret;
}

signed main()
{
	t = read();
	while(t --)
	{
		scanf("%s", s + 1);
		n = strlen(s + 1);
		s[0] = '?';
		init();
		ans = tot = 0;
		rep(i, 1, n)
		{
			rt[0][i] = modify(rt[0][i - 1], 1, n, i + d[0][i]);
			rt[1][i] = modify(rt[1][i - 1], 1, n, i + d[1][i] - 1);
			if(d[0][i] > 1) ans += query(rt[0][i - d[0][i] / 2 - 1], rt[0][i - 1], 1, n, i);
			if(d[0][i] > 0) ans += query(rt[1][i - (d[0][i] - 1) / 2 - 1], rt[1][i], 1, n, i);
		}
		printf("%lld\n", ans);
	}
	return 0;
}