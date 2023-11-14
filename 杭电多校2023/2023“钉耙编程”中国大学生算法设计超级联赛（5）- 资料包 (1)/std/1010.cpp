#include<bits/stdc++.h>
using namespace std;

inline bool getmax(int &a, int b) {return a < b ? a = b, true : false;}

#define N 200007
#define M 6000007

struct trie {
    int tot, son[M][2];
    inline void reset() {
        for (int i = 0; i <= tot; ++i) 
            son[i][0] = son[i][1] = 0;
        tot = 0;
    }
    inline void insert(int x) {
        int nw = 0;
        for (int i = 29; ~i; --i) {
            int cur = ((x >> i) & 1);
            if (!son[nw][cur]) son[nw][cur] = ++tot;
            nw = son[nw][cur];
        }
    }
    inline int xormax(int x) {
        int nw = 0;
        int ans = 0;
        for (int i = 29; ~i; --i) {
            int cur = ((x >> i) & 1);
            if (son[nw][cur ^ 1]) {
                ans |= (1 << i);
                nw = son[nw][cur ^ 1];
            } else nw = son[nw][cur];
        }
        return ans;
    }
} t;

int fa[N], l[N];

int a[N], ans[N], res[N], sz[N], mxs[N], tmpans;

vector<int> e[N];

void dfs(int u, int fat) {
    fa[u] = fat;
    sz[u] = 1;
    for (auto v : e[u])
        if (v != fat) {
            dfs(v, u); sz[u] += sz[v];
            if (sz[v] > sz[mxs[u]]) mxs[u] = v;
        }
}

void add(int u) {
    t.insert(a[u]); 
    getmax(tmpans, t.xormax(a[u]));
    for (auto v : e[u]) if (v != fa[u]) add(v);
}


void dsu(int u) {
    for (auto v : e[u])
        if (v != fa[u] && v != mxs[u]) {
            dsu(v); tmpans = 0; t.reset();
        }
    if (mxs[u]) dsu(mxs[u]);
    for (auto v : e[u])
        if (v != fa[u] && v != mxs[u]) add(v);
    t.insert(a[u]); 
    getmax(tmpans, t.xormax(a[u]));
    res[u] = tmpans;
}

void solve(int nw) {
    tmpans = 0; t.reset(); 
    for(l[0] = 0; nw; nw = fa[nw]) l[++l[0]] = nw;
    for (int i = l[0]; i; --i) {
        ans[l[i]] = tmpans;
        if (i == 1) continue;
        t.insert(a[l[i]]); 
		getmax(tmpans, t.xormax(a[l[i]]));
        for (auto v : e[l[i]]) 
            if (v != l[i - 1] && v != fa[l[i]]) add(v);
    }
}

void work() {
    int n; 
    cin >> n;
    int x = 0, y = 0, tar = 0;
    for (int i = 1; i <= n; ++i) {
        cin >> a[i]; 
        t.insert(a[i]);
        if (getmax(tar, t.xormax(a[i]))) x = i;
    }
    for (int i = 1; i <= n; ++i) 
        if ((a[i] ^ a[x]) == tar) y = i;
    for (int i = 1; i < n; ++i) {
        int u, v; 
        cin >> u >> v;
        e[u].push_back(v); 
        e[v].push_back(u);
    }

    dfs(1, 0); 
    t.reset(); dsu(1);
    memset(ans, -1, sizeof(ans)); 
    solve(x); solve(y);

    int sco = INT_MAX;
    for (int i = 2; i <= n; ++i) {
        if (ans[i] == -1) ans[i] = tar;
        sco = min(sco, abs(ans[i] - res[i]));
    }
    printf("%d\n", sco);
    tmpans = 0; t.reset();
    memset(fa, 0, sizeof(fa));
    memset(sz, 0, sizeof(sz));
    memset(mxs, 0, sizeof(mxs));
    memset(res, 0, sizeof(res));
    for (int i = 1; i <= n; ++i) e[i].clear();
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    int t; cin >> t;
    while (t--) work();
    return 0;
}