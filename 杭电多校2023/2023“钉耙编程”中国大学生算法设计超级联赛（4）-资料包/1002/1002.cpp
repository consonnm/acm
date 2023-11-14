#include<bits/stdc++.h>
using namespace std;
#define MAXN 1000010
#define LOGV 20
#define mid ((l+r)>>1)

inline char gc() {
    static char s[1<<20|1], *p1=s, *p2=s;
    return (p1==p2)&&(p2=(p1=s)+fread(s, 1, 1<<20, stdin), p1==p2)?EOF:*(p1++);
}
inline int read() {
    int x=0;
    char c=0;
    bool neg=0;
    while(!isdigit(c)) neg^=(c=='-'), c=gc();
    while(isdigit(c)) x=x*10+c-'0', c=gc();
    return neg?-x:x;
}

int T,n,V,a[MAXN],ans[MAXN];
int lim[MAXN],cnt[MAXN];
struct edge{int to,idx;};
vector<edge>G[MAXN];
vector<int>vec;

void adde(int u,int v,int idx){
    G[u].push_back((edge){v,idx});
    return;
}

const int MAX=MAXN*(LOGV+2);
int Max[MAX],lc[MAX],rc[MAX];
int root[MAXN],tot,root0;

int pushup(int x,int y){
    return Max[x]>Max[y]?Max[x]:Max[y];
}

void insert(int p,int &x,int node,bool rev,int l=1,int r=V){
    if(!x)x = ++tot;
    if(l==r){Max[x] = (rev?cnt[l]-1:1);return;}
    if(p<=mid)insert(p, lc[x], lc[node], rev, l, mid);
    else insert(p, rc[x], rc[node], rev, mid+1, r);
    Max[x] = pushup(lc[x]?lc[x]:lc[node], rc[x]?rc[x]:rc[node]);
    return;
}

void merge(int &x,int y,int node,bool rev,int l=1,int r=V){
    if(!x||!y){x |= y;return;}
    if(l==r){Max[x] += Max[y]-(rev?cnt[l]:0);return;}
    merge(lc[x], lc[y], lc[node], rev, l, mid);
    merge(rc[x], rc[y], rc[node], rev, mid+1, r);
    Max[x] = pushup(lc[x]?lc[x]:lc[node], rc[x]?rc[x]:rc[node]);
    return;
}

int query(int val,int x,int node,int l=1,int r=V){
    if(Max[x?x:node]<val)return 0;
    if(l==r)return vec[l-1];
    int res=query(val, rc[x], rc[node], mid+1, r);
    if(!res)res = query(val, lc[x], lc[node], l, mid);
    return res;
}

void build(int &node,int l=1,int r=V){
    if(!node)node = ++tot;
    if(l==r){Max[node] = cnt[l];return;}
    build(lc[node], l, mid);
    build(rc[node], mid+1, r);
    Max[node] = pushup(lc[node], rc[node]);
    return;
}

void clear(void){
    for(int i=1;i<=tot;i++)
        Max[i] = lc[i] = rc[i] = 0;
    fill(root+1, root+n+1, 0);
    root0 = tot = 0;
    return;
}

void dfs(int v,int f,bool rev){
    insert(a[v], root[v], root0, rev);
    for(edge e:G[v])
        if(e.to!=f){
            dfs(e.to, v, rev);
            ans[e.idx] = max(ans[e.idx], query(lim[e.idx], root[e.to], root0));
            merge(root[v], root[e.to], root0, rev);
        }
    return;
}

void solve(void){
    n = read();
    for(int i=1;i<=n;i++)
        vec.push_back(a[i] = read());
    for(int i=1,u,v;i<n;i++){
        u = read();v = read();lim[i] = read();
        adde(u, v, i);adde(v, u, i);
    }
    sort(vec.begin(), vec.end());
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
    V = vec.size();
    for(int i=1;i<=n;i++)
        cnt[a[i] = (int)(lower_bound(vec.begin(), vec.end(), a[i])-vec.begin()+1)]++;
    dfs(1, 0, false);clear();
    build(root0);
    dfs(1, 0, true);clear();
    for(int i=1;i<n;i++)
        printf("%d\n",ans[i]);
    for(int i=1;i<=n;i++)G[i].clear();
    fill(cnt+1, cnt+vec.size()+1, 0);
    fill(ans+1, ans+n, 0);vec.clear();
    return;
}

int main(){
    T = read();
    while(T--)solve();
    return 0;
}
