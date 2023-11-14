#include<bits/stdc++.h>
using namespace std;
#define MAXN 500010
#define MAX 1000010
typedef long long ll;
 
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
 
int T,n,m,q;
 
struct Tree{
    
    vector<int>G[MAX];
    int fa[MAX],dep[MAX];
    int siz[MAX],son[MAX];
    int top[MAX],siz0[MAX];
    int dfn[MAX],low[MAX];
    int tot,dfs_clock;
    ll dp[MAX];
    
    void adde(int u,int v){
        G[u].push_back(v);
        return;
    }
    
    void dfs0(int v,int f){
        siz0[v] = (v<=n?1:0);
        dp[v] = 0;
        for(int u:G[v])
            if(u!=f){
                dfs0(u, v);
                siz0[v] += siz0[u];
                dp[v] += (ll)siz0[u]*(siz0[u]-1)>>1;
            }
        dp[v] += (ll)(n-siz0[v])*(n-siz0[v]-1)>>1;
        return;
    }
    
    void dfs1(int v,int f){
        fa[v] = f;siz[v] = 1;
        dep[v] = dep[f]+1;
        son[v] = 0;int Max=0;
        for(int u:G[v])
            if(u!=f){
                dfs1(u, v);
                siz[v] += siz[u];
                if(siz[u]>Max){
                    Max = siz[u];
                    son[v] = u;
                }
            }
        return;
    }
    
    void dfs2(int v,int tp){
        top[v] = tp;
        dfn[v] = ++dfs_clock;
        if(son[v])dfs2(son[v], tp);
        for(int u:G[v])
            if(u!=fa[v]&&u!=son[v])
                dfs2(u, u);
        low[v] = dfs_clock;
        return;
    }
    
    int lca(int u,int v){
        if(dep[top[u]]>dep[top[v]])swap(u, v);
        return top[u]==top[v]?(dep[u]>dep[v]?v:u):lca(u, fa[top[v]]);
    }
    
    int dis(int u,int v){
        int w=lca(u, v);
        return dep[u]+dep[v]-(dep[w]<<1);
    }
    
    int find(int v,int tp){
        for(;top[v]!=top[tp];v=fa[top[v]])
            if(fa[top[v]]==tp)
                return top[v];
        return son[tp];
    }
    
    ll query(void){
        int k;vector<int>vec;
        k = read();
        for(int i=1;i<=k;i++)
            vec.push_back(read());
        if(k==1)
            return ((ll)n*(n-1)>>1)-dp[vec.front()]+1;
        sort(vec.begin(), vec.end(), [=](const int &u,const int &v){
            return this->dfn[u]<this->dfn[v];
        });
        int u=vec.back();vec.pop_back();
        int v=vec.back();vec.pop_back();
        if(dep[u]>dep[v])swap(u, v);
        int dis0=dis(u, v);
        for(int it:vec)
            if(dis(u, it)+dis(it, v)!=dis0){
                if(dfn[v]<=dfn[it]&&low[it]<=low[v]){
                    v = it;dis0 = dis(u, v);
                }
                else if(dfn[u]<=dfn[it]&&low[it]<=low[u]){
                    u = it;dis0 = dis(u, v);
                }
                else if(dfn[u]<=dfn[v]&&low[v]<=low[u]){
                    u = it;dis0 = dis(u, v);
                }
                else return 0;
            }
        if(dfn[u]>dfn[v])swap(u, v);
        if(dfn[u]<=dfn[v]&&low[v]<=low[u]){
            int w=find(v, u);
            return (ll)(n-siz0[w])*siz0[v];
        }
        else return (ll)siz0[u]*siz0[v];
    }
    
    void init(int V){
        for(int i=1;i<=tot;i++)
            G[i].clear();
        tot = V;dfs_clock = 0;
        return;
    }
    
}tree;

int dfn[MAXN],low[MAXN];
int sta[MAXN],top;
int dfs_clock;
vector<int>G[MAXN];
 
void adde(int u,int v){
    G[u].push_back(v);
    return;
}
 
void tarjan(int v,int f){
    dfn[v] = low[v] = ++dfs_clock;
    sta[++top] = v;
    for(int u:G[v]){
        if(u==f)continue;
        if(!dfn[u]){
            tarjan(u, v);
            if(low[u]>=dfn[v]){
                int k=0;++tree.tot;
                tree.adde(v, tree.tot);
                tree.adde(tree.tot, v);
                do{
                    k = sta[top--];
                    tree.adde(tree.tot, k);
                    tree.adde(k, tree.tot);
                }while(k!=u);
            }
            else low[v] = min(low[v], low[u]);
        }
        else low[v] = min(low[v], dfn[u]);
    }
    return;
}
 
void init(int V){
    for(int i=1;i<=V;i++)
        G[i].clear();
    fill(dfn+1, dfn+V+1, 0);
    dfs_clock= 0;
    return;
}
 
void solve(void){
    n = read();m = read();q = read();
    tree.init(n);init(n);
    for(int i=1,u,v;i<=m;i++){
        u = read();v = read();
        adde(u, v);adde(v, u);
    }
    for(int i=1;i<=n;i++)
        if(!dfn[i])tarjan(i, 0);
    tree.dfs0(1, 0);
    tree.dfs1(1, 0);
    tree.dfs2(1, 1);
    for(int i=1;i<=q;i++)
        printf("%lld\n",tree.query());
    return;
}
 
int main(){
	int size(512<<20); // 512M
	__asm__ ( "movq %0, %%rsp\n"::"r"((char*)malloc(size)+size));     
    T = read();
    while(T--)solve();
    exit(0);
}
