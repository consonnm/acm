#include <bits/stdc++.h>
#define rep(i,h,t) for (int i=h;i<=t;i++)
#define dep(i,t,h) for (int i=t;i>=h;i--)
#define ll long long
using namespace std;
string s;
const int N=2e5+10;
const ll m1=1e9+7;
const ll m2=998244353;
const ll m3=-m1-m2;
bool vis[N];
int f[N],sum,son[N],rt;
ll a[N];
vector<int> pq[N];
void gr(int x,int y)
{
  son[x]=1;f[x]=0;
  for (auto v:pq[x])
    if (vis[v]&&v!=y)
    {
      gr(v,x);
      son[x]+=son[v];
      f[x]=max(f[x],son[v]);
    }
  f[x]=max(f[x],sum-son[x]);
  if (f[x]<f[rt]) rt=x;
}
vector<ll> ve;
unordered_map<ll,int> M;
ll ans=0;
void gd(int x,int fa,ll dis,ll o)
{
  dis+=a[x];
  ve.push_back(dis);
  if (M.find(-(dis+o))!=M.end()) ans+=M[-(dis+o)];
  if (dis+o==0) ans++;
  for (auto v:pq[x])
    if (vis[v]&&v!=fa)
    {
      gd(v,x,dis,o);
    }
}
void solve(int x,int y)
{
  vis[x]=0;
  M.clear();
  for (auto v:pq[x])
      if (vis[v]&&v!=x)
      {
          gd(v,x,0,a[x]);
          for (auto v:ve) M[v]++;
          ve.clear();
      }
  for (auto v:pq[x])
    if (vis[v])
    {
      rt=0; sum=son[v];
      gr(v,x);
      solve(rt,y+1);
    }
}
int main()
{
    ios::sync_with_stdio(false);
    int n;
    cin>>n;
    cin>>s;
    f[0]=1e9;
    rep(i,1,n) vis[i]=1;
    rep(i,1,n-1)
    {
        int x,y;
        cin>>x>>y;
        pq[x].push_back(y);
        pq[y].push_back(x);
    }
    for (int i=1;i<=n;i++)
    {
      if (s[i-1]=='a') a[i]=m1;
      if (s[i-1]=='b') a[i]=m2;
      if (s[i-1]=='c') a[i]=m3;
    }
    sum=n;
    gr(1,0);
    solve(rt,0);
    cout<<ans<<endl;
    return 0;
}
/*
6
abcabc
1 2
2 4
1 3
3 5
5 6
*/ 


