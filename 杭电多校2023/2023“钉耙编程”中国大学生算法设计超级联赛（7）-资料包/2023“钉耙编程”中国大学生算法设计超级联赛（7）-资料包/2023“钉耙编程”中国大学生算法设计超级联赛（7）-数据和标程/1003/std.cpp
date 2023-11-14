#include<bits/stdc++.h>

#define pb push_back
#define fi first
#define se second
#define mp make_pair

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<int,int> pii;
typedef pair<ll,ll> pll;
typedef long double ld;

template <typename T> bool chkmin(T &x,T y){return x>y?x=y,1:0;}
template <typename T> bool chkmax(T &x,T y){return x<y?x=y,1:0;}

int readint(){
	int x=0,f=1; char ch=getchar();
	while(ch<'0'||ch>'9'){if(ch=='-')f=-1;ch=getchar();}
	while(ch>='0'&&ch<='9'){x=x*10+ch-'0';ch=getchar();}
	return x*f;
}

const int N=100000,cys=998244353;
int n;
int a[100005],ch[100005][2],f[100005],siz[100005],num[100005],hv[100005],cur[10];
ll fac[100005],inv[100005],ds[100005];
bool can[100005][2];
vector<int> vec[100005];
map<vector<int>,int> res,vis;

ll qpow(ll x,ll p){
	ll ret=1;
	for(;p;p>>=1,x=x*x%cys) if(p&1) ret=ret*x%cys;
	return ret;
}

int getf(int x){return x==f[x]?x:f[x]=getf(f[x]);}

void merge(int x,int y){
	int fx=getf(x),fy=getf(y);
	if(fx==fy) return;
	f[fx]=fy;
	siz[fy]+=siz[fx];
}

void work(int u,int x){
	merge(u,ch[u][0]);
	merge(u,ch[u][1]);
	merge(u,ch[ch[u][x]][x^1]);
}

void solve(){
	n=readint();
	for(int i=1;i<=n;i++) a[i]=readint();
	for(int i=1;i<=n;i++) ch[i][0]=readint(),ch[i][1]=readint();
	for(int i=1;i<=n;i++) f[i]=i,siz[i]=1;
	for(int i=1;i<=n;i++) can[i][0]=can[i][1]=0;
	for(int i=1;i<=n;i++) if(ch[i][0]&&ch[i][1]&&ch[ch[i][0]][1]) can[i][0]=1,work(i,0);
	for(int i=1;i<=n;i++) if(ch[i][0]&&ch[i][1]&&ch[ch[i][1]][0]) can[i][1]=1,work(i,1);
	for(int i=1;i<=n;i++) vec[i].clear();
	for(int i=1;i<=n;i++) vec[getf(i)].pb(i);
	ll ans=1;
	for(int i=1;i<=n;i++){
		if(getf(i)==i){
			if(siz[i]==4){
				for(auto r:vec[i]){
					if(can[r][0]){
						if(a[r]==a[ch[ch[r][0]][1]]&&a[ch[r][0]]==a[ch[r][1]]){
							if(a[r]!=a[ch[r][0]]) ans=ans*2%cys;
						}
						else ans=ans*4%cys;
					}
					else if(can[r][1]){
						if(a[r]==a[ch[ch[r][1]][0]]&&a[ch[r][0]]==a[ch[r][1]]){
							if(a[r]!=a[ch[r][0]]) ans=ans*2%cys;
						}
						else ans=ans*4%cys;
					}
				}
			}
			else if(siz[i]==6){
				vector<int> tmp;
				int tcnt=0;
				for(auto r:vec[i]){
					if(can[r][0]&&(can[ch[r][0]][0]||can[ch[r][0]][1])){
						int u0=ch[r][0],u1=ch[r][1],v1=ch[u0][1],v0=ch[u0][0],w0=can[u0][0]?ch[v0][1]:ch[v1][0];
						if(!hv[a[r]]) hv[a[r]]=++tcnt;
						if(!hv[a[u1]]) hv[a[u1]]=++tcnt;
						if(!hv[a[v1]]) hv[a[v1]]=++tcnt;
						if(!hv[a[u0]]) hv[a[u0]]=++tcnt;
						if(!hv[a[v0]]) hv[a[v0]]=++tcnt;
						if(!hv[a[w0]]) hv[a[w0]]=++tcnt;
						tmp.pb(hv[a[r]]);
						tmp.pb(hv[a[u1]]);
						tmp.pb(hv[a[v1]]);
						tmp.pb(hv[a[u0]]);
						tmp.pb(hv[a[v0]]);
						tmp.pb(hv[a[w0]]);
						ans=ans*res[tmp]%cys;
						// cerr<<"test ";
						// for(auto r:tmp) cerr<<r<<' ';
						// cerr<<res[tmp]<<endl;
					}
					else if(can[r][1]&&(can[ch[r][1]][0]||can[ch[r][1]][1])){
						int u0=ch[r][0],u1=ch[r][1],v0=ch[u1][0],v1=ch[u1][1],w0=can[u1][0]?ch[v0][1]:ch[v1][0];
						if(!hv[a[r]]) hv[a[r]]=++tcnt;
						if(!hv[a[u0]]) hv[a[u0]]=++tcnt;
						if(!hv[a[v0]]) hv[a[v0]]=++tcnt;
						if(!hv[a[u1]]) hv[a[u1]]=++tcnt;
						if(!hv[a[v1]]) hv[a[v1]]=++tcnt;
						if(!hv[a[w0]]) hv[a[w0]]=++tcnt;
						tmp.pb(hv[a[r]]);
						tmp.pb(hv[a[u0]]);
						tmp.pb(hv[a[v0]]);
						tmp.pb(hv[a[u1]]);
						tmp.pb(hv[a[v1]]);
						tmp.pb(hv[a[w0]]);
						ans=ans*res[tmp]%cys;
					}
				}
				for(auto r:vec[i]) hv[a[r]]=0;
			}
			else{
				ans=ans*fac[vec[i].size()]%cys;
				for(auto r:vec[i]) ans=ans*ds[++num[a[r]]]%cys;
				for(auto r:vec[i]) num[a[r]]=0;
			}
		}
	}
	printf("%lld\n",ans);
}

int calc(vector<int> vc){
	map<vector<int>,int> bk;
	int ret=0;
	for(auto r:vis){
		vector<int> tmp(6);
		for(int i=0;i<6;i++) tmp[r.fi[i]-1]=vc[i];
		if(!bk.count(tmp)) bk[tmp]=1,ret++;
	}
	return ret;
}

void dfs(int u,int mx){
	if(u==6){
		vector<int> vc;
		for(int i=0;i<6;i++) vc.pb(cur[i]);
		// cerr<<"init ";
		// for(auto r:vc) cerr<<r<<' ';
		res[vc]=calc(vc);
		// cerr<<res[vc]<<endl;
		return;
	}
	cur[u]=mx+1;
	dfs(u+1,mx+1);
	for(int i=1;i<=mx;i++) cur[u]=i,dfs(u+1,mx);
}

void init(){
	fac[0]=inv[0]=1;
	for(int i=1;i<=N;i++) fac[i]=fac[i-1]*i%cys;
	inv[N]=qpow(fac[N],cys-2);
	for(int i=N-1;i>=1;i--) inv[i]=inv[i+1]*(i+1)%cys;
	for(int i=1;i<=N;i++) ds[i]=fac[i-1]*inv[i]%cys;
	queue<vector<int> > q;
	q.push(vector<int>{1,2,3,4,5,6});
	vis[vector<int>{1,2,3,4,5,6}]=1;
	while(!q.empty()){
		vector<int> t=q.front(); q.pop();
		vector<int> to=t;
		int tmp=t[3];
		for(int i=3;i>=1;i--) to[i]=to[i-1];
		to[0]=tmp;
		if(!vis.count(to)) vis[to]=1,q.push(to);
		to=t;
		tmp=t[5];
		for(int i=5;i>=3;i--) to[i]=to[i-1];
		to[2]=tmp;
		if(!vis.count(to)) vis[to]=1,q.push(to);
	}
	dfs(0,0);
}

int main(){
	init();
	int T=readint();
	while(T--) solve();
	return 0;
}