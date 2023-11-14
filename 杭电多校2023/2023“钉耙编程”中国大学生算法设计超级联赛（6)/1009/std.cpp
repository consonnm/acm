#include<bits/stdc++.h>
#define N 23009
#define inf 2e9+5
using namespace std;
typedef long long ll;
queue<int>q; 
int tot=1,head[N],dis[N],n,du[N],pre[N],fl[N],ans;
int deep[N],cur[N],m,f[N],c[N],d[N],lim[N],st[N],top,L[N],R[N],aa[N];
bool vis[N];
set<int>s[N];
inline ll rd(){
	ll x=0;char c=getchar();bool f=0;
	while(!isdigit(c)){if(c=='-')f=1;c=getchar();}
	while(isdigit(c)){x=(x<<1)+(x<<3)+(c^48);c=getchar();}
	return f?-x:x;
} 
struct edge{int n,to,l,f;}e[600009];
inline void add(int u,int v,int l,int f=0){
	e[++tot].n=head[u];e[tot].to=v;head[u]=tot;e[tot].l=l;e[tot].f=f;
	e[++tot].n=head[v];e[tot].to=u;head[v]=tot;e[tot].l=0;e[tot].f=-f;
}
inline bool spfa(int s,int t){
    memset(dis,0x3f,sizeof(dis));
    dis[s]=0;q.push(s);fl[s]=inf;
	while(!q.empty()){
		int u=q.front();q.pop();vis[u]=0;
		for(int i=head[u];i;i=e[i].n){
			int v=e[i].to;
			if(dis[v]>dis[u]+e[i].f&&e[i].l){
				dis[v]=dis[u]+e[i].f;pre[v]=i;fl[v]=min(fl[u],e[i].l);
				if(!vis[v]){vis[v]=1;q.push(v);}
			}
		}
	} 
	return dis[t]!=0x3f3f3f3f;
}
inline void calc(int s,int t){
    int x=t;
    while(x!=s){
    	int i=pre[x];
    	e[i].l-=fl[t];e[i^1].l+=fl[t];x=e[i^1].to;
	}
	ans+=dis[t]*fl[t];
}
int find(int x){
	return f[x]=f[x]==x?x:find(f[x]);
}
struct nd{
	int u,v,w,id;
	inline bool operator <(const nd &b)const{
		return w<b.w;
	}
}b[N];
void sol(){
	cin>>n>>m;
	for(int i=1;i<=10000;++i)f[i]=i;
	int haha=0;
	for(int i=1;i<=n;++i){
		cin>>c[i]>>d[i]>>lim[i];
		haha+=min(c[i],d[i]);
		ans+=d[i];
	}
	for(int i=1;i<=m;++i){
		cin>>b[i].u>>b[i].v>>b[i].w;
		b[i].id=i;
	}
	for(int i=1;i<=m;++i){
		cin>>L[i];
	}
	for(int i=1;i<=m;++i){
		cin>>R[i];
	}
	sort(b+1,b+m+1);
	for(int i=1;i<=n;++i){
		if(c[i]>=d[i])
			add(i,i+n,1,c[i]-d[i]);
		else{
			add(i+n,i,1,d[i]-c[i]);
			ans+=c[i]-d[i];
			du[i]++;
			du[i+n]--;
		}
		f[i]=i+n;
		s[i+n].insert(i);
	}
	int now=n*2;
	int g=0;
	for(int i=1;i<=m;++i){
		int xx=find(b[i].u),yy=find(b[i].v);
		if(xx!=yy){
			++now;
			add(xx,now,inf,0);
			add(yy,now,inf,0);
			f[xx]=now;f[yy]=now;
			if(s[xx].size()<s[yy].size())swap(s[xx],s[yy]);
			s[now]=s[xx];
			for(auto x:s[yy]){
				s[now].insert(x);
			}
			for(auto x:s[now]){
				if(lim[x]<b[i].w){
					add(now,x,1);
					st[++top]=x;
				}
			}
			while(top)s[now].erase(st[top--]);
			int r=min((int)s[now].size(),L[b[i].id]),l=max(0,(int)s[now].size()-R[b[i].id]);
			if(l>r)g=1;
			if(L[b[i].id]+R[b[i].id]<s[now].size())g=1;
			++now;
			du[now-1]+=l;
			du[now]-=l;
			add(now-1,now,r-l,0);
			f[now-1]=now;
			s[now]=s[now-1];
		}
	}
	for(auto x:s[now]){
		add(now,x,1);
	}
	int S=now+1,T=now+2;
	for(int i=1;i<=now;++i){
		if(du[i]<0){
			add(S,i,-du[i]);
		}
		if(du[i]>0){
			add(i,T,du[i]);
		}
	}
	while(spfa(S,T))calc(S,T);
	for(int i=head[S];i;i=e[i].n){
		if(e[i].l)g=1;
	}
	memset(aa,-1,sizeof(aa));
	for(int u=1;u<=n;++u){
		for(int i=head[u];i;i=e[i].n){
			int v=e[i].to;
			if(v==u+n){
				if(e[i].l)aa[u]=1;
				else aa[u]=0;
			}
		}
	}
	cout<<ans<<endl; 
	for(int i=1;i<=n;++i){
		if(aa[i]==0)ans-=c[i];
		else ans-=d[i];
	}
	if(ans!=0)while(1);
	for(int i=1;i<=n;++i)f[i]=i,f[i+n]=i+n,s[i].insert(i);
	int nm=n;
	for(int i=1;i<=m;++i){
		int xx=find(f[b[i].u]),yy=find(f[b[i].v]);
		if(xx!=yy){
			++nm;
			if(s[xx].size()<s[yy].size())swap(xx,yy);
			s[nm]=s[xx];
			for(auto x:s[yy]){
				s[nm].insert(x);	
			}
			f[xx]=nm;f[yy]=nm;
			for(auto x:s[nm]){
				if(lim[x]<b[i].w){
					st[++top]=x;
				}
			}
			while(top)s[nm].erase(st[top--]);
			for(auto x:s[nm]){
				if(aa[x]==0)L[b[i].id]--;
				else R[b[i].id]--;
			}
		}
	}
	for(int i=1;i<=m;++i){
		if(L[i]<0||R[i]<0)while(1);
	}
	tot=1;ans=top=0;
	memset(head,0,sizeof(head));
	memset(dis,0,sizeof(dis));
	memset(du,0,sizeof(du));
	memset(pre,0,sizeof(pre));
	memset(fl,0,sizeof(fl));
	memset(vis,0,sizeof(vis));
	memset(deep,0,sizeof(deep));
	memset(cur,0,sizeof(cur));
	memset(f,0,sizeof(f));
	memset(c,0,sizeof(c));
	memset(d,0,sizeof(d));
	memset(lim,0,sizeof(lim));
	memset(st,0,sizeof(st));
	memset(L,0,sizeof(L));
	memset(R,0,sizeof(R));
	memset(aa,0,sizeof(aa));
	for(int i=1;i<=10000;++i){
		s[i].clear();
	}
    return;
}
int main(){
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	int T;
	cin>>T;
	while(T--){
		sol();
	}
	return 0;
}

