#include<bits/stdc++.h>
#define N 200009
#define int ll 
using namespace std;
typedef long long ll;
const int mod=998244353;
int ch[N][26],l[N],fa[N];
int cnt,last,tag[N],n,tong[N],rnk[N];
int T[N],top,tot,head[N],siz[N],q,dep[N];
ll ans;
char s[N],t[N];
ll tr1[N*30],tr2[N*30],g[N],f1[N],f2[N];
int L[N*30],R[N*30],ls[N*30],rs[N*30];
inline ll rd(){
	ll x=0;char c=getchar();bool f=0;
	while(!isdigit(c)){if(c=='-')f=1;c=getchar();}
	while(isdigit(c)){x=(x<<1)+(x<<3)+(c^48);c=getchar();}
	return f?-x:x;
}
inline ll power(ll x,ll y){
	ll ans=1;x=x%mod;
	while(y){
		if(y&1)ans=ans*x%mod;
		x=x*x%mod;
		y>>=1;
	}
	return ans;
}
inline ll ni(ll x){
	return power(x,mod-2);
}
struct edge{int n,to;}e[N];
inline void MOD(ll &x){x=x>=mod?x-mod:x;}
inline void add(int u,int v){
	e[++tot].n=head[u];
	e[tot].to=v;
	head[u]=tot;
}
void upd(int &cnt,int l,int r,int x){
	if(!cnt)cnt=++top; 
	if(l==r){
		tr1[cnt]=1ll*x*x%mod;
		ls[cnt]=rs[cnt]=x;
		return;
	}
	int mid=(l+r)>>1;
	if(mid>=x)upd(L[cnt],l,mid,x);
	else upd(R[cnt],mid+1,r,x);
	MOD(tr1[cnt]=tr1[L[cnt]]+tr1[R[cnt]]);
	ls[cnt]=ls[L[cnt]]?ls[L[cnt]]:ls[R[cnt]];
	rs[cnt]=rs[R[cnt]]?rs[R[cnt]]:rs[L[cnt]];
}
void merge(int &u,int v,int l,int r){
	if(!u||!v){
		u=u+v;return; 
	} 
	int mid=(l+r)>>1;
	merge(L[u],L[v],l,mid);
	merge(R[u],R[v],mid+1,r);
	MOD(tr1[u]=tr1[L[u]]+tr1[R[u]]);
	tr2[u]=(tr2[L[u]]+tr2[R[u]]+1ll*rs[L[u]]*ls[R[u]])%mod;
	ls[u]=ls[L[u]]?ls[L[u]]:ls[R[u]];
	rs[u]=rs[R[u]]?rs[R[u]]:rs[L[u]];
}
inline ll C(ll l,ll r){
	if(l>r)return 0;
	return (l+r)*(r-l+1)/2;
}
inline void ins(int x){
	int p=last,np=++cnt;last=np;l[np]=l[p]+1;
	for(;p&&!ch[p][x];p=fa[p])ch[p][x]=np;
	if(!p)fa[np]=1;
    else{
    	int q=ch[p][x];
    	if(l[q]==l[p]+1)fa[np]=q;
    	else{
    		int nq=++cnt;l[nq]=l[p]+1;
    		memcpy(ch[nq],ch[q],sizeof(ch[q]));
    		fa[nq]=fa[q];fa[q]=fa[np]=nq;
    		for(;ch[p][x]==q;p=fa[p])ch[p][x]=nq;
    	}
    }
}
int pp=0;
void dfs(int u){
	if(tag[u]){upd(T[u],1,n,tag[u]);siz[u]++;}
	for(int i=head[u];i;i=e[i].n){
		int v=e[i].to;
		dfs(v);
		merge(T[u],T[v],1,n);
		siz[u]+=siz[v];
	} 

	f1[u]=(-tr1[T[u]]+tr2[T[u]]+1ll*(n+1)*rs[T[u]])%mod;
	f2[u]=(n+1-ls[T[u]])%mod;
	f1[u]=(f1[u]%mod+mod)%mod;
	f2[u]=(f2[u]%mod+mod)%mod;
	g[u]=(l[u]-l[fa[u]])*f1[u]-C(l[fa[u]],l[u]-1)%mod*f2[u];
	g[u]=(g[u]%mod+mod)%mod;
	pp+=siz[u];
}
void dfs2(int u){
	for(int i=head[u];i;i=e[i].n){
		int v=e[i].to;
		dep[v]=dep[u]+1;
		MOD(g[v]+=g[u]);
		dfs2(v);
	}
}
signed main(){
	int TT;
	scanf("%d",&TT);
	int aa=0;
	return 0;
	while(TT--){
		n=rd();q=rd(); 
		scanf("%s",s+1);
		last=cnt=1;
		for(int i=1;i<=n;++i){
			ins(s[i]-'a');
			tag[last]=i;
		}
		for(int i=2;i<=cnt;++i)add(fa[i],i); 
		dfs(1); 
		dfs2(1);
		while(q--){
			scanf("%s",t+1);
			int m=strlen(t+1);
			int now=1,len=0;ans=0;
			for(int j=1;j<=m;++j){
				while(now!=1&&!ch[now][t[j]-'a']){
					now=fa[now];len=l[now];
				}
				if(ch[now][t[j]-'a']){
					now=ch[now][t[j]-'a'];len++;
				}
				if(now!=1)
					ans+=g[fa[now]]+f1[now]*(len-l[fa[now]])%mod-C(l[fa[now]],len-1)%mod*f2[now]%mod;
				ans=(ans%mod+mod)%mod;
			}
			ans=(ans+mod)%mod;
			ans=ans*2*ni(m)%mod*ni(m+1)%mod;
			printf("%lld\n",ans);
		}
		memset(ch,0,sizeof(ch));
		memset(l,0,sizeof(l));
		memset(fa,0,sizeof(fa)); 
		cnt=last=ans=tot=top=0;
		memset(tag,0,sizeof(tag));
		memset(rnk,0,sizeof(rnk));  
		memset(tong,0,sizeof(tong)); 
		memset(head,0,sizeof(head));   
		memset(siz,0,sizeof(siz));
		memset(T,0,sizeof(T));    
		memset(dep,0,sizeof(dep));  
		memset(tr1,0,sizeof(tr1));  
		memset(tr2,0,sizeof(tr2));  
		memset(g,0,sizeof(g));
		memset(f1,0,sizeof(f1));    
		memset(f2,0,sizeof(f2));  
		memset(L,0,sizeof(L));  
		memset(R,0,sizeof(R));  
		memset(ls,0,sizeof(ls));  
		memset(rs,0,sizeof(rs));  
		aa++;
	}
    return 0;
}

