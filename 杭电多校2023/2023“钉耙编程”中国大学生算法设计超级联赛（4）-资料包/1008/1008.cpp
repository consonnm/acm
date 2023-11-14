#include<bits/stdc++.h>
#define rep(i,s,t) for(int i=s;i<=t;++i)
#define pii pair<int,int>
#define fi first
#define se second
#define pb push_back
#define mp make_pair
using namespace std;
mt19937 rnd(time(0));
const int N=1e6+11,mod=998244353;
vector<int>d[N];
int n,tot,pr0,ans=1; 
char s[N];
int fail[N];
int nxt[N],las[N],to[N],phi[N],pr[N];
bool vis[N];
int S[N];
void add(int x,int y){
	nxt[++tot]=las[x];
	las[x]=tot;
	to[tot]=y;
}
void shai_fa(int n){
	for(int i=1;i<=n;++i)
		for(int j=i;j<=n;j+=i)
			d[j].pb(i);
	phi[1]=1;
	pr0=0;
	rep(i,2,n){
		if(!vis[i]){
			pr[++pr0]=i;
			phi[i]=i-1;
		}
		for(int j=1;j<=pr0&&pr[j]*i<=n;++j){
			vis[pr[j]*i]=1;
			if(i%pr[j]==0){
				phi[i*pr[j]]=phi[i]*pr[j];
				break;
			}
			phi[i*pr[j]]=phi[i]*(pr[j]-1);
		}
	}
}
void dfs(int x){
	int sum=1;
	for(auto y:d[x]){
		sum=(sum+1ll*S[y]*phi[y])%mod;
		++S[y];
	}
	ans=1ll*ans*sum%mod;
	for(int e=las[x];e;e=nxt[e])
		dfs(to[e]);
	for(auto y:d[x]){
		--S[y];
	}
}
void solve(){
	scanf("%s",s+1);
	n=strlen(s+1);
	int j=0;
	rep(i,2,n){
		while(j&&s[i]!=s[j+1])
			j=fail[j];
		if(s[i]==s[j+1])
			++j;
		fail[i]=j;
	}
	rep(i,1,n)
		add(fail[i],i);
	dfs(0);
	printf("%d\n",ans);
	ans=1;
	tot=0;
	rep(i,0,n)
		fail[i]=las[i]=0;
}
int main(){
	int size(512<<20); // 512M
	__asm__ ( "movq %0, %%rsp\n"::"r"((char*)malloc(size)+size));
	shai_fa(1e6);
	int T;
	scanf("%d",&T);
	while(T--)
		solve();
	exit(0);
}
