#include<cmath>
#include<cstdio>
#include<unordered_map>
using namespace std;
typedef long long LL;
const int maxs=1000001,maxk=10,MOD=998244353;

int te,K;LL n;
int u[maxs+5],p[maxs+5];bool pri[maxs+5];
int C[maxk+5][maxk+5],val[maxk+5];
int lim,G[maxs+5];
unordered_map<LL,int> f;
int ans;

inline int ADD(int x,int y) {return x+y>=MOD?x+y-MOD:x+y;}
inline int MUL(int x,int y) {return (LL)x*y%MOD;}
int Pow(int w,int b) {int s;for (s=1;b;b>>=1,w=MUL(w,w)) if (b&1) s=MUL(s,w);return s;}
void Make(){
	for (int i=0;i<=maxk;i++){
		C[i][0]=1;
		for (int j=1;j<=i;j++)
			C[i][j]=ADD(C[i-1][j-1],C[i-1][j]);
	}
	for (int i=1;i<=maxk;i++) val[i]=MUL(1<<i,Pow((1<<i)-1,MOD-2));
	u[1]=1;
	for (int i=2;i<=maxs;i++){
		if (!pri[i]) p[++p[0]]=i,u[i]=MOD-1;
		for (int j=1,t;j<=p[0] && (t=i*p[j])<=maxs;j++){
			pri[t]=true;
			if (i%p[j]) u[t]=(MOD-u[i])%MOD;
			else {u[t]=0;break;}
		}
	}
}
int SumF(LL n){
	int ans=MUL(K&1?MOD-1:1,n%MOD);
	int BA=Pow(2,n%(MOD-1)),pw=BA;
	for (int j=1;j<=K;j++){
		int now=MUL(C[K][j],K-j&1?MOD-1:1);
		now=MUL(now,MUL(val[j],ADD(pw,MOD-1)));
		ans=ADD(ans,now);
		pw=MUL(pw,BA);
	}
	return ans;
}
int Sum(LL n){
	if (n<=lim) return G[n];
	if (f.count(n)) return f[n];
	int ans=SumF(n);
	for (LL l=2,r;l<=n;l=r+1){
		r=n/(n/l);
		ans=ADD(ans,MOD-MUL(Sum(n/l),(r-l+1)%MOD));
	}
	return f[n]=ans;
}
int main(){
	Make();
	for (scanf("%d",&te);te;te--){
		scanf("%lld%d",&n,&K);
		lim=pow(n,2.0/3)+1;
		for (int i=1;i<=lim;i++) G[i]=0;
		for (int i=1,pw=2;i<=lim;i++,pw=MUL(pw,2)){
			int F=Pow(pw-1,K);
			for (int j=i,k=1;j<=lim;j+=i,k++)
				G[j]=ADD(G[j],MUL(u[k],F));
		}
		for (int i=1;i<=lim;i++) G[i]=ADD(G[i],G[i-1]);
		f.clear();ans=0;
		for (LL l=1,r;l<=n;l=r+1){
			r=n/(n/l);
			ans=ADD(ans,MUL(ADD(Sum(r),MOD-Sum(l-1)),MUL(n/l%MOD,n/l%MOD)));
		}
		printf("%d\n",ans);
	}
	return 0;
}