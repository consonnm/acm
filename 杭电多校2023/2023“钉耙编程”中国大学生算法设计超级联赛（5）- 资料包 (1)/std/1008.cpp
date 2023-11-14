#include<bits/stdc++.h>
using namespace std;
const int N=5e6+7,p=998244353;
int n,m,a,b,k,lim,rev[N]; long long wn[N],iwn[N];
inline int read(){
	int num=0; char g=getchar(); while(g<48||57<g) g=getchar();
	while(47<g&&g<58) num=(num<<1)+(num<<3)+g-48,g=getchar(); return num;
}
inline long long pows(long long u,int v){
	long long ans=1; while(v>0) {if(v&1) ans=ans*u%p; u=u*u%p,v=v>>1;} return ans;
}
inline void init(){
	for(int d=1;d<4e6;d=d<<1){
		long long w=pows(3,(p-1)/(d<<1)),e=pows(332748118,(p-1)/(d<<1)),wf=1,we=1;
		for(int j=0;j<d;j++,wf=wf*w%p,we=we*e%p) wn[j+d]=wf,iwn[j+d]=we;
	}
}
inline int inc(int a,int b){
	if(a+b>=p) return a+b-p; return a+b;
}
inline void ntt(vector<int>&f){
	for(int i=0;i<lim;i++) if(i<rev[i]) swap(f[i],f[rev[i]]);
	for(int d=1;d<lim;d=d<<1)
		for(int i=0,r=d<<1;i<lim;i+=r){
			for(int j=0;j<d;j++){
				long long a=f[i+j],b=f[i+j+d]*wn[j+d]%p;
				f[i+j]=inc(a,b),f[i+j+d]=inc(a,p-b);
			}
		}
}
inline void intt(vector<int>&f){
	for(int i=0;i<lim;i++) if(i<rev[i]) swap(f[i],f[rev[i]]);
	for(int d=1;d<lim;d=d<<1)
		for(int i=0,r=d<<1;i<lim;i+=r){
			for(int j=0;j<d;j++){
				long long a=f[i+j],b=f[i+j+d]*iwn[j+d]%p;
				f[i+j]=inc(a,b),f[i+j+d]=inc(a,p-b);
			}
		}
	long long inv=pows(lim,p-2);
	for(int i=0;i<lim;i++) f[i]=f[i]*inv%p;
}
inline vector<int> getinv(vector<int>f){
	int l=f.size();
	if(l==1){
		vector<int>e; e.push_back(pows(f[0],p-2)); return e;
	}
	vector<int>h; h.resize((l+1)/2);
	for(int i=0;i<(l+1)/2;i++) h[i]=f[i]; h=getinv(h); 
	lim=1,k=0;  while(lim<=2*l) lim=lim<<1,k++; h.resize(lim),f.resize(lim);
	for(int i=0;i<lim;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));
	ntt(h),ntt(f);
	for(int i=0;i<lim;i++) f[i]=((2*h[i]-1ll*h[i]*h[i]%p*f[i])%p+p)%p;
	intt(f),f.resize(l); return f;
}
inline vector<int> mul(vector<int>a,vector<int>b){
	int L=a.size()+b.size(); int A=a.size(),B=b.size();
	lim=1,k=0;  while(lim<=L) lim=lim<<1,k++; a.resize(lim),b.resize(lim);
	for(int i=0;i<lim;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));
	for(int i=A;i<lim;i++) a[i]=0;
	for(int i=B;i<lim;i++) b[i]=0;
	ntt(a),ntt(b);
	for(int i=0;i<lim;i++) a[i]=1ll*a[i]*b[i]%p;
	intt(a),a.resize(L-1); return a;
}
inline vector<int> inc(vector<int>a,vector<int>b){
	if(a.size()<b.size()) swap(a,b);
	for(int i=0;i<b.size();i++){
		a[i]=a[i]+b[i]; if(a[i]>=p) a[i]-=p; 
	}
	return a;
}
inline vector<int> sub(vector<int>a,vector<int>b){
	for(int i=0;i<b.size();i++) b[i]=p-b[i];
	if(a.size()<b.size()) swap(a,b);
	for(int i=0;i<b.size();i++){
		a[i]=a[i]+b[i]; if(a[i]>=p) a[i]-=p; 
	}
	return a;
}
inline vector<int> clr(int L){
	vector<int>w; w.resize(L);
	for(int i=0;i<L;i++) w[i]=0; return w;
}
int main(){
	init(); int T; cin>>T;
	while(T--){
		n=read(),m=read(),a=read(),b=read(),a=pows(b,p-2)*a%p;
		vector<int>f=clr(m+3),w=clr(m+3);
		for(int i=0;i<=m+2;i++) f[i]=pows(i,m); w[0]=1;
		for(int i=1;i<=m+2;i++) f[i]=inc(f[i],f[i-1]); int t=1;
		for(int i=1;i<=m+2;i++){
			t=t*pows(i,p-2)%p,f[i]=1ll*f[i]*t%p;
			if(i&1) w[i]=p-t; else w[i]=t;
		}
		f=mul(f,w); int ans=0,sum=1;
		for(int i=1;i<=m+2;i++){
			sum=1ll*sum*n%p*a%p,n--;
			ans=(ans+1ll*sum*f[i])%p;
		}
		cout<<ans<<endl;
	}
	return 0;
}