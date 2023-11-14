#include <bits/stdc++.h>
using namespace std;
#define rep(i,h,t) for (int i=h;i<=t;i++)
#define dep(i,t,h) for (int i=t;i>=h;i--)
#define ll long long
#define me(x) memset(x,0,sizeof(x))
#define IL inline
#define rint register int
inline ll rd(){
    ll x=0;char c=getchar();bool f=0;
    while(!isdigit(c)){if(c=='-')f=1;c=getchar();}
    while(isdigit(c)){x=(x<<1)+(x<<3)+(c^48);c=getchar();}
    return f?-x:x;
}
char ss[1<<24],*A=ss,*B=ss;
IL char gc()
{
    return A==B&&(B=(A=ss)+fread(ss,1,1<<24,stdin),A==B)?EOF:*A++;
}
template<class T>void maxa(T &x,T y)
{
    if (y>x) x=y;
}
template<class T>void mina(T &x,T y)
{
    if (y<x) x=y;
}
template<class T>void read(T &x)
{
    int f=1,c; while (c=gc(),c<48||c>57) if (c=='-') f=-1; x=(c^48);
    while(c=gc(),c>47&&c<58) x=x*10+(c^48); x*=f;
}
const int mo=998244353;
ll fsp(int x,int y)
{
    if (y==1) return x;
    ll ans=fsp(x,y/2);
    ans=ans*ans%mo;
    if (y%2==1) ans=ans*x%mo;
    return ans;
}
struct cp {
    ll x,y;
    cp operator +(cp B)
    {
        return (cp){x+B.x,y+B.y};
    }
    cp operator -(cp B)
    {
        return (cp){x-B.x,y-B.y};
    }
    ll operator *(cp B)
    {
        return x*B.y-y*B.x;
    }
    int half() { return y < 0 || (y == 0 && x < 0); }
};
const int N=1e6;
const int G=3;
int f[N],g[N],n;
struct fft{
  int l,n,m;
  int a[N],b[N],inv[N];
  int C[N],D[N];
  fft()
  {
    inv[0]=inv[1]=1;
    rep(i,2,N-1) inv[i]=(1ll*inv[mo%i]*(mo-(mo/i)))%mo; 
  }
  IL void clear()
  {
      rep(i,0,n) a[i]=b[i]=0;
  }
  int ppr[N];
  inline void ntt_init(){
    int lg=0,x;
    for (x=1;x<=m;x*=2) lg++;
    n=x;
    ppr[0]=1;ppr[x]=fsp(31,1<<(21-lg));
    for(int i=x>>1;i;i>>=1) ppr[i]=1ll*ppr[i<<1]*ppr[i<<1]%mo;
    for(int i=1;i<x;i++) ppr[i]=1ll*ppr[i&(i-1)]*ppr[i&-i]%mo;
  }
  inline int del(const int x){
    return x>=mo?x-mo:x;
  }
  inline void DIF(int *f,const int x){
    int len,hl,uni,*s,*i,*w;
    for(len=x,hl=x>>1;hl;len=hl,hl>>=1){
        for(s=f,w=ppr;s!=f+x;s+=len,w++){
            for(i=s;i<s+hl;i++){
                uni=1ll**(i+hl)**w%mo;
                *(i+hl)=del(*i+mo-uni);
                *i=del(*i+uni);
            }
        }
    }
  }
  inline void DIT(int *f,const int x){
    int len,hl,uni,*s,*i,*w;
    for(len=2,hl=1;len<=x;hl=len,len<<=1){
        for(s=f,w=ppr;s!=f+x;s+=len,w++){
            for(i=s;i!=s+hl;i++){
                uni=*i;
                *i=del(uni+*(i+hl));
                *(i+hl)=1ll*(uni+mo-*(i+hl))**w%mo;
            }
        }
    }
    reverse(f+1,f+x);int invx=mo-(mo-1)/x;
    for(i=f;i!=f+x;i++) *i=1ll**i*invx%mo;
  }
  __int128 c[200];
  IL void getcj(int *A,int *B,int len)
  {
  	if (len<=64)
    {
    	rep(i,0,len)
    	  rep(j,0,len)
    	    c[i+j]+=1ll*A[i]*B[j];
		rep(i,0,2*len) B[i]=c[i]%mo,c[i]=0;
		return; 
    }
    m=len*2; ntt_init();
    for (int i=0;i<=len;i++) a[i]=(A[i]%mo+mo)%mo,b[i]=(B[i]%mo+mo)%mo;
    DIF(a,n); DIF(b,n);
    for(int i=0;i<n;i++) a[i]=1ll*a[i]*b[i]%mo;
    DIT(a,n);
    for (int i=0;i<=m;i++) B[i]=a[i];
    clear();
  }
  IL void getinv(int *A,int *B,int len)
  {
    if (len==1) { B[0]=fsp(A[0],mo-2); return; }
    getinv(A,B,(len+1)>>1);
    m=len*2; ntt_init();
    for (int i=0;i<=len;i++) a[i]=A[i],b[i]=B[i];
    DIF(a,n); DIF(b,n);
    for (int i=0;i<n;i++) a[i]=1ll*a[i]*b[i]%mo*b[i]%mo;
    DIT(a,n);
    for (int i=0;i<=len;i++) B[i]=((2*B[i]-a[i])%mo+mo)%mo; 
    clear();
  }
}F;
int sum[N],a[N],b[N];
int gg=0;
vector<int> operator *(const vector<int> &t1,const vector<int> &t2)
{
	int n1=t1.size(),n2=t2.size();
	int n=max(n1,n2); 
	gg+=n; 
	rep(i,0,n1-1) a[i]=t1[i];
	rep(i,0,n2-1) b[i]=t2[i];
	F.getcj(a,b,n);
	int gg=0;
	dep(i,2*n+2,0)
	{
	  gg=i;
	  if (b[i]!=0) break;
    }
    vector<int> c; 
    dep(i,gg,0) c.push_back(b[i]);
    reverse(c.begin(),c.end());
	rep(i,0,2*n) a[i]=b[i]=0;
	return c;
}
struct re{
	vector<int> a[3][3];
}M[N];
void operator += (vector<int> &a,const vector<int> &b){
	if (b.size()>a.size()) a.resize(b.size());
	for (int i=0,si=b.size();i<si;i++) a[i]=(a[i]+b[i])%mo;
}
#define mid ((h+t)/2)
re solve1(int h,int t)
{	
	if (h==t) return M[h];
	re z;
	re x=solve1(h,mid),y=solve1(mid+1,t);
	rep(i,0,2)
	  rep(j,0,2)
	    rep(k,0,2)
	    {
	      z.a[i][k]+=x.a[i][j]*y.a[j][k];
	    }
	return z;
}
re solve3(int h,int t)
{
	if (h==t) return M[h];
	re z;
	re x=solve3(h,mid),y=solve3(mid+1,t);
	rep(i,0,2)
	  rep(j,0,2)
	    rep(k,0,2)
	    {
	      z.a[i][k]+=y.a[i][j]*x.a[j][k];
	    }
	return z;
}
re ans;
int main()
{ 
//    freopen("1.out","w",stdout);
	ios::sync_with_stdio(false);
	int x;
	cin>>n>>x; n++;
	rep(i,0,2*n+x)
	{
	  M[i].a[0][1]={1}; M[i].a[1][2]={1}; M[i].a[2][0]={0,i-x}; M[i].a[2][2]={1};
    }
    auto ans1=solve1(0,2*n);
	auto ans2=solve1(2*n+1,2*n+x);
	rep(i,0,2)
	  rep(j,0,2)
	    rep(k,0,2)
		  ans.a[i][k]+=ans1.a[i][j]*ans2.a[j][k];
	vector<int> po;
    po+=ans.a[2][0];
	po+=ans.a[2][1]; 
	po+=ans.a[2][2];
	sum[0]=-1;
	rep(i,1,n+1) if (i<po.size()) sum[i]=po[i];
	rep(i,1,n+1) sum[i]=(-sum[i]%mo+mo)%mo;
	F.getinv(sum,f,n+3);
	rep(i,0,n+1) if (i%2==0) f[i]=-f[i];
	rep(i,0,2*n)
	{
		M[i].a[0][1]={1}; M[i].a[1][2]={1}; M[i].a[2][0]={0,max(i,1)}; M[i].a[2][2]={1};
	}
	ans=solve3(0,2*n);po.clear();
	po+=ans.a[2][0];
    memset(sum,0,sizeof(sum));
    rep(i,1,n+1) if (i<po.size()) sum[i]=po[i];
    rep(i,1,n+1) if (i%2==0) sum[i]=(-sum[i]%mo+mo)%mo;
    F.getcj(sum,f,n+3);
	po.clear();
    po+=ans1.a[2][0];
    memset(sum,0,sizeof(sum));
    rep(i,1,n+1) if (i<po.size()) sum[i]=po[i];
    rep(i,1,n+1) if (i%2==0) sum[i]=(-sum[i]%mo+mo)%mo;
    F.getcj(sum,f,n+3);
    cout<<(f[n+1]%mo*fsp(n*2-x,mo-2)%mo+mo)%mo<<endl;
	return 0;
}
/*
10000 10000
*/
