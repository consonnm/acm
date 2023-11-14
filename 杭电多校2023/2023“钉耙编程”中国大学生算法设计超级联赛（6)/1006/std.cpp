#include <bits/stdc++.h>
using namespace std;
#define rep(i,h,t) for (int i=h;i<=t;i++)
#define dep(i,t,h) for (int i=t;i>=h;i--)
#define ll long long
#define mp make_pair 
const int N=2010;
const int M=400;
int n,a[N],cnt1[N],cnt2[N];
ll sum[N],f[N];
inline bool pd(int x)
{
    int y=sqrt(x);
    if (y*y==x) return 1;
    return 0;
}
void gg(int x,int y,int z)
{
    int now=sum[y]-sum[x-1];
    rep(i,-300,300) 
      if (pd(now+i)) f[M+i]+=z;
}
int solve()
{
    rep(i,1,n) sum[i]=sum[i-1]+a[i];
    int ans2=0;
    rep(i,1,n)
      rep(j,i,n)
        if (pd(sum[j]-sum[i-1])) ans2++;
    return ans2;
}
int main()
{
    ios::sync_with_stdio(false);
    int T;
    cin>>T;
    while (T--)
    {
        memset(cnt1,0,sizeof(cnt1));
        memset(cnt2,0,sizeof(cnt2));
        memset(sum,0,sizeof(sum));
        memset(f,0,sizeof(f));
	    cin>>n;
	    rep(i,1,n) cin>>a[i];
	    rep(i,1,n) sum[i]=sum[i-1]+a[i];
	    rep(i,1,n)
	    {
	        rep(j,0,i-1)
	          if (pd(sum[i]-sum[j])) cnt1[i]++;
	        cnt1[i]+=cnt1[i-1];
	    }
	    dep(i,n,1)
	    {
	      rep(j,i,n)
	        if (pd(sum[j]-sum[i-1])) cnt2[i]++;
	      cnt2[i]+=cnt2[i+1];
	    }
	    ll ans=0;
	    rep(i,1,n)
	    {
	        rep(j,i,n) gg(i,j,1);
	        rep(j,1,i-1) gg(j,i-1,-1);
	        ll tt=0;
	        rep(j,1,300)
	          tt=max(tt,f[j-a[i]+M]);
	        ans=max(ans,cnt1[i-1]+cnt2[i+1]+tt);
	    }
	    cout<<ans<<endl;
    }
    return 0;
}
