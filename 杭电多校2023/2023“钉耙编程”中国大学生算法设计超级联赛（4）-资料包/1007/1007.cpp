#include<bits/stdc++.h>
#define rep(i,a,b) for(int i=(a);i<=(b);i++)
#define fastio ios_base::sync_with_stdio(false);cin.tie(0);cout.tie(0)
typedef long long ll;
using namespace std;
const int MOD=998244353; 
ll ksm(ll a,ll b,ll p){
	ll ret=1;
	for(;b;b>>=1,a=(__int128)a*a%p)
		if(b&1)ret=(__int128)ret*a%p;
	return ret;
}
bool Miller_Rabin(ll p){
	if(p<2)return 0;
	if(p==2||p==3)return 1;
	ll d=p-1,r=0;
	while(!(d&1))++r,d>>=1;
	for(ll k=0;k<10;++k){
		ll a=rand()%(p-2)+2;
		ll x=ksm(a,d,p);
		if(x==1||x==p-1)continue;
		for(int i=0;i<r-1;++i){
			x=(__int128)x*x%p;
			if(x==p-1)break;
		}
		if(x!=p-1)return 0;
	}
	return 1;
}
int tot;
ll decom[100];
ll Pollard_Rho(ll x){
	ll s=0,t=0;
	ll c=(ll)rand()%(x-1)+1;
	int step=0,goal=1;
	ll val=1;
	for(goal=1;;goal<<=1,s=t,val=1){
		for(step=1;step<=goal;step++){
			t=((__int128)t*t+c)%x;
			val=(__int128)val*abs(t-s)%x;
			if(step%127==0){
				ll d=__gcd(val,x);
				if(d>1)return d;
			}
		}
		ll d=__gcd(val,x);
		if(d>1)return d;
	}
}
void solve(ll x){
	if(x<2)return;
	if(Miller_Rabin(x)){
		decom[++tot]=x;
		return;
	}
	ll p=x;
	while(p==x)p=Pollard_Rho(x);
	solve(x/p);solve(p);
}
int main(){
	fastio;
	srand(time(0));
	int T;cin>>T;
	while(T--){
		tot=0;
		ll n;cin>>n;
		assert(1<=n&&n<=(ll)(1e18)); 
		solve(n);
		sort(decom+1,decom+tot+1);
		tot=unique(decom+1,decom+tot+1)-(decom+1);
		if(tot==1)cout<<decom[1]%MOD;
		else cout<<1;
		if(T)cout<<' ';
	}
	return 0;
}
