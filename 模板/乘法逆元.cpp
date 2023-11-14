#include <bits/stdc++.h>
#define endl '\n'
using namespace std;
typedef long long ll;
#define end '\n' 
ll a[300005];
int ex_gcd(ll a,ll b,ll &x,ll &y){
	if(b==0){
		x=1,y=0;
		return a;
	}
	int d=ex_gcd(b,a%b,x,y);
	ll x0=x,y0=y;
	x=y0;
	y=x0-(a/b)*y0;
	return d;
}
ll qpow(ll a,ll b,ll mod){
	ll result=1;
	while(b){
		if(b&1) result=(result*a)%mod;
		b>>=1;
		a=(a*a)%mod;
	}
	return result;
}
int main(){
	ll n,p;
	 p=998244353;
	cout<<(5*(qpow(9,p-2,p)))%p;
	/*cin>>n>>p;
	for(int i=1;i<=n;i++){
		//ll x,y;
		//ex_gcd(i,p,x,y);
		//x = (x % p + p) % p;
		cout<<qpow(i,p-2,p)<<endl;
	}*/
}
