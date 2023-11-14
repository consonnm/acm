#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 666623333
ll exgcd(ll a,ll b,ll &x,ll &y){
    if(b==0){
        x=1,y=0;
        return a;
    }
    ll d=exgcd(b,a%b,x,y);
    ll temp=x;
    x=y;
    y=temp-a/b*x;
    return d;
}
int main(){
	int T;
    cin>>T;
    while(T--){
        ll a,b,c;
        cin>>a>>b>>c;
        ll x,y;
        ll g=exgcd(a,b,x,y);
        if(c%g!=0) cout<<-1<<endl;
        else{
            x*=c/g,y*=c/g;
            //x=x+n*b/g,y=y+n*a/g;
            ll p=b/g,q=a/g;
            ll x1=(x%p+p)%p;
            ll y1=(y%q+q)%q;
            if(!x1) x1+=p;
            if(!y1) y1+=q;
            if(c-x1*a<=0&&c-y1*b<=0){
                cout<<x1<<" "<<y1<<endl;
            }
            else{
            	ll xmx=-1,xmi=1e9,ymx=-1,ymi=1e9;
            	if(c-x1*a>0) xmx=xmi=x1,ymx=ymi=(c-x1*a)/b;
            	if(c-y1*b>0){
            		xmx=max(xmx,(c-y1*b)/a);
					xmi=min(xmi,(c-y1*b)/a);
					ymx=max(ymx,y1);
					ymi=min(ymi,y1);
				} 
				cout<<(xmx-xmi)/p+1<<" "<<xmi<<" "<<ymi<<" "<<xmx<<" "<<ymx<<endl;
			}
        }
    }
}
