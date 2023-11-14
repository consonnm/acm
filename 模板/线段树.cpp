#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
ll n,m,q,a[100005];
struct tree{
    ll l,r,sum,lza,lzm;
}t[100005*4];
void build(ll p,ll l,ll r){
    t[p].l=l,t[p].r=r;
    t[p].lzm=1,t[p].lza=0;
    if(l==r){
        t[p].sum=a[l];
        return;
    }
    ll mid=l+r>>1;
    build(p<<1,l,mid);
    build(p<<1|1,mid+1,r);
    t[p].sum=(t[p<<1].sum+t[p<<1|1].sum)%q;
}
void pushdown(ll p){
	t[p<<1].sum=(t[p<<1].sum*t[p].lzm+t[p].lza*(t[p<<1].r-t[p<<1].l+1))%q;
	t[p<<1|1].sum=(t[p<<1|1].sum*t[p].lzm+t[p].lza*(t[p<<1|1].r-t[p<<1|1].l+1))%q;
    t[p<<1].lzm=(t[p].lzm*t[p<<1].lzm)%q;
    t[p<<1].lza=(t[p].lzm*t[p<<1].lza+t[p].lza)%q;
    t[p<<1|1].lzm=(t[p].lzm*t[p<<1|1].lzm)%q;
    t[p<<1|1].lza=(t[p].lzm*t[p<<1|1].lza+t[p].lza)%q;
    t[p].lza=0;
    t[p].lzm=1;
}
void add(ll p,ll l,ll r,ll k){
	if(t[p].r<l||t[p].l>r)
		return;
	if(t[p].l>=l&&t[p].r<=r){
		t[p].sum=((t[p].r-t[p].l+1)*k+t[p].sum)%q;
		t[p].lza=(t[p].lza+k)%q;
		return;
	}
	pushdown(p);
	if(t[p<<1].r>=l)
		add(p<<1,l,r,k);
	if(t[p<<1|1].l<=r)
		add(p<<1|1,l,r,k);
	t[p].sum=(t[p<<1].sum+t[p<<1|1].sum)%q;
}
void mul(ll p,ll l,ll r,ll k){
	if(t[p].r<l||t[p].l>r)
		return;
	if(t[p].l>=l&&t[p].r<=r){
		t[p].sum=(t[p].sum*k)%q;
		t[p].lza=(t[p].lza*k)%q;
		t[p].lzm=(t[p].lzm*k)%q;
		return;
	}
	pushdown(p);
	if(t[p<<1].r>=l)
		mul(p<<1,l,r,k);
	if(t[p<<1|1].l<=r)
		mul(p<<1|1,l,r,k);
	t[p].sum=(t[p<<1].sum+t[p<<1|1].sum)%q;
}
ll search(ll p,ll l,ll r){
	if(t[p].r<l||t[p].l>r)
		return 0;
	if(t[p].r<=r&&t[p].l>=l)
		return t[p].sum;
	pushdown(p);
	ll sum=0;
	if(t[p<<1].r>=l)
		sum+=search(p<<1,l,r);
	if(t[p<<1|1].l<=r)
		sum+=search(p<<1|1,l,r);
	return sum;
	
}
int main(){
    cin>>n>>m>>q;
    for(int i=1;i<=n;i++)
        cin>>a[i];
    build(1,1,n);
    for(int i=1;i<=m;i++){
        int b,x,y,k;
        cin>>b;
        if(b==1){
        	cin>>x>>y>>k;
        	mul(1,x,y,k);
		}
		else if(b==2){
			cin>>x>>y>>k;
			add(1,x,y,k);
		}
		else{
			cin>>x>>y;
			cout<<search(1,x,y)%q<<endl;
		}
    }
}
