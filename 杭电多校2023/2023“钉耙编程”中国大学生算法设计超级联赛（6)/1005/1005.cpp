#include<bits/stdc++.h>
#define N 1009
using namespace std;
typedef long long ll;
const ll mod=1e9+7;
int n,m,T;
ll sum[N][N],tag[N][N];
ll a[N][N],b[N][N];
void MOD(ll &x){
	x=x>=mod?x-mod:x;
}
struct nd{
	ll x,y;
	inline nd operator +(const nd &b)const{
		return nd{(x+b.x+mod)%mod,y+b.y};
	}
};
nd op[N][N],op1[N][N],op2[N][N];
int main(){
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	cin>>T;
	while(T--){
		cin>>n>>m;
		for(int i=1;i<=n;++i)
			for(int j=1;j<=m;++j){
				cin>>a[i][j]; 
			}
		for(int i=1;i<=n;++i)
			for(int j=1;j<=m;++j){
				cin>>b[i][j]; 
			}
		for(int i=1;i<=n;++i){
			for(int j=1;j<=m;++j)if(a[i][j]){
				tag[i][j]=min(tag[i-1][j-1],min(tag[i][j-1],tag[i-1][j]))+1;
				op[i-tag[i][j]][j-tag[i][j]]=op[i-tag[i][j]][j-tag[i][j]]+nd{tag[i][j],1};
				op1[i-tag[i][j]][j]=op1[i-tag[i][j]][j]+nd{tag[i][j],1};
				op2[i][j-tag[i][j]]=op2[i][j-tag[i][j]]+nd{tag[i][j],1};
				sum[i][j]+=1ll*(1+tag[i][j])*tag[i][j]/2;
				sum[i][j]%=mod;
				op[i][j]=op[i][j]+nd{0,-1};
				op1[i][j]=op1[i][j]+nd{0,-1};
				op2[i][j]=op2[i][j]+nd{0,-1};
			}
		}
		ll ans=0;
		for(int i=0;i<=n+1;++i){
			for(int j=0;j<=m+1;++j){
				if(i&&j)MOD(op[i-1][j-1].x+=mod-op[i-1][j-1].y);
				if(i)MOD(op1[i-1][j].x+=mod-op1[i-1][j].y);
				if(j)MOD(op2[i][j-1].x+=mod-op2[i][j-1].y);
				if(i&&j)op[i][j]=op[i][j]+op[i-1][j-1];
				if(i)op1[i][j]=op1[i][j]+op1[i-1][j];
				if(j)op2[i][j]=op2[i][j]+op2[i][j-1];
				(sum[i][j]+=op[i][j].x-op1[i][j].x-op2[i][j].x+2ll*mod)%=mod;
			}
		}
		for(int i=n;i>=1;--i){
			for(int j=m;j>=1;--j)sum[i][j]=sum[i-1][j-1];
		}
		for(int i=0;i<=m;++i)sum[0][i]=0;
		for(int i=0;i<=n;++i)sum[i][0]=0;
		for(int i=1;i<=n;++i){
			for(int j=1;j<=m;++j){
				sum[i][j]+=sum[i-1][j]+sum[i][j-1]-sum[i-1][j-1];
				sum[i][j]=(sum[i][j]%mod+mod)%mod;
				(ans+=sum[i][j]*b[i][j]%mod)%=mod; 
			}
		}
		for(int i=0;i<=n+1;++i){
			for(int j=0;j<=m+1;++j){
				tag[i][j]=sum[i][j]=0;
				op[i][j].x=op[i][j].y=0;
				op1[i][j].x=op1[i][j].y=0;
				op2[i][j].x=op2[i][j].y=0;
			}
		}
		cout<<ans<<endl;
	} 
    return 0;
}

