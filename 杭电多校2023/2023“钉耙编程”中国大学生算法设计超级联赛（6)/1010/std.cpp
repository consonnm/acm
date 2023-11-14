#include<bits/stdc++.h>
#define N 100009
using namespace std;
typedef long long ll;
const int mod=1e9+7;
ll k[N][31],b[N][31];
int p[N][31];
int n,q;
int main(){
	ios::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	int T;
	cin>>T;
	while(T--){
		cin>>n>>q;
		for(int i=1;i<=n;++i){
			cin>>k[i][0];
		}
		for(int i=1;i<=n;++i){
			cin>>b[i][0];
		}
		for(int i=1;i<=n;++i){
			cin>>p[i][0];
		}
		for(int i=1;i<=30;++i)
			for(int j=1;j<=n;++j){
				p[j][i]=p[p[j][i-1]][i-1];
				k[j][i]=k[j][i-1]*k[p[j][i-1]][i-1]%mod;
				b[j][i]=(b[j][i-1]*k[p[j][i-1]][i-1]+b[p[j][i-1]][i-1])%mod;
			}
		int x,l,y;
		ll pp=0;
		while(q--){
			cin>>x>>l>>y;
			x=p[x][0];
			for(int i=30;i>=0;--i)if((1ll<<i)&l){
				y=(y*k[x][i]+b[x][i])%mod;
				x=p[x][i];
			}
			cout<<y<<endl;
			pp+=l;
		}
	} 
    return 0;
}

