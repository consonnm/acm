#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 666623333
string a;
ll b[100005];
string Necklace(string a,int len){
	int I=1,J=2;
	string re=a;
	a=" "+a+a;
	while(J<=len){
		for(int k=0;k<len;k++){
			if(a[I+k]<a[J+k]){
				J+=k;
				break;
			}
			if(a[I+k]>a[J+k]){
				int temp=I;
				I=J;
				J=max(J,temp+k);
				break;
			}
		}
		J++;
	}
	for(int t=0;t<len;t++){
		re[t]=a[I+t];
	}
	return re;
}
int main(){
	std::ios::sync_with_stdio(false);
	cin.tie(0); 
	cout.tie(0);
	int T;
	cin>>T;
	while(T--){
		int n,m;
		cin>>n>>m;
		for(int i=1;i<=n;i++){
			cin>>a;
			a=Necklace(a,m);
			b[i]=0;
			for(int t=0;t<m;t++){
				b[i]=(28ll*b[i]+a[t]-'a'+1)%mod;
			}

		}
		int q;
		cin>>q;
		for(int i=1;i<=q;i++){
			int x,y;
			cin>>x>>y;
			if(b[x]==b[y]) cout<<"Yes"<<endl;
			else cout<<"No"<<endl;
		}
		
		
	} 
}
