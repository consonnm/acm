#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 998244353
#define fastio  std::ios::sync_with_stdio(false),cin.tie(0);cout.tie(0);
int t[3000005][65],cnt[3000005],flag,n,q;
int getnum(char a){
	if(a>='a'&&a<='z') return a-'a';
	else if(a>='A'&&a<='Z') return a-'A'+26;
	else return a-'0'+26;
}
void insert(string a){
	int p=0;
	for(int i=0;i<a.size();i++){
		int h=getnum(a[i]);
		if(!t[p][h]) t[p][h]=++flag;
		p=t[p][h];
		cnt[p]++;
	}
}
int find(string a){
	int p=0;
	for(int i=0;i<a.size();i++){
		int h=getnum(a[i]);
		p=t[p][h];
		if(p==0) return 0;
		if(i==a.size()-1) return cnt[p];
	}
}
int main(){
	int T;
	cin>>T;
	while(T--){
		cin>>n>>q;
		for(int i=0;i<=flag;i++){
			for(int j=0;j<64;j++){
				t[i][j]=0;
			}
			cnt[i]=0;
		}
			
		for(int i=1;i<=n;i++){
			string a;
			cin>>a;
			insert(a);
		}
		for(int i=1;i<=q;i++){
			string a;
			cin>>a;
			cout<<find(a)<<endl;
		}
	}
	
}
