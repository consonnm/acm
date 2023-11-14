#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int main(){
	string a1,a="&#";
	getline(cin,a1);
	int ans=0,p[10005];
	for(int i=0;i<a1.size();i++){
		a+=a1[i];
		a+="#";
	}
	int m=0,r=0;
	for(int i=1;i<a.size();i++){
		p[i]=r>i?min(p[2*m-1],r-i):1;
		while(a[i+p[i]]==a[i-p[i]])
			++p[i];
		if(r<i+p[i]){
			r=i+p[i];
			m=i;
		}
		ans=max(ans,p[i]);
	}
	cout<<ans-1;
}
