#include <bits/stdc++.h>
using namespace std;
long long a[100005],c[100005];
long long d,e,f,money;
int main(){
	long long n,m,i,t;
	cin>>n>>m;
	for(i=1;i<=m;i++){
		cin>>a[i];
		if(i==1)
			continue;
		c[min(a[i],a[i-1])]++;
		c[max(a[i],a[i-1])]--;
	}
	for(i=1;i<=n;i++)
		c[i]+=c[i-1];
	for(i=1;i<n;i++){
		cin>>d>>e>>f;
		money+=min(d*c[i],e*c[i]+f);
	}
	n=0;
	cout<<money;
	return 0;
}
