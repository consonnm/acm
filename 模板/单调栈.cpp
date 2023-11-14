#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
ll a[300005];
int lmx[300005],rmx[300005],lmi[300005],rmi[300005];
stack<int> s;
int main(){
	int n;
	cin>>n;
	for(int i=1;i<=n;i++)
		cin>>a[i];
	for(int i=1;i<=n;i++){
		while(!s.empty()&&a[i]>=a[s.top()]) s.pop();
		if(s.empty()) lmx[i]=0;
		else lmx[i]=s.top();
		s.push(i);
	}
	while(!s.empty()) s.pop();
	for(int i=1;i<=n;i++){
		while(!s.empty()&&a[i]<=a[s.top()]) s.pop();
		if(s.empty()) lmi[i]=0;
		else lmi[i]=s.top();
		s.push(i);
	}
	while(!s.empty()) s.pop();
	for(int i=n;i>=1;i--){
		while(!s.empty()&&a[i]>a[s.top()]) s.pop();
		if(s.empty()) rmx[i]=n+1;
		else rmx[i]=s.top();
		s.push(i);
	}
	while(!s.empty()) s.pop();
	for(int i=n;i>=1;i--){
		while(!s.empty()&&a[i]<a[s.top()]) s.pop();
		if(s.empty()) rmi[i]=n+1;
		else rmi[i]=s.top();
		s.push(i);
	}
	ll ans=0;
	for(int i=1;i<=n;i++){
		ans+=a[i]*(i-lmx[i])*(rmx[i]-i);
		ans-=a[i]*(i-lmi[i])*(rmi[i]-i);
	}
	cout<<ans;
}