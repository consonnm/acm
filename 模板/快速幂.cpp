#include <bits/stdc++.h>
using namespace std;
typedef long long ll; 
ll int qpow(ll  base,ll  power,ll q){
	ll result=1;
	while(power>0){
		if(power&1)
			result=(result*base)%q;
		power>>=1;
		base=(base*base)%q;
	}
	return result;
	
}
int main(){
	long long a,b,q;
	cin>>a>>b>>q;
	cout<<a<<"^"<<b<<" mod "<<q<<"="<<qpow(a,b,q)<<endl;
}
