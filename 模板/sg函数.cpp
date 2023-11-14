#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define mod 998244353
int sg[1000005];
int main(){
	int n,k;
	sg[0]=0;
	n=1000,k=7; 
	for(int i=1;i<=k;i++){
		sg[i]=1;
	}
	for(int i=1+k;i<=n;i++){
		set<int> s;
		for(int l=1;l<i;l++){
			int r=i-k-l;
			if(l>0&&r>0) s.insert(sg[l]^sg[r]);
		}
		int mex=0;
		while(s.count(mex)) mex++;
		sg[i]=mex;
	}
	cout<<n<<" "<<k<<endl;
	for(int i=1;i<=n;i++){
		cout<<sg[i];
		if(sg[i]==0){
			cout<<endl;
		}
	}
	
}
