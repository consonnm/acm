#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define mod 998244353
int pre[100000005];
int prime[100000005];
int flag=0;
void getPrime(int n){
	pre[1] = 1;
	for(int i = 2; i <= n; i++){
		if(pre[i]==0) prime[++flag]=i;
		for(int j=1;j<=flag && i*prime[j] <= n;j++) {
			pre[i*prime[j]] =1;
			if(i%prime[j]==0) break;
		}
	}
}
int main(){
	int n,q; 
	cin>>n>>q;
	getPrime(n);
	for(int i=1;i<=q;i++){
		int k;
		scanf("%d",&k);
		printf("%d\n",prime[k]);
	}
}
