#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
const int N=200005; 
int main(){
    int T;
    cin>>T;
    while(T--){
    	ll n;
    	cin>>n;
    	ll ans=0;
    	ll temp=n;
    	for(int i=2;i<=sqrt(n);i++){
    		if(temp%i==0){
    			ll flag=0;
    			while(temp%i==0){
    				temp/=i;
    				flag++;
				}
				ans+=flag*i;
			}
		}
		cout<<ans<<endl;
	}
}
