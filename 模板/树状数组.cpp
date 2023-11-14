#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int n,m,t[500005];
int lowbit(int x){
	return x&-x;
}
int add(int x,int k){
	while(x<=n){
		t[x]+=k;
		x+=lowbit(x);
	}
}
int sum(int x){
	int ans=0;
	while(x>0){
		ans+=t[x];
		x-=lowbit(x);
	}
	return ans;
}
int main(){
    cin>>n>>m;
    int h,a[200005];
    a[0]=0;
    for(int i=1;i<=n;i++){
    cin>>a[i];
        add(i,a[i]-a[i-1]);
	}
    for(int i=1;i<=m;i++){
        int a,b,c,d;
        cin>>a;
        if(a==1){
            cin>>b>>c>>d;
            add(b,d);
            add(c+1,-d);
		}
        else{
            cin>>b;
            cout<<sum(b)<<endl;
		}
	}
    
}
