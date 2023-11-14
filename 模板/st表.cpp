#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
const int N=200005;
int m,;
int n,pre[N][30],a[N];
int find(int l,int r){
    int k=log2(r-l+1);
    return max(pre[l][k],pre[r-(1<<k)+1][k]);
}
void init(){
	for(int i=1;i<=n;i++) pre[i][0]=a[i];
	for(int i=1;i<=20;i++)
        for(int t=1;t+(1<<i)-1<=n;t++)
            pre[t][i]=max(pre[t][i-1],pre[t+(1<<(i-1))][i-1]);
} 
int main(){
    cin>>n>>m;
    for(int i=1;i<=n;i++)
    	scanf("%d",&a[i]);
    for(int i=1;i<=m;i++){
        int l,r;
        scanf("%d%d",&l,&r);
        cout<<find(l,r)<<endl;
    }
}
