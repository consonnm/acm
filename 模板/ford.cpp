#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
struct node{
    int x,y,v;
}p[200005];
int n,m,d[200005];
int ford(){
    for(int i=1;i<=n-1;i++)
        for(int t=1;t<=m;t++)
            d[p[t].x]=min(d[p[t].x],d[p[t].y]+p[t].v); 
    for(int t=1;t<=m;t++){
        if(d[p[t].x]>d[p[t].y]+p[t].v)
            return 0;
    }
    return 1;
}
int main(){
    cin>>n>>m;
    for(int i=1;i<=m;i++)
        cin>>p[i].x>>p[i].y>>p[i].v;
    for(int i=1;i<=n;i++)
        d[i]=1e9;
    if(!ford())
        cout<<"NO"<<endl;
    else{
        for(int i=1;i<=n;i++)
            cout<<d[i]<<" ";
    }
}
