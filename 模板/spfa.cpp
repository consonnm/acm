#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int flag,head[2000005],n,m,x,a,b,via[2000005],dis[2000005],via1[200005];
struct node{
    int to,next,v;
}p[2000005];
void add(int x,int y,int v){
    p[++flag].to=y;
    p[flag].next=head[x];
    head[x]=flag;
    p[flag].v=v;
}
int spfa(){
    queue<int> q;
    via1[0]=1;
    for(int i=1;i<=n;i++)
        dis[i]=-INF;
    dis[0]=0,via[0]=1;
    q.push(0);
    while(!q.empty()){
        int x=q.front();
        q.pop();
        via1[x]=0;
        for(int i=head[x];i;i=p[i].next){
            if(dis[p[i].to]<dis[x]+p[i].v){
                dis[p[i].to]=dis[x]+p[i].v;
                via[p[i].to]++;
                if(via[p[i].to]>=n)
                    return 1;
                if(!via1[p[i].to]){
                    q.push(p[i].to);
                    via1[p[i].to]=1;
                }

            }
        }
    }
    return 0;
}
int main(){
    cin>>n>>m;
    for(int i=1;i<=m;i++){
        cin>>x>>a>>b;
        if(x==1) add(a,b,0),add(b,a,0);
        if(x==2){
        	if(a == b) {printf("-1") ; return 0 ;}
        	add(a,b,1);
		} 
        if(x==3) add(b,a,0);
        if(x==4){
        	if(a == b) {printf("-1") ; return 0 ;}
        	add(b,a,1);
		} 
        if(x==5) add(a,b,0);
    }
    for(int i=n;i>=1;i--){
        add(0,i,1);
    }
    if(spfa()){
        cout<<"-1";
        return 0;
    }
    ll ans=0;
    for(int i=1;i<=n;i++)
        ans+=dis[i];
    cout<<ans;
}