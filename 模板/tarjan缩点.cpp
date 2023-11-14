#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int flag,n,m,x,y,head[200005],dfn[200005],low[200005],via[200005],ind,sum,ans[200005];
struct node{
    int to,next;
}p[200005*2];
void add(int x,int y){
    p[++flag].to=y;
    p[flag].next=head[x];
    head[x]=flag;
}
void tarjan(int x,int fa){
    via[x]=1,dfn[x]=low[x]=++ind;
    int flag1=0;
    for(int i=head[x];i;i=p[i].next){
        int v=p[i].to;
        if(!via[v]){
            if(x==fa)
                flag++;
            tarjan(v,fa);
            low[x]=min(low[x],low[v]);
            if(low[v]>=dfn[x]&&x!=fa&&!flag1)
                ans[++sum]=x,flag1=1;
        }
        else{
            low[x]=min(low[x],dfn[v]);
        }
    }
    if(flag>=2&&x==fa)
        ans[++sum]=fa;
}
int main(){
    cin>>n>>m;
    for(int i=1;i<=m;i++){
        cin>>x>>y;
        add(x,y);
        add(y,x);
    }
    for(int i=1;i<=n;i++){
    	flag=0;
    	if(!via[i])
    		tarjan(i,i);
	}  
    cout<<sum<<endl;
    sort(ans+1,ans+1+sum);
    for(int i=1;i<=sum;i++)
        cout<<ans[i]<<" ";
}
