#include <bits/stdc++.h>
using namespace std;
struct edge{
	int to,next;
}p[200005];
int head[200005],deep[100005],pre[100005][25],flag,n,q;
void add(int x,int y){
	p[++flag].to=y;
	p[flag].next=head[x];
	head[x]=flag;
}
void dfs(int x,int fa){
	deep[x]=deep[fa]+1;
	pre[x][0]=fa;
	for(int i=1;(1<<i)<deep[x];i++)
		pre[x][i]=pre[pre[x][i-1]][i-1];
	for(int i=head[x];i;i=p[i].next){
		if(p[i].to!=fa)
			dfs(p[i].to,x);
	}
}
int lca(int x,int y){
	if(deep[x]<deep[y])
		swap(x,y);
	for(int i=20;i>=0;i--){
		if(deep[pre[x][i]]>=deep[y])
			x=pre[x][i];
	}
	if(x==y)
		return x;
	for(int i=20;i>=0;i--){
		if(pre[x][i]!=pre[y][i]){
			x=pre[x][i];
			y=pre[y][i];
		}
	}
	return pre[x][0];
}
int main(){
	int x,y;
	cin>>n>>q;
	for(int i=1;i<n;i++){
		cin>>x>>y;
		add(x,y);
		add(y,x);
	}
	dfs(1,0);
	int a,b,c,d;
	for(int i=1;i<=q;i++){
		cin>>a>>b>>c>>d;
		int ans1=lca(a,b),ans2=lca(c,d);
		if(deep[ans1]>=deep[ans2]&&(lca(ans1,c)==ans1||lca(ans1,d)==ans1))
			cout<<"Y"<<endl;
		else if(deep[ans2]>=deep[ans1]&&(lca(ans2,a)==ans2||lca(ans2,b)==ans2))
			cout<<"Y"<<endl;
		else
			cout<<"N"<<endl;
	}

}
