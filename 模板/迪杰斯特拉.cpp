#include <bits/stdc++.h>
using namespace std;
int head[100005],vis[100005],flag=0,flag1[100005];
struct node{
	int to,next,dis;
}p[500005];
struct node1{
	int x,y;
	friend bool operator < (node1 a,node1 b){
		return a.y>b.y;
	}
}; 
void add(int a,int b,int c){
	p[++flag].to=b;
	p[flag].dis=c;
	p[flag].next=head[a];
	head[a]=flag;
}
priority_queue<node1> q;
int main(){
	int n,m,a,b,c,s;
	cin>>n>>m>>s;
	for(int i=1;i<=m;i++){
		cin>>a>>b>>c;
		add(a,b,c);
	}
	for(int i=1;i<=n;i++)
		vis[i]=1e9;	
	vis[s]=0;
	q.push({s,0});
	while(!q.empty()){
		int a1=q.top().x;
		q.pop();
		if(flag1[a1])
			continue;
		flag1[a1]=1;
		for(int i=head[a1];i;i=p[i].next)
			if(vis[p[i].to]>vis[a1]+p[i].dis){
				vis[p[i].to]=vis[a1]+p[i].dis;
				q.push({p[i].to,vis[p[i].to]});
			}
	}
	for(int i=1;i<=n;i++){
		cout<<vis[i]<<" ";
	}

}