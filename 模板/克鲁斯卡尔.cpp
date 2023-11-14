#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
int pre[2000006];
int find(int x){
	if(pre[x]==x)
		return x;
	return pre[x]=find(pre[x]);
}
void add(int x,int y){
	int x1=find(x),y1=find(y);
	pre[x1]=y1;
}
struct node{
	int x,y,z;
}p[2000005];
int cmp(node a,node b){
	return a.z<b.z;
}
int main(){
	int n,m;
	cin>>n>>m;
	for(int i=1;i<=m;i++){
		cin>>p[i].x>>p[i].y>>p[i].z;
	}
	for(int i=1;i<=n;i++)
		pre[i]=i;
	int flag1=1,ans=0;
	sort(p+1,p+1+m,cmp);
	for(int i=1;i<=m;i++){
		if(find(p[i].x)!=find(p[i].y)){
			add(p[i].x,p[i].y);
			ans+=p[i].z;
			flag1++;
		}
	}
	if(flag1==n)
		cout<<ans;
	else
		cout<<"orz";
}