#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 666623333
const int N=200005;
struct node{
	int mx,mi;
	vector<int> v;
	int flag;
}p[200005];
struct node1{
	int i;
	int mx;
}id[200005];
int cmp(node1 a,node1 b){
	return a.mx<b.mx;
}
int dp[200005],a[200005];
int n;
int  t[N]; 
void add_max(int x,int k){
	for(;x<=n;x+=-x&x) t[x]=max(t[x],k);
}
int find_max(int x,int y){
	int re=0;
	while(x<=y){
		re=max(re,a[y]);
		y--;
		for(;y-(-y&y)>=x;y-=-y&y){
			re=max(re,t[y]);
		}
	}
	return re;
}
void add_sum(int x,int k){
	for(;x<=n;x+=-x&x) t[x]+=k;
}
int find_sum(int x,int y){
	int re=0;
	for(;x;x-=-x&x) re+=t[x];
	return re;
}
void clear(int x){
	for(;x<=n;x+=-x&x) t[x]=0;;
}
int main(){
	int T;
	cin>>T;
	while(T--){
		scanf("%d",&n);
		for(int i=1;i<=n;i++){
			dp[i]=0;
			p[i].v.clear();
			int k;
			scanf("%d",&k);
			p[i].flag=-1;
			for(int t=1;t<=k;t++){
				int x;
				scanf("%d",&x);
				if(p[i].flag==-1||p[i].v[p[i].flag]<x){
					p[i].v.push_back(x);
					p[i].flag++;
				}
			}
			p[i].mx=p[i].v[p[i].flag];
			id[i].i=i;
			id[i].mx=p[i].mx;
		}
		sort(id+1,id+1+n,cmp);
		int ans=0;
		for(int j=1;j<=n;j++){
			int i=id[j].i;
			for(int t=0;t<=p[i].flag;t++){
				dp[i]=max(dp[i],find_max(1,p[i].v[t]-1)+p[i].flag-t+1);
			}
			a[p[i].mx]=max(a[p[i].mx],dp[i]);
			ans=max(ans,dp[i]);
			add_max(p[i].mx,dp[i]);
		}
		for(int i=1;i<=n;i++){
			a[p[i].mx]=0;
		}
		for(int i=1;i<=n;i++) clear(p[i].mx);
		cout<<ans<<endl;
	}
}
