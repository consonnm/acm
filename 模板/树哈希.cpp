#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 1000000007
#define fastio std::ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
int m;
const int N=200005;
int head[N],flag,base=114514,size[N],d[N],mx=1e9,n;
struct node1{
	int to,next;
}p[N];
void add(int x,int y){
	p[++flag].to=y;
	p[flag].next=head[x];
	head[x]=flag;
}
//ll Has(int x,int fa){
//	ll i,y,res=2333;
//	vector<ll > t;
//	for(i=head[x];i;i=p[i].next){
//		y=p[i].to;if(y==fa)continue;
//		t.push_back(Has(y,x));
//	}
//	sort(t.begin(),t.end());
//	for(i=0;i<t.size();++i)res=((res*base)^t[i])%mod;
//	return res;
//}
ll hsh[N];
mt19937 rnd(random_device{}());
uniform_int_distribution<ll> dist(0,ULLONG_MAX);
const ll s=dist(rnd);
ll xorshift(ll x) {
	x^=x<<13;
	x^=x>>7;
	x^=x<<17;
	return x;
}
ll Has(int x,int fa) {
	ll res=s;
	for (int i=head[x];i;i=p[i].next)
	{
		int y=p[i].to;
		if (y==fa) continue;
		res+=xorshift(Has(y,x));
	}
	return res;
}
void dfs_barycenter(int x,int fa){
	size[x]=1;
	int res=0;
	for(int i=head[x];i;i=p[i].next){
		if(p[i].to==fa) continue;
		dfs_barycenter(p[i].to,x);
		size[x]+=size[p[i].to];
		res=max(res,size[p[i].to]);
	}
	res=max(res,n-size[x]);
	d[x]=res;
	mx=min(mx,res);
}
pair<int,int> find_barycenter(int x,int fa){
	pair<int,int> pos={-1,-1};
	dfs_barycenter(1,0);
	for(int j=1;j<=n;j++){
		if(d[j]==mx){
			if(pos.first==-1) pos.first=j;
			else pos.second=j;
		} 
	}
	return pos;
}
int ans[N];
int main() {
	//fastio;
	//freopen("int.txt","r",stdin);
	//freopen("1.txt","w",stdout);
	cin>>m;
	for(int i=1;i<=m;i++){
		cin>>n;
		for(int j=1;j<=n;j++){
			int y;
			cin>>y;
			if(y) add(j,y),add(y,j);
		}
		pair<int,int> pos=find_barycenter(1,0);
		int flag1=pos.first,flag2=pos.second;
		if(flag2!=-1){
			ans[i]=min(Has(flag1,0),Has(flag2,0));
		}
		else  ans[i]=Has(flag1,0);
		flag=0,mx=1e9;
		for(int j=0;j<=n+100;j++) head[j]=0,size[j]=0,d[j]=0;
	}
	for(int i=1;i<=m;i++){
		int flag=-1;
		for(int t=1;t<=m;t++){
			if(ans[i]==ans[t]){
				cout<<t<<endl;
				flag=1;
				break;
			}
		}
	}
}
