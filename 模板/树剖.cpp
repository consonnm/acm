#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int mp[200005],n,q,r,x,y,z,k,a[2000005],m,deep[200005],f[2000005],head[200005],flag,son[200005],size[2000005],top[200005],dfn[2000005],tim;
struct node{
	int to,next;
}p1[2000005];
struct node1{
	int l,r,sum,lz;
}t[2000005];
void add(int x,int y){
	p1[++flag].to=y;
	p1[flag].next=head[x];
	head[x]=flag;
}
void dfs1(int x,int fa){
	size[x]=1;
	deep[x]=deep[fa]+1;
	f[x]=fa;
	for(int i=head[x];i;i=p1[i].next){
		if(p1[i].to==fa)
			continue;	
		dfs1(p1[i].to,x);
		size[x]+=size[p1[i].to];
		if(size[p1[i].to]>=size[son[x]])
			son[x]=p1[i].to;
	}
}
void dfs2(int x,int fa){
	top[x]=fa;
	dfn[x]=++tim;
	mp[dfn[x]]=a[x];
	if(!son[x])
		return;
	dfs2(son[x],fa);
	for(int i=head[x];i;i=p1[i].next){
		if(p1[i].to==son[x]||p1[i].to==f[x])
			continue;
		dfs2(p1[i].to,p1[i].to);
	}
}
void build(int p,int l,int r){
	t[p].l=l,t[p].r=r,t[p].lz=0;
	if(l==r){
		t[p].sum=mp[l]%m;
		return;
	}
	int mid=(l+r)>>1;
	build(p<<1,l,mid);
	build(p<<1|1,mid+1,r);
	t[p].sum=(t[p<<1].sum+t[p<<1|1].sum)%m;
}
void pushdown(int p){
	t[p<<1].sum=(t[p<<1].sum+(t[p<<1].r-t[p<<1].l+1)*t[p].lz)%m;
	t[p<<1|1].sum=(t[p<<1|1].sum+(t[p<<1|1].r-t[p<<1|1].l+1)*t[p].lz)%m;
	t[p<<1].lz+=t[p].lz;
	t[p<<1|1].lz+=t[p].lz;
	t[p].lz=0;
}
void add(int p,int l,int r,int v){
	if(t[p].l>r||t[p].r<l)
		return;
	if(t[p].r<=r&&t[p].l>=l){
		t[p].sum=(t[p].sum+v*(t[p].r-t[p].l+1))%m;
		t[p].lz+=v;
		return;
	}
	pushdown(p);
	if(t[p<<1].r>=l)
		add(p<<1,l,r,v);
	if(t[p<<1|1].l<=r)
		add(p<<1|1,l,r,v);
	t[p].sum=(t[p<<1].sum+t[p<<1|1].sum)%m;
}
int query(int p,int l,int r){
	if(t[p].r<l||t[p].l>r)
		return 0;
	if(t[p].l>=l&&t[p].r<=r)
		return t[p].sum;
	pushdown(p);
	int re=0;
	if(t[p<<1].r>=l)
		re=(re+query(p<<1,l,r))%m;
	if(t[p<<1|1].l<=r)
		re=(re+query(p<<1|1,l,r))%m;
	return re;
}
void tadd(int x,int y,int v){
	while(top[x]!=top[y]){
		if(deep[top[x]]<deep[top[y]])
			swap(x,y);
		add(1,dfn[top[x]],dfn[x],v);
		x=f[top[x]];	
	}
	if(deep[x]<deep[y])
		swap(x,y);
	add(1,dfn[y],dfn[x],v);
}
void tsum(int x,int y){
	int ans=0;
	while(top[x]!=top[y]){
		if(deep[top[x]]<deep[top[y]])
			swap(x,y);
		ans=(ans+query(1,dfn[top[x]],dfn[x]))%m;
		x=f[top[x]];	
	}
	if(deep[x]<deep[y])
		swap(x,y);
	ans=(ans+query(1,dfn[y],dfn[x]))%m;
	cout<<ans<<endl;
}
int main(){
	cin>>n>>q>>r>>m;
	for(int i=1;i<=n;i++){
		cin>>a[i];
	}
	for(int i=1;i<n;i++){
		cin>>x>>y;
		add(x,y);
		add(y,x);
	}
	dfs1(r,0);
	dfs2(r,r);
	build(1,1,n);
	for(int i=1;i<=q;i++){
		cin>>k;
		if(k==1){
			cin>>x>>y>>z;
			tadd(x,y,z);
		}
		else if(k==2){
			cin>>x>>y;
			tsum(x,y);
		}
		else if(k==3){
			cin>>x>>z;
			add(1,dfn[x],dfn[x]+size[x]-1,z);
		}
		else{
			cin>>x;
			cout<<query(1,dfn[x],dfn[x]+size[x]-1)%m<<endl;
		}

	}
}
