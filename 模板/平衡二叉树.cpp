#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int n,flag,root,a;
struct node{
	int fa,l,r,x,bf;
}p[100005];
void l(int x,int fa){
	p[fa].r=p[x].l,p[p[x].l].fa=fa,p[x].fa=p[fa].fa,p[fa].fa=x,p[x].l=fa,p[x].bf=0,p[fa].bf=0;
	if(p[x].x<p[p[x].fa].x) p[p[x].fa].l=x;
	else p[p[x].fa].r=x;
	if(root==fa) root=x;
}
void r(int x,int fa){
	p[fa].l=p[x].r,p[p[x].r].fa=fa,p[x].fa=p[fa].fa,p[fa].fa=x,p[x].r=fa,p[x].bf=0,p[fa].bf=0;
	if(p[x].x<p[p[x].fa].x) p[p[x].fa].l=x;
	else p[p[x].fa].r=x;
	if(root==fa) root=x;
	
}
void update(int x){
	int fa=p[x].fa;
	if(p[x].x>p[fa].x){
		p[fa].bf++;
	}
	else{
		p[fa].bf--;
	}
	if(p[fa].bf==0) return;
	else if(p[fa].bf==-2||p[fa].bf==2){
		if(p[fa].bf==2&&p[x].bf==1){
			l(x,fa);
		}
		else if(p[fa].bf==-2&&p[x].bf==-1){
			r(x,fa);
		}
		else if(p[fa].bf==2&&p[x].bf==-1){
			int l1=p[x].l;
			int flag2=p[l1].bf;
			r(l1,x);
			l(l1,fa);
			if(flag2==-1) p[x].bf=1,p[l1].bf=0,p[fa].bf=0;
			else if(flag2==0)	p[x].bf=0,p[l1].bf=0,p[fa].bf=0;
			else p[l1].bf=0,p[x].bf=0,p[fa].bf=-1;
		}
		else{
			int r1=p[x].r;
			int flag2=p[r1].bf;
			l(r1,x);
			r(r1,fa);
			if(flag2==-1) p[x].bf=0,p[r1].bf=0,p[fa].bf=1;
			else if(flag2==0)	p[x].bf=0,p[r1].bf=0,p[fa].bf=0;
			else p[r1].bf=0,p[x].bf=-1,p[fa].bf=0;
		}
	}
	else{
		if(fa==root) return ;
		update(fa);
	}
}
void insert(int x,int fa){
	if(flag==0){
		p[++flag].x=x;
		root=flag;
		return ;
	}
	if(x>p[fa].x){
		if(p[fa].r==0) p[++flag].x=x,p[flag].fa=fa,p[fa].r=flag,update(flag);
		else insert(x,p[fa].r);
	}
	else{
		if(p[fa].l==0) p[++flag].x=x,p[flag].fa=fa,p[fa].l=flag,update(flag);
		else insert(x,p[fa].l);
	}
}
int main(){
	cin>>n;
	for(int i=1;i<=n;i++){
		cin>>a;
		insert(a,root);
	}
	cout<<p[root].x;
}
