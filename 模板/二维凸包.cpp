#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int n,top;
struct node{
	double x,y;
}p[2000006],stac[20000005];
double check(node a,node b,node c,node d){
	return (b.x-a.x)*(d.y-c.y)-(b.y-a.y)*(d.x-c.x);
}
double dis(node a,node b){
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)); 
}
int cmp(node a,node b){
	double temp=check(p[1],a,p[1],b);
	if(temp>0)
		return 1;
	if(temp==0&&dis(p[1],a)<dis(p[1],b))
		return 1;
	return 0;
}
int main(){
	cin>>n;
	for(int i=1;i<=n;i++){
		cin>>p[i].x>>p[i].y;
		if(i!=1&&p[i].y<p[1].y) swap(p[i],p[1]);
		else if(i!=1&&p[1].y==p[i].y&&p[i].x<p[1].x) swap(p[i],p[1]);
	}
	sort(p+2,p+1+n,cmp);
	stac[++top]=p[1];
	for(int i=2;i<=n;i++){
		while(top>1&&check(stac[top-1],stac[top],stac[top],p[i])<=0) top--;
		stac[++top]=p[i];
	}
	stac[++top]=p[1];
	double ans=0;
	for(int i=2;i<=top;i++)
		ans+=dis(stac[i],stac[i-1]);
	printf("%.2lf\n",ans);
}
