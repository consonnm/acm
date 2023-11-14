#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int n,m,dp[2000005],fa[200005],flag,head[200005],in[2000005],ans,r,s;
struct node{
	int to,next,x;
}p[100005]; 
void add(int x,int y,int k){
	p[++flag].to=y;
	p[flag].next=head[x];
	p[flag].x=k;
	head[x]=flag;
	in[y]++;
}
void tuopu(){
	for(int i=0;i<n;i++) dp[i]=-1e9;
	dp[r]=0;
	queue<int> q;
	q.push(r);
	while(!q.empty()){
		int x=q.front();
		q.pop();
		for(int i=head[x];i;i=p[i].next){
			if(dp[x]+p[i].x>dp[p[i].to]){
				dp[p[i].to]=dp[x]+p[i].x;
				fa[p[i].to]=x;
			}
			else if(dp[x]+p[i].x==dp[p[i].to]){
				int pa1[55],pa2[55],si1=0,si2=0,h1=fa[p[i].to],h2=x;
				while(h1!=r){
					pa1[++si1]=h1;
					h1=fa[h1];
				}
				while(h2!=r){
					pa2[++si2]=h2;
					h2=fa[h2];
				}
				int flag3=0;
				for(int j=si1,t=si2;j>=1;j--,t--){
					if(pa1[j]>pa2[t]||t==0){
						fa[p[i].to]=x;
						break;
					}  
				}
			}
			in[p[i].to]--;
			if(!in[p[i].to]&&s!=x) q.push(p[i].to);
		}
	}	
}
int main(){
	cin>>n>>m;
	for(int i=1;i<=m;i++){
		int a,b,c;
		cin>>a>>b>>c;
		add(a,b,c);
	}
	cin>>r>>s;
	tuopu();
	if(fa[s]==0) cout<<"none";
	else{
		int pa[55],si=0;
		while(s!=r){
			pa[++si]=s;
			s=fa[s];
		}
		pa[++si]=r;
		for(int i=si;i>=1;i--){
			if(i!=si) cout<<"->";
			cout<<pa[i];
		}
		
	}
} 
