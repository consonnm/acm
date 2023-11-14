#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
int p[505][505],vis[505],mat[505],n,m,e,a,b,ans;
int dfs(int x){
	for(int i=1;i<=m;i++){
		if(!vis[i]&&p[x][i]){
			vis[i]=1;
			if(!mat[i]||dfs(mat[i])){
				mat[i]=x;
				return 1;
			}
		}
	}
	return 0;
}
int main(){
    cin>>n>>m>>e;
    for(int i=1;i<=e;i++){
    	cin>>a>>b;
    	p[a][b]=1;
	}
	for(int i=1;i<=n;i++){
		ans+=dfs(i);
		memset(vis,0,sizeof(vis));
	}
    cout<<ans<<endl;
}
