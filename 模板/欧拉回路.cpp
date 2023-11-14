#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
int in[2000005],out[2000005],n,m,x,y,flag,flag1[2],pos=1,del[2000005];
vector<int> p[2000005];
stack<int> s;
void dfs(int x){
    for(int i=del[x];i<p[x].size();i=del[x]){
        del[x]++;
        dfs(p[x][i]);
    }
    s.push(x);

}
int main(){
    cin>>n>>m;
    for(int i=1;i<=m;i++){
        cin>>x>>y;
        p[x].push_back(y);
        in[y]++;
        out[x]++;
    }
    for(int i=1;i<=n;i++)
    	sort(p[i].begin(),p[i].end());
    for(int i=1;i<=n;i++){
        if(in[i]!=out[i]) flag=1;
        if(out[i]-in[i]==1) flag1[0]++,pos=i;
        if(in[i]-out[i]==1) flag1[1]++;
    }
    if((flag1[0]==flag1[1]&&flag1[0]==1)||flag==0){
        dfs(pos);
        while(!s.empty()){
            cout<<s.top()<<" ";
            s.pop();
        }
    }
    else
        cout<<"NO";
}
