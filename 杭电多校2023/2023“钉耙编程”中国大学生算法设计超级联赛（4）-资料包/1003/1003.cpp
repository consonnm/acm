#include<bits/stdc++.h>
using namespace std;
typedef pair<int, int> P;
#define MAXK 1000010
#define endl '\n'
#define fastio ios_base::sync_with_stdio(false);cin.tie(0);cout.tie(0)

int t,k,buc[MAXK];
vector<P>vec;
void solve(int ans=2e9,int cnt=0){
    cin>>k;
    for(int i=1,siz;i<=k;i++){
        cin>>siz;
        for(int j=1,x;j<=siz;j++){
            cin>>x;
            vec.push_back(P(x,i));
        }
    }
    sort(vec.begin(), vec.end());
    for(int i=0,j=0;i<vec.size();i++){
        cnt += (++buc[vec[i].second]==1);
        for(;j<i&&buc[vec[j].second]>1;j++)
            buc[vec[j].second]--;
        if(cnt==k)ans = min(ans, vec[i].first-vec[j].first);
    }
    cout<<ans<<endl;
    vec.clear();
    fill(buc+1, buc+k+1, 0);
    return;
}

int main(){
    //fastio;
    cin>>t;
    while(t--)solve();
    return 0;
} 
