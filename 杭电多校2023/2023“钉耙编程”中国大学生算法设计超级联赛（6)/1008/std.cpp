#include <bits/stdc++.h>
using namespace std ;
const int N=3e3+7; long long s[N];
int T,n,a[N],dp[N][N],f3[N],f4[N],f5[N],f6[N],d1[N],d2[N];
pair<int,int> f1[N][N],f2[N][N];
inline int getr(int x,int r){
    while(f3[r]<=f5[r]&&f1[r][f3[r]].first>=x) f3[r]++;
    return f1[r][f3[r]].second;
}
inline int getl(int x,int l){
    while(f2[l][f4[l]].first<=x) f4[l]++;
    return f2[l][f4[l]].second;
}
inline void addr(int x,int r,int c){
    while(f3[r]<=f5[r]&&f1[r][f5[r]].second>c) f5[r]--;
    f1[r][++f5[r]]=make_pair(x,c);
}
inline void addl(int x,int l,int c){
    while(f4[l]<=f6[l]&&f2[l][f6[l]].second>c) f6[l]--;
    f2[l][++f6[l]]=make_pair(x,c);
}
int main(){
    cin>>T;
    while(T--){
        scanf("%d",&n);
        for(int i=1;i<=n;i++){
            scanf("%d",&a[i]),dp[i][i]=-(a[i]);
            f1[i][1]=f2[i][1]=make_pair(i,-(a[i]));
            f3[i]=f4[i]=f5[i]=f6[i]=1;
            s[i]=s[i-1]+a[i],d1[i]=i-1,d2[i]=i+1;
        }
        for(int L=2;L<=n;L++){
            for(int l=1;l+L-1<=n;l++){
                int r=l+L-1;
                while(s[d1[l]+1]-s[l-1]<s[r]-s[d1[l]+1]) d1[l]++;
                int x=1e9;
                if(d1[l]!=r-1) x=min(x,getl(d1[l],l));
                if(d1[l]!=l-1) x=min(x,getr(d1[l]+2,r));
                addl(r,l,-x);
                addr(l,r,-x);
                dp[l][r]=-x;
            }
        }
        if(dp[1][n]>=0){
            printf("Alice %d\n",dp[1][n]);
        }
        else{
            printf("Bob %d\n",-dp[1][n]);
        }
    }
    return 0;
}
