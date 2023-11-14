#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
double dist(double x1,double y1,double x2,double y2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
double nl,p[200005][2],dis[50005],via[200005];
double ans=0;
int n; 
double prim(){
    for(int i=1;i<=n;i++){
    	double minn=INF;
		int pos;
        for(int t=1;t<=n;t++){
            if((!via[t])&&dis[t]<minn)
                minn=dis[t],pos=t;
        }
        ans+=minn;via[pos]=1;
        for(int t=1;t<=n;t++)
            dis[t]=min(dis[t],dist(p[pos][0],p[pos][1],p[t][0],p[t][1]));
    }
}
int main(){
    cin>>n;
    for(int i=1;i<=n;i++){
        cin>>p[i][0]>>p[i][1];
        if(i!=1)
        	dis[i]=INF;
    }
    prim();
    printf("%.2lf",ans);
    
    
}
