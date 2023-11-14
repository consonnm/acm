#include<bits/stdc++.h>
using namespace std;
int n;
double p;
void solve(){
    scanf("%d",&n);
    printf("%.9lf ",2-2.0/n);
    if(n==2)printf("1.000000000\n");
    else printf("2.000000000\n");
}
int main(){
    int T;
    scanf("%d",&T);
    while(T--)
        solve();
    return 0;
}
