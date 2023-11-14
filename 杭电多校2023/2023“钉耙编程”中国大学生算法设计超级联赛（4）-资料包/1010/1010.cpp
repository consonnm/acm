#include<stdio.h>
#include<algorithm>

int ans[5][5]=
{
    {1,1,2,2,3},
    {1,1,2,1,1},
    {2,2,2,2,2},
    {2,1,2,1,1},
    {3,1,2,1,1}
};

int sol(int n,int m)
{
    if(n>m)std::swap(n,m);
    if(n==1)return m+1>>1;
    while(n>5)n-=3;
    while(m>5)m-=3;
    return ans[n-1][m-1];
}

int main()
{
    int T;
    scanf("%d",&T);
    while(T--)
    {
        int n,m;
        scanf("%d%d",&n,&m);
        printf("%d\n",sol(n,m));
    }
}
