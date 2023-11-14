#include<stdio.h>
#include<algorithm>

struct Z
{
    int a,b;
    int operator<(Z o) const
    {
        return a+b>o.a+o.b;
    }
};
Z z[100001];

int main()
{
    int T;
    scanf("%d",&T);
    while(T--)
    {
        int n;
        scanf("%d",&n);
        for(int i=0;i<n;i++)
            scanf("%d%d",&z[i].a,&z[i].b);
        std::sort(z,z+n);
        long long sum=0;
        for(int i=0;i<n;i++)
        {
            if(i&1)sum-=z[i].b;
            else sum+=z[i].a;
        }
        printf("%lld\n",sum);
    }
    
}
