#include<stdio.h>
#include<string.h>
#define U unsigned
#define LL long long
#define UL U LL

#define mod 998244353u


U f[41][41][41];
U g[41][41];
U h[41][41][41];
U fact[201];
U ifact[201];

U pow(U a,U b)
{
	U ans=1;
	while(b)
	{
		if(b&1)ans=(UL)ans*a%mod;
		a=(UL)a*a%mod;
		b>>=1;
	}
	return ans;
}

U C(int n,int m)
{
	return (UL)fact[n]*ifact[m]%mod*ifact[n-m]%mod;
}

void sol()
{
	memset(f,0,sizeof f);
	memset(g,0,sizeof g);
	memset(h,0,sizeof h);

	int n,k;
	
	scanf("%d%d",&n,&k);
	U w=pow(2,k);
	fact[0]=1;
	for(int i=1;i<=200;i++)
		fact[i]=(UL)fact[i-1]*i%mod;
	ifact[200]=pow(fact[200],mod-2);
	for(int i=200;i;i--)
		ifact[i-1]=(UL)ifact[i]*i%mod;


	f[0][0][1]=f[0][1][0]=1;


	for(int i=2;i<=n;i++)
	{
		f[0][0][i]=f[0][i][0]=mod-(UL)(i-1)*f[0][i-1][0]%mod;
	}

	for(int c=1;c<=n;c++)
	{
		for(int a=0;a<=n;a++)
		{
			for(int b=0;b<=n;b++)
			{
				if(a>b)
				{
					f[c][a][b]=f[c][b][a];
					continue;
				}

				UL t=0;
				if(b>0)t+=(UL)f[c][a][b-1]*b;
				t+=(UL)f[c-1][a+1][b]*c;
				if(a>0)t+=(UL)f[c][a-1][b]*a;
				t+=(UL)(c-1)*f[c-1][a][b+1];
				if(a>0&&b>0)t+=(UL)a*b*f[c][a-1][b-1];
				t+=(UL)(a*c+b*(c-1))*f[c-1][a][b];
				if(c>=2)t+=(UL)(c-1)*(c-1)*f[c-2][a+1][b+1];

				f[c][a][b]=(mod-t%mod)%mod;
			}
		}
	}

	g[0][0]=1;

	for(int i=0;i<n;i++)
	{
		for(int j=0;j<=i;j++)
		{
			g[i+1][j+1]=(g[i+1][j+1]+(UL)(w-i-j+mod)*(w-i-j-1+mod)%mod*g[i][j])%mod;
			g[i+1][j]=(g[i+1][j]+(UL)(w-i-j+mod)*j*2%mod*g[i][j])%mod;
			if(j)g[i+1][j-1]=(g[i+1][j-1]+(UL)j*j*g[i][j])%mod;
		}
	}

	h[0][0][0]=1;

	for(int c=0;c<=n;c++)
	{
		for(int a=0;a+c<=n;a++)
		{
			for(int b=0;a+b+c<=n;b++)
			if((a+b)%2==0)
			{
				U t=0;
				for(int C=0;C<=c;C++)
					for(int A=0;A<=a;A++)
						for(int B=0;B<=b;B++)
							for(int k1=0;C+k1<=c;k1++)
								for(int k2=0;C+k1+k2<=c;k2++)
								{
									if((A+B+k1+k2)%2==0)
									{
										if(c)
										{
											if(C+k1==0)continue;
											U t2;
											t2=(UL)::C(c-1,C+k1-1)*::C(c-C-k1,k2)%mod*::C(C+k1,C)%mod*::C(a,A)%mod*::C(b,B)%mod*h[c-C-k1-k2][a-A+k2][b-B+k1]%mod*w%mod*f[C][A+k1][B+k2]%mod;
											t=(t+t2)%mod;
										}
										else if(a)
										{
											if(A==0)continue;
											t=(t+(UL)::C(a-1,A-1)*::C(b,B)%mod*h[c-C][a-A][b-B]%mod*w%mod*f[C][A][B])%mod;
										}
										else if(b)
										{
											if(B==0)continue;
											t=(t+(UL)::C(b-1,B-1)*h[c-C][a-A][b-B]%mod*w%mod*f[C][A][B])%mod;
										}
									}
									
								}
							
				if(c||a||b)
				h[c][a][b]=t;
			}
		}
	}

	U G=0;
	for(int i=0;i<=n;i++)
		G=(G+g[n][i])%mod;

	U ans=((UL)(G-h[n][0][0]+mod)*pow(w,mod-2)+h[n][0][0])%mod;

	printf("%u\n",ans);
}

int main()
{
	int T;
	scanf("%d",&T);
	while(T--)sol();
	
}
