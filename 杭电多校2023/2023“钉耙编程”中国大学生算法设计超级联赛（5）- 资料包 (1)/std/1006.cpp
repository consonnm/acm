#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
const int maxn=100000;

int te,n;char s[maxn+5];
int f[maxn+5][20];

#define ID(a,b) (((a)<<2)|(b))
int main(){
	memset(f[0],192,sizeof(f[0]));
	f[0][0]=0;
	for (scanf("%d",&te);te;te--){
		scanf("%s",s+1);n=strlen(s+1);
		for (int i=1;i<=n;i++){
			int c=(s[i]=='R'?1:(s[i]=='G'?2:3));
			for (int a=0;a<=3;a++)
				for (int b=0;b<=3;b++)
					f[i][ID(a,b)]=f[i-1][ID(a,b)];
			for (int a=0;a<=3;a++)
				for (int b=0;b<=3;b++)
					if (a==b && b==c){
						int F=f[i-1][ID(a,b)]+1;
						for (int d=1;d<=3;d++) f[i][d]=max(f[i][d],F);
					} else if (a && b && c && a!=b && a!=c && b!=c){
						int F=f[i-1][ID(a,b)];
						for (int d=1;d<=3;d++)
							for (int e=1;e<=3;e++)
								f[i][ID(d,e)]=max(f[i][ID(d,e)],F);
					} else f[i][ID(b,c)]=max(f[i][ID(b,c)],f[i-1][ID(a,b)]);
		}
		int ans=0;
		for (int a=0;a<=3;a++)
			for (int b=0;b<=3;b++)
				ans=max(ans,f[n][ID(a,b)]);
		printf("%d\n",ans);
	}
	return 0;
}