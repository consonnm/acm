#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
string a,b;
int a1[200005],b1[2000005],c1[200005];
void add(string a,string b){
	for(int i=1;i<=a.size();i++)
		a1[i]=a[a.size()-i]-'0';
	for(int i=1;i<=b.size();i++)
		b1[i]=b[b.size()-i]-'0';
	int maxx=max(a.size(),b.size());
	for(int i=1;i<=maxx;i++){
		a1[i]+=b1[i];
		a1[i+1]+=a1[i]/10;
		a1[i]%=10;
	}
	if(a1[maxx+1]>0)
		maxx++;
	for(int i=maxx;i>=1;i--)
		cout<<a1[i];
}
void mul(string a,string b){
	for(int i=1;i<=a.size();i++)
		a1[i]=a[a.size()-i]-'0';
	for(int i=1;i<=b.size();i++)
		b1[i]=b[b.size()-i]-'0';
	for(int i=1;i<=a.size();i++)
		for(int t=1;t<=b.size();t++){
			c1[i+t-1]+=b1[t]*a1[i];
		}
	int maxx=a.size()+b.size();
	for(int i=1;i<=maxx;i++){
		c1[i+1]+=c1[i]/10;
		c1[i]%=10;
	}
	while(!c1[maxx]&&maxx>1)
		maxx--;
	for(int i=maxx;i>=1;i--)
		cout<<c1[i];
	
}
int main(){
	cin>>a>>b;
	mul(a,b);

}
