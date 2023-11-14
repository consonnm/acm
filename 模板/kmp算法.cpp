#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 998244353
int kmp_next[200005];
int main(){
	string s1,s2;
	cin>>s1>>s2;//s1 text string s2 pattern string
	int n1=s1.size(),n2=s2.size();
	s1=" "+s1,s2=" "+s2;
	int j=0;
	for(int i=2;i<=n2;i++){
		while(j&&s2[i]!=s2[j+1]) j=kmp_next[j];
		if(s2[j+1]==s2[i]) j++;
		kmp_next[i]=j;
	}
	j=0;
	for(int i=1;i<=n1;i++){
		while(j&&s2[j+1]!=s1[i]) j=kmp_next[j];
		if(s2[j+1]==s1[i]) j++;
		if(j==n2){
			cout<<i-n2+1<<endl;
			j=kmp_next[j];
		} 
	}
	for(int i=1;i<=n2;i++) cout<<kmp_next[i]<<" ";
	
}
