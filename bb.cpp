#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define fastio std::ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
#define endl '\n'
const ll mod = 998244353;
const int N = 2e6+10;
int pi[N],n;
string s,a;
void prefix_function() {
  for (int i = 1; i < s.size(); i++) {
    int j = pi[i - 1];
    while (j > 0 && s[i] != s[j]) j = pi[j - 1];
    if (s[i] == s[j]) j++;
    pi[i] = j;
  }
}
int main() {
	while(cin>>a){
		s=a;
		reverse(s.begin(),s.end());
		prefix_function();
		int flag=0,ans=0;
		for(int i=1;i<s.size();i++){
			flag=max(pi[i],flag);
		}
		ans+=flag;
        if(s.size()>pi[pi[s.size()-1]-1]*2)
		ans+=s.size()-pi[pi[s.size()-1]-1]*2;
		cout<<ans<<endl;
	}
	return 0;
}
