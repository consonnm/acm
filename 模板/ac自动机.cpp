#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
#define mod 998244353
int flag=1,n;
string a[200005];
struct node{
	int fail,num;
	int ch[30];
}tr[200005];
void build(){
	for(int i=1;i<=n;i++){
		int p=1;
		for(int t=0;t<a[i].size();t++){
			if(!tr[p].ch[a[i][t]-'a']) tr[p].ch[a[i][t]-'a']=++flag;
			p=tr[p].ch[a[i][t]-'a'];
		}
		tr[p].num++;
	}
}
void getFail(){
	for(int i=0;i<26;i++) tr[0].ch[i]=1;
	queue<int> q;
	q.push(1);
	tr[1].fail=0;
	while(!q.empty()){
		int p=q.front();q.pop();
		for(int i=0;i<26;i++){
			int v=tr[p].ch[i];
			int fail=tr[p].fail;
			if(!v){
				tr[p].ch[i]=tr[fail].ch[i];
				continue;
			}
			tr[v].fail=tr[fail].ch[i];
			q.push(v);
		}
	}
}
int find(string a){
	int p=1,ans=0,len=a.size();
	for(int i=0;i<len;i++){
		int v=a[i]-'a';
		int k=tr[p].ch[v];
		while(k>1&&tr[k].num!=-1){
			ans+=tr[k].num;
			tr[k].num=-1;
			k=tr[k].fail;
		}
		p=tr[p].ch[v];
	}
	return ans;
}
int main(){
	cin>>n;
	for(int i=1;i<=n;i++) cin>>a[i];
	build();
	string b;
	cin>>b;
	getFail();
	cout<<find(b);
	
}