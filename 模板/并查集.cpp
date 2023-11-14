#include <bits/stdc++.h>
using namespace std;
#define mod 1000000007
#define endl '\n'
#define INF 0x3f3f3f3f
const int M=1e9+7;
typedef long long ll;
struct node{
    string x;
    int y;
    friend bool operator <(node a,node b){
        return a.y>b.y;
    }
};
int n,c;
string a,b;
priority_queue <node>q;
int main(){
	cin>>n;
	for(int i=1;i<=n;i++){
		cin>>a;
		if(a=="PUT"){
			cin>>b>>c;
			q.push({b,c});
		}
		else{
			if(q.empty()) cout<<"EMPTY QUEUE!"<<endl;
			else{
				cout<<q.top().x<<endl;
                q.pop();
			}
		}
	}
} 
