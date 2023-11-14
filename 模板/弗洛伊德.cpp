#include <bits/stdc++.h>
using namespace std;
int main(){
	int n,a[105][105];
	cin>>n; 
	for(int i=1;i<=n;i++)
		for(int t=1;t<=n;t++)
			cin>>a[i][t];
	for(int i=1;i<=n;i++)
		for(int t=1;t<=n;t++)
			for(int j=1;j<=n;j++){
				if(a[t][i]&&a[i][j])
					a[t][j]=1;
			}
	for(int i=1;i<=n;i++){
		for(int t=1;t<=n;t++)
			cout<<a[i][t]<<" ";
		cout<<endl;
	}			
}