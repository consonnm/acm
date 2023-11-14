#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define endl '\n'
struct node{
    int l,r,v;
}t[1000005*24];
int flag,root[1000005*24],n,m,a[2000005],b,c,d,va;
int build(int l,int r){
    int id=++flag;
    if(l==r){
        t[id].v=a[l];
        return id;
    }
    int mid=(l+r)>>1;
    t[id].l=build(l,mid);
    t[id].r=build(mid+1,r);
    return id;
}
int clone(int x){
	flag++;
    t[flag]=t[x];
    return flag;
}
int update(int p,int l,int r,int ind,int v){
    p=clone(p);
    if(l==r){
        t[p].v=v;
        return p;
    }
    int mid=(l+r)>>1;
    if(ind<=mid)
        t[p].l=update(t[p].l,l,mid,ind,v);
    else
        t[p].r=update(t[p].r,mid+1,r,ind,v);
    return p;
}
int query(int p,int l,int r,int ind){
    if(l==r)
        return t[p].v;
    int mid=(l+r)>>1;
    if(ind<=mid)
        return query(t[p].l,l,mid,ind);
    else
        return query(t[p].r,mid+1,r,ind);
}
int main(){
    cin>>n>>m;
    ios::sync_with_stdio(false);
    for(int i=1;i<=n;i++)
        cin>>a[i];
    build(1,n);
    int flag1=0;
    root[0]=1;
    for(int i=1;i<=m;i++){
        cin>>va>>b;
        if(b==1){
            cin>>c>>d;
            flag1++;
            root[flag1]=update(root[va],1,n,c,d);
        }
        else{
            cin>>c;
            root[++flag1]=root[va];
            cout<<query(root[va],1,n,c)<<endl;
        }
    }
}
