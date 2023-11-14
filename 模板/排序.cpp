#include <stdio.h>
void sqort(int a[]);
void guibing(int a[],int n);
void choose(int a[],int n);
int main(void){
	int a[14] = {1,4,5,3,5,6,2,5,7,8,23,34,45,76};
	int i=0;
	choose(a,14);
	for(i=0;i<14;i++)
		printf("%d ",a[i]);
	return 0;
}
void fast(int a[],int begin,int end){
	int i=begin,j=end,t;
	int mid=a[(i+j)/2];
	if(begin<end){
	do{
		while(a[i]<mid)
			i++;
		while(a[j]>mid)
			j--;
		if(i<=j){
			t=a[i];
			a[i]=a[j];
			a[j]=t;
			i++;
			j--;
		}
	}while(i<=j);
	fast(a,begin,j);
	fast(a,i,end);
	}
	
}
void guibing(int a[],int n){
	int i,t,temp;
	for(i=1;i<n;i++){
		if(a[i-1]>a[i]){
			temp=a[i];
			for(t=i-1;t>=0&&a[t]>temp;t--)
				a[t+1]=a[t];
			a[t+1]=temp;
		}
	}
}
void choose(int a[],int n)
{
	int i,t,temp,min;
	for(i=0;i<n;i++){
		min=i;
		for(t=i;t<n;t++)
			if(a[min]>a[t])
				min=t;
		temp=a[min];
		a[min]=a[i];
		a[i]=temp;
	}
}
