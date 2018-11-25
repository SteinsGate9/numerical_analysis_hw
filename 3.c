#include <stdio.h>
#include <stdlib.h>
#define Max_size 10000 /* max number of dishes */
#define V(x,m,n) (*(*((x)+(m)-1)+(n)-1))
void Price( int n, double p[] );
double ** Binary(int n);
int main()
{
    int n, i;
    double p[Max_size];



    scanf("%d", &n);
    for (i=0; i<n; i++)
        scanf("%lf", &p[i]);
    Price(n, p);
    for (i=0; i<n; i++)
        printf("%.2f ", p[i]);
    printf("\n");

    return 0;
}
void Price( int n, double p[] )
{
  double **a = Binary(n);  for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i==j)*(*(a+i)+j)=2;else if(i-j==1||i-j==-1)*(*(a+i)+j)=0.5;else *(*(a+i)+j)=0;*(*a+n-1) = 0.5;**(a+n-1)= 0.5;}
  double **l = Binary(n);  for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i==j)*(*(l+i)+j)=1;else *(*(l+i)+j)=0;}
  double **u = Binary(n);

  double *x = (double*)malloc(sizeof(double)*n);for(int i = 0;i<n;i++)x[i]=0;
  double *y = (double*)malloc(sizeof(double)*n);for(int i = 0;i<n;i++)y[i]=0;





  for(int i = 1;i <= n; i++)V(u,1,i) =  (double)(a[i-1][0]/V(l,1,1));
  for(int i = 2;i <= n; i++)V(l,i,1) =  (double)(V(a,i,1)/V(u,1,1));


  for(int i = 2;i <= n-1; i++)
   {  double sum = 0;
     for(int j = i; j <= n; j++)
    {sum = 0;for(int k=1;k<=i-1;k++)sum -= (double)(V(l,i,k)*V(u,k,j));V(u,i,j)=(double)((V(a,i,j)+sum));}

     sum = 0;
     for(int j = i+1; j <= n; j++)
    {sum = 0;for(int k=1;k<=i-1;k++)sum -= (double)(V(l,j,k)*V(u,k,i));V(l,j,i)=(double)((V(a,j,i)+sum)/V(u,i,i));}

  } double sum = 0;


  for(int i = 0;i<n;i++)for(int j =0;j<n;j++){printf("%.4f ",V(l,i+1,j+1));if(j==n-1)printf("\n");}

   for(int i = 1;i<= n-1 ;i++)
   {


       sum -= (double)(V(l,n,i)*V(u,i,n)); }
       V(u,n,n)= (V(a,n,n)+ sum);














  for(int i = 1;i <= n;i++)
  {  double sum =0 ;
      if (i == 1)y[0] = (double)(p[0]/V(l,1,1));
      else
      {for(int k = 1;k<= i-1;k++)sum -= V(l,i,k)*y[k-1];y[i-1] = (double)((p[i-1]+sum));}
  }

  for(int i = n;i >= 1;i--)
  {  double sum =0;
      if (i == n)x[n-1] = (double)(y[n-1]/V(u,n,n));
      else
      {for(int k = i+1;k<= n;k++)sum -= V(u,i,k)*x[k-1];x[i-1] = (double)((y[i-1]+sum)/V(u,i,i));}
  }


  for(int i = 0;i < n;i++) p[i] = x[i];

  return;
}







double ** Binary(int n)
{
    double **a = (double **)malloc(sizeof(double *)*n);for(int i = 0;i<n;i++)*(a+i) = (double*)malloc(sizeof(double)*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){*(*(a+i)+j)=0;}
    return a;
}
