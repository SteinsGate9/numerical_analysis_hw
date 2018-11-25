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
//  double **a = Binary(n);  for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i==j)*(*(a+i)+j)=2;else if(i-j==1||i-j==-1)*(*(a+i)+j)=0.5;else *(*(a+i)+j)=0;*(*a+n-1);**(a+n-1);}
  double *l1 = (double*)malloc(sizeof(double)*(n-1));for(int i = 0;i<n-1;i++)l1[i]=0;
  double *l2 = (double*)malloc(sizeof(double)*n);for(int i = 0;i<n;i++)l2[i]=0;
  double *u1 = (double*)malloc(sizeof(double)*(n-1));for(int i = 0;i<n-1;i++)u1[i]=0;
  double *u2 = (double*)malloc(sizeof(double)*(n-1));for(int i = 0;i<n-1;i++)u2[i]=0;
  double *x = (double*)malloc(sizeof(double)*n);for(int i = 0;i<n;i++)x[i]=0;
  double *y = (double*)malloc(sizeof(double)*n);for(int i = 0;i<n;i++)y[i]=0;


l1[0] = 0.5;
u2[0] = 0.25;
l2[0] =  (double)(2);
u1[0] =  (double)(0.5/l2[0]);

  for(int i = 2;i <= n-1; i++)
   {
     l1[i-2] = (double)(0.5);
     l2[i-1] = 2-l1[i-2]*u1[i-2];
     u1[i-1] = (double)(0.5/l2[i-1]);
  }


for( int i = 2 ;i <= n-2;i++)
{l1[i-1]  = -(u1[i-2])*l1[i-2];
 u2[i-1] = -(0.5)*u2[i-2]/l2[i-1];
}
l1[n-2] = 0.5-(u1[n-3]*l1[n-3]);
u2[n-2] = (0.5-(0.5)*u2[n-3])/l2[n-2];
double sum =0;
 for(int i = 1;i<= n-1 ;i++)
 {
 sum -= l1[i-1]*u2[i-1];
 }
 l2[n-1] = 2.0+sum;




for(int i = 1;i <= n;i++)
  { sum = 0;
      if (i == 1)y[0] = (double)(p[0]/l2[0]);
      else if (i == n){for(int i=0;i<=n-2;i++)sum -= l1[i]*y[i];y[n-1] = (p[n-1]+sum)/l2[n-1];}
      else { sum -= 0.5*y[i-2];y[i-1] = (double)((p[i-1]+sum)/l2[i-1]);}
  }

  for(int i = n;i >= 1;i--)
  { sum = 0;
      if (i == n)x[n-1] = (double)(y[n-1]);
      else if  (i == n-1) x[n-2] = (double)(y[n-2]-x[n-1]*u2[n-2]);
      else { sum =- (u1[i-1]*x[i]+u2[i-1]*x[n-1]);x[i-1] = (double)((y[i-1]+sum));
  }
  }









  for(int i = 0;i < n;i++) p[i] = x[i];
  return;
}







double ** Binary(int n)
{
    double **a = (double **)malloc(sizeof(double *)*n);for(int i = 0;i<n;i++)*(a+i) = (double*)malloc(sizeof(double)*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){*(*(a+i)+j)=0;}
    return a;
}

