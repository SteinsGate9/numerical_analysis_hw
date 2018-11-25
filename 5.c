
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define MAX_SIZE 10
#define V(x,m,n) (*(*((x)+(m)-1)+(n)-1))

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);
double ** Binary(int n);

int main()
{
    int n, MAXN, m, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
    double lambda, TOL;

    scanf("%d", &n);
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            scanf("%lf", &a[i][j]);
    scanf("%lf %d", &TOL, &MAXN);
    scanf("%d", &m);
    for (i=0; i<m; i++) {
        scanf("%lf", &lambda);
        for (j=0; j<n; j++)
            scanf("%lf", &v[j]);
        switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
            case -1:
                printf("%12.8f is an eigenvalue.\n", lambda );
                break;
            case 0:
                printf("Maximum number of iterations exceeded.\n");
                break;
            case 1:
                printf("%12.8f\n", lambda );
                for (k=0; k<n; k++)
                    printf("%12.8f ", v[k]);
                printf("\n");
                break;
        }
    }

    return 0;
}
/*思路:1 A-pI
    2 计算y =(A-PI)-1X
    3 y/= ymax
    4 y = (A-pi)-1x
while(tol <= error)
    */

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN)
{
        double maxx = -100;
        double maxv = -100;
        int times = 0;
    double a2[MAX_SIZE][MAX_SIZE];
    for(int i = 0; i<= n-1; i++) for(int j=0; j <= n-1; j++) a2[i][j] = a[i][j];
    for(int i = 0; i<= n-1; i++) a2[i][i] = a[i][i]-*lambda;


   double l[MAX_SIZE][MAX_SIZE];for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i==j) l[i][j]=1;else l[i][j]=0;}
   double u[MAX_SIZE][MAX_SIZE];

   double x[MAX_SIZE];for(int i = 0;i<n;i++)x[i]=0;
   double y[MAX_SIZE];for(int i = 0;i<n;i++)y[i]=0;
   double z[MAX_SIZE];for(int i = 0;i<n;i++)z[i]=0;







    //求分解 ok
    {


  for(int i = 1;i <= n; i++)V(u,1,i) =  (double)(V(a2,1,i)/V(l,1,1)); if(V(u,1,1) == 0)return -1;
  for(int i = 2;i <= n; i++)V(l,i,1) =  (double)(V(a2,i,1)/V(u,1,1));


  for(int i = 2;i <= n-1; i++)
   {  double sum = 0;
     for(int j = i; j <= n; j++)
    {sum = 0;for(int k=1;k<=i-1;k++)sum -= (double)(V(l,i,k)*V(u,k,j));V(u,i,j)=(double)((V(a2,i,j)+sum));}
if(V(u,i,i) == 0) return -1;
     sum = 0;
     for(int j = i+1; j <= n; j++)
    {sum = 0;for(int k=1;k<=i-1;k++)sum -= (double)(V(l,j,k)*V(u,k,i));V(l,j,i)=(double)((V(a2,j,i)+sum)/V(u,i,i));}

   }
    double sum = 0;

   for(int i = 1;i<= n-1 ;i++)
   {
       sum -= (double)(V(l,n,i)*V(u,i,n)); }

       V(u,n,n)= (V(a2,n,n)+ sum);
       if(V(u,n,n) == 0)return -1;

//  for(int i = 0;i<n;i++)for(int j =0;j<n;j++){printf("%.6f ",V(l,i+1,j+1));if(j==n-1)printf("\n");}
//for(int i = 0;i<n;i++)for(int j =0;j<n;j++){printf("%.6f ",V(u,i+1,j+1));if(j==n-1)printf("\n");}
    }



















 int count = 0;

 double  error = -100;
    //迭代

    do{
    int flag = 1;
    int flag2 = 1;

       maxx = -100;
       maxv = -100;


   for(int i = 0;i <= n-1; i++){if(maxv < fabs(v[i]) ) {maxv= fabs(v[i]);if(v[i] < 0)flag2 = -1;else flag2 = 1;}}
   if (flag2 == -1) maxv = -maxv;
   for(int i = 0;i <= n-1; i++)
        {

            v[i] /= maxv;

        }










   for(int i = 1;i <= n;i++)
  {  double sum =0 ;
      if (i == 1)y[0] = (double)(v[0]/V(l,1,1));
      else
      {for(int k = 1;k<= i-1;k++)sum -= V(l,i,k)*y[k-1];y[i-1] = (double)((v[i-1]+sum));}
  }

   for(int i = n;i >= 1;i--)
  {  double sum =0;
      if (i == n)x[n-1] = (double)(y[n-1]/V(u,n,n));
      else
      {for(int k = i+1;k<= n;k++)sum -= V(u,i,k)*x[k-1];x[i-1] = (double)((y[i-1]+sum)/V(u,i,i));}
  }












        for(int i = 0;i <= n-1; i++){if(maxx < fabs(x[i]) ) {maxx= fabs(x[i]);if(x[i] < 0)flag = -1;else flag = 1;}}
            if(flag == -1)maxx = -maxx;
            if(count == 2)
            {

//             printf("maxx = %f\n",maxx);for(int i = 0;i <= n-1; i++){printf("x = %f\n",x[i]);printf("v = %f\n",v[i]);}
            }

                error= -1;
        for(int i = 0;i <= n-1; i++)
            {
                if(error < fabs(v[i]-x[i]/maxx))
                {
                    error = fabs(v[i]-x[i]/maxx);
//                    printf("this error = %f\n",error);
                }
            }


//        printf("this error = %f\n",error);
        for(int i = 0; i<= n-1; i++){ v[i] = x[i]/maxx; }

//        printf("\n");
        times ++; count ++;
        if(times > MAXN)return 0;
        if(error < TOL ){*lambda = 1/maxx + *lambda; return 1;}
    } while(error >= TOL);


















        return 1;



}

double ** Binary(int n)
{
    double **a = (double **)malloc(sizeof(double *)*n);for(int i = 0;i<n;i++)*(a+i) = (double*)malloc(sizeof(double)*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){*(*(a+i)+j)=0;}
    return a;
}
