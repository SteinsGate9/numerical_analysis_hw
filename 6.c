#include <stdio.h>

#include<malloc.h>
#include<math.h>
#define V(x,m,n) (*(*((x)+(m)-1)+(n)-1))

#define MAX_N 100
void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]);
//int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );
double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] );
double ** Binary(int n);
int I;
int main()
{
    int n, Type, m, i;
    double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
    double s0, sn, Fmax, t0, tm, h, t;

    scanf("%d", &n);
    for (i=0; i<=n; i++)
        scanf("%lf", &x[i]);

    for (i=0; i<=n; i++)
        {scanf("%lf", &f[i]);
//         printf("x[0] = %f\n",x[0]);
         }
    scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);

    Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
    for (i=1; i<=n; i++)
        printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

    scanf("%lf %lf %d", &t0, &tm, &m);
    h = (tm-t0)/(double)m;
    for (i=0; i<=m; i++) {
        t = t0+h*(double)i;
        printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
    }

    return 0;
}

 void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[])
{
//    double juzhen[MAX_N][MAX_N] ;
    double **juzhen = Binary(4*n);

    for(int i = 0; i<4*n; i++)for(int j = 0;j <4*n;j++) juzhen [i][j]=0;
   double *jie = (double* )malloc(sizeof(double)* 4*n);
   for(int i= 0; i<4*n ; i++)jie[i] = 0;

//for(int i = 0;i <= n;i++)
//printf("x[%d] = %f \n",i,x[i]);

// printf("?");
    for(int i = 0; i<= 2*n-1; i++)
    {
        if((i%2) == 0)
        {

        V(juzhen,i+1,i/2+1) = 1;
        jie[i] = f[i/2];
        }
        else
        {
        V(juzhen,i+1,i/2+1) = 1;
        V(juzhen,i+1,i/2+n+1) = x[i/2+1]-x[i/2];
        V(juzhen,i+1,i/2+n*2+1) = pow(x[i/2+1]-x[i/2],2);
        V(juzhen,i+1,i/2+n*3+1) = pow(x[i/2+1]-x[i/2],3);
        jie[i] = f[i/2+1];
        }
    }
    for(int i = 0;i <= n-2; i++)
    {
        V(juzhen,2*n+2*i+1,i+n+1) = 1;
        V(juzhen,2*n+2*i+1,i+n*2+1) = 2*(x[i+1]-x[i]);
        V(juzhen,2*n+2*i+1,i+n*3+1) = 3*pow(x[i+1]-x[i],2);
        V(juzhen,2*n+2*i+1,i+n+1+1) = -1;
//        juzhen[2*n+2*i][i+n*2+1] = -2*x[i+1];
//        juzhen[2*n+2*i][i+n*3+1] = -3*x[i+1]*x[i+1];
        V(juzhen,2*n+2*i+1+1,i+n*2+1) = 2;
        V(juzhen,2*n+2*i+1+1,i+n*3+1) = 6*(x[i+1]-x[i]);
        V(juzhen,2*n+2*i+1+1,i+n*2+1+1) = -2;
//        juzhen[2*n+2*i+1][i+n*3+1] = -6*x[i+1];
    }
    if(Type == 1)
    {
        V(juzhen,4*n-2+1,n+1) = 1;
//        juzhen[4*n-2][n*2] = 2*x[0];
//        juzhen[4*n-2][n*3] = 3*x[0]*x[0];
        jie[4*n-2] = s0;
        V(juzhen,4*n-1+1,n-1+n+1) = 1;
//        juzhen[4*n-1][n-1+n*2] = 2*x[n];
//        juzhen[4*n-1][n-1+n*3] = 3*x[n]*x[n];
        jie[4*n-1] = sn;
    }
    else if(Type ==2)
    {
        jie[4*n-2] = s0;
        jie[4*n-1] = sn;
        V(juzhen,4*n-2+1,n*2+1) = 2;
//        juzhen[4*n-2][n*3] = 6*x[0];
        V(juzhen,4*n-1+1,n-1+n*2+1) = 2;
    V(juzhen,4*n-1+1,n-1+n*3+1) = 6*(x[n]-x[n-1]);
//        juzhen[4*n-3][n*3] = 6*x[n];

    }


//    for(int j =0;j<4*n;j++){printf("%.6f ",jie[j]);if(j==4*n-1)printf("\n");}
//    printf("\n");

    n = 4*n;
  // double **a2 = Binary(n);  for(int i=0;i<n;i++)for(int j=0;j<n;j++){a2[i][j]= juzhen[i][j];}
//  for(int i = 0;i<n;i++)for(int j =0;j<n;j++){printf("%.6f ",V(a2,i+1,j+1));if(j==n-1)printf("\n");}
//  printf("\n");




// partial scaled
{

    for(int i = 1; i<= n; i++)
    {
      double max = -100;
      int maxrow = i;
      for(int j = i; j<= n; j++)
            if(fabs(V(juzhen,j,i)) > max)
            {max = fabs(V(juzhen,j,i));
             maxrow = j;}
      if(max == 0) return;
      if(maxrow != i)
      {
//          printf("maxrow[%d] = %d\n",i,maxrow);
          double *t = (double *)malloc(sizeof(double)*n);
          double temp ;
          for(int i2 = 1; i2<= n; i2++)t[i2-1] = V(juzhen,i,i2);
          for(int i2 = 1; i2<= n; i2++)V(juzhen,i,i2) = V(juzhen,maxrow,i2);
          for(int i2 = 1; i2<= n; i2++)V(juzhen,maxrow,i2) = t[i2-1];
          temp = jie[i-1];
          jie[i-1] = jie[maxrow-1];
          jie[maxrow-1] = temp;
      }

      for(int i3 = i+1; i3 <= n; i3++)
      {
          double m = V(juzhen,i3,i)/V(juzhen,i,i);
          for(int i4 = i; i4<= n; i4++)
            V(juzhen,i3,i4) = V(juzhen,i3,i4) - m*V(juzhen,i,i4);
          jie[i3-1] = jie[i3-1] - m*jie[i-1];
      }
//      for(int i = 0;i<n;i++)for(int j =0;j<n;j++){printf("%.6f ",V(juzhen,i+1,j+1));if(j==n-1)printf("\n");}
//  printf("\n");

    }

    if(V(juzhen,n,n) == 0)return;

    jie[n-1] = jie[n-1]/V(juzhen,n,n);
    for(int i = n-1; i>= 1; i--)
    {
        double sum = 0;
        for(int c = n;c > i; c--)
        {
          sum += V(juzhen,i,c)*jie[c-1];
        }
        jie[i-1] = (jie[i-1]-sum)/V(juzhen,i,i);
    }
    free(juzhen);

    for(int i = 0; i<= n-1; i++){if(jie[i] == -0)jie[i] = 0;}
    n = n/4;
   for(int i = 0; i<= n-1; i++) a[i+1] = jie[i];
   for(int i = n,j=1; i<= 2*n-1; j++,i++){b[j] = jie[i];}
   for(int i = 2*n,j=1; i<= 3*n-1;j++, i++)c[j] = jie[i];
   for(int i = 3*n,j=1; i<= 4*n-1;j++, i++)d[j] = jie[i];

}








}
double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] )
{
  if(t < x[0] || t > x[n] ){
//        printf("t = %f x[0] = %f x[n] = %f\n",t,x[0],x[n]);
  return Fmax;}
  for(int i = 0;i< n; i++)
    if(x[i] <= t && x[i+1] >= t)
   {
      I = i;
//      printf("I = %d\n",I);
      break;

   }
    return ( a[I+1]+b[I+1]*(t-x[I])+c[I+1]*pow((t-x[I]),2)+d[I+1]*pow((t-x[I]),3));


}
double ** Binary(int n)
{
    double **a = (double **)malloc(sizeof(double *)*n);for(int i = 0;i<n;i++)*(a+i) = (double*)malloc(sizeof(double)*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){*(*(a+i)+j)=0;}
    return a;
}
