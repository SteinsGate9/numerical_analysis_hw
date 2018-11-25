#include <stdio.h>
#include <math.h>

#define MAX_SIZE 15
#define bound pow(2, 127)
#define ZERO 1e-9 /* X is considered to be 0 if |X|<ZERO */

int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );

int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );

int main()
{
    int n, MAXN, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE], x[MAX_SIZE];
    double TOL;

    scanf("%d", &n);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++)
            scanf("%lf", &a[i][j]);
        scanf("%lf", &b[i]);
    }
    scanf("%lf %d", &TOL, &MAXN);

    printf("Result of Jacobi method:\n");
    for ( i=0; i<n; i++ )
        x[i] = 0.0;
    k = Jacobi( n, a, b, x, TOL, MAXN );
    switch ( k ) {
        case -2:
            printf("No convergence.\n");
            break;
        case -1:
            printf("Matrix has a zero column.  No unique solution exists.\n");
            break;
        case 0:
            printf("Maximum number of iterations exceeded.\n");
            break;
        default:
            printf("no_iteration = %d\n", k);
            for ( j=0; j<n; j++ )
                printf("%.8f\n", x[j]);
            break;
    }
    printf("Result of Gauss-Seidel method:\n");
    for ( i=0; i<n; i++ )
        x[i] = 0.0;
    k = Gauss_Seidel( n, a, b, x, TOL, MAXN );
    switch ( k ) {
        case -2:
            printf("No convergence.\n");
            break;
        case -1:
            printf("Matrix has a zero column.  No unique solution exists.\n");
            break;
        case 0:
            printf("Maximum number of iterations exceeded.\n");
            break;
        default:
            printf("no_iteration = %d\n", k);
            for ( j=0; j<n; j++ )
                printf("%.8f\n", x[j]);
            break;
    }






    return 0;
}



int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{


 int interation = 0;


        for(int i = 0; i <= n-1; i++)

        {   double MAX = -190000l;int markrow = -100;
            for(int j = i; j <= n-1; j++)if(fabs(a[j][i]) > MAX){MAX = fabs(a[j][i]);markrow = j;}
            if(MAX != 0&&markrow != i){double temp;for(int k =0; k <= n-1; k++){temp = a[i][k];a[i][k] =  a[markrow][k]; a[markrow][k] = temp;}temp = b[i];b[i] = b[markrow];b[markrow] = temp;}
            else if(MAX == 0){
             for(int j = i-1; j >= 0; j--)if(fabs(a[j][i]) > MAX){MAX = fabs(a[j][i]);markrow = j;}
             if(MAX != 0) {for(int k =0; k <= n-1; k++){a[i][k] += a[markrow][k];}b[i]+=b[markrow];}
             else if(MAX == 0)return -1;
            }
        }


        double x0[MAX_SIZE] = {0};
        while(interation <= MAXN )
  {
     double realmax = -100l;
     double max = -100;
     for(int i = 0;i<= n-1;i++)x0[i] = x[i];
    for(int i = 0;i <= n-1; i++)
    {
        double sum = 0;


            for(int j = 0; j <= i-1; j++)
            sum += x0[j]*a[i][j];
            for(int j = i+1; j <= n-1; j++)
            sum += x0[j]*a[i][j];

            x[i] = (b[i]-sum)/a[i][i];
//            printf("x[%d] = %f i = %d\n",i,x[i],interation);

            if (fabs(x[i]-x0[i])>max) max = fabs(x[i]-x0[i]);
            if (fabs(x[i])>realmax)realmax = fabs(x[i]);
    }

       interation ++;
       if(interation > MAXN) {for(int i =0;i<=n-1;i++)printf("%d\n",x[i] );return 0;}
       if(realmax >= bound) return -2;
       if(max <= TOL)return interation;
  }





}

int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{
       int interation = 0;

        int markrow = -100;
        for(int i = 0; i <= n-1; i++)
        {   double MAX = -190000l;
            for(int j = i; j <= n-1; j++)if(fabs(a[i][j]) > MAX){MAX = fabs(a[i][j]);markrow = j;}
            if(MAX != 0&&markrow != i){double temp;for(int k =0; k <= n-1; k++){temp = a[i][k];a[i][k] =  a[markrow][k]; a[markrow][k] = temp;}temp = b[i];b[i] = b[markrow];b[markrow] = temp;}
            else if(MAX == 0){
             for(int j = i-1; j >= 0; j--)if(fabs(a[i][j]) > MAX){MAX = fabs(a[i][j]);markrow = j;}
             if(MAX != 0) {for(int k =0; k <= n-1; k++){a[i][k] += a[markrow][k];}b[i]+=b[markrow];}
             else if(MAX == 0)return -1;
            }
        }















        while(interation <= MAXN )
  {
     double realmax = -100l;
     double max = -100;
    for(int i = 0;i <= n-1; i++)
    {
        double sum = 0;
        double temp=x[i];

            for(int j = 0; j <= i-1; j++)
            sum += x[j]*a[i][j];
            for(int j = i+1; j <= n-1; j++)
            sum += x[j]*a[i][j];

            x[i] = (b[i]-sum)/a[i][i];
//            printf("x[%d] = %f i = %d\n",i,x[i],interation);

            if (fabs(x[i]-temp)>max) max = fabs(x[i]-temp);
            if (fabs(x[i])>realmax)realmax = fabs(x[i]);
    }

       interation ++;
       if(interation > MAXN) return 0;
       if(realmax >= bound) return -2;
       if(max <= TOL)return interation;
  }



}




