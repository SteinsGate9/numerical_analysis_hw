
#include <stdio.h>
#include <math.h>

#define Max_m 200
#define Max_n 5

double f1(double x)
{
    return sin(x);
}

double f2(double x)
{
    return exp(x);
}

int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps );

void print_results( int n, double c[], double eps)
{
    int i;

    printf("%d\n", n);
    for (i=0; i<=n; i++)
        printf("%12.4e ", c[i]);
    printf("\n");
    printf("error = %9.2e\n", eps);
    printf("\n");
}

int main()
{
    int m, i, n;
    double x[Max_m], w[Max_m], c[Max_n+1], eps;

    m = 90;
    for (i=0; i<m; i++) {
        x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f1, m, x, w, c, &eps);
    print_results(n, c, eps);

    m = 200;
    for (i=0; i<m; i++) {
        x[i] = 0.01*(double)i;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f2, m, x, w, c, &eps);
    print_results(n, c, eps);

    return 0;
}


















int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps )
{
   double  fai0[Max_n+1]={0},fai1[Max_n+1]={0},fai2[Max_n+1]={0};
   double  B=0,a=0,b=0,error=0,temp=0,temp2=0;
   int times=1;
   for(int i = 0; i<= Max_n; i++) {c[i] = 0;fai0[i]=0;fai1[i] =0;fai2[i] = 0;}
   fai0[0] = 1;
   for(int i = 0; i< m ; i++)
   {
       a += (*f)(x[i]);
       temp2 += ((*f)(x[i]))*((*f)(x[i]));
       B += x[i];
       b += 1;
   }

//    printf("%f %f",a,temp2);


      c[0] += a/b;
      error = temp2 - a*a/b;
      B /= b;

    fai1[1] = 1;
    fai1[0] = -B;

   a = 0;
   b = 0;
//    printf("error[0] = %f\n",error);
   for(int i = 0;i < m; i++)
  {
        double q =0 ;

   for(int j = 0; j<= Max_n; j++)
   {
     if(fai1[j] != 0)
     {
        a += fai1[j]*pow(x[i],j)*(*f)(x[i]);
        q += fai1[j]*pow(x[i],j);
     }

   }
        b += q*q;
  }
for(int i = 0; i <= Max_n; i++)
 {
//     printf("c[%d] = %f\n",i,c[i]);
}

      c[1] += a/b*fai1[1];
      c[0] += a/b*fai1[0];

     error -= a*a/b;
// printf("error[1] = %f\n",error);
 for(int i = 0; i <= Max_n; i++)
 {
//  printf("fai1[%d] = %.9f\n",i,fai1[i]);
 }
//for(int i = 0; i <= Max_n; i++)
// {printf("c[%d] = %f\n",i,c[i]);
//}
//










while(fabs(error) >= *eps && times < Max_n)
     {
         double bk = 0,ck = 0,bkd = 0,ckd = 0,ak = 0,akf = 0;

//  bk
//for(int j=  0 ; j<= Max_n; j+C+)
//    if(fai1[j] != 0)
//printf("fai1[%d] = %f ",j,fai1[j]); printf("\n");
    for(int i = 0;i < m; i++)
  {
        double q1=0;

   for(int j = 0; j<= Max_n; j++)
   {
     if(fai1[j] != 0)
     {
        q1 += fai1[j]*pow(x[i],j);

     }


   }
        bk  += q1*q1*x[i];
        bkd += q1*q1;
  }

//  printf("bk = %f bkd = %f\n",bk,bkd);


// ck

    for(int i = 0;i < m; i++)
  {
        double qd =0 ,q1=0 ,q2 = 0;

   for(int j = 0; j<= Max_n; j++)
   {
     if(fai1[j] != 0)
     {
        q1 += fai1[j]*pow(x[i],j);
     }
     if(fai0[j] != 0)
        q2 += fai0[j]*pow(x[i],j);

   }
        ck  += q1*q2*x[i];
        ckd += q2*q2;
  }
//  printf("ck = %.14f ckd = %f\n",ck,ckd);

// ¼õ·¨
   for(int i = 0; i<= Max_n-1; i++)
  {
      if(fai1[i] != 0)
      {fai2[i+1] += fai1[i];
       fai2[i] -= (bk/bkd)*fai1[i];}
      if(fai0[i] != 0)
      {
       fai2[i] -= (ck/ckd)*fai0[i];
      }

  }
//  for(int i = 0; i<= Max_n; i++)
//    printf("fai2[%d] = %f\n",i,fai2[i]);
//   for(int i = 0; i<= Max_n; i++)
//    printf("fai1[%d] = %f\n",i,fai1[i]);
// ak
    for(int i = 0;i < m; i++)
  {
        double q =0 ;

   for(int j = 0; j<= Max_n; j++)
   {
     if(fai2[j] != 0)
     {

        q += fai2[j]*pow(x[i],j);
     }

   }
        akf += q*(*f)(x[i]);
        ak += q*q;
  }
//  printf("ak = %f akf = %f\n",ak,akf);
// ¼Ó·¨
   for(int i = 0; i <= Max_n; i++)
   {
       if(fai2[i]!= 0)
      c[i] += akf/ak*fai2[i];
   }
// err
   error -= akf*akf/ak;

  for(int i = 0; i<=Max_n ; i++)
   {
    fai0[i] = fai1[i];
    fai1[i] = fai2[i];
    fai2[i] = 0;
   }
// printf("error[%d] = %f\n",times+1,error);
         times ++;
}


   *eps = error;
   return times;
}

/* Your function will be put here */
