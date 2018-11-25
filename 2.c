#include <stdio.h>
#include <math.h>

#define ZERO 1e-13 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11






double Polynomial_Root(int n, double c[], double a, double b, double EPS);
double Polynomial_Value(int n, double c[],double x);
double Polynomial_Daoshu(int n, double c[],double x);
double Polynomial_Daoshu2(int n, double c[], double x);
double Findroot(int n,double c[],double a,double b,double xo,double *p,double *err,double EPS);
double shenme(int n, double c[],double x, double EPS);
double jdz(double x);
int sgn(double x);
int fuck(double eps,double a,double b);
int shouldtime;
int time2;





double jdz(double x)
{
  if(x >= 0)
    return x;

 else
    return -x;
}


int fuck(double eps,double a,double b)
{

       double interval = b-a;
       return (int)(-1)*log10(eps/interval)/log10(2);

}


int sgn(double x)
{
  if(x >= ZERO)
    return 1;
 else if(jdz(x) < ZERO  )
    return 0;
 else
    return -1;
}


double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{


   if(a == b)
{
    // printf("really?");
     return a;
}


   else if(a != b)
{
     if(a > b)
   {
       double temp;
       temp = b;
       b = a;
       a = temp;
   }



    double x ;//chuzhi
    double shit = 0;
    double x0;
    if(sgn(Polynomial_Value(n,c,a))*sgn(Polynomial_Value(n,c,b)) == 0)
    {
      //  printf("really?\n");;
      //  printf("a = %lf.b = %lf",Polynomial_Value(n,c,a),Polynomial_Value(n,c,b));
      return (Polynomial_Value(n,c,a) == 0)?a:b;
    }



     double result = 0;
     double min = 100;
     double temp;
     double error;
     double flag = 0;






          if(sgn(Polynomial_Value(n,c,a))*sgn(Polynomial_Value(n,c,b)) == -1)
          {
          while(sgn(Polynomial_Value(n,c,a))*sgn(Polynomial_Value(n,c,b)) == -1 && (b-a) >= ZERO)
          {
          x0 = a + (b-a)/2;
          if(sgn(Polynomial_Value(n,c,x0))*sgn(Polynomial_Value(n,c,b)) == -1)a = x0;
          else b = x0;
          }
          if(jdz(Polynomial_Value(n,c,x0)) <= min)
          return x0;
          }


            else
            for( int step = 0; step <= 30; step++)
         {
          {
              x = a + step*1.0/30*(b-a);
              // printf("what??\n");
          }

            temp = Findroot(n,c,a,b,x,&flag,&error,EPS);
            if( flag == 1 && jdz(error) < min)
            {
                min = jdz(error);
                result = temp;
            }
         }
           // printf("step = %d,error = %f,temp = %f,min = %lf\n",step,error,temp,min);
           // printf("flag = %lf\n",Findroot(n,c,a,b,x,&temp,&error,EPS));
         return result;
  }
}






double Findroot(int n,double c[],double a,double b,double x0,double *flag,double *err,double EPS)
{


     shouldtime = 1002;
     int times = 0;
     double x = x0;
     double p0;

      do{


        double value = Polynomial_Value(n,c,x);
        double daoshu = Polynomial_Daoshu(n,c,x);
        double daoshu2 = Polynomial_Daoshu2(n,c,x);


           //运算
        if(sgn(jdz(value)) == 0 )
        {
         // printf("条件1\n");
          *flag = 1;
          break;
        }
        if(sgn(jdz(daoshu*daoshu-value*daoshu2)) == 0 )
        {
          *flag = 0;
          break;
        }
        p0  = x - value*daoshu/(daoshu*daoshu-value*daoshu2);
        if(sgn(jdz(Polynomial_Value(n,c,p0))) == 0 )
        {
          //printf("条件2\n");
          *flag = 1;
          break;
        }
        if(jdz(p0-x) < EPS)
        {
          //printf("条件3\n");
          *flag = 1;
          break;
        }
        if(p0 < a || p0 > b)
        {
          *flag = 0;
          break;
        }
        x = p0;
        times ++;
      }while(times < shouldtime);

       *err = Polynomial_Value(n,c,p0);
       return p0 ;
}

double Polynomial_Value(int n, double c[],double x)
{

   double value = 0.0 ;
    for(int i = n; i >= 0; i--)
    {

     value = value * x + c[i];

    }

      return value;




}

double Polynomial_Daoshu(int n, double c[],double x)
{

    double value = 0 ;
    for(int i = n; i >= 1; i--)
    {
        value = value * x + i * c[i];
    }

      return value;


}

double Polynomial_Daoshu2(int n, double c[], double x)
    {
        double value = 0;
        for(int i = n; i >= 2; i--)
        {
            value = value * x + i * (i - 1) * c[i];
        }
        return value;
    }
