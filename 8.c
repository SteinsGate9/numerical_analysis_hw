#include <stdio.h>
#include <math.h>

double f0( double x, double l, double t )
{
    return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}

double Integral(double a, double b, double (*f)(double x, double y, double z), double eps, double l, double t);

int main()
{
    double a=0.0, b, eps=0.005, l, t;

    scanf("%lf %lf %lf", &l, &b, &t);
    printf("%.2f\n", Integral(a, b, f0, eps, l, t));

    return 0;
}
double Integral(double a, double b, double (*f)(double x, double y, double z), double eps, double l, double t)
{
    double leng = b-a;
    double interval ;
    double h;
    double sum=0;
    double x1[100000];
    int shit = 20;
    double r[100][100] ;

   for(int i =0;i<=99;i++)
    for(int j=0;j<=99;j++)
    {r[i][j] = 0;
//      printf("n = %d\n",shit);
    }

    h = leng;
//    interval = pow(12/leng*eps*1/((l*t*t*t-0.25*l*t*t)/pow(2,1.5)),0.25);
//    double fuck = leng/interval;
//    shit = (int)fuck;
//    printf("n = %d\n",shit);

//    for(int i = 0; i<=shit; i++)
//        x1[i] = (*f)(a+h*i,l,t);
//    for(int i = 1; i<=shit-1; i++)
//        sum+= 2*x1[i];
//        sum += (x1[0]+x1[shit]);

    r[1][1]=(h)/2*((*f)(a,l,t)+(*f)(b,l,t));
    int flag = 0;
//    printf("%f\n",50/2*((*f)(a,l,t)+(*f)(b,l,t)+2*(*f)(a+50,l,t)));
//    printf("%f\n",r[1][1]);
    int i=2;
    for(i = 2; i<=shit&&flag !=1; i++)
    {
       r[2][1]= 0.5*(r[1][1]);


     for(int k =1; k<= pow(2.0,i*1.0-2.0);k++)  r[2][1] += (h*1.0/2.0)*(*f)(a+(k*1.0-0.5)*h,l,t);
     if(fabs(r[2][1]-r[1][1])<eps/2)
        flag =1 ;
//      printf("[2][1] %f \n",r[2][1]);
     double before = 0;
     for(int j =2; j<= i; j++)
     {
     r[2][j] = r[2][j-1]+(r[2][j-1]-r[1][j-1])/(pow(4,j-1)-1);

//     printf("[2][%d] %f \n",j,r[2][j]);
     }
     if(flag ==1 )return r[2][i]/100;
     h/=2.0;
     for(int j =1;j<=i;j++)r[1][j] = r[2][j];

    }

    return r[2][shit]/100;
}
