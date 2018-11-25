#include <stdio.h>
#include <math.h>

void Series_Sum( double sum[] );
double SigDig( double x);

int main()
{
    int i;
    double x, sum[3001];

    Series_Sum( sum );

    x = 0.0;
    for (i=0; i <= 1609; i++)
        printf("%6.2f %16.12f\n", x + (double)i * 0.10, sum[i]);

    return 0;
}




void Series_Sum( double sum[] ){
	int i,k;
	double a;
	double temp=0;
	for( i=0; i<10; i++ ){
		a=0.1*i;
		temp=0;
		for( k=1; k<=pow((1-a)*(2-a)*pow(10,10)/3,1.0/2.1); k++ ){
			temp+=(2.0-a)/(k*(k+1.0)*(k+2.0)*(k+a));
		}
		sum[i]=(1.0-a)*(temp+0.25)+1.0;
	}
	for( i=10; i<3001; i++ ){
		a=0.1*i;
		sum[i]=(a*(a-1)*sum[i-10]+1)/(a*a);
	}
}





