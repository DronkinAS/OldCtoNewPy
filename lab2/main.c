#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 4000


int main (int argc, char** argv)
{
double step, *t, *akf, temp_sum=0, srednee, koeff, dz, mod_akf;  
double *x, *ogib, *time_ogib, temp, temp1, temp2, temp3, temp4,sum1=0,sum2=0;
double Integral;
int i=0,k=0, N=0, m, ind_ogib=0;
FILE *in_date;
FILE *out_date;
FILE *out_ogib_date;
FILE *in_ogib_date;

  in_date = fopen ("realization.dat","r");

  while (!feof(in_date))
  {
     fscanf (in_date, "%lf %lf", &temp, &temp2);
     N=N+1 ;
  }
  fclose (in_date);
  N=N-1;

  x = (double *)calloc(N, sizeof(double));
  t = (double *)calloc(2, sizeof(double));
  akf = (double *)calloc(M+1, sizeof(double));

  in_date = fopen ("realization.dat","r");

  for (i=0; i<2; i++)
  {
   fscanf (in_date, "%lf %lf", &t[i], &temp);
  }
  step=t[1]-t[0];

  fclose (in_date);
  
  in_date = fopen ("realization.dat","r");
  while (!feof(in_date))
  {
     fscanf (in_date, "%lf %lf", &temp, &x[k]);
     temp_sum=temp_sum+x[k];
     k=k+1 ;
  }
  fclose (in_date);
  
  srednee=temp_sum/N;
  
  out_date = fopen ("akf.dat","w");
  for (m=0; m<=M; m++)
  { 
    temp_sum=0;
    k=N-m;
    for(i=0; i<=k; i++)
        {
          koeff = (x[i]-srednee)*(x[i+m]-srednee);
          temp_sum=temp_sum+koeff;
        }
    akf[m]=temp_sum/k;
    
    fprintf (out_date,"%f\t%f\n", m*step,akf[m]/akf[0]);
  }
  
  fclose (out_date);
  
  for (i=2; i<=M; i=i+2)
  {  
    sum1=sum1+fabs(akf[i-1]);
  }
  
  for (i=2; i<=M-2; i=i+2)
  {  
    sum2=sum2+fabs(akf[i]);
  }
  
  Integral=step*(fabs(akf[0])+4*sum1+2*sum2+fabs(akf[M]))/3;
  printf("correlation time =%f\n", Integral/akf[0]);

  
  free(x);
  free(t);
  free(akf);

return(0);
} 
