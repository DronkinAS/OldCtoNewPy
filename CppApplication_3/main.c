/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: Dr.Alex
 *
 * Created on 6 октября 2016 г., 18:14
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
float ak(int m, int n, float l)
{
 return l*sqrtf(4*n*n+m*m)/4;
}

float E1(int m, int n, float a, float b, float x, float y)
{
 return ((m*M_PI)/a)*cos((m*x*M_PI)/a)*sin((n*y*M_PI)/b);
}

float E2(int m, int n, float a, float b, float x, float y)
{
 return ((n*M_PI)/b)*sin((m*x*M_PI)/a)*cos((n*y*M_PI)/b);
}

float H1(int m, int n, float a, float b, float x, float y)
{
 return ((n*M_PI)/b)*sin((m*x*M_PI)/a)*cos((n*y*M_PI)/b);
}

float H2(int m, int n, float a, float b, float x, float y)
{
 return ((m*M_PI)/a)*cos((m*x*M_PI)/a)*sin((n*y*M_PI)/b);
}
int main(int argc, char** argv) {
FILE* f;
 float  a,b,l,x,y,s,E[2],H[2];
 int n,m;
 
 f=fopen("out.dat","w");
 
 x=0; y=0;
 s=0.0001;
 n=1; m=2;
 l=0.037;
 
 if (n>m)
 { 
  a=ak(m,n,l)+0.5*(ak(m+1,n,l)-ak(m,n,l));
 }
 else
 {
  a=ak(m,n,l)+0.5*(ak(m,n+1,l)-ak(m,n,l));
 }
 
 b=2*a;
 printf("%f %f\n",a,b);
 for (x=0;x<=b;x=x+s)
 {
  for (y=0;y<=a;y=y+s)
 {
  E[1]=E2(m,n,a,b,x,y); H[1]=H2(m,n,a,b,x,y);   
  E[0]=E1(m,n,a,b,x,y); H[0]=H1(m,n,a,b,x,y);
  fprintf(f,"%f %f %f %f\n",x,y,sqrtf(E[0]*E[0]+E[1]*E[1]),sqrtf(H[0]*H[0]+H[1]*H[1]));
  printf("%f %f \n", a, b);
  }
  fprintf(f,"\n");
 }
fclose(f);
    return (EXIT_SUCCESS);
}

