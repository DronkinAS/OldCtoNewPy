/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 3 мая 2016 г., 17:52
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
int main(int argc, char** argv) {
    int N,n,i,j,*w,s;
    float max,min,h,*x,*a,k,t,s1,s2;
    FILE *f1= fopen("inp.dat","r"), *f2= fopen("out1.dat","w"), *f3= fopen("out2.dat","w");
    N=50; n=7; s1=0.0; s2=0.0;
    x = (float *)malloc(N*sizeof(float));
    w = (int *)malloc(n*sizeof(int));
    a = (float *)malloc(n*sizeof(float));
    fscanf(f1,"%f\n",&x[0]);
    //printf ("%.3f\n",x[0]);
    min=x[0]; max=x[0];
    for (i=1;i<=(N-1);i++)
    {
        fscanf(f1,"%f\n",&x[i]);
        //printf ("%.3f\n",x[i]);
        if (x[i]<min) min=x[i];
        if (x[i]>max) max=x[i];
    }
    //printf ("%.3f %.3f\n",min,max);
    h=((max-min)/n)*1000;
    s = (int) h; s = s + 1; h =(float) s/1000;
    //printf ("%.3f\n",h);
    a[0]=min; w[0]=0;
    for (i=1;i<=(n-1);i++)
    {
        a[i]=a[i-1]+h;
        w[i]=0;
    }
    for (i=0;i<=(n-1);i++)
    {
        for (j=0;j<=(N-1);j++)
        {
            if ((x[j]<a[i]+h)&&(x[j]>=a[i])) 
                w[i] = w[i] + 1;
        }
        //printf("%.3f %.3f %d\n",a[i],a[i]+h,w[i]);
    }
    for (i=0;i<=(n-1);i++) 
    {
        fprintf(f2,"%.3f %.3f %.3f %d\n",a[i],a[i]+h,a[i]+h/2,w[i]);
        t = (float) w[i]/h;
        s1=s1+(a[i]+h/2)*w[i];
       s2=s2+(a[i]+h/2)*(a[i]+h/2)*w[i];
        for (k=a[i];k<a[i]+h;k=k+(h/100)) 
        {
            fprintf(f3,"%.3f %.3f %.3f\n",k,a[i]+h,t);
        }    
    }
    printf("%.5f %f %f %.5f",s1/N,(s2/N)-(s1/N)*(s1/N),((s2/N)-(s1/N)*(s1/N))*(N/(N-1)),sqrtf((s2/N)-(s1/N)*(s1/N)));
    fclose (f1); fclose(f2);
    return (EXIT_SUCCESS);
}

