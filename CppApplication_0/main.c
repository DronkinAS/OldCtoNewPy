/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: Dr.Alex
 *
 * Created on 26 октября 2016 г., 15:39
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
float fun(float x, float y, float t, float f, float f0, float a, float U)
{
    return ((-a*y)+(x*(2*M_PI*f0)*(2*M_PI*f0)/(1-0.1*U*sin(2*M_PI*f*t))));
}
int main(int argc, char** argv)
{
    FILE  *f=fopen("out.dat","w");
    float x, y, t, h, T, R = 33, L = 0.018, C = 0.00000001, a, f1, f0, U;
    float s1[4], s2[4];
    a = R/L; f0 = 1/(2*M_PI*sqrtf(L*C)); U = R*sqrtf(C/L) + 0.5;
    x=0; y=0.1; h=0.01; T=100;
    for (t=0;t<=T;t=t+h)
    {
        fprintf(f,"%f ",t);
        for (f1=f0;f1<=f0;f1=f1+1)
        {
            s1[0]=y;
            s2[0]=fun(x,y,t,f1,f0,a,U);
            s1[1]=y+(0.5*h)*s2[0];
            s2[1]=fun(x+(0.5*h)*s1[0],y+(0.5*h)*s2[0],t+(0.5*h),f1,f0,a,U);
            s1[2]=y+(0.5*h)*s2[1];
            s2[2]=fun(x+(0.5*h)*s1[1],y+(0.5*h)*s2[1],t+(0.5*h),f1,f0,a,U);
            s1[3]=y+(h)*s2[2];
            s2[3]=fun(x+(h)*s1[2],y+(h)*s2[2],t+(h),f1,f0,a,U);
            x=x+h*(s1[0]+2*s1[1]+2*s1[2]+s1[3])/6;
            y=y+h*(s2[0]+2*s2[1]+2*s2[2]+s2[3])/6;
            fprintf(f,"%f %f %f ",f1,x,y);
        }
        fprintf(f,"\n");
    }
    return (EXIT_SUCCESS);
}

