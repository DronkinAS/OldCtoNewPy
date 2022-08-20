/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 23 апреля 2016 г., 21:21
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float rel (float x, float v)
{
    float w=1, b=0.15;
    return (b*(1-v*v)*v)-(w*w*x);
}

int main(int argc, char** argv) {

    float t,dt,l,g,gtild,w;
    float *x, *k1, *k2, *k3, *k4;
    int n;
    FILE *fp = fopen("s1.dat", "w");
    n = 2;
    x = (float *)malloc(n*sizeof(float));
    k1 = (float *)malloc(n*sizeof(float));
    k2 = (float *)malloc(n*sizeof(float));
    k3 = (float *)malloc(n*sizeof(float));
    k4 = (float *)malloc(n*sizeof(float));
    x[0] = 0.0; x[1] = 1.0; 
    dt = 0.001;
    for (t = 0.0;t <= 500.0;t = (t + dt))
    {
        //printf ("%f %f %f %f %f %f\n",t,x[0],x[1],fabs(x[1]-x[0]),x[2],x[3]);
        //if (t > 100)
        fprintf (fp, "%.4f %.4f %.4f \n",t,x[0],x[1]);
                
                g=rel(x[0],x[1]);
                
                k1[0] = x[1];
                k1[1] = g;
                
                g=rel(x[0]+k1[0]*(0.5*dt),x[1]+k1[1]*(0.5*dt));
                
                k2[0] = g*(0.5*dt);
                k2[1] = g;
                
                g=rel(x[0]+k2[0]*(0.5*dt),x[1]+k2[1]*(0.5*dt));
                                           
                k3[0] = g*(0.5*dt);
                k3[1] = g;
                
                g=rel(x[0]+k3[0]*(dt),x[1]+k3[1]*(dt));
                
                k4[0] = g*(dt);
                k4[1] = g;
                
                x[0]= x[0] + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
                x[1]= x[1] + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
    
          
    }
    return (EXIT_SUCCESS);
}

