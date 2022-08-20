/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: Dr.Alex
 *
 * Created on 5 сентября 2016 г., 22:00
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
float f1 (float y, float z)
{
    return -y-z;
}
float f2 (float x, float y)
{
    float a=0.15;
    return x+a*y;
}
float f3 (float x, float z)
{
    float b=0.2, c=8.0;
    return b+z*(x-c);
}
int main(int argc, char** argv) {
    FILE *f=fopen("out.dat","w");
    float k1[3], k2[3], k3[3], k4[3], x, y, z, t, h=0.001;
    t=0; x=0; y=0; z=0;
    for (t=0; t<=100; t=t+h)
    {
        k1[0]=h*f1(y,z);
        k1[1]=h*f2(x,y);
        k1[2]=h*f3(x,z);
        
        k2[0]=h*f1(y+0.5*k1[1],z+0.5*k1[2]);
        k2[1]=h*f2(x+0.5*k1[0],y+0.5*k1[1]);
        k2[2]=h*f3(x+0.5*k1[0],z+0.5*k1[2]);
        
        k3[0]=h*f1(y+0.5*k2[1],z+0.5*k2[2]);
        k3[1]=h*f2(x+0.5*k2[0],y+0.5*k2[1]);
        k3[2]=h*f3(x+0.5*k2[0],z+0.5*k2[2]);
        
        k4[0]=h*f1(y+k3[1],z+k3[2]);
        k4[1]=h*f2(x+k3[0],y+k3[1]);
        k4[2]=h*f3(x+k3[0],z+k3[2]);
        
        x=x+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
        y=y+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        z=z+(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        
        fprintf(f,"%f %f\n",y,z);
    }
    fclose(f);
    return (EXIT_SUCCESS);
}

