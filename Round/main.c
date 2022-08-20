/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: Dr.Alex
 * Created on 10 ноября 2016 г., 13:30
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
float gamm(float z)
{
    int n,i;
    float G;
    n=1000;
    G=pow(n,z)/z;
    for (i=1;i<(n-1);i++)
    {
        G=G*i/(z+i);
    }
    return G;
}
float bess1(float x, int m)
{
    float J; int i;
    J=(m/2)*(powf(-1,0)/gamm(m+1))*pow(0.5*x,m-1);
    for (i=1; i<=100; i++)
    {
        J=J+((2*i+m)/2)*(powf(-1,i)/gamm(m+1+i))*pow(0.5*x,2*i+m-1);
    }
    return J;
}
float bess2(float x, int m)
{
    float J; int i;
    J=(powf(-1,0)/gamm(m+1))*pow(0.5*x,m);
    for (i=1; i<=100; i++)
    {
        J=J+(powf(-1,i)/gamm(m+1+i))*pow(0.5*x,2*i+m);
    }
    return J;
}
main()
{
    FILE* f=fopen("out1.dat","w"), *g=fopen("out2.dat","w");
    int m,n; float r, R, lk, l, E, H1, H2, w, k, q, z;
    m=0; n=1;
    lk=0.037; l=lk/2;
    k=1/sqrtf(powf(1/l,2)-powf(1/lk,2));
    w=300000000*k;
    R=lk/1.64;
    
    for (q=0; q<=2*M_PI; q=q+M_PI/360)
    {
    for (r=0.0001; r<=R; r=r+0.0001)
    {
        E=bess1(3.828*r/R, m);
        H2=bess2(3.828*r/R, m);
        fprintf(f,"%f %f %f %f\n",r*cos(q),r*sin(q),E,H2);
    }
    fprintf(f,"\n");
    }
    for (z=0; z<=k/2; z=z+k/2000)
    {
    for (r=-R; r<=R; r=r+0.0001)
    {
        E=bess1(3.828*r/R, m)*sin(-(2*M_PI/k)*z);
        H2=bess2(3.828*r/R, m)*cos(-(2*M_PI/k)*z);
        fprintf(g,"%f %f %f %f\n",z,r,E,H2);
    }
    fprintf(g,"\n");
    }
    fclose(f); fclose(g);
}


