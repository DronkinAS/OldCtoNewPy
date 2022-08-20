/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.c
 * Author: dronk
 *
 * Created on 10 апреля 2016 г., 22:52
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float fric (float v)
{
    float v0=1.0,j,g=1.0;
    j=-1*g*(1-(v*v)/(v0*v0));
    return j;
}

float pot (float x1, float x2, float L)
{
  double D=5.0,q=1.0,U,a=5.0,z,f=2000.0;
    //a=2*3.1415926535897932384626433832795*f*sqrtf(m/(2*D));
    if (L > 0)
        z=L - fabs(x1-x2);
    else
        z = fabs(x1-x2);
    U=exp(-a*(z-q));
    //U=2*D*a*(1-U)*U;
    return 2*D*a*(1-U)*U;
}
float sum (float u, float t)
{
    return u+t;
}
/*
 *
 */
int main(int argc, char** argv) {
    float z,dt,g1,g2,gtild2,gtild1,t,m,dU,tr1,tr2,L = 20000;
    float *x, *k1, *k2, *k3, *k4;
    int n;
    FILE *fp = fopen("k0.dat", "w");
    n = 4;
    x = (float *)malloc(n*sizeof(float));
    k1 = (float *)malloc(n*sizeof(float));
    k2 = (float *)malloc(n*sizeof(float));
    k3 = (float *)malloc(n*sizeof(float));
    k4 = (float *)malloc(n*sizeof(float));
    x[0] = 0.0; x[1] = 1; x[2] = 1; x[3] = -1;
    z=fabs(x[0]-x[1]);
    dt = 0.0001; m = 1.0;
    for (t = 0.0;t <= 100.0;t = (t + dt))
    {
        //printf ("%f %f %f %f %f %f\n",t,x[0],x[1],fabs(x[1]-x[0]),x[2],x[3]);
        if (t > 50)
        {fprintf (fp, "%.4f %.4f %.4f %.4f %.4f %.4f\n",t,x[0],x[1],fabs(x[1]-x[0]),x[2],x[3]);}
        if (fabs(x[0]-x[1])<=0.000001)
        {
            if (x[2]>=x[3])
            {
                x[3]=x[3]*(x[3]*x[2])/fabs(x[3]*x[2]);
                x[2]=-x[2];
            }        
            else 
            {
                x[2]=x[2]*(x[3]*x[2])/fabs(x[3]*x[2]);
                x[3]=-x[3];
            }          
        }
                dU = -1*(1/m)*(pot(x[0],x[1],-1) - pot(x[0],x[1],L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr1 = -1*fric(x[2])*x[2];
                tr2 = -1*fric(x[3])*x[3];
                g1 = dU + tr1;
                g2 = -dU + tr2;
                
                k1[0] = x[2];
                k1[1] = x[3];
                k1[2] = g1;
                k1[3] = g2;
                
                dU = -1*(1/m)*(pot(x[0]+k1[0]*(0.5*dt),x[1]+k1[1]*(0.5*dt),-1) - pot(x[0]+k1[0]*(0.5*dt),x[1]+k1[1]*(0.5*dt),L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr1 = -1*fric(x[2]+k1[2]*(0.5*dt))*(x[2]+k1[2]*(0.5*dt));
                tr2 = -1*fric(x[3]+k1[3]*(0.5*dt))*(x[3]+k1[3]*(0.5*dt));
                g1 = dU + tr1;
                g2 = -dU + tr2; 
                
                k2[0] = g1*(0.5*dt);
                k2[1] = g2*(0.5*dt);
                k2[2] = g1;
                k2[3] = g2;
                
                dU = -1*(1/m)*(pot(x[0]+k2[0]*(0.5*dt),x[1]+k2[1]*(0.5*dt),-1) - pot(x[0]+k2[0]*(0.5*dt),x[1]+k2[1]*(0.5*dt),L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr1 = -1*fric(x[2]+k2[2]*(0.5*dt))*(x[2]+k2[2]*(0.5*dt));
                tr2 = -1*fric(x[3]+k2[3]*(0.5*dt))*(x[3]+k2[3]*(0.5*dt));
                g1 = dU + tr1;
                g2 = -dU + tr2; 
                
                k3[0] = g1*(0.5*dt);
                k3[1] = g2*(0.5*dt);
                k3[2] = g1;
                k3[3] = g2;
                
                dU = -1*(1/m)*(pot(x[0]+k3[0]*(dt),x[1]+k3[1]*(dt),-1) - pot(x[0]+k3[0]*(dt),x[1]+k3[1]*(dt),L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr1 = -1*fric(x[2]+k3[2]*(dt))*(x[2]+k3[2]*(dt));
                tr2 = -1*fric(x[3]+k3[3]*(dt))*(x[3]+k3[3]*(dt));
                g1 = dU + tr1;
                g2 = -dU + tr2; 
                
                k4[0] = g1*(dt);
                k4[1] = g2*(dt);
                k4[2] = g1;
                k4[3] = g2;
                
                x[0]= x[0] + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
                x[1]= x[1] + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
                x[2]= x[2] + (dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
                x[3]= x[3] + (dt/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);
    }
    return (EXIT_SUCCESS);
}

