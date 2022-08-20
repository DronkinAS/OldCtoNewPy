/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.c
 * Author: dronk
 *
 * Created on 28 апреля 2016 г., 9:55
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
    float D=1.0,q=1.0,U,a=5.0,z,f=2000.0;
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
    float z,dt,g,t,m,dU,tr,i,j,xx,meanvel, fulltime = 100, L = 20000;
    float *x, *xtild, *k1, *k2, *k3, *k4;
    int n;
    FILE *fp = fopen("out.dat", "w");
    n = 4;
    x = (float *)malloc(n*sizeof(float));
    k1 = (float *)malloc(n*sizeof(float));
    k2 = (float *)malloc(n*sizeof(float));
    k3 = (float *)malloc(n*sizeof(float));
    k4 = (float *)malloc(n*sizeof(float));
    x[0] = 0.0; x[2] = 0.0;
    z=fabs(x[0]-x[1]);
    dt = 0.005; m = 1.0; k1 = 0.0; k2 = 0.0;
    for (i = -5.0; i <= 5.0; i = (i + 0.1))
    {
        for (j = -2.0; j <= 2.0; j = (j + 0.05))
        {
            x[1] = i; x[3] = j;
            x[0]=0; x[2]=1.0;
            meanvel = 0;
            for (t = 0.0;t <= fulltime;t = (t + dt))
            {
                if (fabs(x[0]-x[1])<=1e-04)
                {
                    //xx = x[2];
                    x[2] = -x[2];
                    x[3] = -x[3];
                }
                dU = -1*(1/m)*(pot(x[0],x[1],-1) - pot(x[0],x[1],L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr = -1*fric(x[2])*x[2];
                g =  dU + tr;
                k1[0] = x[2];
                k1[1] = x[3];
                k1[2] = g;
                k1[3] = -g;
                
                dU = -1*(1/m)*(pot(xtild[0],xtild[1],-1) - pot(xtild[0],xtild[1],L))*((x[0] - x[1])/fabs(x[0] - x[1]));
                tr = -1*fric(xtild[2])*xtild[2];
                g = dU + tr;
                
                k2[1] = x[1] + 0.5*dt*(x[3] + xtild[3]);
                k2[0] = x[0] + 0.5*dt*(x[2] + xtild[2]);
                k2[2] = x[2] + 0.5*dt*(g + gtild);
                k2[3] = x[3] - 0.5*dt*(g + gtild);
                
                if ((t>19.9)&&(t<20.0))
                {
                    k1=x[1]; k2=x[0];
                }
                if (t > 20)meanvel += x[2]+x[3];
            }
            meanvel = meanvel/(2*(fulltime - 20)/dt);
            if ((fabs(x[0]-x[1])<=5) && (fabs(meanvel)<=0.5))
            {
                fprintf(fp,"%.4f %.4f\n",i,j);
            }
            printf("%.4f %.4f %.4f %.4f\n",i,j,pot(x[0],x[1],-1),x[1]);
            k1=0.0; k2=0.0;
            xx=0.0;
        }
        
        
    }
    fclose(fp);
    return (EXIT_SUCCESS);
}

