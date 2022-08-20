/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 16 мая 2016 г., 17:42
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    float t,v,d,r,x,y,c=300000000.0, lc, l, dl, h, I;
    int i,j;
    FILE *f=fopen("out.dat","w");
    l=0.0000006; dl=0.00000001;
    lc=l*l/(dl); d=0;
    printf("%.4f\n",lc);
    h=0.00001; v=0.000000005;
    for (t=0;t<=1;t=(t+h))
    
    {   d=d+v*t;
        I=(M_PI)*d/lc;
        I=sin(I)/I;
        I=fabs(I);
        I=1+I*cos(2*(M_PI)*d/l);
        I=2*I;
        fprintf(f,"%.10f %.10f\n",t,I);
        printf("%.10f %.10f %.10f\n",x,y,I);
    }
    fclose(f);
    return (EXIT_SUCCESS);
}

