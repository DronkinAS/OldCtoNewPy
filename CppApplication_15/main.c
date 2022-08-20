/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 16 мая 2016 г., 19:55
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
int main(int argc, char** argv) {
    float T1, T2, T3, R1, R2, R3, df,n,d,l;
    FILE *f=fopen("out.dat","w");
    R1=0.99; R2=0.8; R3=0.4;
    l=0.6; d=6; n=1.5;
    for (l=0.38;l<=0.74;l=l+0.001)
    {
        df=4*M_PI*n*d/l;
        T1=(1-R1)*(1-R1);
        T1=T1/(1+R1*R1-2*R1*cos(df));
        T2=(1-R2)*(1-R2);
        T2=T2/(1+R2*R2-2*R2*cos(df));
        T3=(1-R3)*(1-R3);
        T3=T3/(1+R3*R3-2*R3*cos(df));
        fprintf(f,"%0.4f %0.4f %0.4f %0.4f\n",l,T1,T2,T3);
    }
    fclose(f);
    return (EXIT_SUCCESS);
}

