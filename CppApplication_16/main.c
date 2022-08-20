/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 26 мая 2016 г., 0:16
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
int main(int argc, char** argv) {
    double E=400.0, t=0.005, f0=5300000000000.0, f, U, k;
    int i;
    f=5299999999000; U=0;
    printf("%lf\n",f);
    FILE *g=fopen("out.dat","w");
    for (i =1 ; i <= 2000; i++)
    {
        k=f-f0;
        U=(sin(3.1415926535897932384626433832795*k*t));
        if (f==f0) U=1.0; else 
        U=U/(3.1415926535897932384626433832795*k*t);
        fprintf(g,"%.10f %.10f\n",f/10000000000,U);
        printf("%.10f %.10f\n",f/10000000000,U);
        f=f+1.0;
    }
    printf("%lf\n",f);
    fclose(g);
    return (EXIT_SUCCESS);
}

