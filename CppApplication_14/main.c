/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 16 мая 2016 г., 22:55
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    float f1,f2,b,l,sf,d,df,PI;
    float N=10;
    FILE *f=fopen("out.dat","w");
    l=0.55; b=1; d=4;
    PI=3.1415926535897932384626433832795;
    for (sf=-1;sf<=1;sf=sf+0.001)
    {
        f1=(PI*b*sf)/l;
        f1=pow((sin(f1)/f1) ,2);
        f2=(PI*d*sf)/l;
        f2=pow((sin(N*f2)/sin(f2)) ,2);
        f2=f1*f2;
        f1=N*N*f1;
        fprintf(f,"%.4f %.4f %.4f\n",sf,f1,f2);
    }
    fclose(f);
    return (EXIT_SUCCESS);
}

