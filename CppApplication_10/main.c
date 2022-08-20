/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 1 мая 2016 г., 16:01
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    float a,a0,n1,n2,b1,b2,b,m,n;
    FILE *fp = fopen("out.dat","w");
    n2=1; n1=1.6;
    a0 = asin(n2/n1);
    for (a = a0; a<=M_PI/2; a=(a + M_PI/360))
    {
       m=sin(a);
       n=cos(a);
       b1=n1*n1*m*m;
       b1=b1-n2*n2;
       b1=sqrtf(b1);
       b1=b1/(n1*n);
       b1=2*atan(b1);
       b2=n1*n1*m*m;
       b2=b2-n2*n2;
       b2=sqrtf(b2);
       b2=n1*b2/(n2*n2*n);
       b2=2*atan(b2);
       b=n1*n1*m*m;
       b=b-n2*n2;
       b=sqrtf(b);
       b=b*n/(n1*m*m);
       b=2*atan(b);
       fprintf(fp,"%.4f %.4f %.4f %.4f\n",a,b1,b2,b);
    }
    fclose(fp);
    return (EXIT_SUCCESS);
}

