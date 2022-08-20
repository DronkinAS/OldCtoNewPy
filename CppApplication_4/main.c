/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 5 марта 2016 г., 14:45
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    float t,y1,y2,m,n,h;
    t=0.0; h=0;
    while (t<=((M_PI*0.5)+M_PI/180))
    {
        m=cos(t);
        m=m*1.5;
        n=sin(t)*sin(t);
        n=n/2.25;
        n=sqrtf(1-(n));
        n=n;
        y1=(m-n)/(m+n);
        y1=y1*y1*cos(h)*cos(h);
        m=cos(t);
        m=m;
        n=sin(t)*sin(t);
        n=n/2.25;
        n=sqrtf(1-(n));
        n=1.5*n;
        y2=(m-n)/(m+n);
        y2=y2*y2*sin(h)*sin(h)+y1;
       printf("%f %f\n",t,y2);
        t=t+(M_PI/180);
    }
    return (EXIT_SUCCESS);
}

