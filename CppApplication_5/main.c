/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 2 апреля 2016 г., 17:37
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    float t,y1,y2,y,m,n,h;
    t=0.0; h=M_PI/2;
    while (t<=(M_PI/2)+(M_PI/180))
    {
        m=cos(t);
        n=sin(t)*sin(t);
        n=n/2.25;
        n=sqrtf(1-(n));
        y=m*n*6;
        m=1.5*m;
        y1=1/(m+n);
        y1=y1*y1*0.5;
        m=m/1.5;
        n=1.5*n;
        y2=1/(m+n);
        y2=y2*y2*0.5;
        y=y*(y1+y2);
       printf("%f %f\n",t,y);
        t=t+(M_PI/180.0);
    }
    return (EXIT_SUCCESS);
}

