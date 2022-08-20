/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 19 апреля 2016 г., 21:09
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
int main(int argc, char** argv) {
    int k;
    float q,x;
    q=2;
    printf("%f %f\n",0.0,4.0);
    for (k=1;k<=20;k++)
    {
        x=(16/q)*fabs(sin(k*M_PI/q)/(k*M_PI/q));
        printf ("%d %f\n",k,x);
    }
    return (EXIT_SUCCESS);
}

