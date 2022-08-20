/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.c
 * Author: dronk
 *
 * Created on 26 мая 2016 г., 19:20
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * 
 */
int main(int argc, char** argv) {
    float I,x,y,b=4;
    int i,j;
    FILE *f=fopen("out.dat","w");
    for (x=-10.0;x<=10.0;x=x+0.1)
    {
        for (y=-10;y<=10;y=y+0.1)
        {
            I=exp(-1*(x*x+y*y)/(b*b));
            fprintf(f,"%.10f ",I);
        }
        fprintf(f,"\n");
        
    }
    printf("%f %f",x,y);
    fclose(f);
    return (EXIT_SUCCESS);
}

