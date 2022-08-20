
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
float spec(float w)
{
    return (0.06)/(powf(w,2)*powf(0.007,2)+powf((powf(1.5,2)-powf(w,2)),2));
}
int main(int argc, char** argv) {
    float x,h,f,w,y;
    FILE* g=fopen("out1.dat","w");
    x = 0;
    w = spec(1.5);
    h = 0.001;
    for (x=0 ; x<4.2 ; x=x+h)
    {
        y=(spec(x))/w;
        f=10*(log10f(y));
        fprintf(g,"%f %f\n",x/(2*M_PI),f);
    }
    return (EXIT_SUCCESS);
}

