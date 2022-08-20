

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
    float x,h,f;
    FILE* g=fopen("out.dat","w");
    h=0.01;
    for (x=0;x<=3200;x=x+h)
    {
        f=exp(-1*0.0035*x)*cos(1.5*x);
        fprintf(g,"%f %f\n",x,f);
    }
    }
    return (EXIT_SUCCESS);
}

