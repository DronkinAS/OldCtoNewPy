#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * 
 */
float Ex (float a, float x, float y, float z, int m, int n, int p)
{
    return cos(m*M_PI*x/a)*sin(n*M_PI*y/a)*sin(p*M_PI*z/a);
}
float Ey (float a, float x, float y, float z, int m, int n, int p)
{
    return -1*sin(m*M_PI*x/a)*cos(n*M_PI*y/a)*sin(p*M_PI*z/a);
}
float Hx (float a, float x, float y, float z, int m, int n, int p)
{
    return -1*sin(m*M_PI*x/a)*cos(n*M_PI*y/a)*cos(p*M_PI*z/a);
}
float Hy (float a, float x, float y, float z, int m, int n, int p)
{
    return -1*cos(m*M_PI*x/a)*sin(n*M_PI*y/a)*cos(p*M_PI*z/a);
}
float Hz (float a, float x, float y, float z, int m, int n, int p)
{
    return cos(m*M_PI*x/a)*cos(n*M_PI*y/a)*sin(p*M_PI*z/a);
}
int main(int argc, char** argv) {
    FILE *f1=fopen("out1.dat","w"), *f2=fopen("out2.dat","w"),
         *f3=fopen("out3.dat","w"), *f4=fopen("out4.dat","w"),
         *f5=fopen("out5.dat","w"), *f6=fopen("out6.dat","w");
    int m=1, n=1, p=1;
    float x, y, z, E[3], H[3], a1, a2, f, k1, k2, e;
    f=300000000; e=4;
    a1=sqrt(3)*150000000/f; a2=sqrt(3)*150000000/(f*sqrtf(e));
    k1=(M_PI/a1)*sqrtf(m*m+n*n);
    k2=(M_PI/a2)*sqrtf(m*m+n*n);
    printf("%f %f",k1,k2);
    for (x=0;x<=a1;x=x+0.01*a1)
    {
        for (y=0;y<=a1;y=y+0.01*a1)
        {
            for (z=0;z<=a1;z=z+0.01*a1)
            {
                E[0]=Ex(a1,x,y,z,m,n,p);
                E[1]=Ey(a1,x,y,z,m,n,p);
                H[0]=Hx(a1,x,y,z,m,n,p);
                H[1]=Hy(a1,x,y,z,m,n,p);
                H[2]=Hz(a1,x,y,z,m,n,p);
                if ((z<=(0.251*a1))&&(z>=(0.249*a1))) {fprintf(f1,"%f %f %f %f %f %f %f\n",x,y,E[0],E[1],H[0],H[1],H[2]);}
                if ((x<=(0.251*a1))&&(x>=(0.249*a1))) {fprintf(f2,"%f %f %f %f %f %f %f\n",y,z,E[0],E[1],H[0],H[1],H[2]);}
                if ((y<=(0.251*a1))&&(y>=(0.249*a1))) {fprintf(f3,"%f %f %f %f %f %f %f\n",x,z,E[0],E[1],H[0],H[1],H[2]);}
            }
            fprintf(f2,"\n"); 
        }
        fprintf(f1,"\n"); fprintf(f3,"\n");
    }
    for (x=0;x<=a2;x=x+0.01*a2)
    {
        for (y=0;y<=a2;y=y+0.01*a2)
        {
            for (z=0;z<=a2;z=z+0.01*a2)
            {
                E[0]=Ex(a2,x,y,z,m,n,p);
                E[1]=Ey(a2,x,y,z,m,n,p);
                H[0]=Hx(a2,x,y,z,m,n,p);
                H[1]=Hy(a2,x,y,z,m,n,p);
                H[2]=Hz(a2,x,y,z,m,n,p);
                if ((z<=(0.251*a2))&&(z>=(0.249*a2))) {fprintf(f4,"%f %f %f %f %f %f %f\n",x,y,E[0],E[1],H[0],H[1],H[2]);}
                if ((x<=(0.251*a2))&&(x>=(0.249*a2))) {fprintf(f5,"%f %f %f %f %f %f %f\n",y,z,E[0],E[1],H[0],H[1],H[2]);}
                if ((y<=(0.251*a2))&&(y>=(0.249*a2))) {fprintf(f6,"%f %f %f %f %f %f %f\n",x,z,E[0],E[1],H[0],H[1],H[2]);}
            }
            fprintf(f5,"\n"); 
        }
        fprintf(f4,"\n"); fprintf(f6,"\n");
    }
    fclose(f1); fclose(f2); fclose(f3);
    fclose(f4); fclose(f5); fclose(f6);
    return (EXIT_SUCCESS);
}

