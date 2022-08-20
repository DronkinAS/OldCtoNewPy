

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    double F1(double x1, double x2, double y1, double y2)
    {
        double s = 0.8, r, F;
        r = pow((x2-x1),2) + pow((y2-y1),2);
        r = sqrt(r) - 1;
        if ((r < 2)&&(r > -1))
        {
            F = (1*s*r);
        } 
        else
        {
            F = 0;
        }
        return F;
    }
    
    double F2(double x, double y, double t)
    {
        double s = 1, r, A = 0.5, f = 1, F;
        y = y - A*sin(2*M_PI*f*t);
        r = pow((1-y),2) + pow(x,2);
        r = sqrt(r) - 1;
        if ((r < 2)&&(r > -1))
        {
            F = (1*s*r);
        } 
        else
        {
            F = 0;
        }
        return F;
    }
    
    double F3(double x1, double x2, double y1, double y2)
    {
        double s = 0.8, r, F;
        r = pow((x2 - x1),2) + pow((y2 - y1),2);
        r = sqrt(r) - 1;
        if ((r < 2)&&(r > -1))
        {
            F = (-1*s*r);
        } 
        else
        {
            F = 0;
        }
        return F;
    }
    
    double F4(double x, double y, double t)
    {
        double s = 1, r, A = 0.5, f = 1, F;
        y = y - sin(2*M_PI*f*t + M_PI);
        r = pow((1-y),2) + pow((x-1),2);
        r = sqrt(r) - 1;
        if ((r < 2)&&(r > -1))
        {
            F = (1*s*r);
        } 
        else
        {
            F = 0;
        }
        return F;
    }
    
    int i = 0;
    double t = 0.0, h = 0.001;
    double x1 = 0, y1 = 0, x2 = 1, y2 = 0;
    double a1 = 0, a2 = 0, b1 = 0, b2 = 0;
    double k[8][4];
    
    FILE *out = fopen("out.dat","w"), *out1 = fopen("out1.dat","w"), *out2 = fopen("out2.dat","w"); 
        
    for (t = 0; t <= 10000; t = t + h)
    {
        
        
        k[0][0] = F1(x1,x2,y1,y2);
        k[1][0] = F2(x1,y1,t);
        k[2][0] = a1;
        k[3][0] = b1;
        k[4][0] = F3(x1,x2,y1,y2);
        k[5][0] = F4(x2,y2,t);
        k[6][0] = a2;
        k[7][0] = b2;
        
        k[0][1] = F1(x1 + (h/2)*k[2][0],x2 + (h/2)*k[6][0],y1 + (h/2)*k[3][0],y2 + (h/2)*k[7][0]);
        k[1][1] = F2(x1 + (h/2)*k[2][0],y1 + (h/2)*k[3][0],t + (h/2));
        k[2][1] = a1 + (h/2)*k[0][0];
        k[3][1] = b1 + (h/2)*k[1][0];;
        k[4][1] = F3(x1 + (h/2)*k[2][0],x2 + (h/2)*k[6][0],y1 + (h/2)*k[3][0],y2 + (h/2)*k[7][0]);
        k[5][1] = F4(x2 + (h/2)*k[6][0],y2 + (h/2)*k[7][0],t + (h/2));
        k[6][1] = a2 + (h/2)*k[4][0];
        k[7][1] = b2 + (h/2)*k[5][0];
        
        k[0][2] = F1(x1 + (h/2)*k[2][1],x2 + (h/2)*k[6][1],y1 + (h/2)*k[3][1],y2 + (h/2)*k[7][1]);
        k[1][2] = F2(x1 + (h/2)*k[2][1],y1 + (h/2)*k[3][1],t + (h/2));
        k[2][2] = a1 + (h/2)*k[0][1];
        k[3][2] = b1 + (h/2)*k[1][1];
        k[4][2] = F3(x1 + (h/2)*k[2][1],x2 + (h/2)*k[6][1],y1 + (h/2)*k[3][1],y2 + (h/2)*k[7][1]);
        k[5][2] = F4(x2 + (h/2)*k[6][1],y2 + (h/2)*k[7][1],t + (h/2));
        k[6][2] = a2 + (h/2)*k[4][1];
        k[7][2] = b2 + (h/2)*k[5][1];
        
        k[0][3] = F1(x1 + (h)*k[2][2],x2 + (h)*k[6][2],y1 + (h)*k[3][2],y2 + (h)*k[7][2]);
        k[1][3] = F2(x1 + (h)*k[2][2],y1 + (h)*k[3][2],t + (h));
        k[2][3] = a1 + (h)*k[0][2];
        k[3][3] = b1 + (h)*k[1][2];;
        k[4][3] = F3(x1 + (h)*k[2][2],x2 + (h)*k[6][2],y1 + (h)*k[3][2],y2 + (h)*k[7][2]);
        k[5][3] = F4(x2 + (h)*k[6][2],y2 + (h)*k[7][2],t + (h));
        k[6][3] = a2 + (h)*k[4][2];
        k[7][3] = b2 + (h)*k[5][2];
        
        a1 = a1 + (h/6)*(k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3]);
        b1 = b1 + (h/6)*(k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3]);
        x1 = x1 + (h/6)*(k[2][0] + 2*k[2][1] + 2*k[2][2] + k[2][3]);
        y1 = y1 + (h/6)*(k[3][0] + 2*k[3][1] + 2*k[3][2] + k[3][3]);
        a2 = a2 + (h/6)*(k[4][0] + 2*k[4][1] + 2*k[4][2] + k[4][3]);
        b2 = b2 + (h/6)*(k[5][0] + 2*k[5][1] + 2*k[5][2] + k[5][3]);
        x2 = x2 + (h/6)*(k[6][0] + 2*k[6][1] + 2*k[6][2] + k[6][3]);
        y2 = y2 + (h/6)*(k[7][0] + 2*k[7][1] + 2*k[7][2] + k[7][3]);
        
        if (t>=9600)   
        {
            fprintf(out,"%f %f %f\n",t,sqrt((x1*x1+y1*y1)),sqrt((x2*x2+y2*y2)));
            fprintf(out1,"%f %f %f %f %f\n",t,x1,a1,y1,b1);
            fprintf(out2,"%f %f %f %f %f\n",t,x2,a2,y2,b2);
        }
        i = i + 1;
    }
    
    fclose(out);
    fclose(out1);
    fclose(out2);
    
    return (EXIT_SUCCESS);
}

