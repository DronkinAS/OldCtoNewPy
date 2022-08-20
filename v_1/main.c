#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "control_linkage_backlash.h"
#include "control_linkage_backlash.c"

double* read_matrix(char file_name[], int N, int M) 
{
    int i,j;
    double *A;
    FILE *f;
    
    f = fopen(file_name,"r");
    
    A = (double*)malloc(N * M * sizeof(double));   
    
    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            fscanf(f,"%lf",(A + i*M + j));
        }
    }
    
    fclose(f);
    
    return A;
}

double* straight(double* A, double* B, double *b, int N, int M)
{
    int i, j;
    double *C, s = 0.0;
    
    C = (double*)malloc(N * sizeof(double));
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
        *(C + i) = 0.0;
    }

    for (i = 0 ; i <= N-1 ; i++)
    {
        for (j = 0 ; j <= M-1 ; j++)
        {
            *(C + i) += (*(A + j)) * (*(B + i*M + j));
        }
        s = (*(C + i)) + (*(b + i));
        *(C + i) = 2/(1 + exp(-1 * s)) - 1;
    }
    return C;
}

double* delta_out(double out, double error)
{
    int i;
    double *x;
    
    x = (double*)malloc(sizeof(double));
    
    if (( fabs(out) > 0.99 ) && (fabs(error) >= 0.0001))
    
    {
        out = out * 0.99;
    }

    //printf("%lf\n",out);
    *x = 0.5*error*(1 - pow(out,2));

    return x;
}

double* back(double* W, double* delta,
             double *out,
             int N, int M)
{
    int i,j;
    double *x, s;
    
    x = (double*)malloc(M * sizeof(double));
    
    for ( i = 0 ; i <= M-1 ; i++)
    {
        *(x + i) = 0.0;
    }
    
    for ( i = 0 ; i <= M-1 ; i++)
    {
        s = 0.0;
        
        for ( j = 0 ; j <= N-1 ; j++)
        {
            s += (*(delta + j)) * (*(W + i + j*M));
        }
        
        if (( fabs(*(out + i)) > 0.95 ) && ( fabs(s) >= 0.0001))
        {
            *(out + i) = *(out + i) * 0.95;
        }
        
        *(x + i) = 0.5*s*(1 - pow(*(out + i),2));
    }
        
    return x;
}

void correction_w(double E, double A,
        char file_name1[], char file_name2[],
        double* w, double* dw, double* d, double* out,
        int N, int M)
{
    int i, j;
    double grad_1, grad_2, grad;
    FILE *f, *g;
    
    f = fopen(file_name1,"w");
    g = fopen(file_name2,"w");
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
        for ( j = 0 ; j <= M-1 ; j++)
        {
            grad_1 = *(out + j);
            grad_2 = *(d + i);
            grad = grad_1 * grad_2;
            
            *(dw + i*M + j) = (E * grad) + (A * (*(dw + i*M + j)));
            *(w + i*M + j) += *(dw + i*M + j);
            
            fprintf(f,"%lf ",*(w + i*M + j));
            fprintf(g,"%lf ",*(dw + i*M + j));
            
        }
        fprintf(f,"\n");
        fprintf(g,"\n");
    }
    fclose(f);
    fclose(g);
}

void correction_b(double E, double A,
        char file_name1[], char file_name2[],
        double* b, double* db, double* d, int N)
{
    int i;
    double grad;
    FILE *f, *g;
    
    f = fopen(file_name1,"w");
    g = fopen(file_name2,"w");
    
    for ( i = 0 ; i <= N-1 ; i++)
    {
            grad = (*(d + i));
            *(db + i) = E*grad;
            *(db + i) += A * (*(db + i));
            *(b + i) += *(db + i);
 
            fprintf(f,"%lf ",*(b + i));
            fprintf(g,"%lf ",*(db +i));
    }
    
    fclose(f);
    fclose(g);
}

void converter_U(double A, double *U, double angle)
{
    double Ua, Ub, Uc;

    Ua = A * sin(angle);
    Ub = A * sin(angle + 4*M_PI/3);
    Uc = A * sin(angle + 2*M_PI/3);

    *U = Ua*cos(angle) + ((Ub - Uc)/sqrt(3.0))*sin(angle);
    *(U + 1) = ((Ub - Uc)/sqrt(3.0))*cos(angle) - Ua*sin(angle);
}

double calculation_Id(double L, double R, double p,
                      double U, double *I, double w)
{
    double Id, Iq, out;

    Id = *I;
    Iq = *(I + 1);
    out = (-1*R*Id)/L + p*w*Iq + U/L;

    return out;
}


double calculation_Iq(double L, double R, double p, double F,
                      double U, double *I, double w)
{
    double Id, Iq, out;

    Id = *I;
    Iq = *(I + 1);
    out = (-1*R*Iq)/L - p*w*Id + U/L - (p*w*F)/L;

    return out;
}

double calculation_Mn(double Mn_0, double Mn_1, double Mn_2,
                      double w, double kr, double nr)
{
    double Mn;

    Mn = Mn_1 + Mn_2;
    Mn = Mn/(kr*nr);
    Mn = Mn + Mn_0;

    if (w < 0)
    {
        Mn = -1*Mn;
    }

    return Mn;
}

double calculation_J(double J_0, double J_1, double J_2,
                     double kr, double nr)
{
    double J;

    J = J_1 + J_2;
    J = J/(kr*kr*nr);
    J = J + J_0;

    return J;
}

double PID(double *error, double dt)
{
    double Td = 0.04, Kp = 50.0, Ti = 0.1, signal;
    
    signal = *(error + 1) * Kp;
    signal += Td * ( *(error + 1) - *(error) ) / dt;
    signal += ( *error + *(error + 1)* dt ) / Ti;
    
    return signal;
}

int main(int argc, char** argv) {

    int i,j;

    double *in, *out, *W1, *W2, *W3, *W4, *W5, *b1, *b2, *b3, *b4, *b5,
           *L1, *L2, *L3, *L4, *delta_0, *delta_1, *delta_2, *delta_3, *delta_4,
           *dW1, *dW2, *dW3, *dW4, *dW5, *db1, *db2, *db3, *db4, *db5,
           t, dt = 0.0001, angle, rotation, delta_rotation,
           angle_out, rotation_out,
           angle_l, rotation_l,
           M, Mn, J, *U, *I, Id, Iq,
           ideal, *error, control_signal, n_out, k = 0.0, T = 0.5,
           Mn_0 = 0.00639, Mn_1 = 0.0, Mn_2 = 0.0,
           J_0 = 0.0000020, J_1 = 0.0000005, J_2 = 0.0,
           L = 0.000226, R = 1.39, F = 0.02052, p = 1.0,
           kr = 103.0, nr = 0.7;

    control_linkage_backlash_DATA lu;

    FILE *g;

    control_linkage_backlash_init (&lu);
    lu.freq = 10000.0;

    W1 = read_matrix("W1.dat",60,7);
    W2 = read_matrix("W2.dat",100,60);
    W3 = read_matrix("W3.dat",100,100);
    W4 = read_matrix("W4.dat",60,100);
    W5 = read_matrix("W5.dat",1,60);

    b1 = read_matrix("b1.dat",1,60);
    b2 = read_matrix("b2.dat",1,100);
    b3 = read_matrix("b3.dat",1,100);
    b4 = read_matrix("b4.dat",1,60);
    b5 = read_matrix("b5.dat",1,1);

    dW1 = read_matrix("dW1.dat",60,7);
    dW2 = read_matrix("dW2.dat",100,60);
    dW3 = read_matrix("dW3.dat",100,100);
    dW4 = read_matrix("dW4.dat",60,100);
    dW5 = read_matrix("dW5.dat",1,60);

    db1 = read_matrix("db1.dat",1,60);
    db2 = read_matrix("db2.dat",1,100);
    db3 = read_matrix("db3.dat",1,100);
    db4 = read_matrix("db4.dat",1,60);
    db5 = read_matrix("db5.dat",1,1);
    
    g = fopen("out.dat","w");

    in = (double*)malloc(7 * sizeof(double));
    I = (double*)malloc(2 * sizeof(double));
    U = (double*)malloc(2 * sizeof(double));
    error = (double*)malloc(2 * sizeof(double));

    for ( i = 0 ; i <= 6 ; i++)
    {
        *(in + i) = 0.0;
    }

    n_out = 0.0;

    angle = 0.0; rotation = 0.0; M = 0.0;
    angle_out = 0.0; rotation_out = 0.0;
    *I = 0.0; *(I + 1) = 0.0;
    *U = 0.0; *(U + 1) = 0.0;

    Mn = calculation_Mn(Mn_0,Mn_1,Mn_2,rotation,kr,nr);
    J = calculation_J(J_0,J_1,J_2,kr,nr);

    control_signal = 0.0;
    angle_l = 0.0;

    for ( t = 0.0 ; t <= 25.0 ; t = t + dt )
    {
        ideal = M_PI/6;
        if (t >= 2.5) ideal = M_PI/2;
        if ((t >= 5.0 + k*T) && (t <= 5.0 + (k+0.5)*T)) ideal = M_PI; 
        if ((t > 5.0 + (k+0.5)*T) && (t <= 5.0 + (k+1)*T)) ideal = 0.0; 
        if (t >= 15.0) ideal = sin(2*M_PI*t);
        if (t >= 20.0) ideal = 2*sin(M_PI*t);
        if (t >= 5.0 + (k+1)*T) k++;
        //ideal = 2*sin(M_PI*t);
        
        converter_U(n_out,U,angle);

        Id = calculation_Id(L,R,p,*U,I,rotation);
        Iq = calculation_Iq(L,R,p,F,*(U + 1),I,rotation);

        *I += Id*dt;
        *(I + 1) += Iq*dt;

        M = p*F * *(I+1);

        Mn = calculation_Mn(Mn_0,Mn_1,Mn_2,rotation,kr,nr);

        delta_rotation = ((M - Mn)*dt)/J;
        rotation += delta_rotation;

        angle += rotation*dt;

        if (angle >= (2.0 * M_PI))
          {
              angle += -2.0 * M_PI;
          }

        if (angle <= (-2.0 * M_PI))
          {
              angle += 2.0 * M_PI;
          }

        rotation_out = rotation/kr;
        angle_out += rotation_out * dt;

        if (angle_out >= (2.0 * M_PI))
          {
              angle_out += -2.0 * M_PI;
          }

        if (angle_out <= (-2.0 * M_PI))
          {
              angle_out += 2.0 * M_PI;
          }

        rotation_l = angle_l;

        lu.Vhod = angle_out;
        control_linkage_backlash_tick(&lu);

        angle_l = lu.Vyhod;
        rotation_l += angle_l * dt;

        if (angle_l >= (2.0 * M_PI))
        {
           angle_l += -2.0 * M_PI;
        }

        if (angle_l <= (-2.0 * M_PI))
        {
           angle_l += 2.0 * M_PI;
        }

        *in = control_signal;
        *(in + 1) = *(in + 2);
        *(in + 2) = *(in + 3);
        *(in + 3) = angle_l;
        *(in + 4) = *(in + 5);
        *(in + 5) = rotation_l;
        *(in + 6) = dt;

        L1 = straight(in,W1,b1,60,7);
        L2 = straight(L1,W2,b2,100,60);
        L3 = straight(L2,W3,b3,100,100);
        L4 = straight(L3,W4,b4,60,100);
        out = straight(L4,W5,b5,1,60);

        n_out = *out * (-24);

        *error = *(error + 1);
        *(error + 1) = ideal - angle_l;
          
        control_signal = PID(error,dt);
        control_signal = -1 * control_signal;
          
        if ( control_signal >= 24.0 )
        {
            control_signal = 24.0;
        }
          
        if ( control_signal <= -24.0 )
        {
            control_signal = -24.0;
        }

        if ( fabs(*(error + 1)) > 0.0001 )
        {
           delta_0 = delta_out(*out,*(error+1));
           printf("%lf %lf %lf %lf\n",t,*(out),*delta_0,*(error+1));
           delta_1 = back(W5,delta_0,L4,1,60);
           delta_2 = back(W4,delta_1,L3,60,100);
           delta_3 = back(W3,delta_2,L2,100,100);
           delta_4 = back(W2,delta_3,L1,100,60);

           correction_w (0.0005,0.05,"W1.dat","dW1.dat",W1,dW1,delta_4,in,60,7);
           correction_w (0.0005,0.05,"W2.dat","dW2.dat",W2,dW2,delta_3,L1,100,60);
           correction_w (0.0005,0.05,"W3.dat","dW3.dat",W3,dW3,delta_2,L2,100,100);
           correction_w (0.0005,0.05,"W4.dat","dW4.dat",W4,dW4,delta_1,L3,60,100);
           correction_w (0.0005,0.05,"W5.dat","dW5.dat",W5,dW5,delta_0,L4,1,60);
           
           correction_b (0.0005,0.05,"b1.dat","db1.dat",b1,db1,delta_4,60);
           correction_b (0.0005,0.05,"b2.dat","db2.dat",b2,db2,delta_3,100);
           correction_b (0.0005,0.05,"b3.dat","db3.dat",b3,db3,delta_2,100);
           correction_b (0.0005,0.05,"b4.dat","db4.dat",b4,db4,delta_1,60);
           correction_b (0.0005,0.05,"b5.dat","db5.dat",b5,db5,delta_0,1);

           free(W1);
           free(W2);
           free(W3);
           free(W4);
           free(W5);

           free(b1);
           free(b2);
           free(b3);
           free(b4);
           free(b5);

           free(delta_0);
           free(delta_1);
           free(delta_2);
           free(delta_3);
           free(delta_4);

           free(dW1);
           free(dW2);
           free(dW3);
           free(dW4);
           free(dW5);

           free(db1);
           free(db2);
           free(db3);
           free(db4);
           free(db5);

           W1 = read_matrix("W1.dat",60,7);
           W2 = read_matrix("W2.dat",100,60);
           W3 = read_matrix("W3.dat",100,100);
           W4 = read_matrix("W4.dat",60,100);
           W5 = read_matrix("W5.dat",1,60);

           b1 = read_matrix("b1.dat",1,60);
           b2 = read_matrix("b2.dat",1,100);
           b3 = read_matrix("b3.dat",1,100);
           b4 = read_matrix("b4.dat",1,60);
           b5 = read_matrix("b5.dat",1,1);

           dW1 = read_matrix("dW1.dat",60,7);
           dW2 = read_matrix("dW2.dat",100,60);
           dW3 = read_matrix("dW3.dat",100,100);
           dW4 = read_matrix("dW4.dat",60,100);
           dW5 = read_matrix("dW5.dat",1,60);

           db1 = read_matrix("db1.dat",1,60);
           db2 = read_matrix("db2.dat",1,100);
           db3 = read_matrix("db3.dat",1,100);
           db4 = read_matrix("db4.dat",1,60);
           db5 = read_matrix("db5.dat",1,1);

        }

        fprintf(g,"%lf %lf %lf %lf %lf %lf\n",t,angle_l,ideal,n_out,*(error+1),n_out);

        free(out);
        free(L1);
        free(L2);
        free(L3);
        free(L4);

    }

    free(I);
    free(U);
    free(in);
    fclose(g);

    return (EXIT_SUCCESS);
}

