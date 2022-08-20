#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "control_linkage_backlash.h"
#include "control_linkage_backlash.c"
#include "neural_network_function.h"
#include "neural_network_function.c"
#include "electric_drive_model.h"
#include "electric_drive_model.c"


int main(int argc, char** argv) {

    int i,j,z = 0;
    
    double *in, *out, *W1, *W2, *W3, *W4, *W5, *b1, *b2, *b3, *b4, *b5,
           *L1, *L2, *L3, *L4, *delta_0, *delta_1, *delta_2, *delta_3, *delta_4,
           *dW1, *dW2, *dW3, *dW4, *dW5, *db1, *db2, *db3, *db4, *db5,
           t, dt = 0.0001, angle, rotation, delta_rotation,
           angle_out, rotation_out,
           angle_l, rotation_l,
           M, Mn, J, *U, *I, Id, Iq,
           ideal, *error, control_signal, n_out, k = 0.0, T = 8.0,
           Mn_0 = 0.00639, Mn_1 = 0.0, Mn_2 = 0.0,
           J_0 = 0.0000020, J_1 = 0.0000005, J_2 = 0.0,
           L = 0.000226, R = 1.39, F = 0.02052, p = 1.0,
           kr = 103.0, nr = 0.7,
           max = 0.0, all = 0.0;

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
    angle_l = 0.0; rotation_l = 0.0;
    
    *I = 0.0; *(I + 1) = 0.0;
    *U = 0.0; *(U + 1) = 0.0;

    Mn = calculation_Mn(Mn_0,Mn_1,Mn_2,rotation,kr,nr);
    J = calculation_J(J_0,J_1,J_2,kr,nr);

    control_signal = 0.0;
    
    z = 0;

    for ( t = 0.0 ; t <= 4.0 ; t = t + dt )
    {
        //ideal = M_PI/60;
        /*if (t >= 5.0) ideal = M_PI/90;
        if ((t >= 10.0 + k*T) && (t < 10.0 + (k+0.5)*T)) ideal = M_PI/90; 
        if ((t >= 10.0 + (k+0.5)*T) && (t < 10.0 + (k+1)*T)) ideal = 0.0; 
        if (t >= 20.0) ideal = (M_PI/6)*sin(2*M_PI*t);
        if (t >= 30.0) ideal = (M_PI/18)*sin(M_PI*t);
        if (t >= 10.0 + (k+1)*T) k++;*/
        //if (t >= (k+1)*T) k++;
        //if ((t >= k*T) && (t < (k+0.5)*T)) ideal = (M_PI/30)*sin(t - k*T); 
        //if ((t >= (k+0.5)*T) && (t < (k+1)*T)) ideal = M_PI/30;
        ideal = (M_PI/30)*sin(M_PI*t);
        
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
        
        *error = *(error + 1);
        *(error + 1) = ideal - angle_l;
        
        if (fabs(*(error+1)) >= max) max = fabs(*(error+1));
        all += fabs(*(error + 1));
        
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
        out = nnout(L4,W5,b5,1,60);

        n_out = *out * (-24);
        printf("%lf %lf %lf\n",t,*(out),*(error+1));

        if ( (z % 1 == 0) && (fabs(*(error + 1)) > 10.0) )
        {
           
           delta_0 = delta_out(*out,*(error+1));
           delta_1 = back(W5,delta_0,L4,1,60);
           delta_2 = back(W4,delta_1,L3,60,100);
           delta_3 = back(W3,delta_2,L2,100,100);
           delta_4 = back(W2,delta_3,L1,100,60);

           correction_w (0.01,0.2,W1,dW1,delta_4,in,60,7);
           correction_w (0.01,0.2,W2,dW2,delta_3,L1,100,60);
           correction_w (0.01,0.2,W3,dW3,delta_2,L2,100,100);
           correction_w (0.01,0.2,W4,dW4,delta_1,L3,60,100);
           correction_w (0.01,0.2,W5,dW5,delta_0,L4,1,60);
           
           correction_b (0.01,0.2,b1,db1,delta_4,60);
           correction_b (0.01,0.2,b2,db2,delta_3,100);
           correction_b (0.01,0.2,b3,db3,delta_2,100);
           correction_b (0.01,0.2,b4,db4,delta_1,60);
           correction_b (0.01,0.2,b5,db5,delta_0,1);
      
        }

        if (z % 10 == 0) fprintf(g,"%lf %lf %lf %lf %lf %lf\n",t,angle_l,ideal,n_out,*(error+1),n_out);
        z++;
        
        free(out);
        free(L1);
        free(L2);
        free(L3);
        free(L4);

    }
    
    all = all/z;
    printf("%lf %lf",max,all);
    
    write_matrix("W1.dat",W1,60,7);
    write_matrix("W2.dat",W2,100,60);
    write_matrix("W3.dat",W3,100,100);
    write_matrix("W4.dat",W4,60,100);
    write_matrix("W5.dat",W5,1,60);

    write_matrix("b1.dat",b1,1,60);
    write_matrix("b2.dat",b2,1,100);
    write_matrix("b3.dat",b3,1,100);
    write_matrix("b4.dat",b4,1,60);
    write_matrix("b5.dat",b5,1,1);

    write_matrix("dW1.dat",dW1,60,7);
    write_matrix("dW2.dat",dW2,100,60);
    write_matrix("dW3.dat",dW3,100,100);
    write_matrix("dW4.dat",dW4,60,100);
    write_matrix("dW5.dat",dW5,1,60);

    write_matrix("db1.dat",db1,1,60);
    write_matrix("db2.dat",db2,1,100);
    write_matrix("db3.dat",db3,1,100);
    write_matrix("db4.dat",db4,1,60);
    write_matrix("db5.dat",db5,1,1);
    
    free(I);
    free(U);
    free(in);
    fclose(g);
    
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

    return (EXIT_SUCCESS);
}

