#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "control_linkage_backlash.h"
#include "control_linkage_backlash.c"

void converter_U(double U0, double *U, double angle)
{
    double Ua, Ub, Uc;

    Ua = U0 * sin(angle);
    Ub = U0 * sin(angle + 4*M_PI/3);
    Uc = U0 * sin(angle + 2*M_PI/3);

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

double PID(double *error, double *integral, double dt)
{
    double Td = 0.1, Kp = 100.0, Ti = 0.4, signal;
    
    signal = *(error + 1);
    signal += Td * ( *(error + 1) - *(error) ) / dt; 
    *integral += *(error + 1) * dt;
    signal += *integral / Ti;
    signal = Kp * signal;
    
    return signal;
}

int main(int argc, char** argv) {

    int z = 0;
    double t, dt = 0.0001, i, angle, rotation, delta_rotation,
           angle_out, rotation_out, angle_l,
           M, Mn, J, *U, *I, Id, Iq,
           ideal, *error, integral, control_signal,
           Mn_0 = 0.00639, Mn_1 = 0.0, Mn_2 = 0.0,
           J_0 = 0.0000020, J_1 = 0.0000005, J_2 = 0.0,
           L = 0.000226, R = 1.39, F = 0.02052, p = 1.0,
           kr = 103.0, nr = 0.7,
           k = 0.0, T = 4.0,
           all = 0.0, max = 0.0;

    control_linkage_backlash_DATA lu;

    FILE *g;

    control_linkage_backlash_init (&lu);
    lu.freq = 10000.0;

    g = fopen("out_motor.dat","w");

    U = (double*)malloc(2 * sizeof(double));
    I = (double*)malloc(2 * sizeof(double));
    error = (double*)malloc(2 * sizeof(double));

    angle = 0.0; rotation = 0.0; M = 0.0;
    *I = 0.0; *(I + 1) = 0.0;
    angle_out = 0.0; rotation_out = 0.0;
    angle_l = 0.0;
    

    Mn = calculation_Mn(Mn_0,Mn_1,Mn_2,rotation,kr,nr);
    J = calculation_J(J_0,J_1,J_2,kr,nr);

    control_signal = 0.0;
    *error = 0.0; *(error + 1) = 0.0;
    integral = 0.0;
    
    for ( t = 0.0 ; t <= 0.4 ; t+=dt )
    {
          //ideal = M_PI/30;
          //if (t >= (k+1)*T) k++;
          //if ((t >= k*T) && (t < (k+0.5)*T)) ideal = 0; 
          //if ((t >= (k+0.5)*T) && (t < (k+1)*T)) ideal = (M_PI/60);
          ideal = (M_PI/60)*sin(5*M_PI*t);
        
          converter_U(control_signal,U,angle);

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
          
          /*angle_l[1] += rotation_out*dt;
          
          if (angle_l[1] >= (2.0 * M_PI))
          {
              angle_l[1] += -2.0 * M_PI;
          }

          if (angle_l[1] <= (-2.0 * M_PI))
          {
              angle_l[1] += 2.0 * M_PI;
          }
          
          if ( l[0] * l[1] < 0)
          {
             if ( fabs(angle_l[0] - angle_l[1]) >= luft )
             {
                 angle_l[1] = angle_l[0];
                 l[0] = l[1];
             }
          }
          else
          {
             angle_l[0] = angle_l[1];
             l[0] = l[1];

             if ((rotation_out) >= 0)
             {
                l[1] = 1.0;
             }
             else
             {
                l[1] = -1.0;
             }
          }*/

          lu.Vhod = angle_out;
          control_linkage_backlash_tick(&lu);

          angle_l = lu.Vyhod ;

          if (angle_l >= (2.0 * M_PI))
          {
              angle_l += -2.0 * M_PI;
          }

          if (angle_l <= (-2.0 * M_PI))
          {
              angle_l += 2.0 * M_PI;
          }

          printf("%lf %lf \n",t,angle_l);

          *error = *(error + 1);
          *(error + 1) = ideal - angle_l;
          
          if (fabs(*(error+1)) >= max) max = fabs(*(error+1));
          all += fabs(*(error + 1));
          
          control_signal = PID(error,&integral,dt);
          control_signal = -24 * control_signal;
          
          if ( control_signal >= 24.0 )
          {
              control_signal = 24.0;
          }
          
          if ( control_signal <= -24.0 )
          {
              control_signal = -24.0;
          }

          if (z % 100 == 0) fprintf(g,"%lf %lf %lf %lf %lf\n",t,angle_out,ideal,angle_l,control_signal);
          z++;
    }
    
    all = all/z;
    printf("%lf %lf",max,all);

    free (I);
    free (U);

    fclose(g);

    return (EXIT_SUCCESS);
}
